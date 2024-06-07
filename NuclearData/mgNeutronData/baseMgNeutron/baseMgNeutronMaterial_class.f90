module baseMgNeutronMaterial_class

  use numPrecision
  use endfConstants
  use genericProcedures, only : fatalError, numToChar
  use RNG_class,         only : RNG
  use dictionary_class,  only : dictionary
  use dictDeck_class,    only : dictDeck

  ! Nuclear Data Interfaces
  use materialHandle_inter,    only : materialHandle
  use mgNeutronMaterial_inter, only : mgNeutronMaterial, kill_super => kill
  use neutronXSPackages_class, only : neutronMacroXSs

  ! Reaction objects
  use reactionMG_inter,        only : reactionMG
  use fissionMG_class,         only : fissionMG
  use multiScatterMG_class,    only : multiScatterMG
  use multiScatterP1MG_class,  only : multiScatterP1MG

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: baseMgNeutronMaterial_TptrCast
  public :: baseMgNeutronMaterial_CptrCast

  ! Public data location parameters
  ! Use them if accessing data entries directly
  integer(shortInt), parameter, public :: TOTAL_XS      = 1
  integer(shortInt), parameter, public :: IESCATTER_XS  = 2
  integer(shortInt), parameter, public :: CAPTURE_XS    = 3
  integer(shortInt), parameter, public :: FISSION_XS    = 4
  integer(shortInt), parameter, public :: NU_FISSION    = 5
  integer(shortInt), parameter, public :: NU            = 6

  !!
  !! Basic type of MG material data
  !!
  !! Stores MG data in a table.
  !! Fission is treated as a seperate reaction
  !! All other scattering reactions are lumped into single multiplicative scattering,
  !! which is stored as INELASTIC scatering in macroXSs package! After all it is inelastic in
  !! the sense that outgoing group can change. Diffrent types of multiplicative scattering can be
  !! build. See doc of "init" procedure for details.
  !!
  !! Public members:
  !!   data -> Rank 2 array with all XSs data
  !!
  !! Interface:
  !!   materialHandle interface
  !!   mgNeutronMaterial interface
  !!   init -> initialise Basic MG Material from dictionary and config keyword
  !!   nGroups -> returns number of energy groups
  !!
  !! Note:
  !!   Order of "data" array is: data(XS_type, Group #)
  !!   Dictionary with data must contain following entries:
  !!     -> numberOfGroups
  !!     -> capture [nGx1]
  !!     -> scatteringMultiplicity [nGxnG]
  !!     -> P0 [nGxnG]
  !!   Optional entries:
  !!     -> fission [nGx1]
  !!     -> nu [nGx1]
  !!     -> chi [nGx1]
  !!     -> P# [nGxnG]
  !!
  type, public, extends(mgNeutronMaterial) :: baseMgNeutronMaterial
    real(defReal),dimension(:,:), allocatable :: data
    real(defReal),dimension(:,:), allocatable :: data_sigma
    class(multiScatterMG), allocatable        :: scatter
    class(multiScatterMG), allocatable        :: scatter_sigma
    type(fissionMG), allocatable              :: fission

  contains
    ! Superclass procedures
    procedure :: kill
    procedure :: getMacroXSs_byG
    procedure :: getTotalXS_det
    procedure :: getMacroXSs_uncertain
    procedure :: getTotalXS_uncertain

    ! Local procedures
    procedure :: init
    procedure :: nGroups

  end type baseMgNeutronMaterial

contains

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(baseMgNeutronMaterial), intent(inout) :: self

    ! Call superclass procedure
    call kill_super(self)

    ! Kill local content
    if(allocated(self % data))    deallocate(self % data)
    if(allocated(self % data_sigma))    deallocate(self % data_sigma)
    if(allocated(self % scatter)) deallocate(self % scatter)
    if(allocated(self % fission)) deallocate(self % fission)

  end subroutine kill

  !!
  !! Load Macroscopic XSs into the provided package for a given group index G
  !!
  !! See mgNeutronMaterial documentation for more details
  !!
  subroutine getMacroXSs_byG(self, xss, G, rand)
    class(baseMgNeutronMaterial), intent(in) :: self
    type(neutronMacroXSs), intent(out)       :: xss
    integer(shortInt), intent(in)            :: G
    class(RNG), intent(inout)                :: rand
    character(100), parameter :: Here = ' getMacroXSs (baseMgNeutronMaterial_class.f90)'

    ! Verify bounds
    if(G < 1 .or. self % nGroups() < G) then
      call fatalError(Here,'Invalid group number: '//numToChar(G)// &
                           ' Data has only: ' // numToChar(self % nGroups()))
    end if

    ! Get XSs
    xss % total            = self % data(TOTAL_XS, G)
    xss % elasticScatter   = ZERO
    xss % inelasticScatter = self % data(IESCATTER_XS, G)
    xss % capture          = self % data(CAPTURE_XS, G)

    if(self % isFissile()) then
      xss % fission        = self % data(FISSION_XS, G)
      xss % nuFission      = self % data(NU_FISSION, G)
    else
      xss % fission        = ZERO
      xss % nuFission      = ZERO
    end if

  end subroutine getMacroXSs_byG

  !!
  !! Load uncertain Macroscopic XSs into the provided package for a given group index G
  !!
  !! See mgNeutronMaterial documentation for more details
  !!
  subroutine getMacroXSs_uncertain(self, xss, G, rand, X)
    class(baseMgNeutronMaterial), intent(in)    :: self
    type(neutronMacroXSs), intent(out)      :: xss
    integer(shortInt), intent(in)           :: G
    class(RNG), intent(inout)               :: rand
    real(defReal), dimension(:,:), intent(in) :: X
    character(100), parameter :: Here = ' getMacroXSs_uncertain (baseMgNeutronMaterial_class.f90)'


    ! Verify bounds
    if(G < 1 .or. self % nGroups() < G) then
      call fatalError(Here,'Invalid group number: '//numToChar(G)// &
                           ' Data has only: ' // numToChar(self % nGroups()))
    end if

    ! Get XSs
    xss % elasticScatter   = ZERO
    ! xss % inelasticScatter = self % data(IESCATTER_XS, G) + X(1) * self % data_sigma(IESCATTER_XS, G)
    xss % inelasticScatter = self % data(IESCATTER_XS, G)
    xss % capture          = self % data(CAPTURE_XS, G) + X(G, 1) * self % data_sigma(CAPTURE_XS, G)

    ! Careful about nu : do not use fissionMG % release when using uncertain nu
    if(self % isFissile()) then
      xss % fission        = self % data(FISSION_XS, G) + X(G, 2) * self % data_sigma(FISSION_XS, G)
      xss % nu             = self % data(NU, G) + X(G, 3) * self % data_sigma(NU, G)
      xss % nu             = self % data(NU, G)
      xss % nuFission      = xss % fission * xss % nu
    else
      xss % fission        = ZERO
      xss % nu             = ZERO
      xss % nuFission      = ZERO
    end if

    xss % total            = xss % elasticScatter + xss % inelasticScatter + xss % capture + xss % fission
  end subroutine getMacroXSs_uncertain

    !!
  !! Return uncertain Total XSs for energy group G
  !!
  !! See mgNeutronMaterial documentationfor details
  !!
  function getTotalXS_uncertain(self, G, rand, X) result(xs)
    class(baseMgNeutronMaterial), intent(in) :: self
    integer(shortInt), intent(in)            :: G
    class(RNG), intent(inout)                :: rand
    real(defReal), dimension(:,:), intent(in)  :: X
    real(defReal)                            :: xs
    character(100), parameter :: Here = ' getTotalXs_uncertain (baseMgNeutronMaterial_class.f90)'

    ! Verify bounds
    if (G < 1 .or. self % nGroups() < G) then
      call fatalError(Here,'Invalid group number: '//numToChar(G)// &
                           ' Data has only: ' // numToChar(self % nGroups()))
      xs = ZERO ! Avoid warning
    end if
    
    xs = ZERO
    ! xs = xs + self % data(IESCATTER_XS, G) + X(1) * self % data_sigma(IESCATTER_XS, G)
    xs = xs + self % data(IESCATTER_XS, G)
    xs = xs + self % data(CAPTURE_XS, G) + X(G, 2) * self % data_sigma(CAPTURE_XS, G)
    if (self % isFissile()) xs = xs + self % data(FISSION_XS, G) + X(G, 3) * self % data_sigma(FISSION_XS, G)

  end function getTotalXS_uncertain

  !!
  !! Return Total XSs for energy group G
  !!
  !! See mgNeutronMaterial documentationfor details
  !!
  function getTotalXS_det(self, G, rand) result(xs)
    class(baseMgNeutronMaterial), intent(in) :: self
    integer(shortInt), intent(in)            :: G
    class(RNG), intent(inout)                :: rand
    real(defReal)                            :: xs
    character(100), parameter :: Here = ' getTotalXS_det (baseMgNeutronMaterial_class.f90)'

    ! Verify bounds
    if (G < 1 .or. self % nGroups() < G) then
      call fatalError(Here,'Invalid group number: '//numToChar(G)// &
                           ' Data has only: ' // numToChar(self % nGroups()))
      xs = ZERO ! Avoid warning
    end if
    xs = self % data(TOTAL_XS, G)

  end function getTotalXS_det


  !!
  !! Initialise Base MG Neutron Material fromdictionary
  !!
  !! Args:
  !!   dict       [in] -> Input dictionary with all required XSs
  !!   scatterKey [in] -> String with keyword to choose approperiate multiplicative scatering
  !!                        type
  !! Errors:
  !!   FatalError if scatteKey is invalid
  !!   FatalError if data in dictionary is invalid (inconsistant # of groups;
  !!     -ve entries in P0 XSs)
  !!
  !! Note:
  !!   Some time in the future scattering MG reaction objects will have factory. For now
  !!   the factory is hardcoded into this procedure. Not the best solution but is fine at this
  !!   stage. The following scatterKey are supported:
  !!     -> P0
  !!     -> P1
  !!
  subroutine init(self, dict, scatterKey)
    class(baseMgNeutronMaterial), intent(inout) :: self
    class(dictionary),target, intent(in)        :: dict
    character(nameLen), intent(in)              :: scatterKey
    integer(shortInt)                           :: nG, N, i
    real(defReal), dimension(:), allocatable    :: temp
    type(dictDeck)                              :: deck
    character(100), parameter :: Here = 'init (baseMgNeutronMaterial_class.f90)'


    ! Read number of groups
    call dict % get(nG, 'numberOfGroups')
    if(nG < 1) call fatalError(Here,'Number of groups is invalid' // numToChar(nG))

    ! Set fissile flag
    call self % set(fissile = dict % isPresent('fission'))

    ! Build scattering reaction
    ! Prepare input deck
    deck % dict => dict

    ! Choose Scattering type
    select case(scatterKey)
      case ('P0')
        allocate( multiScatterMG :: self % scatter)

      case ('P1')
        allocate( multiScatterP1MG :: self % scatter)

      case default
        call fatalError(Here,'scatterKey: '//trim(scatterKey)//'is wrong. Must be P0 or P1')

    end select

    ! Initialise
    call self % scatter % init(deck, macroAllScatter)

    ! Deal with fission
    if(self % isFissile()) allocate(self % fission)
    if(self % isFissile()) call self % fission % init(deck, macroFission)

    ! Allocate space for data
    if(self % isFissile()) then
      N = 6
    else
      N = 3
    end if

    allocate(self % data(N, nG))
    allocate(self % data_sigma(N, nG))

    ! Load cross sections
    call dict % get(temp, 'capture')
    if(size(temp) /= nG) then
      call fatalError(Here,'Capture XSs have wong size. Must be: ' &
                          // numToChar(nG)//' is '//numToChar(size(temp)))
    end if
    self % data(CAPTURE_XS,:) = temp

    ! Extract values of scattering XS
    if(size(self % scatter % scatterXSs) /= nG) then
      call fatalError(Here, 'Somthing went wrong. Inconsistant # of groups in material and reaction&
                            &. Clearly programming error.')
    end if
    self % data(IESCATTER_XS,:) = self % scatter % scatterXSs

    ! Load Fission-data
    if( self % isFissile()) then
      ! Load Fission
      call dict % get(temp, 'fission')
      if(size(temp) /= nG) then
        call fatalError(Here,'Fission XSs have wong size. Must be: ' &
                            // numToChar(nG)//' is '//numToChar(size(temp)))
      end if
      self % data(FISSION_XS,:) = temp

      ! Calculate nuFission
      call dict % get(temp, 'nu')
      if(size(temp) /= nG) then
        call fatalError(Here,'Nu vector has wong size. Must be: ' &
                            // numToChar(nG)//' is '//numToChar(size(temp)))
      end if
      self % data(NU_FISSION,:) = temp * self % data(FISSION_XS,:)
      self % data(NU, :) = temp 
    end if

    ! Calculate total XS
    do i =1,nG
      self % data(TOTAL_XS, i) = self % data(IESCATTER_XS, i) + self % data(CAPTURE_XS, i)
      if(self % isFissile()) then
        self % data(TOTAL_XS, i) = self % data(TOTAL_XS, i) + self % data(FISSION_XS, i)
      end if
    end do

    ! Load uncertainty on cross sections (Careful, in practice ignore IESCATTERING and NU for now)
    if (dict % isPresent('capture_sigma')) then
      call dict % get(temp, 'capture_sigma')
      if(size(temp) /= nG) then
        call fatalError(Here,'Capture XSs uncertainties have wong size. Must be: ' &
                          // numToChar(nG)//' is '//numToChar(size(temp)))
      end if
      self % data_sigma(CAPTURE_XS,:) = temp
      self % data_sigma(TOTAL_XS, :) = self % data_sigma(TOTAL_XS, :) + temp
    else
      self % data_sigma(CAPTURE_XS,:) = ZERO
    end if

    ! Get uncertainty on scattering XS (as the column-wise sum of the uncertainties in P0...)
    if (dict % isPresent('scatter_sigma')) then 
      self % data_sigma(IESCATTER_XS,:) = self % scatter % scatterXSs_sigma
      self % data_sigma(TOTAL_XS, :) = self % data_sigma(TOTAL_XS, :) + temp
    else
      self % data_sigma(IESCATTER_XS,:) = ZERO
    end if

    ! Load Fission-data uncertainties
    if( self % isFissile()) then
      if (dict % isPresent('fission_sigma')) then 
        ! Load Fission
        call dict % get(temp, 'fission_sigma')
        if(size(temp) /= nG) then
          call fatalError(Here,'Fission XSs uncertainties have wong size. Must be: ' &
                            // numToChar(nG)//' is '//numToChar(size(temp)))
        end if
        self % data_sigma(FISSION_XS,:) = temp
        self % data_sigma(TOTAL_XS, :) = self % data_sigma(TOTAL_XS, :) + temp
      else
        self % data_sigma(FISSION_XS,:) = ZERO
      end if

      ! Load nu uncertainties
      if (dict % isPresent('nu_sigma')) then 
        call dict % get(temp, 'nu_sigma')
        if(size(temp) /= nG) then
          call fatalError(Here,'Nu uncertainty vector has wong size. Must be: ' &
                            // numToChar(nG)//' is '//numToChar(size(temp)))
        end if
        self % data_sigma(NU,:) = temp
      else
        self % data_sigma(NU,:) = ZERO
      end if
    end if
  end subroutine init

  !!
  !! Return number of energy groups
  !!
  !! Args:
  !!   None
  !!
  !! Errors:
  !!   None
  !!
  pure function nGroups(self) result(nG)
    class(baseMgNeutronMaterial), intent(in) :: self
    integer(shortInt)                        :: nG

    if(allocated(self % data)) then
      nG = size(self % data,2)
    else
      nG = 0
    end if

  end function nGroups

  !!
  !! Cast materialHandle pointer to baseMgNeutronMaterial type pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class materialHandle
  !!
  !! Result:
  !!   Null if source is not of baseMgNeutronMaterial type
  !!   Target points to source if source is baseMgNeutronMaterialtype
  !!
  pure function baseMgNeutronMaterial_TptrCast(source) result(ptr)
    class(materialHandle), pointer, intent(in) :: source
    type(baseMgNeutronMaterial), pointer           :: ptr

    select type(source)
      type is(baseMgNeutronMaterial)
        ptr => source

      class default
        ptr => null()
    end select

  end function baseMgNeutronMaterial_TptrCast

  !!
  !! Cast materialHandle pointer to baseMgNeutronMaterial class pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class materialHandle
  !!
  !! Result:
  !!   Null if source is not of baseMgNeutronMaterial class
  !!   Target points to source if source is baseMgNeutronMaterial class
  !!
  pure function baseMgNeutronMaterial_CptrCast(source) result(ptr)
    class(materialHandle), pointer, intent(in) :: source
    class(baseMgNeutronMaterial), pointer          :: ptr

    select type(source)
      class is(baseMgNeutronMaterial)
        ptr => source

      class default
        ptr => null()
    end select

  end function baseMgNeutronMaterial_CptrCast


end module baseMgNeutronMaterial_class
