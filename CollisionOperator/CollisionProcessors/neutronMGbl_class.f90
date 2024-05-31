module neutronMGbls_class

  use numPrecision
  use endfConstants
  use genericProcedures,             only : fatalError, rotateVector, numToChar
  use dictionary_class,              only : dictionary
  use RNG_class,                     only : RNG

  ! Particle types
  use particle_class,                only : particle, particleState, printType, P_NEUTRON
  use particleDungeon_class,         only : particleDungeon

  ! Abstract interface
  use collisionProcessor_inter,      only : collisionProcessor, collisionData, init_super => init

  ! Nuclear Data Interface
  use nuclearDataReg_mod,            only : ndReg_getNeutronMG => getNeutronMG
  use nuclearDatabase_inter,         only : nuclearDatabase
  use mgNeutronDatabase_inter,       only : mgNeutronDatabase
  use mgNeutronMaterial_inter,       only : mgNeutronMaterial, mgNeutronMaterial_CptrCast
  use reactionHandle_inter,          only : reactionHandle
  use multiScatterMG_class,          only : multiScatterMG, multiScatterMG_CptrCast
  use fissionMG_class,               only : fissionMG, fissionMG_TptrCast

  ! Cross section packages
  use neutronXsPackages_class,       only : neutronMacroXSs


  ! Nuclear Data
  !use nuclearData_inter,              only : nuclearData
  !use perMaterialNuclearDataMG_inter, only : perMaterialNuclearDataMG

  ! Cross-section packages to interface with nuclear data
  !use xsMacroSet_class,               only : xsMacroSet, xsMacroSet_ptr

  implicit none
  private

  !!
  !! Standard (default) scalar collision processor for MG neutrons
  !!   -> Preforms implicit fission site generation
  !!   -> Preforms analog capture
  !!   -> Treats fission as capture (only implicit generation of 2nd-ary neutrons)
  !!   -> Does not create secondary non-neutron projectiles
  !!
  !! Settings:
  !!  NONE
  !!
  !! Sample dictionary input:
  !!   collProcName {
  !!   type            neutronMGbls; 
  !!   }
  !!
  type, public, extends(collisionProcessor):: neutronMGbls
    private
    class(mgNeutronDatabase), pointer, public:: xsData => null()
    class(mgNeutronMaterial), pointer, public:: mat    => null()

    !! Russian roulette and splitting settings
    real(defReal) :: minWgt
    real(defReal) :: maxWgt
    real(defReal) :: avWgt
  contains
    ! Initialisation procedure
    procedure:: init

    ! Implementation of customisable procedures
    procedure:: sampleCollision
    procedure:: implicit
    procedure:: elastic
    procedure:: inelastic
    procedure:: capture
    procedure:: fission
    procedure:: cutoffs

    ! Local procedures
    procedure,private :: russianRoulette
    procedure,private :: split
  end type neutronMGbls

contains

  !!
  !! Initialise from dictionary
  !!
  subroutine init(self, dict)
    class(neutronMGbls), intent(inout):: self
    class(dictionary), intent(in)      :: dict
    character(100), parameter:: Here = 'init (neutronMGbls_class.f90)'

    ! Call superclass
    call init_super(self, dict)

    ! Roulette and splitting settings

    call dict % getOrDefault(self % minWgt, 'minWgt', 0.4_defReal)
    call dict % getOrDefault(self % maxWgt, 'maxWgt', 2.0_defReal)
    call dict % getOrDefault(self % avWgt, 'avWgt', 1.0_defReal)

  end subroutine init

  !!
  !! Samples collision without any implicit treatment
  !!
  subroutine sampleCollision(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGbls), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon), intent(inout):: thisCycle
    class(particleDungeon), intent(inout):: nextCycle
    type(neutronMacroXSs)                :: macroXSs
    real(defReal)                        :: r
    character(100), parameter:: Here =' sampleCollision (neutronMGbls_class.f90)'

    ! Verify that particle is MG neutron
    if( .not. p % isMG .or. p % type /= P_NEUTRON) then
      call fatalError(Here, 'Supports only MG Neutron. Was given CE '//printType(p % type))
    end if

    ! Verify and load nuclear data pointer
    self % xsData => ndReg_getNeutronMG()
    if(.not.associated(self % xsData)) call fatalError(Here, "Failed to get active database for MG Neutron")

    ! Get and verify material pointer
    self % mat => mgNeutronMaterial_CptrCast( self % xsData % getMaterial( p % matIdx()))
    if(.not.associated(self % mat)) call fatalError(Here, "Failed to get MG Neutron Material")

    ! Select Main reaction channel
    call self % mat % getMacroXSs(macroXSs, p % G, p % pRNG)
    r = p % pRNG % get()
    
    ! Call branchless collide procedure instead of regular
    collDat % MT = macroXSs % invert_bl(r)

  end subroutine sampleCollision

  !!
  !! Modify statistical weight of particle for branchless col
  !!
  subroutine implicit(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGbls), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon), intent(inout):: thisCycle
    class(particleDungeon), intent(inout):: nextCycle
    type(neutronMacroXSs)                :: macroXSs
    real(defReal)                        :: effectiveXStot, weight_modifier
    character(100), parameter:: Here = 'implicit (neutronMGbls_class.f90)'

    call self % mat % getMacroXSs(macroXSs, p % G, p % pRNG)
    
    ! Compute effective total cross section for branchless collisions
    effectiveXStot = macroXSs % elasticScatter+macroXSs % inelasticScatter &
      + macroXSs % nuFission
    
    ! Modify the neutron statistical weight
    weight_modifier = effectiveXStot/macroXSs % total
    p % w = p % w * weight_modifier

  end subroutine implicit

  !!
  !! Elastic Scattering
  !!
  subroutine elastic(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGbls), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon), intent(inout):: thisCycle
    class(particleDungeon), intent(inout):: nextCycle

    ! Do nothing. Should not be called

  end subroutine elastic

  !!
  !! Preform scattering
  !!
  subroutine inelastic(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGbls), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon), intent(inout):: thisCycle
    class(particleDungeon), intent(inout):: nextCycle
    class(multiScatterMG), pointer        :: scatter
    integer(shortInt)                    :: G_out   ! Post-collision energy group
    real(defReal)                        :: phi     ! Azimuthal scatter angle
    real(defReal)                        :: w_mul   ! Weight multiplier
    character(100), parameter:: Here = "inelastic (neutronMGbls_class.f90)"

    ! Assign MT number
    collDat % MT = macroIEscatter

    ! Get Scatter object
    scatter => multiScatterMG_CptrCast( self % xsData % getReaction(macroIEscatter, collDat % matIdx))
    if(.not.associated(scatter)) call fatalError(Here, "Failed to get scattering reaction object for MG neutron")

    ! Sample Mu and G_out
    call scatter % sampleOut(collDat % muL, phi, G_out, p % G, p % pRNG)

    ! Read scattering multiplicity
    w_mul = scatter % production(p % G, G_out)

    ! Update neutron state
    p % G = G_out
    p % w = p % w*w_mul
    call p % rotate(collDat % muL, phi)

  end subroutine inelastic

  !!
  !! Preform capture
  !!
  subroutine capture(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGbls), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon), intent(inout):: thisCycle
    class(particleDungeon), intent(inout):: nextCycle

    !! Should never happen with branchless collisions
    
  end subroutine capture

  !!
  !! Preform fission
  !!
  subroutine fission(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGbls), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon), intent(inout):: thisCycle
    class(particleDungeon), intent(inout):: nextCycle
    type(fissionMG), pointer:: fission_reaction
    type(particleState):: pTemp
    real(defReal), dimension(3):: r, dir
    real(defReal):: mu, phi, w0
    integer(shortInt):: fission_yield
    integer(shortInt):: g_out, i
    character(100), parameter:: Here = 'fission (neutronMGbls_class.f90)'

    ! Get Fission reaction object
    fission_reaction => fissionMG_TptrCast( self % xsData % getReaction(macroFission, collDat % matIdx))
    if (.not.associated(fission_reaction)) call fatalError(Here, 'Failed to getrive fissionMG reaction object')

    ! Store new sites in the next cycle dungeon

    w0 = p % preHistory % wgt
    fission_yield = int(abs(p % pRNG % get() + p % w / (p % k_eff * w0)),kind=shortInt)
    do i = 1, fission_yield
      r   = p % rGlobal()
      call fission_reaction % sampleOut(mu, phi, G_out, p % G, p % pRNG)
      dir = rotateVector(p % dirGlobal(), mu, phi)

      ! Copy extra detail from parent particle (i.e. time, flags ect.)
      pTemp       = p

      ! Overwrite position, direction, energy group
      pTemp % r   = r
      pTemp % dir = dir
      pTemp % G   = G_out
      pTemp % wgt = 1.0

      call nextCycle % detain(pTemp)
    end do

    p % isDead = .true.

  end subroutine fission

  !!
  !! Applay cutoffs or post-collision implicit treatment
  !!
  subroutine cutoffs(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGbls), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle


    ! Splitting with fixed threshold
    if ((p % w > self % maxWgt) .and. (p % isDead .eqv. .false.)) then
      call self % split(p, thisCycle, self % maxWgt)
    ! Roulette with fixed threshold and survival weight
    elseif ((p % w < self % minWgt) .and. (p % isDead .eqv. .false.)) then
      call self % russianRoulette(p, self % avWgt)
    end if

  end subroutine cutoffs

    !!
  !! Perform Russian roulette on a particle
  !!
  subroutine russianRoulette(self, p, avWgt)
    class(neutronMGbls), intent(inout):: self
    class(particle), intent(inout)     :: p
    real(defReal), intent(in)          :: avWgt

    if (p % pRNG % get() < (ONE-p % w/avWgt)) then
      p % isDead = .true.
    else
      p % w = avWgt
    end if

  end subroutine russianRoulette

  !!
  !! Split particle which has too large a weight
  !!
  subroutine split(self, p, thisCycle, maxWgt)
    class(neutronMGbls), intent(inout)    :: self
    class(particle), intent(inout)        :: p
    class(particleDungeon), intent(inout):: thisCycle
    real(defReal), intent(in)             :: maxWgt
    integer(shortInt)                     :: mult, i

    ! This value must be at least 2
    mult = ceiling(p % w/maxWgt)
    ! Decrease weight
    p % w = p % w/mult

    ! Add split particle's to the dungeon
    do i = 1, mult-1
      call thisCycle % detain(p)
    end do

  end subroutine split

end module neutronMGbls_class
