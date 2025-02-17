!!
!! This module brakes standard rules
!! It contains a library of XS Packages for Neutron particle type
!!
!!
module neutronXsPackages_class

  use numPrecision
  use endfConstants

  implicit none
  private

  !!
  !! Neutron MACROscopic Reaction XSS
  !!
  !! Public Members:
  !!   total            -> total Cross-Section [1/cm]
  !!   elasticScatter   -> sum of MT = 2 elastic neutron scattering [1/cm]
  !!   inelasticScatter -> sum of all neutron producing reaction that are not elastic scattering
  !!     or fission. [1/cm]
  !!   capture          -> sum of all reactions without secendary neutrons excluding fission [1/cm]
  !!   fission          -> total Fission MT = 18 Cross-section [1/cm]
  !!   nuFission        -> total average neutron production Cross-section [1/cm]
  !!
  !!  Interface:
  !!    clean -> Set all XSs to 0.0
  !!    add   -> Add a nuclide microscopic XSs to macroscopic
  !!    get   -> Return XS by MT number
  !!
  type, public:: neutronMacroXSs
    real(defReal):: total            = ZERO
    real(defReal):: elasticScatter   = ZERO
    real(defReal):: inelasticScatter = ZERO
    real(defReal):: capture          = ZERO
    real(defReal):: fission          = ZERO
    real(defReal):: nuFission        = ZERO
  contains
    procedure:: clean => clean_neutronMacroXSs
    procedure:: add   => add_neutronMacroXSs
    procedure:: get
    procedure:: invert => invert_macroXSs
    procedure:: invert_bl => invert_macroXSs_bl
  end type neutronMacroXSs


  !!
  !! Neutron microscopic Reaction XSS
  !!
  !! Public Members:
  !!   total            -> total Cross-Section [barn]
  !!   elasticScatter   -> MT = 2 elastic neutron scattering [barn]
  !!   inelasticScatter -> all neutron producing reaction that are not elastic scattering
  !!     or fission. [barn]
  !!   capture          -> all reactions without secendary neutrons excluding fission [barn]
  !!   fission          -> total Fission MT = 18 Cross-section [barn]
  !!   nuFission        -> total average neutron production Cross-section [barn]
  !!
  type, public:: neutronMicroXSs
    real(defReal):: total            = ZERO
    real(defReal):: elasticScatter   = ZERO
    real(defReal):: inelasticScatter = ZERO
    real(defReal):: capture          = ZERO
    real(defReal):: fission          = ZERO
    real(defReal):: nuFission        = ZERO
  contains
    procedure:: invert => invert_microXSs
    procedure :: invert_bl => invert_microXSs_bl
  end type neutronMicroXSs

contains

  !!
  !! Clean neutron MacroXSs
  !!
  !! Sets all XSs to 0.0
  !!
  !! Args:
  !!   None
  !!
  !! Errors:
  !!   None
  !!
  elemental subroutine clean_neutronMacroXSs(self)
    class(neutronMacroXSs), intent(inout):: self

    self % total            = ZERO
    self % elasticScatter   = ZERO
    self % inelasticScatter = ZERO
    self % capture          = ZERO
    self % fission          = ZERO
    self % nuFission        = ZERO

  end subroutine clean_neutronMacroXSs

  !!
  !! Add nuclide XSs on Macroscopic XSs
  !!
  !! Takes microscopic XSs*density and adds them to neutronMacroXSs
  !!
  !! Args:
  !!   micro [in] -> microscopic XSs
  !!   dens  [in] -> nuclide density in [1/barn/cm]
  !!
  !! Errors:
  !!   None
  !!
  elemental subroutine add_neutronMacroXSs(self, micro, dens)
    class(neutronMacroXSs), intent(inout):: self
    type(neutronMicroXSs), intent(in)     :: micro
    real(defReal), intent(in)             :: dens

    self % total            = self % total            + dens*micro % total
    self % elasticScatter   = self % elasticScatter   + dens*micro % elasticScatter
    self % inelasticScatter = self % inelasticScatter+dens*micro % inelasticScatter
    self % capture          = self % capture          + dens*micro % capture
    self % fission          = self % fission          + dens*micro % fission
    self % nuFission        = self % nuFission        + dens*micro % nuFission

  end subroutine add_neutronMacroXSs

  !!
  !! Return XSs by MT number
  !!
  !! Args:
  !!   MT [in] -> Requested MT number
  !!
  !! Result:
  !!   Value of the XS
  !!
  !! Errors:
  !!   Returns 0.0 for invalid MT
  !!
  elemental function get(self, MT) result(xs)
    class(neutronMacroXSs), intent(in):: self
    integer(shortInt), intent(in)      :: MT
    real(defReal)                      :: xs

     select case(MT)
      case(macroTotal)
        xs = self % total

      case(macroCapture)
        xs = self % capture

      case(macroEscatter)
        xs = self % elasticScatter

      case(macroFission)
        xs = self % fission

      case(macroNuFission)
        xs = self % nuFission

      case(macroAbsorbtion)
        xs = self % fission+self % capture

      case default
        xs = ZERO

    end select

  end function get

  !! Use a real r in < 0; 1 > to sample reaction for branchless algorithm from Macroscopic XSs
  !!
  !! Args:
  !!    r [in] -> Real number in < 0; 1>
  !!    flag [in] -> flag 'true' for preserving the interface
  !!
  !!  Result:
  !!    One of the macroscopic MT number
  !!      elasticScatter = macroEscatter
  !!      inelasticScatter = macroIEscatter
  !!      fission = macroFission

  elemental function invert_macroXSs_bl(self, r) result(MT)
    class(neutronMacroXSs), intent(in):: self
    real(defReal), intent(in):: r
    integer(shortInt):: MT
    real(defReal):: effectiveXStot
    
    effectiveXStot = self % elasticScatter + self % inelasticScatter &
                    + self % nuFission
    
    if (r*effectiveXStot < self % elasticScatter) then
      MT = N_N_ELASTIC
    elseif (r*effectiveXStot < self % inelasticScatter+self % elasticScatter) then
      MT = N_N_INELASTIC
    elseif (r*effectiveXStot < self % inelasticScatter + self % elasticScatter + self % nuFission) then
      MT = N_FISSION
    else 
      MT = huge(0_shortInt)
    end if
  end function invert_macroXSs_bl
  !!
  !! Use a real r in < 0; 1 > to sample reaction from Macroscopic XSs
  !!
  !! This function might be common thus is type-bound procedure for conveniance
  !!
  !! Args:
  !!   r [in] -> Real number in < 1; 0>
  !!
  !! Result:
  !!   One of the Macroscopic MT numbers
  !!     elasticScatter   = macroEscatter
  !!     inelasticScatter = macroIEscatter
  !!     capture          = macroCapture
  !!     fission          = macroFission
  !!
  !! Errors::
  !!   If r < 0 then returns macroEscatter
  !!   If r > 1 then returns macroFission
  !!
  elemental function invert_macroXSs(self, r) result(MT)
    class(neutronMacroXSs), intent(in):: self
    real(defReal), intent(in)          :: r
    integer(shortInt)                  :: MT
    real(defReal)                      :: xs
    integer(shortInt)                  :: C

    ! Elastic Scattering
    C = 1
    xs = self % total*r - self % elasticScatter
    if (xs > ZERO) C = C+1

    ! Inelastic Scattering
    xs = xs-self % inelasticScatter
    if(xs > ZERO) C = C+1

    ! Capture
    xs = xs-self % capture
    if(xs > ZERO) C = C+1

    ! Choose MT number
    select case(C)
      case(1)
        MT = macroEscatter

      case(2)
        MT = macroIEscatter

      case(3)
        MT = macroCapture

      case(4)
        MT = macroFission

      case default  ! Should never happen -> Avoid compiler error and return nonsense number
        MT = huge(C)

    end select

  end function invert_macroXSs


  !!
  !! Use a real r in < 0; 1 > to sample reaction from Microscopic XSs
  !!
  !! This function involves a bit of code so is written for conviniance
  !!
  !! Args:
  !!   r [in] -> Real number in < 0; 1>
  !!
  !! Result:
  !!   MT number of the reaction:
  !!     elastic scatter   = N_N_elastic
  !!     inelastic scatter = N_N_inelastic
  !!     capture           = N_diasp
  !!     fission           = N_FISSION
  !!
  !! Errors:
  !!   If r < 0 then returns N_N_elastic
  !!   if r > 1 then returns N_FISSION
  !!
  elemental function invert_microXSs(self, r) result(MT)
    class(neutronMicroXSs), intent(in):: self
    real(defReal), intent(in)          :: r
    integer(shortInt)                  :: MT
    real(defReal)                      :: xs
    integer(shortInt)                  :: C

    ! Elastic Scattering
    C = 1
    xs = self % total*r - self % elasticScatter
    if (xs > ZERO) C = C+1

    ! Inelastic Scattering
    xs = xs-self % inelasticScatter
    if(xs > ZERO) C = C+1

    ! Capture
    xs = xs-self % capture
    if(xs > ZERO) C = C+1

    ! Choose MT number
    select case(C)
      case(1)
        MT = N_N_elastic

      case(2)
        MT = N_N_inelastic

      case(3)
        MT = N_disap

      case(4)
        MT = N_fission

      case default  ! Should never happen -> Avoid compiler error and return nonsense number
        MT = huge(C)
    end select

  end function invert_microXSs

elemental function invert_microXSs_bl(self, r) result(MT)
  class(neutronMicroXSs), intent(in) :: self
  real(defReal), intent(in) :: r
  integer(shortInt) :: MT
  real(defReal) :: E_eff


  E_eff = r * (self % elasticScatter + self % inelasticScatter + self % nuFission)

  if (E_eff < self % elasticScatter) then
    MT = N_N_elastic
  elseif (E_eff < self % elasticScatter + self % inelasticScatter) then
    MT = N_N_inelastic
  elseif (E_eff < self % elasticScatter + self % inelasticScatter + self % nuFission) then
    MT = N_fission
  else
    MT = huge(0_shortInt)
  end if
end function invert_microXSs_bl
end module neutronXsPackages_class
