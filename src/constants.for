      MODULE CONSTANTS
      !** speed of light in SI units (m/s)
      real*8,parameter :: CLIGHT_SI = 299792458d0 !(2.99792458d8)
      !** speed of light in CGS units (cm/s)
      real*8,parameter :: CLIGHT_CGS = 299792458d2
      real*8,parameter :: PI =  3.14159265358979323846_8

      real*8, parameter :: kBoltzmann = 0.69503476 ! in cm^(-1)/K

      ! Boltzamnn constant (erg / K)
      REAL*8, PARAMETER :: BoltzmannConstantCGS = 1.3806488D-16

      ! Boltzmann constant (eV / K)
      REAL*8, PARAMETER :: BoltzmannConstantEV = 8.6173324D-5

      ! Planck constant (eV * s)
      REAL*8, PARAMETER :: PlanckConstantEV =    6.58211928D-16

      ! Solar radius (cm)
      REAL*8, PARAMETER :: SolarRadiusKM =       6.955D+5
  
      END MODULE
