/******************************************************************************
 *                                                                            *
 * CONSTANTS.H                                                                *
 *                                                                            *
 * NUMERICAL CONSTANTS IN CGS UNITS                                           *
 *                                                                            *
 ******************************************************************************/

// Fundamental constants
#define EE (4.80320680e-10)       // Electron charge
#define CL (2.99792458e10)        // Speed of light
#define ME (9.1093826e-28)        // Electron mass
#define MP (1.67262171e-24)       // Proton mass
#define MN (1.67492728e-24)       // Neutron mass
#define HPL (6.6260693e-27)       // Planck constant
#define HBAR (HPL / (2. * M_PI))  // Reduced Planck constant
#define KBOL (1.3806505e-16)      // Boltzmann constant
#define GNEWT (6.6742e-8)         // Gravitational constant
#define SIG (5.670400e-5)         // Stefan-Boltzmann constant
#define AR (4 * SIG / CL)         // Radiation constant
#define THOMSON (0.665245873e-24) // Thomson cross section
#define COULOMB_LOG (20.)         // Coulomb logarithm
#define ALPHAFS (0.007299270073)  // Fine structure constant ~ 1./137.
#define GFERM (1.435850814e-49)   // Fermi constant
#define GA (-1.272323)            // Axial-vector coupling
#define GA2 (GA * GA)
#define S2THW (0.222321) // sin^2(Theta_W), Theta_W = Weinberg angle
#define S4THW (S2THW * S2THW)
#define NUSIGMA0 (1.7611737037e-44) // Fundamental neutrino cross section

// Frequency scale of neutrino oscillations
#define ROOT2 (1.4142135623730951)
#define NUFERM ((ROOT2*GFERM) / HBAR)

// Unit conversions
#define EV (1.60217653e-12)   // Electron-volt
#define MEV (1.0e6 * EV)      // Mega-Electron-Volt
#define GEV (1.0e9 * EV)      // Giga-Electron-Volt
#define JY (1.e-23)           // Jansky
#define PC (3.085678e18)      // Parsec
#define AU (1.49597870691e13) // Astronomical unit
#define YEAR (31536000.)
#define DAY (86400.)
#define HOUR (3600.)
#define MSUN (1.989e33) // Solar mass
#define RSUN (6.96e10)  // Solar radius
#define LSUN (3.827e33) // Solar luminosity
