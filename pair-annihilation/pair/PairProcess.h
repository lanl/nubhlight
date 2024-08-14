/* ========================================
 * Neutrino opacities: pair process
 * ========================================
 * Copyright (C) Maitraya Bhattacharyya
 * Implements pair reaction for neutrinos
 * See: Pons et. al. (1998), Bruenn (1985)
 *
 * History:
 * 01/2023: maitraya
 */

#include <tuple>

#ifndef SRC_PAIRPROCESS_H
#define SRC_PAIRPROCESS_H

namespace pairprocess {

    /* ---------------------------------------------------------
     * PairProcess class
     * Contains functions to compute kernels and opacities
     * ---------------------------------------------------------
     */

    class PairProcess {

    public:
        double eta_pr;
        double filt = 0.;
        const double G_sqr = 1.55 * 1E-33;                      // [cm^3 MeV^-2 s^-1], from Bruenn Eqn. (C51)
        //const double G_sqr = 1.0;
        double planck = 1;
        double c_speed = 1;
        const double sin_sqr_thetaw = 0.23120;                  // square of sine of Weinberg angle, approximate value (https://en.wikipedia.org/wiki/Weinberg_angle)
        double alpha_1;
        double alpha_2;
        const double alpha_1_e = 1. + 2. * sin_sqr_thetaw;      // Table 1, Pons et. al. (1998)
        const double alpha_2_e = 2. * sin_sqr_thetaw;           // Table 1, Pons et. al. (1998)
        const double alpha_1_mutau = -1. + 2. * sin_sqr_thetaw; // Table 1, Pons et. al. (1998)
        const double alpha_2_mutau = 2. * sin_sqr_thetaw;       // Table 1, Pons et. al. (1998)

        double a[4][12];                                        // a_ln from Eqn. (11) Pons et. al.
        double c[4][3];                                         // c_ln from Eqn. (11) Pons et. al.
        double d[4][3];                                         // d_ln from Eqn. (11) Pons et. al.

        PairProcess(char nutype, double eta_prime, double filtpar = 0.);

        double T(int l, double alpha, double tolerance = 1e-6); // T_l(alpha) from Appendix B Pons et. al.

        double F(int k, double eta, double x1); // F_k(eta, x1) from Appendix B Pons et. al.

        double G(int n, double a, double b, double eta, double y, double z); // G_n(a,b) from Eqn. (12) Pons et. al.

        double Psi(int l, double y, double z, double eta);  // Psi_l(y,z) from Eqn. (11) Pons et. al.

        double Phi(int l, double omega, double omega_prime, double eta,
                   double T);  // Phi_l(omega,omega') from Eqn. (10) Pons et. al.

        std::tuple<double, double> R_TP(double omega, double omega_prime, double cos_theta, double eta, double T,
                                        int lmax = 3); // (R^p_TP,R^a_TP) from Eqn. (3) & (2) Pons et. al.

        std::tuple<double, double> EtaLambdaInv(double omega, double mu, double eta, double T, int lmax);
        double FermiDirac(double omega_prime, double eta_prime);
    };
}


#endif //SRC_PAIRPROCESS_H
