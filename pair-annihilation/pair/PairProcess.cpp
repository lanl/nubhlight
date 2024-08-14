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

#include <cmath>
#include <iostream>
#include <tuple>
#include <cassert>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include "PairProcess.h"
#include "../integration/QuadratureIntegration.h"

namespace pairprocess {

    /**
     * @brief Constructor for PairProcess class
     *
     * The constructor only allows a character array. Choose 'x' for mu/tau neutrinos and 'e' of electron
     * neutrinos. Sets constants accordingly.
     *
     * @param nutype 'x' or 'e' for the type of neutrino species.
     */
    PairProcess::PairProcess(char nutype, double eta_prime, double filtpar) {

        eta_pr = eta_prime;
        filt = filtpar;

        if (nutype == 'x') {
            alpha_1 = alpha_1_mutau;
            alpha_2 = alpha_2_mutau;
        } else if (nutype == 'e') {
            alpha_1 = alpha_1_e;
            alpha_2 = alpha_2_e;
        } else {
            throw std::invalid_argument("Expected character x or e as input.");
        }
    }

    /**
     * @brief Calculates T_l(alpha) from appendix B, Pons et. al.
     *
     * Calculate:
     *          T_l(alpha) = sum_{n=1}^{inf} (-1)^{n+1} e^{-n alpha} / n^l
     *
     * The computation is truncated when a certain accuracy is reached.
     *
     * @param l integer >= 1
     * @param alpha double >= 0
     * @param tolerance accuracy of computation, set to 1e-6 by default.
     * @return T_l(alpha), double
     */
    double PairProcess::T(int l, double alpha, double tolerance) { //@TODO: check this function after lunch

        assert(alpha >= 0 && l >= 1);

        if (alpha == 0 and l != 1) {

            return pow(2., -l) * (pow(2., l) - 2.) * boost::math::zeta<double>(l);  // Computed in Mathematica

        } else if (alpha == 0 and l == 1) {

            return log(2.0);

        } else {

            double val = abs(42. * tolerance);
            double result{0.};
            int n{1};

            while (tolerance < abs(val)) {
                val = pow(-1., n + 1) * exp(-n * alpha) / pow(n, l);
                result += val;
                n++;
            }

            return result;

        }
    }

    double PairProcess::FermiDirac(double omega_prime, double eta_prime) {
        return 1. / (exp(omega_prime - eta_prime) + 1.);
    }

    /**
     * @brief Calculates F_k(eta,x1) from appendix B, Pons et. al.
     *
     * @param k an integer
     * @param eta a double
     * @param x1 a double
     * @return F_k(eta,x1), a double
     */
    double PairProcess::F(int k, double eta, double x1) {

        double result{0.};

        if (eta < 0.) {

            double sum{0.};

            for (size_t l = 0; l <= k; l++) {
                sum += T(k + 1 - l, x1 - eta) * pow(x1, l) / boost::math::factorial<double>(l);
            }

            result = boost::math::factorial<double>(k) * (T(k + 1., -eta) - sum);

        } else if (0. <= eta && eta <= x1) {

            double sum{0.};
            double sum_2{0.};

            if (k > 0) { // When k = 0, this sum should be equal to zero
                for (size_t l = 0; l <= int((k - 1.) / 2.); l++) {
                    sum += T(2 * l + 2, 0.) * pow(eta, k - 1. - 2. * l) / boost::math::factorial<double>(k - 1 - 2 * l);
                }
            }

            for (size_t l = 0; l <= k; l++) {
                sum_2 += T(k + 1 - l, x1 - eta) * pow(x1, l) / boost::math::factorial<double>(l);
            }

            result = pow(eta, k + 1.) / (k + 1.) +
                     boost::math::factorial<double>(k) * (2.0 * sum + pow(-1., k) * T(k + 1, eta) - sum_2);

        } else {

            double sum{0.};

            for (size_t l = 0; l <= k; l++) {
                sum += pow(-1., k - l) * T(k + 1 - l, eta - x1) * pow(x1, l) / boost::math::factorial<double>(l);
            }

            result = pow(x1, k + 1.) / (k + 1.) +
                     boost::math::factorial<double>(k) * (pow(-1., k) * T(k + 1, eta) - sum);

        }

        return result;
    }

    /**
     * Calculates G_n(a,b) from Eqn. (12) Pons et. al. Special form required for n = 0 case.
     *
     * @param n an integer
     * @param a dimensionless quantity, double
     * @param b dimensionless quantity, double
     * @param eta the electron degeneracy parameter, a double
     * @param y the neutrino energy divided by temperature y = omega / T to make it dimensionless, a double
     * @param z the anti-neutrino enegy, also made dimensionless by dividing with temperature, a double
     * @return a double
     */
    double PairProcess::G(int n, double a, double b, double eta, double y, double z) {

        double result{0.};

        result = F(n, eta, b) - F(n, eta, a) - F(n, eta + y + z, b) + F(n, eta + y + z, a);

        return result;
    }

    /**
     * @brief Calculates Psi_l(y,z) from Eqn. (11) Pons et. al.
     *
     * @param l degree of Legendre polynomial, integer
     * @param y dimensionless quantity, double
     * @param z dimensionless quantity, double
     * @param eta electron degeneracy parameter, double
     * @return Psi_l(y,z), a double
     */
    double PairProcess::Psi(int l, double y, double z, double eta) {

        assert(0 <= l);
        assert(l <= 3);

        // Coefficients for l = 0
        a[0][3] = 8. / (3. * y * y);
        a[0][4] = -4. / (3. * y * y * z);
        a[0][5] = 4. / (15. * y * y * z * z);
        c[0][0] = (4. * y / (z * z)) * (2. * z * z / 3. + y * z + 2. * y * y / 5.);
        c[0][1] = -(4. * y / (3. * z * z)) * (3. * y + 4. * z);
        c[0][2] = 8. * y / (3. * z * z);
        d[0][0] = 4. * z * z * z / (15. * y * y);
        d[0][1] = -4. * z * z / (3. * y * y);
        d[0][2] = 8. * z / (3. * y * y);

        // Coefficients for l = 1
        a[1][3] = 8. / (3. * y * y);
        a[1][4] = -4. * (4. * y + 3. * z) / (3. * y * y * y * z);
        a[1][5] = 4. * (13. * y + 18. * z) / (15. * y * y * y * z * z);
        a[1][6] = -4. * (y + 3. * z) / (5. * y * y * y * z * z * z);
        a[1][7] = 16. / (35. * y * y * y * z * z * z);
        c[1][0] = -(4. * y / (z * z * z)) *
                  (2. * y * y * y / 7. + 4. * y * y * z / 5. + 4. * y * z * z / 5. + z * z * z / 3.);
        c[1][1] = (4. * y / (z * z * z)) * (4. * y * y / 5. + 7. * y * z / 5. + 2. * z * z / 3.);
        c[1][2] = -(4. * y / (z * z * z)) * (3. * y / 5. + z / 3.);
        d[1][0] = -(4. * z * z * z / (105. * y * y * y)) * (14. * y + 9. * z);
        d[1][1] = (4. * z * z / (5. * y * y * y)) * (7. * y / 3. + 2. * z);
        d[1][2] = -(4. * z / (y * y * y)) * (y / 3. + 3. * z / 5.);

        // Coefficients for l = 2
        a[2][3] = 8. / (3. * y * y);
        a[2][4] = -4. * (10. * y + 9. * z) / (3. * y * y * y * z);
        a[2][5] = 4. * (73. * y * y + 126. * y * z + 36. * z * z) / (15. * y * y * y * y * z * z);
        a[2][6] = -12. * (y * y + 3. * y * z + 8. * z * z / 5.) / (y * y * y * y * z * z * z);
        a[2][7] = 48. * (2. * y * y + 13. * y * z + 12. * z * z) / (35. * y * y * y * y * z * z * z * z);
        a[2][8] = -24. * (y + 2. * z) / (7. * y * y * y * y * z * z * z * z);
        a[2][9] = 8. / (7. * y * y * y * y * z * z * z * z);
        c[2][0] = 4. * y * (2. * y * y * y * y / 7. + 6. * y * y * y * z / 7 + 32. * y * y * z * z / 35. +
                            2. * y * z * z * z / 5. + z * z * z * z / 15.) / (z * z * z * z);
        c[2][1] =
                -4. * y * (6. * y * y * y / 7. + 12. * y * y * z / 7. + y * z * z + 2. * z * z * z / 15.) / (z * z * z *
                                                                                                             z);
        c[2][2] = 8. * y * (z * z / 6. + 3 * y * z / 2. + 12 * y * y / 7.) / (5. * z * z * z * z);
        d[2][0] = 4. * z * z * z * (16. * y * y + 27. * y * z + 12. * z * z) / (105. * y * y * y * y);
        d[2][1] = -4. * z * z * (18. * z * z / 35. + 6. * y * z / 7. + y * y / 3.) / (y * y * y * y);
        d[2][2] = 8. * z * (y * y / 6. + 3. * y * z / 2. + 12. * z * z / 7.) / (5. * y * y * y * y);

        // Coefficients for l = 3
        a[3][3] = 8. / (3. * y * y);
        a[3][4] = -4. * (19. * y + 18. * z) / (3. * y * y * y * z);
        a[3][5] = 4. * (253. * y * y + 468. * y * z + 180. * z * z) / (15. * y * y * y * y * z * z);
        a[3][6] = -8. * (149. * y * y * y + 447. * y * y * z + 330. * y * z * z + 50. * z * z * z) /
                  (15. * y * y * y * y * y * z * z * z);
        a[3][7] = 8. * (116. * y * y * y + 2916. * y * y * z / 5. + 696. * y * z * z + 200. * z * z * z) /
                  (21. * y * y * y * y * y * z * z * z * z);
        a[3][8] = -40. * (5. * y * y * y + 54. * y * y * z + 108. * y * z * z + 50. * z * z * z) /
                  (21. * y * y * y * y * y * z * z * z * z * z);
        a[3][9] = 40. * (10. * y * y + 43. * y * z + 100. * z * z / 3.) / (21. * y * y * y * y * y * z * z * z * z * z);
        a[3][10] = -40. * (3. * y + 5. * z) / (9. * y * y * y * y * y * z * z * z * z * z);
        a[3][11] = 320. / (99. * y * y * y * y * y * z * z * z * z * z);
        c[3][0] = -4. * y * y * (10. * y * y * y * y / 33. + 20. * y * y * y * z / 21. + 68. * y * y * z * z / 63. +
                                 18. * y * z * z * z / 35. + 3. * z * z * z * z / 35.) / (z * z * z * z * z);
        c[3][1] = 4. * y * y *
                  (20. * y * y * y / 21. + 130. * y * y * z / 63. + 48. * y * z * z / 35. + 9. * z * z * z / 35.) / (z *
                                                                                                                     z *
                                                                                                                     z *
                                                                                                                     z *
                                                                                                                     z);
        c[3][2] = -4. * y * y * (50. * y * y / 63. + 6. * y * z / 7. + 6. * z * z / 35.) / (z * z * z * z * z);
        d[3][0] = -4. * z * z * z * (9. * y * y * y + 34. * y * y * z + 40. * y * z * z + 500. * z * z * z / 33.) /
                  (105. * y * y * y * y * y);
        d[3][1] = 4. * z * z * (3. * y * y * y + 24. * y * y * z + 130. * y * z * z / 3. + 200. * z * z * z / 9.) /
                  (35. * y * y * y * y * y);
        d[3][2] = -4. * z * z * (50. * z * z / 63. + 6. * z * y / 7. + 6. * y * y / 35.) / (y * y * y * y * y);

        double result{0.};

        for (size_t n = 0; n <= 2; n++) {
            result += c[l][n] * G(n, y, y + z, eta, y, z) + d[l][n] * G(n, z, y + z, eta, y, z);
        }

        for (size_t n = 3; n <= 2 * l + 5; n++) {
            result += a[l][n] * (G(n, 0, std::min(y, z), eta, y, z) - G(n, std::max(y, z), y + z, eta, y, z));
        }

        return result;
    }

    /**
     * @brief Calculates Phi_l(y,z) from Eqn. (10) of Pons et. al.
     *
     * @param l degree of Legendre polynomial, integer
     * @param omega energy of neutrino, double
     * @param omega_prime energy of anti-neutrino, double
     * @param eta electron degeneracy parameter, double
     * @param T temperature in units of energy, double
     * @return Phi_l(y,z), double
     */
    double PairProcess::Phi(int l, double omega, double omega_prime, double eta, double T) {

        double y = omega / T;
        double z = omega_prime / T;

        // double result = G_sqr * T * T * (alpha_1 * Psi(l, y, z, eta) + alpha_2 * Psi(l, z, y, eta)) /
        //(M_PI * (1. - exp(y + z)));
        double result = G_sqr * T * T * (alpha_1 * alpha_1 * Psi(l, y, z, eta) + alpha_2 * alpha_2 * Psi(l, z, y, eta)) /
                        (M_PI * (1. - exp(y + z))); // @TODO: This expression has been modified to add squares to the constants

        return result;
    }

    /**
     * @brief Calculates the production and absorption kernels
     *
     * @param omega energy of neutrino, double
     * @param omega_prime energy of anti-neutrino, double
     * @param cos_theta cosine of the angle between neutrino and anti-neutrino directions, double
     * @param eta electron degeneracy parameter, double
     * @param T the temperature in units of energy, double
     * @param lmax upper limit of Legendre sum
     * @return production and absorption kernels
     */
    std::tuple<double, double>
    PairProcess::R_TP(double omega, double omega_prime, double cos_theta, double eta, double T, int lmax) {

        double R_TP_p{0.};
        double R_TP_a{0.};

        lmax = 3 * (lmax > 3) + lmax * (lmax <= 3); // lmax cannot be greater than 4

        for (size_t l = 0; l <= lmax; l++) {
            double Pl = boost::math::legendre_p<double>(l, cos_theta);
            double Phival = Phi(l, omega, omega_prime, eta, T);

            R_TP_p += (1. / (1. + filt * l * l * (l + 1) * (l + 1))) * (2. * l + 1.) * Phival * Pl / 2.;
        }

        R_TP_a = exp((omega + omega_prime) / T) * R_TP_p;

        return std::make_tuple(R_TP_p, R_TP_a);
    }

    std::tuple<double, double> PairProcess::EtaLambdaInv(double omega, double mu, double eta, double T, int lmax) {
        /* auto IntegrationObject = new integration::QuadratureIntegration();

        lmax = 3 * (lmax > 3) + lmax * (lmax <= 3);

        double emissivity{0.};
        double absorpsivity{0.};

        std::function<double(double, double)> om2Phi0 = [this, &omega, &eta, &T](double omega_prime, double mu) {
            return omega_prime * omega_prime * Phi(0, omega, omega_prime, eta, T);
        };

        emissivity = 0.5 * (2. * M_PI) / (c_speed * pow(planck * c_speed, 3)) * IntegrationObject->QuadIntegrate(om2Phi0);
        absorpsivity = emissivity;

        for (size_t i = 0; i <= lmax; i++) {

            std::function<double(double, double)> om2PhilPlFbar = [this, &i, &omega, &eta, &T](double omega_prime, double mu) {
                return omega_prime * omega_prime * Phi(i, omega, omega_prime, eta, T) * boost::math::legendre_p<double>(i, mu) *
                       FermiDirac(omega_prime, eta_pr);
            };
            std::function<double(double, double)> om2PhilPlFbarFactor = [this, &i, &omega, &eta, &T](double omega_prime, double mu) {
                return omega_prime * omega_prime * Phi(i, omega, omega_prime, eta, T) * boost::math::legendre_p<double>(i, mu) * (1. - exp((omega + omega_prime) / T)) *
                       FermiDirac(omega_prime, eta_pr);
            };

            emissivity -= (2. * M_PI) / (c_speed * pow(planck * c_speed, 3)) * ((2. * +1.) / 2.) * boost::math::legendre_p<double>(i, mu) *
                          (1. / (1. + filt * i * i * (i + 1) * (i + 1))) * IntegrationObject->QuadIntegrate(om2PhilPlFbar);
            absorpsivity -= (2. * M_PI) / (c_speed * pow(planck * c_speed, 3)) * ((2. * i + 1.) / 2.) * boost::math::legendre_p<double>(i, mu) *
                            (1. / (1. + filt * i * i * (i + 1) * (i + 1))) * IntegrationObject->QuadIntegrate(om2PhilPlFbarFactor);

        }

        return std::make_tuple(emissivity, absorpsivity); */

    }
}
