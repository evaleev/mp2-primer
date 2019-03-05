/******************************
 * MP2 Primer                 *
 * An extension to HF Primer  *
 *                            *
 ******************************/
#include <cstddef>
#include "hartree-fock.h"
double compute_2body_ints(double mu, double neu, double lambda, double sigma)
{
    return 1.0;
}

double transform_ints(const Matrix& ao_mat,
                      size_t mu_cf, size_t neu_cf,
                      size_t lambda_cf, size_t sigma_cf)
{
    const auto n = ao_mat.rows();
    size_t mu, neu, lambda, sigma;
    double integral = 1.0;
    for ( mu = 0; mu < n; ++mu ) {
        double int_1 = 1.0;
        for ( neu = 0; neu < n; ++neu ) {
            double int_2 = 1.0;
            for ( lambda = 0; lambda < n; ++lambda ) {
                double int_3 = 1.0;
                for ( sigma = 0; sigma < n; ++sigma ) {
                    int_3 += ao_mat(sigma_cf, sigma)*compute_2body_ints(1,1,1,1);
                }
            }
        }
    }
     
    return 1.0;
}
