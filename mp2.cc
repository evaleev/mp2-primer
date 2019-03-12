/******************************
 * MP2 Primer                 *
 * An extension to HF Primer  *
 *                            *
 ******************************/
#include <cstddef>
#include <iostream>
#include <vector>
#include "hartree-fock.h"

void fill_ao_ints_vec(std::vector<libint2::Shell>& shells, std::vector<double>& ao_ints_vector) { 
     
    /* double mu, double neu, double lambda, double sigma) { */
    // returns 2e integral on ao basis

    using std::cout;
    using std::endl;

    using libint2::Shell;
    using libint2::Engine;
    using libint2::Operator;

    ao_ints_vector.clear();
     
    libint2::initialize();

    // construct the electron repulsion integrals engine
    Engine engine(Operator::coulomb, max_nprim(shells), max_l(shells), 0);

    auto shell2bf = map_shell_to_basis_function(shells);

    // buf[0] points to the target shell set after every call  to engine.compute()
    const auto& buf = engine.results();

    // loop over shell pairs of the Fock matrix, {s1,s2}
    // Fock matrix is symmetric, but skipping it here for simplicity (see compute_2body_fock)
    for(auto s1=0; s1!=shells.size(); ++s1) {

        auto bf1_first = shell2bf[s1]; // first basis function in this shell
        auto n1 = shells[s1].size();

        for(auto s2=0; s2!=shells.size(); ++s2) {

            auto bf2_first = shell2bf[s2];
            auto n2 = shells[s2].size();

            // loop over shell pairs of the density matrix, {s3,s4}
            // again symmetry is not used for simplicity
            for(auto s3=0; s3!=shells.size(); ++s3) {

                auto bf3_first = shell2bf[s3];
                auto n3 = shells[s3].size();

                for(auto s4=0; s4!=shells.size(); ++s4) {

                    auto bf4_first = shell2bf[s4];
                    auto n4 = shells[s4].size();

                    // Coulomb contribution to the Fock matrix is from {s1,s2,s3,s4} integrals
                    engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
                    const auto* buf_1234 = buf[0];
                    if (buf_1234 == nullptr)
                        continue; // if all integrals screened out, skip to next quartet

                    // we don't have an analog of Eigen for tensors (yet ... see github.com/BTAS/BTAS, under development)
                    // hence some manual labor here:
                    // 1) loop over every integral in the shell set (= nested loops over basis functions in each shell)
                    // and 2) add contribution from each integral
                    for(auto f1=0, f1234=0; f1!=n1; ++f1) {
                        const auto bf1 = f1 + bf1_first;
                        for(auto f2=0; f2!=n2; ++f2) {
                            const auto bf2 = f2 + bf2_first;
                            for(auto f3=0; f3!=n3; ++f3) {
                                const auto bf3 = f3 + bf3_first;
                                for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                                    const auto bf4 = f4 + bf4_first;
                                    ao_ints_vector.push_back( buf_1234[f1234] );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

double int_2e_mo(Eigen::MatrixXd& coff_mat, std::vector<double>& ao_ints_vector,
                    size_t p, size_t q, size_t r, size_t s) {
    // specified p, q, r and s: indices of the two pairs of molecular orbitals
    // returns the two-electron repulsion integral
    return 1.0;
}

int main() {
    using std::cout;
    using std::endl;
    const auto filename = "h2o.xyz";
    std::vector<Atom> atoms = read_geometry(filename);
    auto shells = make_sto3g_basis(atoms);
    std::vector<double> ao_ints_vector;
    fill_ao_ints_vec(shells, ao_ints_vector);
    cout << ao_ints_vector.size() << endl;
    return 0;
}
