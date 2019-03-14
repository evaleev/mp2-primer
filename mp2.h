#include <cstddef>
#include <vector>
#include <math.h>
#include "hartree-fock.h"
void fill_ao_ints_vec(std::vector<libint2::Shell>& shells, std::vector<double>& ao_ints_vector);
double int_2e_mo(Matrix& coff_mat, std::vector<double>& ao_ints_vector, size_t p, size_t q, size_t r, size_t s);
