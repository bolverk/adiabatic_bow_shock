#ifndef CALC_INIT_COND_HPP
#define CALC_INIT_COND_HPP 1

#include <vector>
#include "source/newtonian/two_dimensional/computational_cell_2d.hpp"
#include "source/tessellation/tessellation.hpp"

using std::vector;

vector<ComputationalCell> calc_init_cond(const Tessellation& tess);

#endif // CALC_INIT_COND_HPP
