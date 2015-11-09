#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "source/tessellation/geometry.hpp"
#include "source/tessellation/tessellation.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/test_2d/random_pert.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/mpi/MeshPointsMPI.hpp"
#include "source/tessellation/shape_2d.hpp"
#include "source/newtonian/test_2d/piecewise.hpp"
#include "sim_data.hpp"
#include "my_main_loop.hpp"

using namespace std;

int main(void)
{
  SimData sim_data;
  hdsim& sim = sim_data.getSim();

  my_main_loop(sim);

  return 0;
}

