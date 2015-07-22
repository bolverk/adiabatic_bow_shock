#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "source/tessellation/geometry.hpp"
#include "source/tessellation/tessellation.hpp"
#include "source/newtonian/two_dimensional/interpolations/pcm2d.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/newtonian/test_2d/random_pert.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/mpi/mpi_macro.hpp"
#include "source/mpi/MeshPointsMPI.hpp"
#include "source/tessellation/shape_2d.hpp"
#include "source/newtonian/test_2d/piecewise.hpp"
#include "source/newtonian/test_2d/consecutive_snapshots.hpp"
#include "source/newtonian/test_2d/multiple_diagnostics.hpp"
#include "source/misc/vector_initialiser.hpp"
#include "sim_data.hpp"

using namespace std;
using namespace simulation2d;

namespace {

  void my_main_loop(hdsim& sim)
  {
    const double tf = 10;
    SafeTimeTermination term_cond(tf,1e6);
    MultipleDiagnostics diag
      (VectorInitialiser<DiagnosticFunction*>
       (new WriteTime("time.txt"))
       (new ConsecutiveSnapshots
	(new ConstantTimeInterval(tf/1000),
	 new Rubric("snapshot_",".h5")))());
    main_loop(sim,
	      term_cond,
	      &hdsim::TimeAdvance,
	      &diag);
  }
}

int main(void)
{
#ifdef RICH_MPI
  MPI_Init(NULL, NULL);
#endif

  SimData sim_data;
  hdsim& sim = sim_data.getSim();

  my_main_loop(sim);

#ifdef RICH_MPI
  write_snapshot_to_hdf5(sim, "process_"+int2str(get_mpi_rank())+"_final.h5");
#else
  write_snapshot_to_hdf5(sim, "final.h5");
#endif


#ifdef RICH_MPI
  MPI_Finalize();
#endif

  return 0;
}

