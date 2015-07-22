#ifdef RICH_MPI
#include <mpi.h>
#endif
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "source/tessellation/geometry.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/tessellation/tessellation.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/interpolations/pcm2d.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/newtonian/test_2d/random_pert.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/misc/simple_io.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/mpi/mpi_macro.hpp"
#include "source/mpi/MeshPointsMPI.hpp"
#include "source/tessellation/shape_2d.hpp"
#include "source/newtonian/test_2d/piecewise.hpp"
#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"
#include "source/newtonian/test_2d/consecutive_snapshots.hpp"
#include "source/newtonian/test_2d/multiple_diagnostics.hpp"
#include "source/newtonian/two_dimensional/source_terms/cylindrical_complementary.hpp"
#include "source/misc/vector_initialiser.hpp"
#include "calc_init_cond.hpp"
#include "custom_cell_updater.hpp"
#include "custom_flux_calculator.hpp"

using namespace std;
using namespace simulation2d;

namespace {

  class SimData
  {
  public:

    SimData(void):
      pg_(Vector2D(0,0), Vector2D(0,1)),
      outer_(Vector2D(0,-0.2),
	     Vector2D(0.5,1.8)),
      init_points_(cartesian_mesh(100*2,400*2,
				  outer_.getBoundary().first,
				  outer_.getBoundary().second)),
      tess_(init_points_, outer_),
      eos_(5./3.),
      point_motion_(),
      rs_(),
      force_(pg_.getAxis()),
      tsf_(0.3),
      fc_(rs_),
      eu_(),
      cu_(),
      sim_(tess_,
	   outer_,
	   pg_,
	   calc_init_cond(tess_),
	   eos_,
	   point_motion_,
	   force_,
	   tsf_,
	   fc_,
	   eu_,
	   cu_) {}

    hdsim& getSim(void)
    {
      return sim_;
    }

  private:
    const CylindricalSymmetry pg_;
    const SquareBox outer_;
    const vector<Vector2D> init_points_;
    VoronoiMesh tess_;
    PCM2D interp_method_;
    const IdealGas eos_;
    Eulerian point_motion_;
    const Hllc rs_;
    CylindricalComplementary force_;
    const SimpleCFL tsf_;
    //    const SimpleFluxCalculator fc_;
    const CustomFluxCalculator fc_;
    const SimpleExtensiveUpdater eu_;
    //    const SimpleCellUpdater cu_;
    const CustomCellUpdater cu_;
    hdsim sim_;
  };

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

