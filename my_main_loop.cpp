#include "my_main_loop.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/test_2d/multiple_diagnostics.hpp"
#include "source/misc/vector_initialiser.hpp"
#include "source/newtonian/test_2d/consecutive_snapshots.hpp"
#include "volume_appendix.hpp"
#include "temperature_appendix.hpp"
#include "write_wall_time.hpp"
#include <boost/lexical_cast.hpp>

using namespace simulation2d;

void my_main_loop(hdsim& sim)
{
  const boost::mpi::communicator world;
  const string rank_s = boost::lexical_cast<string>(world.rank());
  const double tf = 10;
  SafeTimeTermination term_cond(tf,1e6);
  MultipleDiagnostics diag
    (VectorInitialiser<DiagnosticFunction*>
     (new WriteTime("time_"+rank_s+".txt"))
     (new WriteWallTime("wall_time_"+rank_s+".txt"))
     (new ConsecutiveSnapshots
      (new ConstantTimeInterval(tf/1000),
       new Rubric("snapshot_"+rank_s+"_",".h5"),
       VectorInitialiser<DiagnosticAppendix*>
       (new TemperatureAppendix(1,1))
       (new VolumeAppendix())()
       ))());
  main_loop(sim,
	    term_cond,
	    &hdsim::TimeAdvance,
	    &diag);
  write_snapshot_to_hdf5(sim, "final_"+rank_s+".h5");
}
