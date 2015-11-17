#include "sim_data.hpp"
#include "source/misc/vector_initialiser.hpp"
#include <boost/mpi/collectives.hpp>
#include "source/mpi/MeshPointsMPI.hpp"
#include "source/tessellation/RoundGrid.hpp"

namespace {
  vector<Vector2D> process_positions
  (const SquareBox& boundary)
  {
    const Vector2D lower_left = boundary.getBoundary().first;
    const Vector2D upper_right = boundary.getBoundary().second;
    const boost::mpi::communicator world;
    vector<Vector2D> res(static_cast<size_t>(world.size()));
    assert(res.size()>1);
    if(world.rank()==0)
      res = RandSquare
	(world.size(),
	 lower_left.x,
	 upper_right.x,
	 lower_left.y,
	 upper_right.y);
    boost::mpi::broadcast(world, res, 0);
    return RoundGrid(res,&boundary,20);
  }
}

SimData::SimData(void):
  pg_(Vector2D(0,0), Vector2D(0,1)),
  outer_(Vector2D(0,-0.5),
	 Vector2D(0.5,0.5)),
  proctess_(process_positions(outer_),outer_),
  init_points_
  (SquareMeshM
   (400,
    800,
    proctess_,
    outer_.getBoundary().first,
    outer_.getBoundary().second)),
  tess_(proctess_,init_points_, outer_),
  eos_(5./3.),
  point_motion_(),
  evc_(),
  rs_(),
  geom_source_(pg_.getAxis()),
  wind_source_
  (1e-2*4*M_PI/10/(4.0*M_PI*pow(0.01,3)/3.),
   10,
   1e-3,
   0.01),
  source_
  (VectorInitialiser<SourceTerm*>
   (&geom_source_)
   (&wind_source_)
   ()),
  tsf_(0.3),
  fc_(rs_),
  eu_(),
  cu_(),
  sim_
  (proctess_,
   tess_,
   outer_,
   pg_,
   calc_init_cond(tess_),
   eos_,
   point_motion_,
   evc_,
   source_,
   tsf_,
   fc_,
   eu_,
   cu_) {}

hdsim& SimData::getSim(void)
{
  return sim_;
}
