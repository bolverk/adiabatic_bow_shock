#include "sim_data.hpp"

SimData::SimData(void):
  pg_(Vector2D(0,0), Vector2D(0,1)),
  outer_(Vector2D(0,-0.05),
	 Vector2D(0.1,0.05)),
  init_points_(cartesian_mesh(4*100,4*100,
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

hdsim& SimData::getSim(void)
{
  return sim_;
}
