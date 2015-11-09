#include "sim_data.hpp"

SimData::SimData(void):
  pg_(Vector2D(0,0), Vector2D(0,1)),
  outer_(Vector2D(0,-0.2),
	 Vector2D(0.25,1.8)),
  init_points_(cartesian_mesh(50*2,400*2,
			      outer_.getBoundary().first,
			      outer_.getBoundary().second)),
  tess_(init_points_, outer_),
  eos_(5./3.),
  point_motion_(),
  evc_(),
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
       evc_,
       force_,
       tsf_,
       fc_,
       eu_,
       cu_) {}

hdsim& SimData::getSim(void)
{
  return sim_;
}
