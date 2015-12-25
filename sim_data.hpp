#ifndef SIM_DATA_HPP
#define SIM_DATA_HPP 1

#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/two_dimensional/source_terms/cylindrical_complementary.hpp"
#include "custom_flux_calculator.hpp"
#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"
#include "custom_cell_updater.hpp"
#include "source/misc/mesh_generator.hpp"
#include "calc_init_cond.hpp"
#include "source/newtonian/two_dimensional/stationary_box.hpp"
#include "wind.hpp"
#include "source/newtonian/two_dimensional/source_terms/SeveralSources.hpp"

class SimData
{
public:

  SimData(void);

  hdsim& getSim(void);

private:
  const CylindricalSymmetry pg_;
  const SquareBox outer_;
  const vector<Vector2D> init_points_;
  VoronoiMesh tess_;
  const IdealGas eos_;
  Eulerian point_motion_;
  const StationaryBox evc_;
  const Hllc rs_;
  CylindricalComplementary geom_source_;
  Wind wind_source_;
  SeveralSources source_;
  const SimpleCFL tsf_;
  //  const SimpleFluxCalculator fc_;
  const CustomFluxCalculator fc_;
  const SimpleExtensiveUpdater eu_;
  const SimpleCellUpdater cu_;
  //  const CustomCellUpdater cu_;
  hdsim sim_;
};

#endif // SIM_DATA_HPP
