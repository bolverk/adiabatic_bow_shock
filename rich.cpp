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
#include "source/newtonian/two_dimensional/simple_flux_calculator.hpp"
#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"
#include "source/newtonian/test_2d/consecutive_snapshots.hpp"
#include "source/newtonian/test_2d/multiple_diagnostics.hpp"
#include "source/newtonian/two_dimensional/source_terms/cylindrical_complementary.hpp"
#include "source/misc/vector_initialiser.hpp"

using namespace std;
using namespace simulation2d;

namespace {

#ifdef RICH_MPI
  vector<Vector2D> process_positions(const SquareBox& boundary)
  {
    const Vector2D lower_left = boundary.getBoundary().first;
    const Vector2D upper_right = boundary.getBoundary().second;
    vector<Vector2D> res(get_mpi_size());
    if(get_mpi_rank()==0){
      res = RandSquare(get_mpi_size(),
		       lower_left.x,upper_right.x,
		       lower_left.y,upper_right.y);
    }
    MPI_VectorBcast_Vector2D(res,0,MPI_COMM_WORLD,get_mpi_rank());
    return res;
  }
#endif

  vector<ComputationalCell> calc_init_cond(const Tessellation& tess)
  {
    vector<ComputationalCell> res(static_cast<size_t>(tess.GetPointNo()));
    const Circle obstacle(Vector2D(0,0),0.01);
    for(size_t i=0;i<res.size();++i){
      const Vector2D r = tess.GetMeshPoint(static_cast<int>(i));
      res[i].density = 1;
      res[i].pressure = 1e-9;
      res[i].velocity = Vector2D(0,1);
      res[i].stickers["obstacle"] = false;
      if(obstacle(r)){
	res[i].velocity = Vector2D(0,0);
	res[i].stickers["obstacle"] = true;
      }
    }
    return res;
  }

  class CustomCellUpdater: public CellUpdater
  {
  public:

    CustomCellUpdater(void) {}

    vector<ComputationalCell> operator()
    (const Tessellation& /*tess*/,
     const PhysicalGeometry& /*pg*/,
     const EquationOfState& eos,
     const vector<Extensive>& extensives,
     const vector<ComputationalCell>& old,
     const CacheData& cd) const
    {
      vector<ComputationalCell> res = old;
      for(size_t i=0;i<extensives.size();++i){
	if(old[i].stickers.find("obstacle")->second)
	  continue;
	const double volume = cd.volumes[i];
	res[i].density = extensives[i].mass/volume;
	res[i].velocity = extensives[i].momentum / extensives[i].mass;
	const double total_energy = extensives[i].energy/extensives[i].mass;
	const double kinetic_energy =
	  0.5*ScalarProd(res[i].velocity, res[i].velocity);
	const double energy = total_energy - kinetic_energy;
	res[i].pressure = eos.de2p(res[i].density, energy);
	for(std::map<std::string,double>::const_iterator it =
	      extensives[i].tracers.begin();
	    it!=extensives[i].tracers.end();++it)
	  res[i].tracers[it->first] = it->second/extensives[i].mass;
      }
      return res;
    }
  };

  double calc_tracer_flux(const Edge& edge,
			  const Tessellation& tess,
			  const vector<ComputationalCell>& cells,
			  const string& name,
			  const Conserved& hf)
  {
    if(hf.Mass>0 &&
       edge.neighbors.first>0 &&
       edge.neighbors.first<tess.GetPointNo())
      return hf.Mass*
	cells[static_cast<size_t>(edge.neighbors.first)].tracers.find(name)->second;
    if(hf.Mass<0 &&
       edge.neighbors.second>0 &&
       edge.neighbors.second<tess.GetPointNo())
      return hf.Mass*
	cells[static_cast<size_t>(edge.neighbors.second)].tracers.find(name)->second;
    return 0;	      
  }

  class CustomFluxCalculator: public FluxCalculator
  {
  public:

    CustomFluxCalculator(const RiemannSolver& rs):
      rs_(rs) {}

    vector<Extensive> operator()
    (const Tessellation& tess,
     const vector<Vector2D>& point_velocities,
     const vector<ComputationalCell>& cells,
     const vector<Extensive>& /*extensive*/,
     const EquationOfState& eos,
     const double /*time*/,
     const double /*dt*/) const
    {
      vector<Extensive> res(tess.getAllEdges().size());
      for(size_t i=0;i<tess.getAllEdges().size();++i){
	const Conserved hydro_flux = 
	  calcHydroFlux(tess, point_velocities,
			cells, eos, i);
	res[i].mass = hydro_flux.Mass;
	res[i].momentum = hydro_flux.Momentum;
	res[i].energy = hydro_flux.Energy;
	for(std::map<std::string,double>::const_iterator it =
	      cells.front().tracers.begin();
	    it!=cells.front().tracers.end();++it)
	  res[i].tracers[it->first] = 
	    calc_tracer_flux(tess.getAllEdges()[i],
			     tess,cells,it->first,hydro_flux);
      }
      return res;
    }

  private:
    const RiemannSolver& rs_;

    const Conserved calcHydroFlux
    (const Tessellation& tess,
     const vector<Vector2D>& point_velocities,
     const vector<ComputationalCell>& cells,
     const EquationOfState& eos,
     const size_t i) const
    {
      const Edge& edge = tess.GetEdge(static_cast<int>(i));
      const std::pair<bool,bool> flags 
	(edge.neighbors.first>=0 && edge.neighbors.first<tess.GetPointNo(),
	 edge.neighbors.second>=0 && edge.neighbors.second<tess.GetPointNo());
      assert(flags.first || flags.second);
      if(!flags.first){
	const size_t right_index =
	  static_cast<size_t>(edge.neighbors.second);
	const ComputationalCell& right_cell = cells[right_index];
	if(right_cell.stickers.find("obstacle")->second)
	  return Conserved();
	const Vector2D p = Parallel(edge);
	const Primitive right = convert_to_primitive(right_cell,eos);
	//const Vector2D pos = tess.GetMeshPoint(edge.neighbors.second);
	const Primitive left = 
	  CalcPrimitive(1,
			1e-9,
			Vector2D(0,1),
			eos);
	const Vector2D n = remove_parallel_component
	  (tess.GetMeshPoint(edge.neighbors.second) -
	   edge.vertices.first, p);
	return rotate_solve_rotate_back
	  (rs_, left, right, 0, n, p);
      }
      if(!flags.second){
	const size_t left_index = 
	  static_cast<size_t>(edge.neighbors.first);
	const ComputationalCell& left_cell = cells[left_index];
	if(left_cell.stickers.find("obstacle")->second)
	  return Conserved();
	const Primitive left = convert_to_primitive(left_cell, eos);
	const Vector2D p = Parallel(edge);
	//	const Vector2D pos = tess.GetMeshPoint(edge.neighbors.first);
	const Primitive right = 
	  CalcPrimitive(1,
			1e-9,
			Vector2D(0,1),
			eos);
	const Vector2D n = remove_parallel_component
	  (edge.vertices.second - 
	   tess.GetMeshPoint(edge.neighbors.first),p);
	return rotate_solve_rotate_back
	  (rs_,left,right,0,n,p);
      }
      const size_t left_index =
	static_cast<size_t>(edge.neighbors.first);
      const size_t right_index =
	static_cast<size_t>(edge.neighbors.second);
      const ComputationalCell& left_cell = cells[left_index];
      const ComputationalCell& right_cell = cells[right_index];
      if(left_cell.stickers.find("obstacle")->second &&
	 right_cell.stickers.find("obstacle")->second)
	return Conserved();
      const Vector2D p = Parallel(edge);
      const Vector2D n =
	tess.GetMeshPoint(edge.neighbors.second) -
	tess.GetMeshPoint(edge.neighbors.first);
      const double velocity = Projection
	(tess.CalcFaceVelocity
	 (point_velocities[left_index],
	  point_velocities[right_index],
	  tess.GetCellCM(edge.neighbors.first),
	  tess.GetCellCM(edge.neighbors.second),
	  calc_centroid(edge)),n);
      if(left_cell.stickers.find("obstacle")->second){
	const Primitive right = convert_to_primitive(right_cell,eos);
	const Primitive left = reflect(right,p);
	return rotate_solve_rotate_back
	  (rs_,left,right,velocity,n,p);
      }
      if(right_cell.stickers.find("obstacle")->second){
	const Primitive left = convert_to_primitive(left_cell,eos);
	const Primitive right = reflect(left,p);
	return rotate_solve_rotate_back
	  (rs_,left,right,velocity,n,p);
      }
      const Primitive left =
	convert_to_primitive(left_cell,eos);
      const Primitive right = 
	convert_to_primitive(right_cell,eos);
      return rotate_solve_rotate_back
	(rs_,left,right,velocity,n,p);
    }
  };

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

