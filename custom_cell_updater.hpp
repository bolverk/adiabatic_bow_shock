#ifndef CUSTOM_CELL_UPDATER_HPP
#define CUSTOM_CELL_UPDATER_HPP 1

#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"

class CustomCellUpdater: public CellUpdater
{
public:

  CustomCellUpdater(void);

  vector<ComputationalCell> operator()
  (const Tessellation& /*tess*/,
   const PhysicalGeometry& /*pg*/,
   const EquationOfState& eos,
   const vector<Extensive>& extensives,
   const vector<ComputationalCell>& old,
   const CacheData& cd) const;
};

#endif // CUSTOM_CELL_UPDATER_HPP
