#ifndef VOLUME_APPENDIX_HPP
#define VOLUME_APPENDIX_HPP 1

#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"

class VolumeAppendix: public DiagnosticAppendix
{
public:

  VolumeAppendix(void);

  string getName(void) const;

  vector<double> operator()(const hdsim& sim) const;
};

#endif // VOLUME_APPENDIX_HPP
