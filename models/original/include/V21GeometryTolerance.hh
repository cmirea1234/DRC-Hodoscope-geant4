#ifndef V21_GEOMETRY_TOLERANCE_HH
#define V21_GEOMETRY_TOLERANCE_HH

#include "globals.hh"

namespace hodoscope::v21 {

struct GeometryToleranceSnapshot {
  G4double requestedWorldMaximumExtent = 0.;
  G4double surfaceTolerance = 0.;
  G4double radialTolerance = 0.;
  G4bool configuredBeforeGeometry = false;
};

// This must be the first Geant4 geometry operation in main(): Geant4 only
// accepts SetWorldMaximumExtent() while the solid store is still empty.
void ConfigureGeometryToleranceBeforeGeometry();

// V21 is intentionally SerialOnly. The snapshot records the exact values
// requested and observed by that process before its run manager was created.
const GeometryToleranceSnapshot& GetGeometryToleranceSnapshot();

}  // namespace hodoscope::v21

#endif
