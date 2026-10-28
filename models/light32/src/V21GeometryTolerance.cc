#include "V21GeometryTolerance.hh"

#include "G4Exception.hh"
#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include <iomanip>

namespace hodoscope::v21 {
namespace {
GeometryToleranceSnapshot gSnapshot;
}

void ConfigureGeometryToleranceBeforeGeometry() {
  if (gSnapshot.configuredBeforeGeometry) {
    G4Exception("ConfigureGeometryToleranceBeforeGeometry",
                "HODO_V21_TOLERANCE_REPEAT", FatalException,
                "V21 geometry tolerance configuration was requested twice.");
  }

  // This is a geometry-scale numerical setting, not an empirical efficiency
  // parameter. Geant4 11.4.1 computes Cartesian and radial tolerances as
  // worldExtent*1e-11, so the official API maps 5 m to 50 pm.
  const G4double requestedWorldMaximumExtent = 5.*m;
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(
      requestedWorldMaximumExtent);
  const auto* tolerance = G4GeometryTolerance::GetInstance();

  gSnapshot.requestedWorldMaximumExtent = requestedWorldMaximumExtent;
  gSnapshot.surfaceTolerance = tolerance->GetSurfaceTolerance();
  gSnapshot.radialTolerance = tolerance->GetRadialTolerance();
  gSnapshot.configuredBeforeGeometry = true;

  G4cout << std::setprecision(17)
         << "V21_GEOMETRY_TOLERANCE"
         << " requested_world_extent_mm="
         << gSnapshot.requestedWorldMaximumExtent/mm
         << " surface_tolerance_mm=" << gSnapshot.surfaceTolerance/mm
         << " radial_tolerance_mm=" << gSnapshot.radialTolerance/mm
         << " surface_tolerance_pm="
         << gSnapshot.surfaceTolerance/(1.e-9*mm)
         << " radial_tolerance_pm="
         << gSnapshot.radialTolerance/(1.e-9*mm)
         << G4endl;
}

const GeometryToleranceSnapshot& GetGeometryToleranceSnapshot() {
  if (!gSnapshot.configuredBeforeGeometry) {
    G4Exception("GetGeometryToleranceSnapshot",
                "HODO_V21_TOLERANCE_UNCONFIGURED", FatalException,
                "V21 geometry tolerance was not configured before geometry.");
  }
  return gSnapshot;
}

}  // namespace hodoscope::v21
