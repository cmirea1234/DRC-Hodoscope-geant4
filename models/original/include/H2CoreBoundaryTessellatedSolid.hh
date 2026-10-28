#pragma once

#include "G4TessellatedSolid.hh"
#include "G4GeometryTolerance.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exception.hh"
#include <cmath>
#include <iomanip>
#include <sstream>

// Only resolves the Geant4 tessellated surface double-zero ambiguity.
// Vertices, Inside, normals, safety, extent, voxelization and materials are
// unchanged. A bounded ray query does NOT move a track or change its weight.
class H2CoreBoundaryTessellatedSolid final : public G4TessellatedSolid {
 public:
  explicit H2CoreBoundaryTessellatedSolid(const G4String& name)
      : G4TessellatedSolid(name) {}
  using G4TessellatedSolid::DistanceToIn;
  G4double DistanceToIn(const G4ThreeVector& p,
                       const G4ThreeVector& v) const override {
    const auto din=G4TessellatedSolid::DistanceToIn(p,v);
    if(din!=0.)return din;
    const auto dout=G4TessellatedSolid::DistanceToOut(p,v,false,nullptr,nullptr);
    if(dout!=0.)return din;
    if(Inside(p)!=kSurface)return din;
    const auto probe=Probe(p,v);
    G4double resolved=din;
    if(probe.state==kOutside) {
      const auto next=G4TessellatedSolid::DistanceToIn(probe.point,v);
      resolved=next>=kInfinity ? kInfinity : probe.distance+next;
    }
    if(fReports++<4) {
      std::ostringstream out;
      out << std::setprecision(17) << "H2_CORE_BOUNDARY_RESOLVED solid=" << GetName()
          << " p=" << p << " v=" << v << " inside=" << Inside(p)
          << " normal=" << SurfaceNormal(p) << " dot=" << SurfaceNormal(p).dot(v)
          << " probe_mm=" << probe.distance/mm << " forward=" << probe.state
          << " corrected_in_mm=" << resolved/mm;
      G4cout << out.str() << G4endl;
    }
    return resolved;
  }
  using G4TessellatedSolid::DistanceToOut;
  G4double DistanceToOut(const G4ThreeVector& p,const G4ThreeVector& v,
                        const G4bool calcNorm=false,G4bool* validNorm=nullptr,
                        G4ThreeVector* n=nullptr) const override {
    const auto base=G4TessellatedSolid::DistanceToOut(p,v,calcNorm,validNorm,n);
    if(base!=0. || Inside(p)!=kSurface ||
       G4TessellatedSolid::DistanceToIn(p,v)!=0.)return base;
    const auto probe=Probe(p,v);
    if(probe.state==kOutside)return base;
    // Recompute normal/validity at the *actual future exit*, not at p.
    const auto next=G4TessellatedSolid::DistanceToOut(
        probe.point,v,calcNorm,validNorm,n);
    if(next<=0. || next>=kInfinity || !std::isfinite(next)) {
      G4Exception("H2CoreBoundaryTessellatedSolid::DistanceToOut",
                  "HODO_CORE_BOUNDARY_UNRESOLVED",FatalException,
                  "Interior ray probe has no finite positive exit.");
    }
    return probe.distance+next;
  }
  G4VSolid* Clone() const override {
    return new H2CoreBoundaryTessellatedSolid(*this);
  }
 private:
  struct RayProbe { G4ThreeVector point; G4double distance; EInside state; };
  RayProbe Probe(const G4ThreeVector& p,const G4ThreeVector& v) const {
    const auto tol=G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
    for(int factor=4;factor<=128;factor*=2) {
      const auto distance=factor*tol;
      const auto q=p+distance*v;
      const auto state=G4TessellatedSolid::Inside(q);
      if(state!=kSurface)return {q,distance,state};
    }
    // No fabricated escape step or suppressed warning for an unresolved ray.
    G4ExceptionDescription msg;
    msg << "solid=" << GetName() << " p=" << p << " v=" << v;
    G4Exception("H2CoreBoundaryTessellatedSolid::Probe",
                "HODO_CORE_BOUNDARY_UNRESOLVED",FatalException,msg);
    return {p,0.,kSurface};
  }
  mutable unsigned fReports=0;
};
