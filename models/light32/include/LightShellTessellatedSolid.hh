#pragma once
#include "G4TessellatedSolid.hh"
#include "G4Exception.hh"
#include <array>
#include <cmath>
#include <limits>

// Used by light-model H1 core and cladding shells. Same triangle mesh/material/
// voxelization as G4TessellatedSolid. Refine only
// finite entry distances whose endpoint is not on the native solid surface
// (except an already-inside start with zero entry distance).
// No surface-tolerance widening, normal displacement, or optical scaling.
class LightShellTessellatedSolid final : public G4TessellatedSolid {
 public:
  explicit LightShellTessellatedSolid(const G4String& name):G4TessellatedSolid(name) {}
  using G4TessellatedSolid::DistanceToIn;
  G4VSolid* Clone() const override { return new LightShellTessellatedSolid(*this); }

  G4double DistanceToIn(const G4ThreeVector& p,const G4ThreeVector& d) const override {
    const auto base=G4TessellatedSolid::DistanceToIn(p,d);
    return RefineEntry(*this,p,d,base);
  }

  static G4double RefineEntry(const G4TessellatedSolid& solid,
                             const G4ThreeVector& p,const G4ThreeVector& d,
                             G4double base) {
    if(!(base>=0. && base<kInfinity))return base;
    const auto state=solid.G4TessellatedSolid::Inside(p+base*d);
    if(state==kSurface || (base==0. && state==kInside))return base;
    // Native facet intersections admit tolerance-wide edge extensions. At a
    // shared edge their intersection may still be outside the closed mesh.
    // Recompute the first front-facing triangle hit without edge extension.
    long double closest=std::numeric_limits<long double>::infinity();
    for(int i=0;i<solid.GetNumberOfFacets();++i) {
      const auto* f=solid.GetFacet(i);
      if(f->GetNumberOfVertices()!=3)continue;
      const auto t=TriangleEntry(p,d,f->GetVertex(0),f->GetVertex(1),f->GetVertex(2));
      if(t>=0 && t<closest)closest=t;
    }
    if(!std::isfinite(closest))return kInfinity;
    const auto refined=static_cast<double>(closest);
    if(solid.G4TessellatedSolid::Inside(p+refined*d)!=kSurface) {
      G4ExceptionDescription message;
      message<<"Native and refined mesh entries remain inconsistent: "<<solid.GetName()
             <<" base="<<base<<" refined="<<refined<<" p="<<p<<" d="<<d;
      G4Exception("LightShellTessellatedSolid::DistanceToIn","LIGHT_SHELL_ENTRY_INCONSISTENT",
                  FatalException,message);
    }
    return refined;
  }

  static long double TriangleEntry(const G4ThreeVector& p,const G4ThreeVector& d,
                                  const G4ThreeVector& a,const G4ThreeVector& b,
                                  const G4ThreeVector& c) {
    using Vec=std::array<long double,3>;
    const auto sub=[](const G4ThreeVector& x,const G4ThreeVector& y)->Vec {
      return {static_cast<long double>(x.x())-y.x(),static_cast<long double>(x.y())-y.y(),
              static_cast<long double>(x.z())-y.z()};
    };
    const auto cross=[](const Vec& x,const Vec& y)->Vec {
      return {x[1]*y[2]-x[2]*y[1],x[2]*y[0]-x[0]*y[2],x[0]*y[1]-x[1]*y[0]};
    };
    const auto dot=[](const Vec& x,const Vec& y) {
      return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
    };
    const Vec e1=sub(b,a),e2=sub(c,a),s=sub(p,a),v{d.x(),d.y(),d.z()};
    const auto h=cross(v,e2);const auto det=dot(e1,h);
    const auto miss=std::numeric_limits<long double>::infinity();
    if(det<=0)return miss; // Front-facing entry only; tangent/exit excluded.
    const auto u=dot(s,h)/det;
    if(u<0 || u>1)return miss;
    const auto q=cross(s,e1);const auto baryV=dot(v,q)/det;
    if(baryV<0 || u+baryV>1)return miss;
    const auto t=dot(e2,q)/det;
    return t>=0?t:miss;
  }
};
