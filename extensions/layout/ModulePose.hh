#pragma once
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include <array>
#include <fstream>
#include <cstdlib>
#include <stdexcept>
#include <cmath>
struct ModulePose {
  G4ThreeVector pivot, position, angles;
  G4RotationMatrix rotation;
  bool identity=true;
  G4ThreeVector point(const G4ThreeVector& p) const {return identity?p:rotation*(p-pivot)+position;}
  G4ThreeVector inversePoint(const G4ThreeVector& p) const {return identity?p:rotation.inverse()*(p-position)+pivot;}
  G4ThreeVector direction(const G4ThreeVector& p) const {return identity?p:rotation*p;}
  G4ThreeVector inverseDirection(const G4ThreeVector& p) const {return identity?p:rotation.inverse()*p;}
};
struct ModuleLayout {
  std::array<ModulePose,4> poses;
  G4ThreeVector worldHalf{250*mm,250*mm,250*mm};
  bool configured=false, applied=false;
  ModuleLayout() {
    const std::array<G4ThreeVector,4> anchors{{{-25*mm,0,80*mm},{0,-25*mm,-75*mm},{0,0,20*mm},{0,0,-20*mm}}};
    for(int i=0;i<4;++i)poses[i].pivot=poses[i].position=anchors[i];
    const char* path=std::getenv("KNU_LAYOUT_FILE");if(!path)return;
    std::ifstream input(path);std::string magic;input>>magic;
    if(magic!="KNU_LAYOUT_V1")throw std::runtime_error("Invalid layout file");
    double x,y,z;input>>x>>y>>z;worldHalf={x*mm,y*mm,z*mm};
    for(int i=0;i<4;++i){int id;double rx,ry,rz;input>>id>>x>>y>>z>>rx>>ry>>rz;
      if(!input||id!=i)throw std::runtime_error("Invalid layout module order");
      auto& p=poses[i];p.position={x*mm,y*mm,z*mm};p.angles={rx*deg,ry*deg,rz*deg};
      p.identity=(p.position==p.pivot && rx==0 && ry==0 && rz==0);
      p.rotation.rotateX(rx*deg);p.rotation.rotateY(ry*deg);p.rotation.rotateZ(rz*deg);
    }
    std::string extra;if(input>>extra)throw std::runtime_error("Unexpected layout data");
    for(int i=0;i<3;++i)if(!std::isfinite(worldHalf[i])||worldHalf[i]<=0)throw std::runtime_error("Invalid world size");
    for(const auto& p:poses)for(int i=0;i<3;++i)if(!std::isfinite(p.position[i])||!std::isfinite(p.angles[i]))throw std::runtime_error("Non-finite pose");
    configured=true;
  }
  bool changed()const {for(const auto& p:poses)if(!p.identity)return true;return false;}
};
