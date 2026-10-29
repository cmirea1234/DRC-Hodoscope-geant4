"""Build an original-based layout extension without editing the frozen payload."""
import argparse,shutil,json,hashlib,difflib
from pathlib import Path
from verify_bundle import BASE,verify
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();verify()
root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
shutil.copytree(BASE/'models/original',root/'models/original')
shutil.copy2(BASE/'extensions/layout/ModulePose.hh',root/'models/original/include/ModulePose.hh')
shutil.copy2(BASE/'extensions/layout/ModulePlacement.inc',root/'models/original/src/ModulePlacement.inc')
(root/'CMakeLists.txt').write_text('''cmake_minimum_required(VERSION 3.16)
project(OriginalLayout LANGUAGES CXX)
find_package(Geant4 11.4.1 EXACT REQUIRED)
include(${Geant4_USE_FILE})
add_executable(hodoscope-original-layout models/original/src/main.cc models/original/src/DetectorConstruction.cc models/original/src/ActionInitialization.cc models/original/src/V21GeometryTolerance.cc)
target_include_directories(hodoscope-original-layout PRIVATE models/original/include)
target_compile_features(hodoscope-original-layout PRIVATE cxx_std_17)
target_link_libraries(hodoscope-original-layout PRIVATE ${Geant4_LIBRARIES})
''')
patch=[]
def save(path,old,new):
 path.write_text(new);patch.extend(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=str(path.relative_to(root)),tofile=str(path.relative_to(root))))
h=root/'models/original/include/DetectorConstruction.hh';old=h.read_text();s='#include "ModulePose.hh"\n#include <vector>\n'+old
s=s.replace(' private:\n  void DefineMaterials();',''' public:
  G4double ValidationSurfaceDistance(int surface,const G4ThreeVector& p) const;
  G4ThreeVector ModuleDirection(int module,const G4ThreeVector& d) const {
    return fModuleLayout.applied?fModuleLayout.poses.at(module).direction(d):d;
  }
 private:
  ModuleLayout fModuleLayout;
  std::vector<std::pair<int,G4VPhysicalVolume*>> fModuleVolumes;
  void RegisterModuleVolume(int module,G4VPhysicalVolume* volume);
  void ApplyModuleLayout(G4VPhysicalVolume* worldPV);
  void DefineMaterials();''')
funcs=['PlaneNominalToWorldPosition','PlaneWorldToNominalPosition','PlaneNominalToWorldDirection','PlaneWorldToNominalDirection','TriggerWorldToLocalPosition','TriggerWorldToLocalDirection']
for name in funcs:
 anchor=f'  G4ThreeVector {name}('
 start=s.index(anchor);end=s.index(';',start)+1;s=s[:end]+'\n'+s[start:end].replace(name,'Legacy'+name)+s[end:]
anchor='  G4ThreeVector GetTriggerCentre(G4int trigger) const;';s=s.replace(anchor,anchor+'\n  G4ThreeVector LegacyGetTriggerCentre(G4int trigger) const;')
save(h,old,s)
p=root/'models/original/src/DetectorConstruction.cc';old=p.read_text();s=old
for name in funcs+['GetTriggerCentre']:s=s.replace('DetectorConstruction::'+name+'(','DetectorConstruction::Legacy'+name+'(')
# The legacy trigger inverse uses the legacy centre, not the transformed centre.
s=s.replace('const auto relative=position-GetTriggerCentre(trigger);','const auto relative=position-LegacyGetTriggerCentre(trigger);')
s=s.replace('const auto worldXY = fUseBuilding ? 2.5*m : 250.*mm;','const auto worldXY = fUseBuilding ? 2.5*m : fModuleLayout.worldHalf.x();')
s=s.replace('const auto worldZ = fUseBuilding ? 2.5*m : 250.*mm;','const auto worldZ = fUseBuilding ? 2.5*m : fModuleLayout.worldHalf.z();')
s=s.replace('new G4Box("World", worldXY, worldXY, worldZ)','new G4Box("World", worldXY, fUseBuilding?worldXY:fModuleLayout.worldHalf.y(), worldZ)')
s=s.replace('  return worldPV;','  ApplyModuleLayout(worldPV);\n  return worldPV;',1)
anchor='      if (deterministicSquareFanoutCheck)\n        priorSquareCentrelines.push_back(centres);'
assert s.count(anchor)==1
s=s.replace(anchor,'''      for(auto* pv:{outerPV,innerPV,corePV,static_cast<G4VPhysicalVolume*>(greasePV),static_cast<G4VPhysicalVolume*>(windowPV),static_cast<G4VPhysicalVolume*>(sensorPV)})
        RegisterModuleVolume(2+plane/2,pv);
'''+anchor)
anchor='  new G4PVPlacement(rotation, position, logical, "Trigger"+suffix,\n                    world, false, 0, true);';assert s.count(anchor)==1
s=s.replace(anchor,anchor.replace('  new G4PVPlacement','  auto* triggerPV=new G4PVPlacement')+'\n  if(trigger<2)RegisterModuleVolume(trigger,triggerPV);')
s+='\n#include "ModulePlacement.inc"\n'
for name in funcs:
 plane=name.startswith('Plane');arg='plane' if plane else 'trigger';vec='direction' if name.endswith('Direction') else 'position';idx='2+plane/2' if plane else 'trigger'
 forward=name.startswith('PlaneNominal');method=('direction' if vec=='direction' else 'point') if forward else ('inverseDirection' if vec=='direction' else 'inversePoint')
 body=f'Legacy{name}({arg},{vec})'
 if forward:body=f'fModuleLayout.poses[{idx}].{method}({body})'
 else:body=f'Legacy{name}({arg},fModuleLayout.poses[{idx}].{method}({vec}))'
 s+=f'''\nG4ThreeVector DetectorConstruction::{name}(G4int {arg},const G4ThreeVector& {vec}) const {{
 if(!fModuleLayout.applied || {arg}<0 || {arg}>={4 if plane else 2})return Legacy{name}({arg},{vec});
 return {body};
}}\n'''
s+='''G4ThreeVector DetectorConstruction::GetTriggerCentre(G4int trigger) const {
 auto p=LegacyGetTriggerCentre(trigger);
 return fModuleLayout.applied && trigger>=0 && trigger<2?fModuleLayout.poses[trigger].point(p):p;
}\n'''
save(p,old,s)
p=root/'models/original/src/ActionInitialization.cc';old=p.read_text();s=old
anchor='''      const G4double a=pre->GetPosition().z();
      const G4double b=post->GetPosition().z();
      const G4double delta=b-a;''';assert s.count(anchor)==1;s=s.replace(anchor,'')
anchor='        const G4double z=fDetector->GetValidationSurfaceZ(surface);'
s=s.replace(anchor,'''        const G4double a=fDetector->ValidationSurfaceDistance(surface,pre->GetPosition());
        const G4double b=fDetector->ValidationSurfaceDistance(surface,post->GetPosition());
        const G4double delta=b-a;
        const G4double z=0.;''',1);s=s.replace('const G4ThreeVector axisX(std::cos(angle),std::sin(angle),0.);','const auto axisX=fDetector->ModuleDirection(trigger,G4ThreeVector(std::cos(angle),std::sin(angle),0.));')
s=s.replace('const G4ThreeVector axisY(-std::sin(angle),std::cos(angle),0.);','const auto axisY=fDetector->ModuleDirection(trigger,G4ThreeVector(-std::sin(angle),std::cos(angle),0.));')
s=s.replace('const G4ThreeVector axisZ(0.,0.,1.);','const auto axisZ=fDetector->ModuleDirection(trigger,G4ThreeVector(0.,0.,1.));')
save(p,old,s)
(root/'layout_extension.patch').write_text(''.join(patch))
(root/'source_sha256.json').write_text(json.dumps({str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in root.rglob('*') if p.is_file()},indent=2)+'\n')
print(root)
