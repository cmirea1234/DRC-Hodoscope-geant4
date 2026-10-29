"""Isolated observer-only boundary traces on the existing disk-photon fixture."""
from pathlib import Path
import shutil,difflib,json,hashlib
base=Path(__file__).resolve().parents[1]
root=base/'build-knu/boundary_trace_20261029';root.mkdir(exist_ok=False)
source=root/'source';shutil.copytree(base/'build-knu/photon_disk_20261029/source',source)
observer=r'''
    // Diagnostic observer only: do not change tracks, materials or RNG state.
    const int knuEvent=G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    if(atBoundary && (knuEvent==0 || knuEvent==16 || knuEvent==17 ||
                     knuEvent==20 || knuEvent==25 || knuEvent==53)) {
      G4ThreeVector normal;
      int surfaceState=-1;
      if(prePV) {
        const auto transform=prePoint->GetTouchableHandle()->GetHistory()->GetTopTransform();
        const auto localPoint=transform.TransformPoint(postPoint->GetPosition());
        const auto* solid=prePV->GetLogicalVolume()->GetSolid();
        surfaceState=static_cast<int>(solid->Inside(localPoint));
        if(surfaceState==static_cast<int>(kSurface))
          normal=transform.Inverse().TransformAxis(solid->SurfaceNormal(localPoint));
      }
      const auto rindex=[&](const G4Material* material) {
        if(!material || !material->GetMaterialPropertiesTable())return -1.;
        auto* property=material->GetMaterialPropertiesTable()->GetProperty("RINDEX");
        return property ? property->Value(track->GetTotalEnergy()) : -1.;
      };
      const auto p=postPoint->GetPosition();
      const auto d=prePoint->GetMomentumDirection();
      const auto q=postPoint->GetMomentumDirection();
      std::ostringstream s;
      s << std::setprecision(17) << "KNU_BOUNDARY"
        << " event=" << knuEvent << " step=" << track->GetCurrentStepNumber()
        << " plane=" << audit.plane << " channel=" << audit.originChannel
        << " pre=" << (prePV?prePV->GetName():"none")
        << " post=" << (postPV?postPV->GetName():"none")
        << " pre_region=" << static_cast<int>(preRegion)
        << " candidate_region=" << static_cast<int>(candidatePostRegion)
        << " status=" << static_cast<int>(status)
        << " tir=" << (status==TotalInternalReflection)
        << " reflected=" << reflected << " length_mm=" << length/mm
        << " x=" << p.x()/mm << " y=" << p.y()/mm << " z=" << p.z()/mm
        << " dx=" << d.x() << " dy=" << d.y() << " dz=" << d.z()
        << " out_dx=" << q.x() << " out_dy=" << q.y() << " out_dz=" << q.z()
        << " nx=" << normal.x() << " ny=" << normal.y() << " nz=" << normal.z()
        << " normal_mag=" << normal.mag() << " surface_state=" << surfaceState
        << " n1=" << rindex(prePoint->GetMaterial())
        << " n2=" << rindex(postPoint->GetMaterial());
      G4cout << s.str() << G4endl;
    }
'''
patch=[]
for model in ('original','light32'):
 p=source/'models'/model/'src/ActionInitialization.cc';old=p.read_text()
 s='#include "G4TouchableHistory.hh"\n#include "G4LogicalVolume.hh"\n#include "G4VSolid.hh"\n'+old
 anchor='    const G4bool reflected=atBoundary&&IsReflection(status);'
 assert s.count(anchor)==1
 # Place after length declaration so all required values exist.
 anchor='    const auto length=std::max(0.,step->GetStepLength());'
 assert s.count(anchor)==1;s=s.replace(anchor,anchor+observer)
 p.write_text(s);patch+=list(difflib.unified_diff(old.splitlines(True),s.splitlines(True),fromfile=model+'/ActionInitialization.cc',tofile=model+'/ActionInitialization.cc'))
(root/'observer.patch').write_text(''.join(patch))
(root/'source_sha256.json').write_text(json.dumps({str(p.relative_to(source)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source.rglob('*') if p.is_file()},indent=2)+'\n')
s=(base/'scripts/knu_photon_disk.sh').read_text().replace('disk_${1}_${2}_','boundary_${1}_${2}_').replace('build-knu/photon_disk_20261029','build-knu/boundary_trace_20261029').replace('disk_source.patch','observer.patch')
s+='\n.venv/bin/python scripts/check_knu_boundary_trace.py "$attempt"\n'
(base/'scripts/knu_boundary_trace.sh').write_text(s);(base/'scripts/knu_boundary_trace.sh').chmod(0o755)
s=(base/'condor/knu_photon_disk.sub').read_text().replace('knu_photon_disk.sh','knu_boundary_trace.sh').replace('/disk.$','/boundary.$')
(base/'condor/knu_boundary_trace.sub').write_text(s)
print(root)
