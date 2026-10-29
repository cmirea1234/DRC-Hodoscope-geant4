"""Create isolated timer-only source copies; never change the frozen payload."""
from pathlib import Path
import shutil,json,hashlib,difflib
base=Path(__file__).resolve().parents[1]
root=base/'build-knu/profile_20261029'
root.mkdir(exist_ok=False)
source=root/'source'; source.mkdir()
shutil.copy2(base/'CMakeLists.txt',source/'CMakeLists.txt')
shutil.copytree(base/'models',source/'models')
header=r'''#pragma once
#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <cstdlib>
namespace knuprof {
using Clock=std::chrono::steady_clock;
enum Stage {other,geometry,primary,run_begin,run_end,event_begin,event_end,stacking,nonoptical,optical,skipped_optical,waveform,audit,N};
inline const char* names[N]={"process_other","geometry_construct","primary_sampling","run_begin","run_end_write_close","event_begin","event_end_adc_and_rows","stacking_registration","nonoptical_tracking","optical_tracking","skipped_optical_tracking","waveform_peak","photon_audit"};
inline std::array<double,N> seconds{};
inline std::array<unsigned long long,N> calls{};
inline std::array<int,64> stack{};
inline int depth=0;
inline auto last=Clock::now();
inline void charge(){auto now=Clock::now();if(depth)seconds[stack[depth-1]]+=std::chrono::duration<double>(now-last).count();last=now;}
inline void start(Stage s){charge();if(depth==64)std::abort();stack[depth++]=s;++calls[s];}
inline void stop(){charge();if(!depth)std::abort();--depth;}
struct Scope{explicit Scope(Stage s){start(s);}~Scope(){stop();}};
inline void dump(){std::ofstream f("profile.json");f<<std::setprecision(15)<<"{\n";for(int i=0;i<N;++i){f<<"  \""<<names[i]<<"\": {\"exclusive_wall_seconds\": "<<seconds[i]<<", \"calls\": "<<calls[i]<<"}"<<(i+1<N?",":"")<<"\n";}f<<"}\n";}
}
'''
patch=[]
for model in ('original','light32'):
    folder=source/'models'/model
    (folder/'include/KnuProfile.hh').write_text(header)
    p=folder/'src/ActionInitialization.cc';old=p.read_text();s='#include "KnuProfile.hh"\n'+old
    anchors={
      'void GeneratePrimaries(G4Event* event) override {':'primary',
      'void BeginOfRunAction(const G4Run*) override {':'run_begin',
      'void EndOfRunAction(const G4Run* base) override {':'run_end',
      'void BeginOfEventAction(const G4Event* event) override {':'event_begin',
      'void EndOfEventAction(const G4Event* event) override {':'event_end',
      'void WritePhotonAudit(G4int eventId) {':'audit',
      'G4double WaveformPeak(const std::vector<G4double>& times,G4int channel) const {':'waveform',
      'G4ClassificationOfNewTrack OpticalLastStackingAction::ClassifyNewTrack(\n    const G4Track* track) {':'stacking'}
    for anchor,stage in anchors.items():
        assert s.count(anchor)==1,anchor
        s=s.replace(anchor,anchor+'\n    knuprof::Scope knuTimer(knuprof::'+stage+');')
    anchor='void PreUserTrackingAction(const G4Track* track) override {'
    assert s.count(anchor)==1
    s=s.replace(anchor,anchor+'''\n    const bool knuOptical=track->GetDefinition()==G4OpticalPhoton::Definition();
    const bool knuSkipped=knuOptical && fEvent->Detector()->GetSkipOpticalUntilT1T2() && !fEvent->HasTriggerCoincidence();
    knuprof::start(!knuOptical ? knuprof::nonoptical : (knuSkipped ? knuprof::skipped_optical : knuprof::optical));''')
    anchor='void PostUserTrackingAction(const G4Track* track) override {\n    fEvent->FinalizeOpticalTrack(track);'
    assert s.count(anchor)==1;s=s.replace(anchor,anchor+'\n    knuprof::stop();')
    p.write_text(s);patch+=list(difflib.unified_diff(old.splitlines(True),s.splitlines(True),fromfile=f'{model}/ActionInitialization.cc',tofile=f'{model}/ActionInitialization.cc'))
    p=folder/'src/DetectorConstruction.cc';old=p.read_text();s='#include "KnuProfile.hh"\n'+old
    anchor='G4VPhysicalVolume* DetectorConstruction::Construct() {'
    assert s.count(anchor)==1;s=s.replace(anchor,anchor+'\n  knuprof::Scope knuTimer(knuprof::geometry);')
    p.write_text(s);patch+=list(difflib.unified_diff(old.splitlines(True),s.splitlines(True),fromfile=f'{model}/DetectorConstruction.cc',tofile=f'{model}/DetectorConstruction.cc'))
    p=folder/'src/main.cc';old=p.read_text();s='#include "KnuProfile.hh"\n'+old
    s=s.replace('int main(int argc, char** argv) {','int main(int argc, char** argv) {\n  knuprof::start(knuprof::other);')
    s=s.replace('return commandStatus == 0 ? 0 : 20;','knuprof::stop();\n  knuprof::dump();\n  return commandStatus == 0 ? 0 : 20;')
    p.write_text(s);patch+=list(difflib.unified_diff(old.splitlines(True),s.splitlines(True),fromfile=f'{model}/main.cc',tofile=f'{model}/main.cc'))
(root/'instrumentation.patch').write_text(''.join(patch))
(root/'source_sha256.json').write_text(json.dumps({str(p.relative_to(source)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source.rglob('*') if p.is_file()},indent=2)+'\n')
print(source)
