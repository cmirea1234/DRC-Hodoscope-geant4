#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "HitClassifier.hh"
#include "V21GeometryTolerance.hh"

#include "G4Event.hh"
#include "G4AnalysisManager.hh"
#include "G4AutoLock.hh"
#include "G4IonisParamMat.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalParameters.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4Poisson.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Scintillation.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHandle.hh"
#include "G4PhysicalConstants.hh"
#include "G4Track.hh"
#include "G4UserStackingAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VUserPrimaryParticleInformation.hh"
#include "Randomize.hh"

#include <algorithm>
#include <atomic>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {
std::atomic<G4int> gAcceptedT1T2{0};
std::atomic<G4int> gAcceptedTriggerCoincidence{0};
std::atomic<G4int> gLastProgressBucket{0};
// The executable is SerialOnly and a geometry-probe process is deliberately
// single-run.  Keeping this process-wide makes an accidental second beamOn
// fail before analysis booking can silently shift hard-coded legacy IDs.
std::atomic<G4int> gGeometryProbeRunBegins{0};
std::atomic<G4int> gGenerationProbeRunBegins{0};
std::atomic<G4int> gGenerationProbeLoadedRows{0};
std::atomic<G4double> gGenerationObservedMaxCoreStep{0.};
G4Mutex gAnalysisBookingMutex = G4MUTEX_INITIALIZER;

constexpr G4int kPrimaryTriggerCrossingNtupleId = 12;
constexpr G4int kValidationTriggerMapNtupleId = 13;
constexpr G4int kValidationGeometryProvenanceNtupleId = 14;

// Stage-3 is an isolated evidence lane: Gate-1/2 trees must be disabled, so
// its four trees append directly after the unchanged legacy IDs 0..6.
constexpr G4int kValidationGenerationStepNtupleId = 7;
constexpr G4int kValidationGenerationSourceNtupleId = 8;
constexpr G4int kStage3GenerationEventNtupleId = 9;
constexpr G4int kValidationGenerationProvenanceNtupleId = 10;

void ObserveGenerationCoreStepLength(G4double length) {
  auto observed=gGenerationObservedMaxCoreStep.load(std::memory_order_relaxed);
  while(length>observed &&
        !gGenerationObservedMaxCoreStep.compare_exchange_weak(
            observed,length,std::memory_order_relaxed)) {}
}

G4double DeterministicUnitHash(std::uint64_t value) {
  // SplitMix64 finalizer.  This deliberately does not use Geant4's global RNG:
  // changing an effective-light scan point must not alter subsequent primary
  // generation or optical-transport random streams in same-seed comparisons.
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  value = value ^ (value >> 31);
  return static_cast<G4double>(value >> 11) *
      (1.0/9007199254740992.0);
}

G4bool KeepEffectiveLightPe(
    G4int eventId, G4int channel, std::size_t photonIndex, G4double scale) {
  if(scale>=1.) return true;
  if(scale<=0.) return false;
  const auto key =
      (static_cast<std::uint64_t>(static_cast<std::uint32_t>(eventId)) << 32) ^
      (static_cast<std::uint64_t>(static_cast<std::uint32_t>(channel)) << 16) ^
      static_cast<std::uint64_t>(photonIndex);
  return DeterministicUnitHash(key) < scale;
}

G4double ModifiedGaisserLogEWeight(G4double energyGeV, G4double cosTheta) {
  const G4double safeCos = std::max(1.e-6, cosTheta);
  const G4double lowEnergy =
      1. + 3.64/(energyGeV*std::pow(safeCos, 1.29));
  const G4double shiftedEnergy = energyGeV*lowEnergy;
  const G4double spectrum =
      0.14*std::pow(shiftedEnergy, -2.7)*
      (1./(1.+1.1*energyGeV*safeCos/115.) +
       0.054/(1.+1.1*energyGeV*safeCos/850.));
  // Sampling in log(E): dE = E d(log E).
  return std::max(0., energyGeV*spectrum);
}

// V15 audit codes are intentionally integer-only so the ROOT schema remains
// stable and compact. The audit is diagnostic-only and never draws random
// numbers or changes a track.
enum class AuditRegion : G4int {
  unknown = 0,
  active_core = 1,
  fanout_core = 2,
  active_inner_clad = 3,
  fanout_inner_clad = 4,
  active_outer_clad = 5,
  fanout_outer_clad = 6,
  grease = 7,
  window = 8,
  sipm = 9,
  world = 10,
  other = 11,
  far_end_core = 12,
  far_end_inner_clad = 13,
  far_end_outer_clad = 14
};

enum class TerminalReason : G4int {
  unclassified = 0,
  sipm_detected = 1,
  sipm_pde_rejected = 2,
  bulk_core = 3,
  bulk_inner_clad = 4,
  bulk_outer_clad = 5,
  bulk_coupling = 6,
  bulk_other = 7,
  world_escape = 8,
  degradation_stress = 9,
  degradation_crack = 10,
  degradation_packing = 11,
  degradation_adhesive = 12,
  degradation_h1_cladding = 13,
  degradation_h2_cladding = 14,
  fanout_rejected = 15,
  endface_rejected = 16,
  coupling_rejected = 17,
  boundary_absorption = 18,
  other_terminal = 19
};

enum AuditStage : std::uint32_t {
  stage_active_exit = 1u << 0,
  stage_core_to_grease = 1u << 1,
  stage_clad_to_grease = 1u << 2,
  stage_window_entry = 1u << 3,
  stage_sipm_boundary = 1u << 4
};

// Preserve V15's ancestor lookup for nested geometry, then fall back to the
// current fan-out physical-volume name/copy for V19 sibling Boolean shells.
// The two-pass order is essential: a nested core still has copy 0, while its
// H1_X/H1_Y ancestor carries the real channel. Thus the V18-compatible nested
// case is unchanged and sibling-shell geometry uses the fallback.
void ResolvePlaneChannel(const G4TouchableHandle& touch,
                         G4int& plane,G4int& globalChannel) {
  plane=-1;
  globalChannel=-1;
  if(!touch)return;
  for(G4int depth=0;depth<=touch->GetHistoryDepth();++depth) {
    const auto* volume=touch->GetVolume(depth);
    if(!volume)continue;
    const auto& name=volume->GetName();
    const auto copy=touch->GetCopyNumber(depth);
    if(name=="H1_X") {
      plane=0;globalChannel=copy;return;
    }
    if(name=="H1_Y") {
      plane=1;globalChannel=copy;return;
    }
    if(name=="H2_X") {
      plane=2;globalChannel=64+copy;return;
    }
    if(name=="H2_Y") {
      plane=3;globalChannel=64+copy;return;
    }
  }
  const auto* current=touch->GetVolume(0);
  if(!current)return;
  const auto& name=current->GetName();
  const auto copy=touch->GetCopyNumber(0);
  if(name.rfind("H1XFanout",0)==0) {
    plane=0;globalChannel=copy;return;
  }
  if(name.rfind("H1YFanout",0)==0) {
    plane=1;globalChannel=copy;return;
  }
  if(name.rfind("H2XFanout",0)==0) {
    plane=2;globalChannel=64+copy;return;
  }
  if(name.rfind("H2YFanout",0)==0) {
    plane=3;globalChannel=64+copy;return;
  }
}

struct PhotonBirthContext {
  G4int parentTrack = -1;
  G4int parentStep = -1;
  G4int plane = -1;
  G4int globalChannel = -1;
  G4ThreeVector position{};
};

struct PhotonAudit {
  G4int track = -1;
  G4int parent = -1;
  G4int plane = -1;
  G4int originChannel = -1;
  G4int originPhysicalCopy = -1;
  G4int originLocalChannel = -1;
  G4int terminalGlobalChannel = -1;
  G4int birthContextMatched = 0;
  G4int parentStep = -1;
  G4int creatorProcess = 4;
  AuditRegion originRegion = AuditRegion::unknown;
  AuditRegion lastRegion = AuditRegion::unknown;
  AuditRegion terminalRegion = AuditRegion::unknown;
  TerminalReason terminal = TerminalReason::unclassified;
  std::uint32_t stageMask = 0;
  G4bool everCladding = false;
  G4bool reenteredCore = false;
  G4bool wasInCore = false;
  G4bool wasInCladding = false;
  G4double wavelength = 0.;
  G4double originTime = 0.;
  G4double terminalTime = 0.;
  G4ThreeVector origin{};
  G4ThreeVector initialDirection{};
  G4ThreeVector terminalPosition{};
  G4double activeCorePath = 0.;
  G4double fanoutCorePath = 0.;
  G4double activeInnerCladPath = 0.;
  G4double fanoutInnerCladPath = 0.;
  G4double activeOuterCladPath = 0.;
  G4double fanoutOuterCladPath = 0.;
  G4double greasePath = 0.;
  G4double windowPath = 0.;
  G4double otherPath = 0.;
  G4double totalTrackLength = 0.;
  G4int coreToClad = 0;
  G4int cladToCore = 0;
  G4int coreGreaseAttempts = 0;
  G4int cladGreaseAttempts = 0;
  G4int greaseWindowAttempts = 0;
  G4int sipmAttempts = 0;
  G4int sipmReflections = 0;
  G4int tirCount = 0;
  G4int fresnelReflectionCount = 0;
  G4int otherReflectionCount = 0;
  G4int duplicateTerminalMarks = 0;
  G4int steps = 0;
  G4bool gate4BoundaryTraceRecorded = false;
};

enum class FibreVolumeClass : G4int {
  unknown = 0,
  core = 1,
  inner_cladding = 2,
  outer_cladding = 3
};

struct FibreCrossing {
  G4int track = -1;
  G4int parent = -1;
  G4int pdg = 0;
  G4int crossingIndex = -1;
  G4int plane = -1;
  G4int physicalCopy = -1;
  G4int localChannel = -1;
  G4int globalChannel = -1;
  FibreVolumeClass volumeClass = FibreVolumeClass::unknown;
  G4int firstStep = -1;
  G4int lastStep = -1;
  G4ThreeVector entryWorld{};
  G4ThreeVector exitWorld{};
  G4ThreeVector entryLocal{};
  G4ThreeVector exitLocal{};
  G4ThreeVector entryDirection{};
  G4ThreeVector exitDirection{};
  G4double entryEnergy = 0.;
  G4double exitEnergy = 0.;
  G4double entryTime = 0.;
  G4double exitTime = 0.;
  G4double path = 0.;
  G4double corePath = 0.;
  G4double activeCorePath = 0.;
  G4double farEndCorePath = 0.;
  G4double edep = 0.;
};

struct PrimarySurfaceState {
  G4bool reached = false;
  G4int track = 1;
  G4ThreeVector position{};
  G4ThreeVector direction{};
  G4double kineticEnergy = 0.;
  G4double time = 0.;
};

struct PrimaryTriggerCrossing {
  G4bool reached = false;
  G4int track = 1;
  G4int parent = 0;
  G4int pdg = 0;
  G4int steps = 0;
  G4ThreeVector entryWorld{};
  G4ThreeVector exitWorld{};
  G4ThreeVector entryDirection{};
  G4ThreeVector exitDirection{};
  G4double entryTime = 0.;
  G4double exitTime = 0.;
  G4double path = 0.;
  G4double edep = 0.;
};

struct GenerationProbeRow {
  G4int event=-1;
  G4int cohort=-1;
  G4int targetPlane=-1;
  G4int targetChannel=-1;
  G4int expectedVolumeClass=-1;
  G4int expectedRegion=-1;
  G4int chordBin=-1;
  G4int repeatIndex=-1;
  G4ThreeVector position{};
  G4ThreeVector direction{};
  G4double energy=0.;
  G4double analyticCoreChord=0.;
  G4double analyticTargetPath=0.;
  G4int caseIndex=-1;
};

class GenerationProbePrimaryInfo final
    : public G4VUserPrimaryParticleInformation {
 public:
  explicit GenerationProbePrimaryInfo(const GenerationProbeRow& row)
      : value(row) {}
  void Print() const override {
    G4cout << "[gate3 source] row=" << value.event
           << " cohort=" << value.cohort
           << " plane=" << value.targetPlane
           << " channel=" << value.targetChannel << G4endl;
  }
  GenerationProbeRow value;
};

struct GenerationSecondaryCounts {
  G4int scintillation=0;
  G4int cerenkov=0;
  G4int other=0;
  G4int scintillationActive=0;
  G4int scintillationFanout=0;
  G4int cerenkovActive=0;
  G4int cerenkovFanout=0;
  G4int otherActive=0;
  G4int otherFanout=0;
  // Independent evidence that active/fan-out creator counts were classified
  // from each secondary's actual birth coordinate, not the parent pre-step.
  // -1 is the fail-closed sentinel when the corresponding count is zero.
  G4double activeMaxAbsLongitudinal=-1.;
  G4double fanoutMinAbsLongitudinal=-1.;
  G4int Total() const { return scintillation+cerenkov+other; }
};

struct GenerationPlaneLedger {
  G4double primaryCorePath=0.;
  G4double primaryActiveCorePath=0.;
  G4double primaryCladdingPath=0.;
  G4double primaryActiveCladdingPath=0.;
  G4double primaryCoreEdep=0.;
  G4double primaryCoreNiel=0.;
  G4double primaryCoreVisibleEdep=0.;
  G4double allCoreEdep=0.;
  G4double allCoreNiel=0.;
  G4double primaryCladdingEdep=0.;
  G4double primaryCladdingNiel=0.;
  G4double allCladdingEdep=0.;
  G4double allCladdingNiel=0.;
  G4int primaryCoreNonzeroNielSteps=0;
  G4int primaryScintillation=0;
  G4int secondaryScintillation=0;
  G4int primaryCerenkov=0;
  G4int secondaryCerenkov=0;
  G4int primaryOther=0;
  G4int secondaryOther=0;
  G4int scintillationActive=0;
  G4int scintillationFanout=0;
  G4int cerenkovActive=0;
  G4int cerenkovFanout=0;
  G4int otherActive=0;
  G4int otherFanout=0;
  G4int claddingScintillation=0;
  G4int claddingCerenkov=0;
  G4int claddingOther=0;
  G4int directOpticalTotal=0;
  G4int chargedFibreSteps=0;
  G4int primaryCoreSteps=0;
};

class PrimaryGenerator final : public G4VUserPrimaryGeneratorAction {
 public:
  explicit PrimaryGenerator(const DetectorConstruction* detector) : fDetector(detector) {
    fGun.SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("mu-"));
    fGun.SetParticleEnergy(4.*GeV);
    fGun.SetParticlePosition({0., 0., 150.*mm});
    fGun.SetParticleMomentumDirection({0., 0., -1.});
    if(fDetector->GetGeometryProbeMode())LoadGeometryProbeTable();
    if(fDetector->GetGenerationProbeMode())LoadGenerationProbeTable();
  }
  void GeneratePrimaries(G4Event* event) override {
    // Preserve a single global ID across primary/transport/digitization and
    // all truth tables when replaying a contiguous serial chunk. G4Run still
    // counts the actual number of processed events independently of this ID.
    const auto offset=fDetector->GetReplayEventIdOffset();
    if(offset>0) {
      if(!fDetector->GetCosmicBeam() ||
         fDetector->GetTransportSeedEventOffset()!=0 ||
         event->GetEventID()>std::numeric_limits<G4int>::max()-offset)
        G4Exception("PrimaryGenerator::GeneratePrimaries",
                    "HODO_REPLAY_OFFSET",FatalException,
                    "Global replay offset requires cosmicBeam and transportSeedEventOffset=0 without ID overflow.");
      event->SetEventID(event->GetEventID()+offset);
    }
    const G4double sourceZ = fDetector->GetUseBuilding() ? 2000.*mm : 150.*mm;
    if(fDetector->GetGate6ElectronicsProbeMode()) {
      auto* geantino=G4ParticleTable::GetParticleTable()->FindParticle("geantino");
      if(!geantino)
        G4Exception("PrimaryGenerator::GeneratePrimaries",
                    "HODO_GATE6_GEANTINO",FatalException,
                    "Geantino particle definition is unavailable.");
      fGun.SetParticleDefinition(geantino);
      fGun.SetParticleEnergy(1.*GeV);
      fGun.SetParticlePosition({0.,0.,140.*mm});
      fGun.SetParticleMomentumDirection({0.,0.,1.});
      fGun.GeneratePrimaryVertex(event);
      return;
    } else if(fDetector->GetGeometryProbeMode()) {
      if(fGeometryProbeRows.empty())LoadGeometryProbeTable();
      const auto rowIndex=static_cast<std::size_t>(event->GetEventID())+
          static_cast<std::size_t>(fDetector->GetGeometryProbeRowOffset());
      if(rowIndex>=fGeometryProbeRows.size()) {
        G4ExceptionDescription message;
        message << "Geometry-probe event " << event->GetEventID()
                << " with row offset "
                << fDetector->GetGeometryProbeRowOffset()
                << " exceeds table rows=" << fGeometryProbeRows.size() << ".";
        G4Exception("PrimaryGenerator::GeneratePrimaries",
                    "HODO_GATE2_PROBE_RANGE",FatalException,message);
      }
      const auto& row=fGeometryProbeRows[rowIndex];
      auto* geantino=G4ParticleTable::GetParticleTable()->FindParticle("geantino");
      if(!geantino)
        G4Exception("PrimaryGenerator::GeneratePrimaries",
                    "HODO_GATE2_GEANTINO",FatalException,
                    "Geantino particle definition is unavailable.");
      fGun.SetParticleDefinition(geantino);
      fGun.SetParticleEnergy(1.*GeV);
      fGun.SetParticlePosition(row.position);
      fGun.SetParticleMomentumDirection(row.direction);
      fGun.GeneratePrimaryVertex(event);
      return;
    } else if(fDetector->GetGenerationProbeMode()) {
      if(fGenerationProbeRows.empty())LoadGenerationProbeTable();
      const auto rowIndex=static_cast<std::size_t>(event->GetEventID())+
          static_cast<std::size_t>(fDetector->GetGenerationProbeRowOffset());
      if(rowIndex>=fGenerationProbeRows.size()) {
        G4ExceptionDescription message;
        message << "Generation-probe event " << event->GetEventID()
                << " with row offset "
                << fDetector->GetGenerationProbeRowOffset()
                << " exceeds table rows=" << fGenerationProbeRows.size()
                << ".";
        G4Exception("PrimaryGenerator::GeneratePrimaries",
                    "HODO_GATE3_PROBE_RANGE",FatalException,message);
      }
      const auto& row=fGenerationProbeRows[rowIndex];
      auto* muon=G4ParticleTable::GetParticleTable()->FindParticle("mu-");
      if(!muon)
        G4Exception("PrimaryGenerator::GeneratePrimaries",
                    "HODO_GATE3_MUON",FatalException,
                    "mu- particle definition is unavailable.");
      fGun.SetParticleDefinition(muon);
      fGun.SetParticleEnergy(row.energy);
      fGun.SetParticlePosition(row.position);
      fGun.SetParticleMomentumDirection(row.direction);
      fGun.GeneratePrimaryVertex(event);
      auto* vertex=event->GetPrimaryVertex();
      auto* primary=vertex?vertex->GetPrimary():nullptr;
      if(!primary)
        G4Exception("PrimaryGenerator::GeneratePrimaries",
                    "HODO_GATE3_PRIMARY_INFO",FatalException,
                    "Generated Stage-3 primary is unavailable.");
      primary->SetUserInformation(new GenerationProbePrimaryInfo(row));
      return;
    } else if(fDetector->GetFixedPhotonTest()) {
      const auto uniform = [event, this](std::uint64_t stream) {
        std::uint64_t value=static_cast<std::uint64_t>(
            fDetector->GetCosmicSeed());
        value += 0x9e3779b97f4a7c15ULL*
                 (static_cast<std::uint64_t>(event->GetEventID())+1);
        value += 0xbf58476d1ce4e5b9ULL*(stream+101);
        value=(value^(value>>30))*0xbf58476d1ce4e5b9ULL;
        value=(value^(value>>27))*0x94d049bb133111ebULL;
        value^=value>>31;
        return (value>>11)*(1.0/9007199254740992.0);
      };
      const auto plane=fDetector->GetFixedPhotonPlane();
      const G4bool xAxis=plane=="h1x"||plane=="h2x";
      const G4int planeIndex=plane=="h1x" ? 0 : plane=="h1y" ? 1
                             : plane=="h2x" ? 2 : 3;
      const G4int channel=std::clamp(
          fDetector->GetFixedPhotonChannelForEvent(event->GetEventID()),0,31);
      const G4double transverse=
          fDetector->GetPlaneSourceTransverse(planeIndex,channel);
      const G4double z=fDetector->GetPlaneSourceZ(planeIndex,channel);
      G4double phi=CLHEP::twopi*uniform(0);
      G4double longitudinal=2.*uniform(1)-1.;
      const auto& angularMode=fDetector->GetFixedPhotonAngularMode();
      if(angularMode=="forward")
        longitudinal=uniform(1);
      else if(angularMode=="axis") {
        longitudinal=1.;
        phi=0.;
      } else if(angularMode=="polar") {
        longitudinal=std::cos(fDetector->GetFixedPhotonPolarAngle());
        phi=fDetector->GetFixedPhotonAzimuth();
      }
      const G4double radial=std::sqrt(std::max(0.,1.-longitudinal*longitudinal));
      const G4ThreeVector nominalDirection=xAxis
          ? G4ThreeVector(longitudinal,radial*std::cos(phi),
                          radial*std::sin(phi))
          : G4ThreeVector(radial*std::cos(phi),longitudinal,
                          radial*std::sin(phi));
      G4ThreeVector direction=fDetector->PlaneNominalToWorldDirection(
          planeIndex,nominalDirection);
      fGun.SetParticleDefinition(G4OpticalPhoton::Definition());
      G4double wavelength=fDetector->GetFixedPhotonWavelength();
      if(fDetector->GetFixedPhotonSpectrum()=="bcf20_5point") {
        constexpr std::array<G4double,5> wavelengths{
          521.,490.,459.,431.,400.};
        constexpr std::array<G4double,5> weights{.30,1.00,.42,.07,.01};
        const G4double pick=uniform(2)*
            std::accumulate(weights.begin(),weights.end(),0.);
        G4double cumulative=0.;
        for(std::size_t i=0;i<weights.size();++i) {
          cumulative+=weights[i];
          if(pick<=cumulative) { wavelength=wavelengths[i]*nm;break; }
        }
      }
      fGun.SetParticleEnergy(h_Planck*c_light/wavelength);
      const auto longitudinalOffset=
          fDetector->GetFixedPhotonLongitudinalOffset();
      const G4ThreeVector nominalPosition=xAxis
          ? G4ThreeVector(longitudinalOffset,transverse,z)
          : G4ThreeVector(transverse,longitudinalOffset,z);
      G4ThreeVector position=fDetector->PlaneNominalToWorldPosition(
          planeIndex,nominalPosition);
      if(fDetector->GetFixedPhotonSource()=="terminal") {
        const G4int globalChannel=planeIndex*32+channel;
        if(!fDetector->HasFanoutEndpoint(globalChannel)) {
          G4ExceptionDescription message;
          message << "No constructed fan-out endpoint for global channel "
                  << globalChannel << ".";
          G4Exception("PrimaryGenerator::GeneratePrimaries",
                      "HODO_GATE1_ENDPOINT_SOURCE",FatalException,message);
        }
        direction=fDetector->GetFanoutTerminalDirection(globalChannel).unit();
        // Start well inside the core compared with the configured surface
        // tolerance, while remaining negligible beside the 0.10 mm grease gap.
        position=fDetector->GetFanoutEndpoint(globalChannel)-direction*1.*um;
      }
      G4ThreeVector reference(0.,0.,1.);
      if(std::abs(direction.dot(reference))>.99)reference={1.,0.,0.};
      const auto polarization=direction.cross(reference).unit();
      fGun.SetParticlePosition(position);
      fGun.SetParticleMomentumDirection(direction);
      fGun.SetParticlePolarization(polarization);
    } else if (fDetector->GetCosmicBeam()) {
      const auto uniform = [event, this](std::uint64_t stream) {
        std::uint64_t value=static_cast<std::uint64_t>(
            fDetector->GetCosmicSeed());
        value += 0x9e3779b97f4a7c15ULL*
                 (static_cast<std::uint64_t>(event->GetEventID())+1);
        value += 0xbf58476d1ce4e5b9ULL*(stream+1);
        value=(value^(value>>30))*0xbf58476d1ce4e5b9ULL;
        value=(value^(value>>27))*0x94d049bb133111ebULL;
        value^=value>>31;
        return (value>>11)*(1.0/9007199254740992.0);
      };
      const G4double cmin = std::cos(fDetector->GetMaxZenith());
      G4double zenithPower = 2.;
      const auto& angularModel = fDetector->GetCosmicAngularDistribution();
      G4ThreeVector direction;
      if (angularModel == "modified_gaisser_joint") {
        const G4double eMin =
            std::max(0.01, fDetector->GetCosmicMinEnergy()/GeV);
        const G4double eMax =
            std::max(eMin*1.001, fDetector->GetCosmicMaxEnergy()/GeV);
        EnsureJointGaisserCache(cmin, eMin, eMax);
        const auto& cache = fJointGaisserCache;
        const G4double pick = uniform(0)*cache.total;
        const auto it = std::lower_bound(
            cache.cumulative.begin(), cache.cumulative.end(), pick);
        const G4int flatBin = std::min<G4int>(
            cache.nCos*cache.nEnergy - 1,
            std::distance(cache.cumulative.begin(), it));
        const G4int cosBin = flatBin / cache.nEnergy;
        const G4int energyBin = flatBin % cache.nEnergy;
        const G4double cosTheta =
            cache.cmin + (cosBin + uniform(1))*cache.dCos;
        const G4double theta = std::acos(std::clamp(cosTheta, -1., 1.));
        const G4double phi = 360.*deg*uniform(2);
        direction = G4ThreeVector(std::sin(theta)*std::cos(phi),
                                  std::sin(theta)*std::sin(phi),
                                  -cosTheta);
        const G4double energy =
            std::exp(cache.logMin + (energyBin + uniform(3))*cache.dLog);
        fGun.SetParticleEnergy(energy*GeV);
      } else {
      if (angularModel == "uniform_solid_angle") {
        zenithPower = 0.;
      } else if (angularModel == "flat_plane_cos2") {
        // Horizontal generation plane: an angular intensity I(theta) ∝ cos^2
        // contributes one extra projection factor cos(theta), so the
        // generated plane-crossing event distribution is
        // p(theta) ∝ cos^3(theta) sin(theta), equivalently
        // p(cos(theta)) ∝ cos^3(theta).
        zenithPower = 3.;
      } else if (angularModel == "cos_power") {
        zenithPower = fDetector->GetCosmicZenithPower();
      }
      // Sample p(cos(theta)) proportional to cos^n(theta) over
      // cos(maxZenith)..1.  The default n=2 preserves the predecessor
      // hard-coded cbrt expression; flat_plane_cos2 uses n=3 for the
      // horizontal-plane projection correction.
      const G4double exponent = zenithPower + 1.;
      const G4double cminPow = std::pow(cmin, exponent);
      const G4double cosTheta =
          std::pow(cminPow + uniform(0)*(1.-cminPow), 1./exponent);
      const G4double theta = std::acos(cosTheta);
      const G4double phi = 360.*deg*uniform(1);
      direction = G4ThreeVector(std::sin(theta)*std::cos(phi),
                                std::sin(theta)*std::sin(phi), -cosTheta);
      if (fDetector->GetCosmicEnergySpectrum() == "modified_gaisser") {
        const G4double eMin = std::max(0.01, fDetector->GetCosmicMinEnergy()/GeV);
        const G4double eMax = std::max(eMin*1.001,
                                      fDetector->GetCosmicMaxEnergy()/GeV);
        constexpr G4int bins = 512;
        const G4double logMin = std::log(eMin);
        const G4double dLog = (std::log(eMax)-logMin)/bins;
        std::array<G4double,bins> cumulative{};
        G4double total = 0.;
        for (G4int i=0; i<bins; ++i) {
          const G4double energy = std::exp(logMin+(i+.5)*dLog);
          // Modified Gaisser sea-level differential spectrum.  The extra
          // energy factor is the Jacobian dE = E d(log E).
          total += ModifiedGaisserLogEWeight(energy, cosTheta);
          cumulative[i] = total;
        }
        const G4double pick = uniform(4)*total;
        const auto it = std::lower_bound(cumulative.begin(), cumulative.end(), pick);
        const G4int bin = std::min<G4int>(
            bins-1, std::distance(cumulative.begin(), it));
        const G4double energy =
            std::exp(logMin+(bin+uniform(5))*dLog);
        fGun.SetParticleEnergy(energy*GeV);
      }
      }
      const G4int targetStreamX =
          angularModel == "modified_gaisser_joint" ? 4 : 2;
      const G4int targetStreamY =
          angularModel == "modified_gaisser_joint" ? 5 : 3;
      const G4double targetX =
          (2.*uniform(targetStreamX)-1.)*fDetector->GetCosmicTargetHalfX();
      const G4double targetY =
          (2.*uniform(targetStreamY)-1.)*fDetector->GetCosmicTargetHalfY();
      const G4double scale =
          (sourceZ-fDetector->GetTriggerCentre(0).z())/(-direction.z());
      fGun.SetParticleMomentumDirection(direction);
      if(fDetector->GetCosmicPositionModel()=="uniform_source_plane") {
        // Streams 6/7 are distinct from all direction/energy streams above.
        // No detector-target projection, ray rejection or direction steering.
        fGun.SetParticlePosition({(2.*uniform(6)-1.)*fDetector->GetCosmicSourceHalfX(),
                                 (2.*uniform(7)-1.)*fDetector->GetCosmicSourceHalfY(),sourceZ});
      } else {
        fGun.SetParticlePosition({targetX-scale*direction.x(),
                                  targetY-scale*direction.y(), sourceZ});
      }
    } else if (fDetector->GetUseBuilding()) {
      const auto old = fGun.GetParticlePosition();
      fGun.SetParticlePosition({old.x(), old.y(), sourceZ});
    }
    fGun.GeneratePrimaryVertex(event);
  }
 private:
  struct JointGaisserCache {
    G4bool ready = false;
    G4double cmin = 0.;
    G4double eMin = 0.;
    G4double eMax = 0.;
    G4double logMin = 0.;
    G4double dLog = 0.;
    G4double dCos = 0.;
    G4double total = 0.;
    G4int nCos = 128;
    G4int nEnergy = 512;
    std::vector<G4double> cumulative;
  };

  void EnsureJointGaisserCache(G4double cmin, G4double eMin, G4double eMax) {
    if (fJointGaisserCache.ready &&
        std::abs(fJointGaisserCache.cmin-cmin) < 1.e-12 &&
        std::abs(fJointGaisserCache.eMin-eMin) < 1.e-12 &&
        std::abs(fJointGaisserCache.eMax-eMax) < 1.e-9) {
      return;
    }
    auto& cache = fJointGaisserCache;
    cache.ready = false;
    cache.cmin = cmin;
    cache.eMin = eMin;
    cache.eMax = eMax;
    cache.logMin = std::log(eMin);
    cache.dLog = (std::log(eMax)-cache.logMin)/cache.nEnergy;
    cache.dCos = (1.-cmin)/cache.nCos;
    cache.total = 0.;
    cache.cumulative.assign(
        static_cast<std::size_t>(cache.nCos*cache.nEnergy), 0.);
    for (G4int ix=0; ix<cache.nCos; ++ix) {
      const G4double cosTheta = cmin + (ix+.5)*cache.dCos;
      for (G4int ie=0; ie<cache.nEnergy; ++ie) {
        const G4double energy = std::exp(cache.logMin+(ie+.5)*cache.dLog);
        // Horizontal generation plane flux:
        // p(E, cos(theta)) ∝ dI/dE(E,theta) * cos(theta).
        // The logE Jacobian is already included by ModifiedGaisserLogEWeight.
        cache.total += cosTheta*ModifiedGaisserLogEWeight(energy, cosTheta);
        cache.cumulative[static_cast<std::size_t>(ix*cache.nEnergy+ie)] =
            cache.total;
      }
    }
    if (!(cache.total > 0.) || !std::isfinite(cache.total)) {
      G4Exception("PrimaryGenerator::EnsureJointGaisserCache",
                  "HODO_COSMIC_JOINT_CACHE", FatalException,
                  "modified_gaisser_joint cache has non-positive total.");
    }
    cache.ready = true;
  }

  struct GeometryProbeRow {
    G4int event=-1;
    G4int cohort=-1;
    G4int targetPlane=-1;
    G4int targetChannel=-1;
    G4int expectedClass=-1;
    G4ThreeVector position{};
    G4ThreeVector direction{};
  };
  void LoadGeometryProbeTable() {
    std::ifstream input(fDetector->GetGeometryProbeFile());
    if(!input) {
      G4ExceptionDescription message;
      message << "Cannot open geometry-probe TSV: "
              << fDetector->GetGeometryProbeFile();
      G4Exception("PrimaryGenerator::LoadGeometryProbeTable",
                  "HODO_GATE2_PROBE_OPEN",FatalException,message);
    }
    std::string line;
    G4int lineNumber=0;
    std::array<G4int,8> cohortRows{};
    while(std::getline(input,line)) {
      ++lineNumber;
      const auto first=line.find_first_not_of(" \t\r");
      if(first==std::string::npos||line[first]=='#')continue;
      GeometryProbeRow row;
      G4double x=0.,y=0.,z=0.,dx=0.,dy=0.,dz=0.;
      std::istringstream parser(line);
      if(!(parser>>row.event>>row.cohort>>row.targetPlane>>row.targetChannel
                 >>row.expectedClass>>x>>y>>z>>dx>>dy>>dz)) {
        G4ExceptionDescription message;
        message << "Malformed geometry-probe TSV row at line " << lineNumber
                << ": " << line;
        G4Exception("PrimaryGenerator::LoadGeometryProbeTable",
                    "HODO_GATE2_PROBE_PARSE",FatalException,message);
      }
      std::string extra;
      if(parser>>extra) {
        G4ExceptionDescription message;
        message << "Extra geometry-probe field at line " << lineNumber
                << ": " << extra;
        G4Exception("PrimaryGenerator::LoadGeometryProbeTable",
                    "HODO_GATE2_PROBE_EXTRA",FatalException,message);
      }
      if(row.event!=static_cast<G4int>(fGeometryProbeRows.size())) {
        G4ExceptionDescription message;
        message << "Geometry-probe event IDs must be contiguous from zero; "
                << "line " << lineNumber << " has event=" << row.event
                << " expected=" << fGeometryProbeRows.size() << ".";
        G4Exception("PrimaryGenerator::LoadGeometryProbeTable",
                    "HODO_GATE2_PROBE_EVENT",FatalException,message);
      }
      const auto expectedCohort=[](G4int event) {
        if(event<128)return 0;
        if(event<256)return 1;
        if(event<384)return 2;
        if(event<506)return 3;
        if(event<518)return 4;
        if(event<523)return 5;
        if(event<529)return 6;
        if(event<DetectorConstruction::kGeometryProbeExpectedRows)return 7;
        return -1;
      }(row.event);
      if(row.cohort!=expectedCohort) {
        G4ExceptionDescription message;
        message << "Geometry-probe cohort/order mismatch at line "
                << lineNumber << ": event=" << row.event
                << " cohort=" << row.cohort
                << " expected_cohort=" << expectedCohort << ".";
        G4Exception("PrimaryGenerator::LoadGeometryProbeTable",
                    "HODO_GATE2_PROBE_COHORT",FatalException,message);
      }
      const auto planeInRange=row.targetPlane>=0&&row.targetPlane<4;
      const auto channelInRange=row.targetChannel>=0&&row.targetChannel<32;
      G4bool metadataValid=false;
      switch(row.cohort) {
        case 0:
        case 1:
          metadataValid=planeInRange&&channelInRange&&row.expectedClass==1;
          break;
        case 2:
          metadataValid=planeInRange&&channelInRange&&row.expectedClass==2;
          break;
        case 3: {
          const G4bool leftChannel=planeInRange&&row.targetChannel>=0&&
              row.targetChannel<=30&&
              (row.targetPlane>=2||row.targetChannel!=15);
          metadataValid=leftChannel&&row.expectedClass==0;
          break;
        }
        case 4:
          metadataValid=planeInRange&&row.targetChannel==-1&&
              row.expectedClass==0;
          break;
        case 5:
        case 6:
          metadataValid=row.targetPlane==-1&&row.targetChannel==-1&&
              row.expectedClass==-1;
          break;
        case 7:
          metadataValid=planeInRange&&row.targetChannel==16&&
              row.expectedClass==1;
          break;
        default:
          metadataValid=false;
      }
      if(!metadataValid) {
        G4ExceptionDescription message;
        message << "Geometry-probe metadata is outside the frozen contract at "
                << "line " << lineNumber << ": cohort=" << row.cohort
                << " plane=" << row.targetPlane
                << " channel=" << row.targetChannel
                << " class=" << row.expectedClass << ".";
        G4Exception("PrimaryGenerator::LoadGeometryProbeTable",
                    "HODO_GATE2_PROBE_METADATA",FatalException,message);
      }
      if(row.cohort>=0&&row.cohort<static_cast<G4int>(cohortRows.size()))
        ++cohortRows[row.cohort];
      if(!std::isfinite(x)||!std::isfinite(y)||!std::isfinite(z)||
         !std::isfinite(dx)||!std::isfinite(dy)||!std::isfinite(dz)) {
        G4ExceptionDescription message;
        message << "Non-finite geometry-probe coordinate/direction at line "
                << lineNumber << ".";
        G4Exception("PrimaryGenerator::LoadGeometryProbeTable",
                    "HODO_GATE2_PROBE_FINITE",FatalException,message);
      }
      row.position={x*mm,y*mm,z*mm};
      row.direction={dx,dy,dz};
      if(!std::isfinite(row.direction.mag2())||row.direction.mag2()<=0.) {
        G4ExceptionDescription message;
        message << "Invalid geometry-probe direction at line " << lineNumber;
        G4Exception("PrimaryGenerator::LoadGeometryProbeTable",
                    "HODO_GATE2_PROBE_DIRECTION",FatalException,message);
      }
      row.direction=row.direction.unit();
      fGeometryProbeRows.push_back(row);
    }
    constexpr std::array<G4int,8> expectedCohortRows{{
        128,128,128,122,12,5,6,4}};
    if(fGeometryProbeRows.size()!=
           static_cast<std::size_t>(
               DetectorConstruction::kGeometryProbeExpectedRows)||
       cohortRows!=expectedCohortRows||
       fDetector->GetGeometryProbeRowOffset()>=
           static_cast<G4int>(fGeometryProbeRows.size())) {
      G4ExceptionDescription message;
      message << "Geometry-probe table violates the frozen 533-row contract: "
              << fDetector->GetGeometryProbeFile()
              << " rows=" << fGeometryProbeRows.size()
              << " offset=" << fDetector->GetGeometryProbeRowOffset()
              << " cohort_counts=";
      for(std::size_t cohort=0;cohort<cohortRows.size();++cohort)
        message << (cohort?",":"") << cohortRows[cohort];
      G4Exception("PrimaryGenerator::LoadGeometryProbeTable",
                  "HODO_GATE2_PROBE_ROW_CONTRACT",FatalException,message);
    }
    G4cout << "[gate2 geometry probe] loaded rows="
           << fGeometryProbeRows.size()
           << " offset=" << fDetector->GetGeometryProbeRowOffset()
           << " file=" << fDetector->GetGeometryProbeFile() << G4endl;
  }
  void LoadGenerationProbeTable() {
    std::ifstream input(fDetector->GetGenerationProbeFile());
    if(!input) {
      G4ExceptionDescription message;
      message << "Cannot open Stage-3 generation-probe TSV: "
              << fDetector->GetGenerationProbeFile();
      G4Exception("PrimaryGenerator::LoadGenerationProbeTable",
                  "HODO_GATE3_PROBE_OPEN",FatalException,message);
    }
    std::string line;
    G4int lineNumber=0;
    std::array<G4int,5> cohortRows{};
    while(std::getline(input,line)) {
      ++lineNumber;
      const auto first=line.find_first_not_of(" \t\r");
      if(first==std::string::npos||line[first]=='#')continue;
      GenerationProbeRow row;
      G4double x=0.,y=0.,z=0.,dx=0.,dy=0.,dz=0.,energyGeV=0.;
      G4double analyticCoreChordMm=0.,analyticTargetPathMm=0.;
      std::istringstream parser(line);
      if(!(parser>>row.event>>row.cohort>>row.targetPlane
                 >>row.targetChannel>>row.expectedVolumeClass
                 >>row.expectedRegion>>row.chordBin>>row.repeatIndex
                 >>x>>y>>z>>dx>>dy>>dz>>energyGeV
                 >>analyticCoreChordMm>>analyticTargetPathMm
                 >>row.caseIndex)) {
        G4ExceptionDescription message;
        message << "Malformed Stage-3 18-column TSV row at line "
                << lineNumber << ": " << line;
        G4Exception("PrimaryGenerator::LoadGenerationProbeTable",
                    "HODO_GATE3_PROBE_PARSE",FatalException,message);
      }
      std::string extra;
      if(parser>>extra) {
        G4ExceptionDescription message;
        message << "Extra Stage-3 TSV field at line " << lineNumber
                << ": " << extra;
        G4Exception("PrimaryGenerator::LoadGenerationProbeTable",
                    "HODO_GATE3_PROBE_EXTRA",FatalException,message);
      }
      const auto expectedEvent=static_cast<G4int>(fGenerationProbeRows.size());
      if(row.event!=expectedEvent) {
        G4ExceptionDescription message;
        message << "Stage-3 event IDs must be contiguous from zero; line "
                << lineNumber << " has event=" << row.event
                << " expected=" << expectedEvent << ".";
        G4Exception("PrimaryGenerator::LoadGenerationProbeTable",
                    "HODO_GATE3_PROBE_EVENT",FatalException,message);
      }

      const G4int expectedCohort=row.event<196608?0:
          row.event<196864?1:row.event<196880?2:
          row.event<196896?3:
          row.event<DetectorConstruction::kGenerationProbeExpectedRows?4:-1;
      G4bool metadataValid=false;
      if(expectedCohort==0) {
        const G4int within=row.event;
        metadataValid=row.cohort==0&&row.targetPlane==within/49152&&
            row.targetChannel==16&&row.expectedVolumeClass==1&&
            row.expectedRegion==1&&row.chordBin==(within%49152)/4096&&
            row.repeatIndex==within%4096&&row.caseIndex==row.chordBin;
      } else if(expectedCohort==1) {
        const G4int within=row.event-196608;
        metadataValid=row.cohort==1&&row.targetPlane==within/64&&
            row.targetChannel==16&&row.expectedVolumeClass==2&&
            row.expectedRegion==1&&row.chordBin==-1&&
            row.repeatIndex==within%16&&row.caseIndex==(within%64)/16;
      } else if(expectedCohort==2) {
        const G4int within=row.event-196864;
        metadataValid=row.cohort==2&&row.targetPlane==within/4&&
            row.targetChannel==-1&&row.expectedVolumeClass==0&&
            row.expectedRegion==0&&row.chordBin==-1&&
            row.repeatIndex==-1&&row.caseIndex==within%4;
      } else if(expectedCohort==3) {
        const G4int within=row.event-196880;
        const G4int caseIndex=within%4;
        constexpr std::array<G4int,4> volume{{0,1,1,1}};
        constexpr std::array<G4int,4> region{{0,1,1,2}};
        metadataValid=row.cohort==3&&row.targetPlane==within/4&&
            row.targetChannel==16&&
            row.expectedVolumeClass==volume[caseIndex]&&
            row.expectedRegion==region[caseIndex]&&row.chordBin==-1&&
            row.repeatIndex==-1&&row.caseIndex==caseIndex;
      } else if(expectedCohort==4) {
        const G4int within=row.event-196896;
        metadataValid=row.cohort==4&&row.targetPlane==within/16&&
            row.targetChannel==16&&row.expectedVolumeClass==1&&
            row.expectedRegion==3&&row.chordBin==-1&&
            row.repeatIndex==within%16&&row.caseIndex==0;
      }
      if(!metadataValid) {
        G4ExceptionDescription message;
        message << "Stage-3 metadata/order violates the frozen 196960-row "
                << "contract at line " << lineNumber
                << ": event=" << row.event
                << " cohort=" << row.cohort
                << " plane=" << row.targetPlane
                << " channel=" << row.targetChannel
                << " volume=" << row.expectedVolumeClass
                << " region=" << row.expectedRegion
                << " chord_bin=" << row.chordBin
                << " repeat=" << row.repeatIndex
                << " case=" << row.caseIndex << ".";
        G4Exception("PrimaryGenerator::LoadGenerationProbeTable",
                    "HODO_GATE3_PROBE_METADATA",FatalException,message);
      }
      const std::array<G4double,10> floating{{
          x,y,z,dx,dy,dz,energyGeV,analyticCoreChordMm,
          analyticTargetPathMm,dx*dx+dy*dy+dz*dz}};
      if(!std::all_of(floating.begin(),floating.end(),
                      [](G4double value){return std::isfinite(value);})||
         std::abs(floating.back()-1.)>1.e-12||
         std::abs(energyGeV-4.)>1.e-9||analyticCoreChordMm<0.||
         analyticTargetPathMm<0.||
         analyticTargetPathMm+1.e-12<analyticCoreChordMm) {
        G4ExceptionDescription message;
        message << "Invalid Stage-3 coordinate/direction/energy/path at line "
                << lineNumber << ": |dir|^2=" << floating.back()
                << " energy_GeV=" << energyGeV
                << " core_chord_mm=" << analyticCoreChordMm
                << " target_path_mm=" << analyticTargetPathMm << ".";
        G4Exception("PrimaryGenerator::LoadGenerationProbeTable",
                    "HODO_GATE3_PROBE_NUMERIC",FatalException,message);
      }
      row.position={x*mm,y*mm,z*mm};
      row.direction={dx,dy,dz};
      row.energy=energyGeV*GeV;
      row.analyticCoreChord=analyticCoreChordMm*mm;
      row.analyticTargetPath=analyticTargetPathMm*mm;
      if(row.cohort>=0&&row.cohort<5)++cohortRows[row.cohort];
      fGenerationProbeRows.push_back(row);
    }
    constexpr std::array<G4int,5> expectedCohortRows{{196608,256,16,16,64}};
    if(fGenerationProbeRows.size()!=static_cast<std::size_t>(
           DetectorConstruction::kGenerationProbeExpectedRows)||
       cohortRows!=expectedCohortRows||
       fDetector->GetGenerationProbeRowOffset()<0||
       fDetector->GetGenerationProbeRowOffset()>=
           static_cast<G4int>(fGenerationProbeRows.size())) {
      G4ExceptionDescription message;
      message << "Generation-probe table violates the frozen 196960-row "
              << "contract: " << fDetector->GetGenerationProbeFile()
              << " rows=" << fGenerationProbeRows.size()
              << " offset=" << fDetector->GetGenerationProbeRowOffset()
              << " cohort_counts=" << cohortRows[0] << "," << cohortRows[1]
              << "," << cohortRows[2] << "," << cohortRows[3]
              << "," << cohortRows[4] << ".";
      G4Exception("PrimaryGenerator::LoadGenerationProbeTable",
                  "HODO_GATE3_PROBE_ROW_CONTRACT",FatalException,message);
    }
    gGenerationProbeLoadedRows.store(
        static_cast<G4int>(fGenerationProbeRows.size()),
        std::memory_order_relaxed);
    G4cout << "[gate3 generation probe] loaded rows="
           << fGenerationProbeRows.size()
           << " offset=" << fDetector->GetGenerationProbeRowOffset()
           << " file=" << fDetector->GetGenerationProbeFile() << G4endl;
  }
  G4ParticleGun fGun{1};
  const DetectorConstruction* fDetector;
  JointGaisserCache fJointGaisserCache;
  std::vector<GeometryProbeRow> fGeometryProbeRows;
  std::vector<GenerationProbeRow> fGenerationProbeRows;
};

// Keep all optical tracks in the waiting stack until the complete charged and
// non-optical cascade has been transported. This prevents changes in optical
// survival models from consuming random numbers between delta-electron and
// other secondary-transport steps, which otherwise changes T1/T2 energy
// deposition in an optical-only sensitivity comparison.
class OpticalEventAction;

class OpticalLastStackingAction final : public G4UserStackingAction {
 public:
  explicit OpticalLastStackingAction(OpticalEventAction* event)
      : fEvent(event) {}
  G4ClassificationOfNewTrack ClassifyNewTrack(
      const G4Track* track) override;
 private:
  OpticalEventAction* fEvent;
};

class OpticalRun final : public G4Run {
 public:
  void Count(G4bool t1, G4bool t2, G4bool h1, G4bool h2,
             G4bool t3, G4bool t4,
             G4bool roundX, G4bool roundY,
             G4bool squareX, G4bool squareY,
             G4bool parkRoundX, G4bool parkRoundY,
             G4bool parkSquareX, G4bool parkSquareY,
             G4bool parkRoundTimed, G4bool parkSquareTimed,
             G4bool parkRoundNoise, G4bool parkSquareNoise,
             G4long produced, G4long detected,
             G4long producedH1, G4long producedH2,
             G4long detectedH1, G4long detectedH2,
             G4double edepH1, G4double edepH2,
             const std::array<G4long,4>& stageH1,
             const std::array<G4long,4>& stageH2,
             G4long rejectedH1,G4long rejectedH2,
             const std::array<G4double,4>& primaryCoreTrack,
             const std::array<G4double,4>& primaryActiveCoreTrack) {
    fProduced += produced;
    fDetected += detected;
    fProducedH1 += producedH1;
    fProducedH2 += producedH2;
    fDetectedH1 += detectedH1;
    fDetectedH2 += detectedH2;
    fEdepH1 += edepH1;
    fEdepH2 += edepH2;
    for(G4int i=0;i<4;++i){fStageH1[i]+=stageH1[i];fStageH2[i]+=stageH2[i];}
    fRejectedH1+=rejectedH1;fRejectedH2+=rejectedH2;
    if (t1 && t2) {
      ++fTriggered;
      if (t3) ++fTriggeredT123;
      if (t3 && t4) ++fTriggeredT1234;
      if (h1) ++fH1;
      if (h2) ++fH2;
      if (h1 && h2) ++fBoth;
      if (roundX) ++fRoundX;
      if (roundY) ++fRoundY;
      if (squareX) ++fSquareX;
      if (squareY) ++fSquareY;
      if (roundX && roundY) ++fRoundPair;
      if (squareX && squareY) ++fSquarePair;
      if (roundX && roundY && squareX && squareY) ++fAllPlanes;
      if (parkRoundX) ++fParkRoundX;
      if (parkRoundY) ++fParkRoundY;
      if (parkSquareX) ++fParkSquareX;
      if (parkSquareY) ++fParkSquareY;
      if (parkRoundX && parkRoundY) ++fParkRoundPair;
      if (parkSquareX && parkSquareY) ++fParkSquarePair;
      if (parkRoundTimed) ++fParkRoundTimed;
      if (parkSquareTimed) ++fParkSquareTimed;
      if (parkRoundNoise) ++fParkRoundNoise;
      if (parkSquareNoise) ++fParkSquareNoise;
      CountStructural(primaryCoreTrack,true);
      CountStructural(primaryActiveCoreTrack,false);
    }
    if (t1 && t2 && t3) {
      CountStructuralT123(primaryCoreTrack,true);
      CountStructuralT123(primaryActiveCoreTrack,false);
    }
    if (t1 && t2 && t3 && t4) {
      CountStructuralT1234(primaryCoreTrack,true);
      CountStructuralT1234(primaryActiveCoreTrack,false);
    }
  }
  void CountVirtualRound(G4bool t1, G4bool t2, G4bool t3,
                         G4bool h1r, G4bool h2r, G4bool h3r) {
    if (t1 && t2) {
      if (h1r) ++fVirtualT12H1R;
      if (h2r) ++fVirtualT12H2R;
      if (h3r) ++fVirtualT12H3R;
      if (h1r && h2r) ++fVirtualT12H1H2;
      if (h1r && h2r && h3r) ++fVirtualT12H1H2H3;
    }
    if (t1 && t2 && t3) {
      if (h1r) ++fVirtualT123H1R;
      if (h2r) ++fVirtualT123H2R;
      if (h3r) ++fVirtualT123H3R;
      if (h1r && h2r) ++fVirtualT123H1H2;
      if (h1r && h2r && h3r) ++fVirtualT123H1H2H3;
    }
  }
  void CountStructural(const std::array<G4double,4>& track,
                       G4bool totalCore) {
    constexpr G4double eps=1.e-6*mm;
    const G4bool h1=track[0]>eps||track[1]>eps;
    const G4bool h2=track[2]>eps||track[3]>eps;
    auto& counts=totalCore ? fStructuralTotalCore : fStructuralActiveCore;
    if(!h1&&!h2)++counts[0];
    else if(h1&&!h2)++counts[1];
    else if(!h1&&h2)++counts[2];
    else ++counts[3];
    auto& planeCounts=totalCore ? fStructuralPlaneTotalCore
                                : fStructuralPlaneActiveCore;
    for(G4int plane=0;plane<4;++plane)
      if(track[plane]>eps)++planeCounts[plane];
  }
  void CountStructuralT123(const std::array<G4double,4>& track,
                           G4bool totalCore) {
    CountStructuralInto(track,totalCore,
        fStructuralT123TotalCore,fStructuralT123ActiveCore,
        fStructuralT123PlaneTotalCore,fStructuralT123PlaneActiveCore);
  }
  void CountStructuralT1234(const std::array<G4double,4>& track,
                            G4bool totalCore) {
    CountStructuralInto(track,totalCore,
        fStructuralT1234TotalCore,fStructuralT1234ActiveCore,
        fStructuralT1234PlaneTotalCore,fStructuralT1234PlaneActiveCore);
  }
  static void CountStructuralInto(
      const std::array<G4double,4>& track,
      G4bool totalCore,
      std::array<G4long,4>& totalCounts,
      std::array<G4long,4>& activeCounts,
      std::array<G4long,4>& totalPlaneCounts,
      std::array<G4long,4>& activePlaneCounts) {
    constexpr G4double eps=1.e-6*mm;
    const G4bool h1=track[0]>eps||track[1]>eps;
    const G4bool h2=track[2]>eps||track[3]>eps;
    auto& counts=totalCore ? totalCounts : activeCounts;
    if(!h1&&!h2)++counts[0];
    else if(h1&&!h2)++counts[1];
    else if(!h1&&h2)++counts[2];
    else ++counts[3];
    auto& planeCounts=totalCore ? totalPlaneCounts : activePlaneCounts;
    for(G4int plane=0;plane<4;++plane)
      if(track[plane]>eps)++planeCounts[plane];
  }
  void Merge(const G4Run* other) override {
    const auto* rhs = static_cast<const OpticalRun*>(other);
    fTriggered += rhs->fTriggered;
    fTriggeredT123 += rhs->fTriggeredT123;
    fTriggeredT1234 += rhs->fTriggeredT1234;
    fH1 += rhs->fH1;
    fH2 += rhs->fH2;
    fBoth += rhs->fBoth;
    fRoundX += rhs->fRoundX;
    fRoundY += rhs->fRoundY;
    fSquareX += rhs->fSquareX;
    fSquareY += rhs->fSquareY;
    fRoundPair += rhs->fRoundPair;
    fSquarePair += rhs->fSquarePair;
    fAllPlanes += rhs->fAllPlanes;
    fParkRoundX += rhs->fParkRoundX;
    fParkRoundY += rhs->fParkRoundY;
    fParkSquareX += rhs->fParkSquareX;
    fParkSquareY += rhs->fParkSquareY;
    fParkRoundPair += rhs->fParkRoundPair;
    fParkSquarePair += rhs->fParkSquarePair;
    fParkRoundTimed += rhs->fParkRoundTimed;
    fParkSquareTimed += rhs->fParkSquareTimed;
    fParkRoundNoise += rhs->fParkRoundNoise;
    fParkSquareNoise += rhs->fParkSquareNoise;
    fProduced += rhs->fProduced;
    fDetected += rhs->fDetected;
    fProducedH1 += rhs->fProducedH1;
    fProducedH2 += rhs->fProducedH2;
    fDetectedH1 += rhs->fDetectedH1;
    fDetectedH2 += rhs->fDetectedH2;
    fEdepH1 += rhs->fEdepH1;
    fEdepH2 += rhs->fEdepH2;
    for(G4int i=0;i<4;++i){fStageH1[i]+=rhs->fStageH1[i];fStageH2[i]+=rhs->fStageH2[i];}
    for(G4int i=0;i<4;++i) {
      fStructuralTotalCore[i]+=rhs->fStructuralTotalCore[i];
      fStructuralActiveCore[i]+=rhs->fStructuralActiveCore[i];
      fStructuralPlaneTotalCore[i]+=rhs->fStructuralPlaneTotalCore[i];
      fStructuralPlaneActiveCore[i]+=rhs->fStructuralPlaneActiveCore[i];
      fStructuralT123TotalCore[i]+=rhs->fStructuralT123TotalCore[i];
      fStructuralT123ActiveCore[i]+=rhs->fStructuralT123ActiveCore[i];
      fStructuralT123PlaneTotalCore[i]+=rhs->fStructuralT123PlaneTotalCore[i];
      fStructuralT123PlaneActiveCore[i]+=rhs->fStructuralT123PlaneActiveCore[i];
      fStructuralT1234TotalCore[i]+=rhs->fStructuralT1234TotalCore[i];
      fStructuralT1234ActiveCore[i]+=rhs->fStructuralT1234ActiveCore[i];
      fStructuralT1234PlaneTotalCore[i]+=rhs->fStructuralT1234PlaneTotalCore[i];
      fStructuralT1234PlaneActiveCore[i]+=rhs->fStructuralT1234PlaneActiveCore[i];
    }
    fRejectedH1+=rhs->fRejectedH1;fRejectedH2+=rhs->fRejectedH2;
    fVirtualT12H1R += rhs->fVirtualT12H1R;
    fVirtualT12H2R += rhs->fVirtualT12H2R;
    fVirtualT12H3R += rhs->fVirtualT12H3R;
    fVirtualT12H1H2 += rhs->fVirtualT12H1H2;
    fVirtualT12H1H2H3 += rhs->fVirtualT12H1H2H3;
    fVirtualT123H1R += rhs->fVirtualT123H1R;
    fVirtualT123H2R += rhs->fVirtualT123H2R;
    fVirtualT123H3R += rhs->fVirtualT123H3R;
    fVirtualT123H1H2 += rhs->fVirtualT123H1H2;
    fVirtualT123H1H2H3 += rhs->fVirtualT123H1H2H3;
    G4Run::Merge(other);
  }
  G4long Triggered() const { return fTriggered; }
  G4long TriggeredT123() const { return fTriggeredT123; }
  G4long TriggeredT1234() const { return fTriggeredT1234; }
  G4long H1() const { return fH1; }
  G4long H2() const { return fH2; }
  G4long Both() const { return fBoth; }
  G4long RoundX() const { return fRoundX; }
  G4long RoundY() const { return fRoundY; }
  G4long SquareX() const { return fSquareX; }
  G4long SquareY() const { return fSquareY; }
  G4long RoundPair() const { return fRoundPair; }
  G4long SquarePair() const { return fSquarePair; }
  G4long AllPlanes() const { return fAllPlanes; }
  G4long ParkRoundX() const { return fParkRoundX; }
  G4long ParkRoundY() const { return fParkRoundY; }
  G4long ParkSquareX() const { return fParkSquareX; }
  G4long ParkSquareY() const { return fParkSquareY; }
  G4long ParkRoundPair() const { return fParkRoundPair; }
  G4long ParkSquarePair() const { return fParkSquarePair; }
  G4long ParkRoundTimed() const { return fParkRoundTimed; }
  G4long ParkSquareTimed() const { return fParkSquareTimed; }
  G4long ParkRoundNoise() const { return fParkRoundNoise; }
  G4long ParkSquareNoise() const { return fParkSquareNoise; }
  G4long Produced() const { return fProduced; }
  G4long Detected() const { return fDetected; }
  G4long ProducedH1() const { return fProducedH1; }
  G4long ProducedH2() const { return fProducedH2; }
  G4long DetectedH1() const { return fDetectedH1; }
  G4long DetectedH2() const { return fDetectedH2; }
  G4double EdepH1() const { return fEdepH1; }
  G4double EdepH2() const { return fEdepH2; }
  G4long StageH1(G4int i) const { return fStageH1[i]; }
  G4long StageH2(G4int i) const { return fStageH2[i]; }
  G4long RejectedH1() const { return fRejectedH1; }
  G4long RejectedH2() const { return fRejectedH2; }
  G4long StructuralTotalCore(G4int i) const {
    return i>=0&&i<4 ? fStructuralTotalCore[i] : 0;
  }
  G4long StructuralActiveCore(G4int i) const {
    return i>=0&&i<4 ? fStructuralActiveCore[i] : 0;
  }
  G4long StructuralPlaneTotalCore(G4int i) const {
    return i>=0&&i<4 ? fStructuralPlaneTotalCore[i] : 0;
  }
  G4long StructuralPlaneActiveCore(G4int i) const {
    return i>=0&&i<4 ? fStructuralPlaneActiveCore[i] : 0;
  }
  G4long StructuralT123ActiveCore(G4int i) const {
    return i>=0&&i<4 ? fStructuralT123ActiveCore[i] : 0;
  }
  G4long StructuralT1234ActiveCore(G4int i) const {
    return i>=0&&i<4 ? fStructuralT1234ActiveCore[i] : 0;
  }
  G4long StructuralT123PlaneActiveCore(G4int i) const {
    return i>=0&&i<4 ? fStructuralT123PlaneActiveCore[i] : 0;
  }
  G4long StructuralT1234PlaneActiveCore(G4int i) const {
    return i>=0&&i<4 ? fStructuralT1234PlaneActiveCore[i] : 0;
  }
  G4long VirtualT12H1R() const { return fVirtualT12H1R; }
  G4long VirtualT12H2R() const { return fVirtualT12H2R; }
  G4long VirtualT12H3R() const { return fVirtualT12H3R; }
  G4long VirtualT12H1H2() const { return fVirtualT12H1H2; }
  G4long VirtualT12H1H2H3() const { return fVirtualT12H1H2H3; }
  G4long VirtualT123H1R() const { return fVirtualT123H1R; }
  G4long VirtualT123H2R() const { return fVirtualT123H2R; }
  G4long VirtualT123H3R() const { return fVirtualT123H3R; }
  G4long VirtualT123H1H2() const { return fVirtualT123H1H2; }
  G4long VirtualT123H1H2H3() const { return fVirtualT123H1H2H3; }
 private:
  G4long fTriggered=0, fTriggeredT123=0, fTriggeredT1234=0;
  G4long fH1=0, fH2=0, fBoth=0;
  G4long fRoundX=0, fRoundY=0, fSquareX=0, fSquareY=0;
  G4long fRoundPair=0, fSquarePair=0, fAllPlanes=0;
  G4long fParkRoundX=0, fParkRoundY=0, fParkSquareX=0, fParkSquareY=0;
  G4long fParkRoundPair=0, fParkSquarePair=0;
  G4long fParkRoundTimed=0, fParkSquareTimed=0;
  G4long fParkRoundNoise=0, fParkSquareNoise=0;
  G4long fProduced=0, fDetected=0;
  G4long fProducedH1=0, fProducedH2=0, fDetectedH1=0, fDetectedH2=0;
  G4double fEdepH1=0., fEdepH2=0.;
  std::array<G4long,4> fStageH1{},fStageH2{};
  G4long fRejectedH1=0,fRejectedH2=0;
  // Structural path-only categories for T1&&T2 events.
  // Order: T-only, H1-only, H2-only, H1&H2.
  std::array<G4long,4> fStructuralTotalCore{};
  std::array<G4long,4> fStructuralActiveCore{};
  // Plane order: MID1 round X, MID9 round Y, MID15 square X, MID13 square Y.
  std::array<G4long,4> fStructuralPlaneTotalCore{};
  std::array<G4long,4> fStructuralPlaneActiveCore{};
  std::array<G4long,4> fStructuralT123TotalCore{};
  std::array<G4long,4> fStructuralT123ActiveCore{};
  std::array<G4long,4> fStructuralT123PlaneTotalCore{};
  std::array<G4long,4> fStructuralT123PlaneActiveCore{};
  std::array<G4long,4> fStructuralT1234TotalCore{};
  std::array<G4long,4> fStructuralT1234ActiveCore{};
  std::array<G4long,4> fStructuralT1234PlaneTotalCore{};
  std::array<G4long,4> fStructuralT1234PlaneActiveCore{};
  G4long fVirtualT12H1R=0, fVirtualT12H2R=0, fVirtualT12H3R=0;
  G4long fVirtualT12H1H2=0, fVirtualT12H1H2H3=0;
  G4long fVirtualT123H1R=0, fVirtualT123H2R=0, fVirtualT123H3R=0;
  G4long fVirtualT123H1H2=0, fVirtualT123H1H2H3=0;
};

class OpticalRunAction final : public G4UserRunAction {
 public:
  explicit OpticalRunAction(const DetectorConstruction* detector) : fDetector(detector) {}
  G4Run* GenerateRun() override { return new OpticalRun; }
  void BeginOfRunAction(const G4Run*) override {
    const G4bool generationMode=fDetector->GetGenerationProbeMode();
    const G4bool recordGeneration=fDetector->GetRecordGenerationNtuples();
    if(generationMode!=recordGeneration)
      G4Exception("OpticalRunAction::BeginOfRunAction",
                  "HODO_GATE3_MODE_RECORD_CONTRACT",FatalException,
                  "generationProbeMode and recordGenerationNtuples must be "
                  "enabled or disabled together.");
    if(generationMode) {
      const auto* runManager=G4RunManager::GetRunManager();
      if(!runManager||
         runManager->GetRunManagerType()!=G4RunManager::sequentialRM)
        G4Exception("OpticalRunAction::BeginOfRunAction",
                    "HODO_GATE3_SERIAL_ONLY",FatalException,
                    "Stage-3 generationProbeMode requires the SerialOnly "
                    "kernel.");
      if(fDetector->GetRecordSpatialClosureNtuples()||
         fDetector->GetRecordGeometryControlNtuples()||
         fDetector->GetGeometryProbeMode())
        G4Exception("OpticalRunAction::BeginOfRunAction",
                    "HODO_GATE3_DIAGNOSTIC_CONFLICT",FatalException,
                    "Stage-3 requires Gate-1 spatial and Gate-2 geometry "
                    "diagnostic trees to remain disabled so IDs 7..10 are "
                    "stable.");
      const auto previous=
          gGenerationProbeRunBegins.fetch_add(1,std::memory_order_relaxed);
      if(previous!=0)
        G4Exception("OpticalRunAction::BeginOfRunAction",
                    "HODO_GATE3_ONE_RUN_PER_PROCESS",FatalException,
                    "Stage-3 generationProbeMode permits exactly one run per "
                    "process; start a new process for every lane.");
      const auto eventCount=runManager->GetNumberOfEventsToBeProcessed();
      if(eventCount<=0||fDetector->GetGenerationProbeRowOffset()+eventCount>
             DetectorConstruction::kGenerationProbeExpectedRows) {
        G4ExceptionDescription message;
        message << "Stage-3 beamOn range is outside the frozen source table: "
                << "offset=" << fDetector->GetGenerationProbeRowOffset()
                << " events=" << eventCount
                << " rows="
                << DetectorConstruction::kGenerationProbeExpectedRows << ".";
        G4Exception("OpticalRunAction::BeginOfRunAction",
                    "HODO_GATE3_BEAMON_RANGE",FatalException,message);
      }
      gGenerationObservedMaxCoreStep.store(0.,std::memory_order_relaxed);
      // G4RunMessenger passes /random/setSeeds through a stack-local array.
      // CLHEP MixMaxRng retains that pointer in HepRandomEngine::theSeeds, so
      // G4Random::getTheSeeds() is a dangling-pointer readback after the UI
      // command returns and must never be dereferenced here.  getTheSeed()
      // reads MixMaxRng's owned first seed and is the only safe initialized
      // seed readback exposed by this engine.  The second column is therefore
      // an explicit zero sentinel meaning "unavailable", not an invented
      // copy of the second macro literal.  The frozen macro contract checks
      // the complete /random/setSeeds pair separately.
      fInitialRandomSeed0=G4Random::getTheSeed();
      fInitialRandomSeed1=0;
    }
    if(fDetector->GetGeometryProbeMode()) {
      const auto* runManager=G4RunManager::GetRunManager();
      if(!runManager||
         runManager->GetRunManagerType()!=G4RunManager::sequentialRM)
        G4Exception("OpticalRunAction::BeginOfRunAction",
                    "HODO_GATE2_SERIAL_ONLY",FatalException,
                    "Stage-2 geometryProbeMode requires the SerialOnly kernel.");
      const auto previous=
          gGeometryProbeRunBegins.fetch_add(1,std::memory_order_relaxed);
      if(previous!=0)
        G4Exception("OpticalRunAction::BeginOfRunAction",
                    "HODO_GATE2_ONE_RUN_PER_PROCESS",FatalException,
                    "Stage-2 geometryProbeMode permits exactly one run per "
                    "process; start a new process for every control.");
    }
    // G4AnalysisManager is thread-local, but ROOT ntuple-merging setup touches
    // shared booking/merge state. Serialise the unusually wide (227-column)
    // schema creation and file opening; concurrent worker booking caused an
    // intermittent crash at the first basket flush near event 1000.
    G4AutoLock bookingLock(&gAnalysisBookingMutex);
    if (IsMaster()) {
      gAcceptedT1T2.store(0, std::memory_order_relaxed);
      gAcceptedTriggerCoincidence.store(0, std::memory_order_relaxed);
      gLastProgressBucket.store(0, std::memory_order_relaxed);
    }
    auto* analysis = G4AnalysisManager::Instance();
    analysis->SetVerboseLevel(0);
    const auto* bookingRunManager=G4RunManager::GetRunManager();
    if(!bookingRunManager ||
       bookingRunManager->GetRunManagerType()!=G4RunManager::sequentialRM)
      analysis->SetNtupleMerging(true);
    analysis->CreateNtuple("events", "optical-v21 event diagnostics");
    analysis->CreateNtupleIColumn("event");
    analysis->CreateNtupleDColumn("h1_attenuation_m");
    analysis->CreateNtupleIColumn("h1_single_clad");
    analysis->CreateNtupleIColumn("h1_emission_22point");
    analysis->CreateNtupleIColumn("produced");
    analysis->CreateNtupleIColumn("detected");
    analysis->CreateNtupleIColumn("pde_rejected");
    analysis->CreateNtupleIColumn("escaped");
    analysis->CreateNtupleIColumn("bulk_abs_core");
    analysis->CreateNtupleIColumn("bulk_abs_clad");
    analysis->CreateNtupleIColumn("bulk_abs_coupling");
    analysis->CreateNtupleIColumn("bulk_abs_other");
    analysis->CreateNtupleIColumn("microcell_count");
    analysis->CreateNtupleDColumn("dark_count_rate_khz");
    analysis->CreateNtupleDColumn("pde_scale");
    analysis->CreateNtupleDColumn("coupling_efficiency");
    analysis->CreateNtupleDColumn("bend_transmission");
    analysis->CreateNtupleDColumn("scintillation_yield_scale");
    analysis->CreateNtupleIColumn("t1");
    analysis->CreateNtupleIColumn("t2");
    analysis->CreateNtupleIColumn("h1");
    analysis->CreateNtupleIColumn("h2");
    analysis->CreateNtupleDColumn("t1_edep_mev");
    analysis->CreateNtupleDColumn("t2_edep_mev");
    analysis->CreateNtupleIColumn("accepted_t1t2_cumulative");
    analysis->CreateNtupleIColumn("target_coincidences");
    analysis->CreateNtupleDColumn("primary_x_mm");
    analysis->CreateNtupleDColumn("primary_y_mm");
    analysis->CreateNtupleDColumn("primary_z_mm");
    analysis->CreateNtupleDColumn("primary_dir_x");
    analysis->CreateNtupleDColumn("primary_dir_y");
    analysis->CreateNtupleDColumn("primary_dir_z");
    analysis->CreateNtupleDColumn("primary_zenith_deg");
    analysis->CreateNtupleDColumn("h1_adc_max");
    analysis->CreateNtupleDColumn("h2_adc_max");
    analysis->CreateNtupleIColumn("loss_stress");
    analysis->CreateNtupleIColumn("loss_crack");
    analysis->CreateNtupleIColumn("loss_packing");
    analysis->CreateNtupleIColumn("loss_adhesive");
    analysis->CreateNtupleDColumn("primary_energy_gev");
    analysis->CreateNtupleDColumn("energy_at_t1_gev");
    analysis->CreateNtupleDColumn("edep_concrete_mev");
    analysis->CreateNtupleDColumn("edep_blackout_mev");
    analysis->CreateNtupleIColumn("loss_h1_cladding_damage");
    analysis->CreateNtupleIColumn("loss_h2_cladding_damage");
    analysis->CreateNtupleIColumn("measured_noise_enabled");
    analysis->CreateNtupleDColumn("adc_noise_scale");
    analysis->CreateNtupleDColumn("h1_adc_noiseless_max");
    analysis->CreateNtupleDColumn("h2_adc_noiseless_max");
    analysis->CreateNtupleDColumn("h1_adc_noise_only_max");
    analysis->CreateNtupleDColumn("h2_adc_noise_only_max");
    analysis->CreateNtupleIColumn("h1_noise_false_positive");
    analysis->CreateNtupleIColumn("h2_noise_false_positive");
    analysis->CreateNtupleIColumn("h1_noise_lost");
    analysis->CreateNtupleIColumn("h2_noise_lost");
    analysis->CreateNtupleIColumn("waveform_digitization_enabled");
    analysis->CreateNtupleDColumn("adc_sample_period_ns");
    analysis->CreateNtupleDColumn("single_pe_fwhm_samples");
    analysis->CreateNtupleDColumn("h1_adc_prompt_noiseless_max");
    analysis->CreateNtupleDColumn("h2_adc_prompt_noiseless_max");
    analysis->CreateNtupleDColumn("h1_adc_waveform_noiseless_max");
    analysis->CreateNtupleDColumn("h2_adc_waveform_noiseless_max");
    analysis->CreateNtupleIColumn("h2_single_clad");
    analysis->CreateNtupleDColumn("h1_total_edep_mev");
    analysis->CreateNtupleDColumn("h2_total_edep_mev");
    analysis->CreateNtupleDColumn("h1_core_edep_mev");
    analysis->CreateNtupleDColumn("h2_core_edep_mev");
    analysis->CreateNtupleIColumn("h1_photons_produced");
    analysis->CreateNtupleIColumn("h2_photons_produced");
    analysis->CreateNtupleIColumn("h1_photons_detected");
    analysis->CreateNtupleIColumn("h2_photons_detected");
    analysis->CreateNtupleDColumn("primary_projected_x_z0_mm");
    analysis->CreateNtupleDColumn("primary_projected_y_z0_mm");
    analysis->CreateNtupleDColumn("trigger1_time_ns");
    analysis->CreateNtupleDColumn("adc_peak_window_delay_ns");
    analysis->CreateNtupleDColumn("h1_first_arrival_rel_t1_ns");
    analysis->CreateNtupleDColumn("h2_first_arrival_rel_t1_ns");
    analysis->CreateNtupleDColumn("h1_predicted_peak_sample");
    analysis->CreateNtupleDColumn("h2_predicted_peak_sample");
    analysis->CreateNtupleDColumn("round_x_adc_max");
    analysis->CreateNtupleDColumn("round_y_adc_max");
    analysis->CreateNtupleDColumn("square_x_adc_max");
    analysis->CreateNtupleDColumn("square_y_adc_max");
    analysis->CreateNtupleIColumn("round_x_hit");
    analysis->CreateNtupleIColumn("round_y_hit");
    analysis->CreateNtupleIColumn("square_x_hit");
    analysis->CreateNtupleIColumn("square_y_hit");
    analysis->CreateNtupleIColumn("round_xy_hit");
    analysis->CreateNtupleIColumn("square_xy_hit");
    analysis->CreateNtupleIColumn("all_four_planes_hit");
    analysis->CreateNtupleIColumn("park_round_x_hit");
    analysis->CreateNtupleIColumn("park_round_y_hit");
    analysis->CreateNtupleIColumn("park_square_x_hit");
    analysis->CreateNtupleIColumn("park_square_y_hit");
    analysis->CreateNtupleIColumn("park_round_xy_hit");
    analysis->CreateNtupleIColumn("park_square_xy_hit");
    analysis->CreateNtupleIColumn("park_round_noise_false_hit");
    analysis->CreateNtupleIColumn("park_square_noise_false_hit");
    analysis->CreateNtupleDColumn("park_round_xy_dt_ns");
    analysis->CreateNtupleDColumn("park_square_xy_dt_ns");
    analysis->CreateNtupleIColumn("park_round_xy_timed_hit");
    analysis->CreateNtupleIColumn("park_square_xy_timed_hit");
    analysis->CreateNtupleIColumn("h1_sipm_boundary_photons");
    analysis->CreateNtupleIColumn("h2_sipm_boundary_photons");
    analysis->CreateNtupleIColumn("h1_pde_rejected_photons");
    analysis->CreateNtupleIColumn("h2_pde_rejected_photons");
    analysis->CreateNtupleIColumn("round_x_pde_accepted_pe");
    analysis->CreateNtupleIColumn("round_y_pde_accepted_pe");
    analysis->CreateNtupleIColumn("square_x_pde_accepted_pe");
    analysis->CreateNtupleIColumn("square_y_pde_accepted_pe");
    analysis->CreateNtupleDColumn("round_x_waveform_noiseless_max_adc");
    analysis->CreateNtupleDColumn("round_y_waveform_noiseless_max_adc");
    analysis->CreateNtupleDColumn("square_x_waveform_noiseless_max_adc");
    analysis->CreateNtupleDColumn("square_y_waveform_noiseless_max_adc");
    analysis->CreateNtupleDColumn("mid1_primary_core_track_mm");
    analysis->CreateNtupleDColumn("mid9_primary_core_track_mm");
    analysis->CreateNtupleDColumn("mid15_primary_core_track_mm");
    analysis->CreateNtupleDColumn("mid13_primary_core_track_mm");
    analysis->CreateNtupleIColumn("mid1_primary_unique_cores");
    analysis->CreateNtupleIColumn("mid9_primary_unique_cores");
    analysis->CreateNtupleIColumn("mid15_primary_unique_cores");
    analysis->CreateNtupleIColumn("mid13_primary_unique_cores");
    analysis->CreateNtupleDColumn("mid1_core_edep_mev");
    analysis->CreateNtupleDColumn("mid9_core_edep_mev");
    analysis->CreateNtupleDColumn("mid15_core_edep_mev");
    analysis->CreateNtupleDColumn("mid13_core_edep_mev");
    analysis->CreateNtupleIColumn("mid1_photons_produced");
    analysis->CreateNtupleIColumn("mid9_photons_produced");
    analysis->CreateNtupleIColumn("mid15_photons_produced");
    analysis->CreateNtupleIColumn("mid13_photons_produced");
    analysis->CreateNtupleDColumn("mid1_primary_core_edep_mev");
    analysis->CreateNtupleDColumn("mid9_primary_core_edep_mev");
    analysis->CreateNtupleDColumn("mid15_primary_core_edep_mev");
    analysis->CreateNtupleDColumn("mid13_primary_core_edep_mev");
    analysis->CreateNtupleDColumn("mid1_secondary_core_edep_mev");
    analysis->CreateNtupleDColumn("mid9_secondary_core_edep_mev");
    analysis->CreateNtupleDColumn("mid15_secondary_core_edep_mev");
    analysis->CreateNtupleDColumn("mid13_secondary_core_edep_mev");
    analysis->CreateNtupleIColumn("mid1_primary_photons_produced");
    analysis->CreateNtupleIColumn("mid9_primary_photons_produced");
    analysis->CreateNtupleIColumn("mid15_primary_photons_produced");
    analysis->CreateNtupleIColumn("mid13_primary_photons_produced");
    analysis->CreateNtupleIColumn("mid1_secondary_photons_produced");
    analysis->CreateNtupleIColumn("mid9_secondary_photons_produced");
    analysis->CreateNtupleIColumn("mid15_secondary_photons_produced");
    analysis->CreateNtupleIColumn("mid13_secondary_photons_produced");
    analysis->CreateNtupleIColumn("loss_fanout_exit");
    analysis->CreateNtupleIColumn("loss_end_face");
    analysis->CreateNtupleIColumn("loss_grease_window_coupling");
    analysis->CreateNtupleIColumn("h1_stage_active_exit");
    analysis->CreateNtupleIColumn("h1_stage_fanout_exit");
    analysis->CreateNtupleIColumn("h1_stage_window_entry");
    analysis->CreateNtupleIColumn("h2_stage_active_exit");
    analysis->CreateNtupleIColumn("h2_stage_fanout_exit");
    analysis->CreateNtupleIColumn("h2_stage_window_entry");
    analysis->CreateNtupleIColumn("fixed_photon_test");
    analysis->CreateNtupleIColumn("fixed_photon_channel");
    analysis->CreateNtupleDColumn("fixed_photon_wavelength_nm");
    analysis->CreateNtupleIColumn("fixed_photon_global_channel");
    analysis->CreateNtupleDColumn("fixed_photon_min_radius_mm");
    analysis->CreateNtupleDColumn("fixed_photon_channel_transmission");
    const std::array<G4String,4> planeNames{"h1x","h1y","h2x","h2y"};
    for(const auto& plane:planeNames) {
      analysis->CreateNtupleIColumn(plane+"_park_best_local_channel");
      analysis->CreateNtupleIColumn(plane+"_park_best_pe");
      analysis->CreateNtupleDColumn(plane+"_park_best_noiseless_adc");
      analysis->CreateNtupleDColumn(plane+"_park_best_final_adc");
      analysis->CreateNtupleDColumn(plane+"_park_best_threshold_adc");
      analysis->CreateNtupleDColumn(plane+"_park_best_margin_adc");
      analysis->CreateNtupleDColumn(plane+"_park_best_gain_adc_per_pe");
      analysis->CreateNtupleDColumn(plane+"_park_best_noise_rms_adc");
      analysis->CreateNtupleDColumn(plane+"_park_best_noise_only_adc");
    }
    for(const auto& plane:planeNames) {
      analysis->CreateNtupleIColumn(plane+"_photons_produced_active");
      analysis->CreateNtupleIColumn(plane+"_photons_produced_fanout");
    }
    analysis->CreateNtupleDColumn("h2_attenuation_m");
    analysis->CreateNtupleDColumn("h1_scintillation_yield_per_mev");
    analysis->CreateNtupleDColumn("h2_scintillation_yield_per_mev");
    analysis->CreateNtupleDColumn("h1_birks_mm_per_mev");
    analysis->CreateNtupleDColumn("h2_birks_mm_per_mev");
    analysis->CreateNtupleDColumn("trigger_birks_mm_per_mev");
    const std::array<G4String,4> stageNames{
      "active_exit","fanout_exit","window_entry","sipm_boundary"
    };
    for(const auto& plane:planeNames)
      for(const auto& stage:stageNames)
        analysis->CreateNtupleIColumn(plane+"_stage_"+stage);
    if (fDetector->GetRecordActiveCoreChordBranches()) {
      analysis->CreateNtupleDColumn("mid1_primary_active_core_track_mm");
      analysis->CreateNtupleDColumn("mid9_primary_active_core_track_mm");
      analysis->CreateNtupleDColumn("mid15_primary_active_core_track_mm");
      analysis->CreateNtupleDColumn("mid13_primary_active_core_track_mm");
    }
    analysis->CreateNtupleIColumn("t3");
    analysis->CreateNtupleIColumn("t4");
    analysis->CreateNtupleDColumn("t3_edep_mev");
    analysis->CreateNtupleDColumn("t4_edep_mev");
    analysis->CreateNtupleIColumn("selected_trigger_coincidence");
    analysis->CreateNtupleIColumn("accepted_selected_trigger_cumulative");
    if(fDetector->GetActivePlaneOverlap()||fDetector->GetFibreFarEndExtension()>0.)
      for(const auto& plane:planeNames)
        analysis->CreateNtupleIColumn(plane+"_photons_produced_far_end");
    analysis->FinishNtuple();
    analysis->CreateNtuple("photons", "SiPM boundary photon diagnostics");
    analysis->CreateNtupleIColumn("event");
    analysis->CreateNtupleIColumn("channel");
    analysis->CreateNtupleIColumn("outcome"); // 1=detected, 0=PDE rejected
    analysis->CreateNtupleDColumn("wavelength_nm");
    analysis->CreateNtupleDColumn("arrival_time_ns");
    analysis->FinishNtuple();
    analysis->CreateNtuple("channels",
        "Optional per-event per-channel electronics diagnostics");
    analysis->CreateNtupleIColumn("event");
    analysis->CreateNtupleIColumn("global_channel");
    analysis->CreateNtupleIColumn("plane");
    analysis->CreateNtupleIColumn("local_channel");
    analysis->CreateNtupleDColumn("final_adc");
    analysis->CreateNtupleDColumn("noise_only_adc");
    analysis->CreateNtupleDColumn("noiseless_adc");
    analysis->CreateNtupleDColumn("threshold_adc");
    analysis->CreateNtupleDColumn("gain_adc_per_pe");
    analysis->CreateNtupleDColumn("noise_rms_adc");
    analysis->FinishNtuple();
    analysis->CreateNtuple("photon_audit",
        "V15 fibre-origin photon lineage; diagnostic-only");
    analysis->CreateNtupleIColumn("event");
    analysis->CreateNtupleIColumn("track");
    analysis->CreateNtupleIColumn("parent");
    analysis->CreateNtupleIColumn("plane");
    analysis->CreateNtupleIColumn("origin_global_channel");
    analysis->CreateNtupleIColumn("creator_process");
    analysis->CreateNtupleIColumn("origin_region");
    analysis->CreateNtupleIColumn("terminal_reason");
    analysis->CreateNtupleIColumn("terminal_region");
    analysis->CreateNtupleIColumn("stage_mask");
    analysis->CreateNtupleIColumn("ever_cladding");
    analysis->CreateNtupleIColumn("reentered_core");
    analysis->CreateNtupleDColumn("wavelength_nm");
    analysis->CreateNtupleDColumn("origin_x_mm");
    analysis->CreateNtupleDColumn("origin_y_mm");
    analysis->CreateNtupleDColumn("origin_z_mm");
    analysis->CreateNtupleDColumn("initial_dir_x");
    analysis->CreateNtupleDColumn("initial_dir_y");
    analysis->CreateNtupleDColumn("initial_dir_z");
    analysis->CreateNtupleDColumn("active_core_path_mm");
    analysis->CreateNtupleDColumn("fanout_core_path_mm");
    analysis->CreateNtupleDColumn("active_inner_clad_path_mm");
    analysis->CreateNtupleDColumn("fanout_inner_clad_path_mm");
    analysis->CreateNtupleDColumn("active_outer_clad_path_mm");
    analysis->CreateNtupleDColumn("fanout_outer_clad_path_mm");
    analysis->CreateNtupleDColumn("grease_path_mm");
    analysis->CreateNtupleDColumn("window_path_mm");
    analysis->CreateNtupleDColumn("other_path_mm");
    analysis->CreateNtupleIColumn("core_to_clad");
    analysis->CreateNtupleIColumn("clad_to_core");
    analysis->CreateNtupleIColumn("core_grease_attempts");
    analysis->CreateNtupleIColumn("clad_grease_attempts");
    analysis->CreateNtupleIColumn("grease_window_attempts");
    analysis->CreateNtupleIColumn("sipm_attempts");
    analysis->CreateNtupleIColumn("sipm_reflections");
    analysis->CreateNtupleIColumn("tir_count");
    analysis->CreateNtupleIColumn("fresnel_reflection_count");
    analysis->CreateNtupleIColumn("other_reflection_count");
    analysis->CreateNtupleIColumn("duplicate_terminal_marks");
    analysis->CreateNtupleDColumn("terminal_x_mm");
    analysis->CreateNtupleDColumn("terminal_y_mm");
    analysis->CreateNtupleDColumn("terminal_z_mm");
    analysis->CreateNtupleDColumn("total_track_length_mm");
    analysis->CreateNtupleIColumn("steps");
    analysis->CreateNtupleDColumn("origin_time_ns");
    analysis->CreateNtupleDColumn("terminal_time_ns");
    analysis->FinishNtuple();
    analysis->CreateNtuple("loss_summary",
        "V15 per-event per-plane fibre-photon closure");
    analysis->CreateNtupleIColumn("event");
    analysis->CreateNtupleIColumn("plane");
    analysis->CreateNtupleIColumn("registered");
    analysis->CreateNtupleIColumn("terminal_total");
    analysis->CreateNtupleIColumn("sipm_detected");
    analysis->CreateNtupleIColumn("sipm_pde_rejected");
    analysis->CreateNtupleIColumn("bulk_core");
    analysis->CreateNtupleIColumn("bulk_inner_clad");
    analysis->CreateNtupleIColumn("bulk_outer_clad");
    analysis->CreateNtupleIColumn("bulk_coupling");
    analysis->CreateNtupleIColumn("bulk_other");
    analysis->CreateNtupleIColumn("world_escape");
    analysis->CreateNtupleIColumn("degradation");
    analysis->CreateNtupleIColumn("fanout_rejected");
    analysis->CreateNtupleIColumn("endface_rejected");
    analysis->CreateNtupleIColumn("coupling_rejected");
    analysis->CreateNtupleIColumn("boundary_absorption");
    analysis->CreateNtupleIColumn("other_terminal");
    analysis->CreateNtupleIColumn("unclassified");
    analysis->CreateNtupleIColumn("stage_active_exit");
    analysis->CreateNtupleIColumn("stage_core_to_grease");
    analysis->CreateNtupleIColumn("stage_clad_to_grease");
    analysis->CreateNtupleIColumn("stage_window_entry");
    analysis->CreateNtupleIColumn("stage_sipm_boundary");
    analysis->CreateNtupleIColumn("core_to_clad_tracks");
    analysis->CreateNtupleIColumn("clad_to_core_tracks");
    analysis->CreateNtupleIColumn("duplicate_terminal_marks");
    analysis->CreateNtupleIColumn("terminal_minus_registered");
    analysis->CreateNtupleIColumn("legacy_produced_plane");
    analysis->CreateNtupleIColumn("registered_minus_legacy_produced");
    analysis->FinishNtuple();
    // Keep the canonical 227-column events tree unchanged unless the exact
    // V16 active-chord option is requested. Forward-geometry assumptions live
    // in a separate ntuple and never masquerade as confirmed measurements.
    analysis->CreateNtuple("v21_configuration",
        "V21 geometry provenance including world-scale tolerances");
    analysis->CreateNtupleIColumn("round_tessellated_sides");
    analysis->CreateNtupleIColumn("round_cladding_hierarchy_code");
    for (const auto& plane:planeNames)
      analysis->CreateNtupleIColumn(plane+"_layout_code");
    for (const auto& plane:planeNames)
      analysis->CreateNtupleIColumn(plane+"_channel_order_code");
    for (const auto& plane:planeNames)
      analysis->CreateNtupleDColumn(plane+"_row_z_spacing_mm");
    analysis->CreateNtupleIColumn("mixed_layout_pattern_assumed");
    analysis->CreateNtupleIColumn("active_chord_branches");
    analysis->CreateNtupleIColumn("layout_semantics_version");
    analysis->CreateNtupleIColumn("channel_order_confirmed");
    analysis->CreateNtupleIColumn("row_z_confirmed");
    analysis->CreateNtupleIColumn("fanout_preview_enabled");
    analysis->CreateNtupleDColumn("h1_offset_x_mm");
    analysis->CreateNtupleDColumn("h1_offset_y_mm");
    analysis->CreateNtupleDColumn("h2_offset_x_mm");
    analysis->CreateNtupleDColumn("h2_offset_y_mm");
    analysis->CreateNtupleDColumn("world_maximum_extent_mm");
    analysis->CreateNtupleDColumn("surface_tolerance_mm");
    analysis->CreateNtupleDColumn("radial_tolerance_mm");
    analysis->CreateNtupleIColumn("run_manager_mode_code");
    analysis->CreateNtupleDColumn("round_polygon_phase_rad");
    analysis->CreateNtupleIColumn("round_sweep_representation_code");
    analysis->CreateNtupleDColumn("round_miter_joint_overlap_mm");
    analysis->CreateNtupleIColumn("round_active_facet_refinement");
    analysis->CreateNtupleIColumn("square_cladding_hierarchy_code");
    analysis->CreateNtupleDColumn("fibre_far_end_extension_mm");
    analysis->CreateNtupleIColumn("active_plane_overlap");
    analysis->CreateNtupleIColumn("h2_square_fibre_count");
    analysis->CreateNtupleIColumn("selected_diagnostics_only");
    analysis->FinishNtuple();
    analysis->CreateNtuple("v21_geometry_map",
        "V21 requested per-channel source slot and 4x8 route mapping");
    analysis->CreateNtupleIColumn("plane");
    analysis->CreateNtupleIColumn("local_channel");
    analysis->CreateNtupleDColumn("source_transverse_mm");
    analysis->CreateNtupleDColumn("source_z_mm");
    analysis->CreateNtupleIColumn("fanout_row");
    analysis->CreateNtupleIColumn("fanout_column");
    analysis->CreateNtupleDColumn("active_low_mm");
    analysis->CreateNtupleDColumn("active_high_mm");
    analysis->CreateNtupleDColumn("far_end_mm");
    analysis->CreateNtupleDColumn("effective_outer_size_mm");
    analysis->CreateNtupleIColumn("physical_fibre_present");
    analysis->FinishNtuple();

    // Preserve the predecessor ROOT key set when the passive diagnostic is
    // disabled.  When enabled, the five Gate-1 trees append after legacy IDs
    // 0..6 and therefore receive stable IDs 7..11.
    if(fDetector->GetRecordSpatialClosureNtuples()) {
      analysis->CreateNtuple("primary_surface_truth",
        "Gate-1 propagated primary state at six fixed scoring surfaces");
    analysis->CreateNtupleIColumn("event");
    analysis->CreateNtupleIColumn("primary_track");
    analysis->CreateNtupleIColumn("surface");
    analysis->CreateNtupleIColumn("reached");
    analysis->CreateNtupleDColumn("x_mm");
    analysis->CreateNtupleDColumn("y_mm");
    analysis->CreateNtupleDColumn("z_mm");
    analysis->CreateNtupleDColumn("dir_x");
    analysis->CreateNtupleDColumn("dir_y");
    analysis->CreateNtupleDColumn("dir_z");
    analysis->CreateNtupleDColumn("kinetic_energy_gev");
    analysis->CreateNtupleDColumn("time_ns");
    analysis->CreateNtupleDColumn("surface_z_mm");
      analysis->FinishNtuple();

      analysis->CreateNtuple("muon_fibre_crossing",
        "Gate-1 contiguous charged-particle fibre crossings");
    analysis->CreateNtupleIColumn("event");
    analysis->CreateNtupleIColumn("track");
    analysis->CreateNtupleIColumn("parent");
    analysis->CreateNtupleIColumn("pdg");
    analysis->CreateNtupleIColumn("crossing_index");
    analysis->CreateNtupleIColumn("plane");
    analysis->CreateNtupleIColumn("physical_copy");
    analysis->CreateNtupleIColumn("local_channel");
    analysis->CreateNtupleIColumn("global_channel");
    analysis->CreateNtupleIColumn("volume_class");
    analysis->CreateNtupleIColumn("first_step");
    analysis->CreateNtupleIColumn("last_step");
    analysis->CreateNtupleDColumn("entry_x_mm");
    analysis->CreateNtupleDColumn("entry_y_mm");
    analysis->CreateNtupleDColumn("entry_z_mm");
    analysis->CreateNtupleDColumn("exit_x_mm");
    analysis->CreateNtupleDColumn("exit_y_mm");
    analysis->CreateNtupleDColumn("exit_z_mm");
    analysis->CreateNtupleDColumn("entry_u_mm");
    analysis->CreateNtupleDColumn("entry_v_mm");
    analysis->CreateNtupleDColumn("entry_w_mm");
    analysis->CreateNtupleDColumn("exit_u_mm");
    analysis->CreateNtupleDColumn("exit_v_mm");
    analysis->CreateNtupleDColumn("exit_w_mm");
    analysis->CreateNtupleDColumn("entry_dir_x");
    analysis->CreateNtupleDColumn("entry_dir_y");
    analysis->CreateNtupleDColumn("entry_dir_z");
    analysis->CreateNtupleDColumn("exit_dir_x");
    analysis->CreateNtupleDColumn("exit_dir_y");
    analysis->CreateNtupleDColumn("exit_dir_z");
    analysis->CreateNtupleDColumn("entry_energy_gev");
    analysis->CreateNtupleDColumn("exit_energy_gev");
    analysis->CreateNtupleDColumn("entry_time_ns");
    analysis->CreateNtupleDColumn("exit_time_ns");
    analysis->CreateNtupleDColumn("path_mm");
    analysis->CreateNtupleDColumn("core_path_mm");
    analysis->CreateNtupleDColumn("active_core_path_mm");
    analysis->CreateNtupleDColumn("edep_mev");
    analysis->CreateNtupleDColumn("far_end_core_path_mm");
    analysis->CreateNtupleDColumn("fanout_side_core_path_mm");
      analysis->FinishNtuple();

      analysis->CreateNtuple("photon_spatial_link",
        "Gate-1 photon origin-to-terminal channel relation");
    analysis->CreateNtupleIColumn("event");
    analysis->CreateNtupleIColumn("track");
    analysis->CreateNtupleIColumn("parent");
    analysis->CreateNtupleIColumn("plane");
    analysis->CreateNtupleIColumn("origin_physical_copy");
    analysis->CreateNtupleIColumn("origin_local_channel");
    analysis->CreateNtupleIColumn("origin_global_channel");
    analysis->CreateNtupleIColumn("terminal_global_channel");
    analysis->CreateNtupleIColumn("creator_process");
    analysis->CreateNtupleIColumn("terminal_reason");
    analysis->CreateNtupleIColumn("terminal_region");
    analysis->CreateNtupleDColumn("origin_time_ns");
    analysis->CreateNtupleDColumn("terminal_time_ns");
    analysis->CreateNtupleIColumn("birth_context_matched");
    analysis->CreateNtupleIColumn("parent_step");
      analysis->FinishNtuple();

      analysis->CreateNtuple("event_plane_reco",
        "Gate-1 per-event per-plane truth-to-ADC closure");
    analysis->CreateNtupleIColumn("event");
    analysis->CreateNtupleIColumn("plane");
    analysis->CreateNtupleIColumn("mid");
    analysis->CreateNtupleIColumn("primary_unique_core_channels");
    analysis->CreateNtupleDColumn("primary_core_path_mm");
    analysis->CreateNtupleDColumn("primary_active_core_path_mm");
    analysis->CreateNtupleDColumn("core_edep_mev");
    analysis->CreateNtupleIColumn("photons_produced");
    analysis->CreateNtupleIColumn("sipm_boundary_photons");
    analysis->CreateNtupleIColumn("pde_accepted_pe");
    analysis->CreateNtupleIColumn("pde_rejected_photons");
    analysis->CreateNtupleIColumn("best_global_channel");
    analysis->CreateNtupleIColumn("best_local_channel");
    analysis->CreateNtupleIColumn("best_detected_pe");
    analysis->CreateNtupleDColumn("best_noiseless_adc");
    analysis->CreateNtupleDColumn("best_final_adc");
    analysis->CreateNtupleDColumn("threshold_adc");
    analysis->CreateNtupleIColumn("hit");
      analysis->FinishNtuple();

      analysis->CreateNtuple("validation_geometry_map",
        "Gate-1 extended static channel identity and active-frame map");
    analysis->CreateNtupleIColumn("plane");
    analysis->CreateNtupleIColumn("mid");
    analysis->CreateNtupleIColumn("physical_copy");
    analysis->CreateNtupleIColumn("local_channel");
    analysis->CreateNtupleIColumn("global_channel");
    analysis->CreateNtupleDColumn("centre_x_mm");
    analysis->CreateNtupleDColumn("centre_y_mm");
    analysis->CreateNtupleDColumn("centre_z_mm");
    analysis->CreateNtupleDColumn("axis_x");
    analysis->CreateNtupleDColumn("axis_y");
    analysis->CreateNtupleDColumn("axis_z");
    analysis->CreateNtupleDColumn("transverse_x");
    analysis->CreateNtupleDColumn("transverse_y");
    analysis->CreateNtupleDColumn("transverse_z");
    analysis->CreateNtupleDColumn("normal_x");
    analysis->CreateNtupleDColumn("normal_y");
    analysis->CreateNtupleDColumn("normal_z");
    analysis->CreateNtupleDColumn("core_half_extent_mm");
    analysis->CreateNtupleDColumn("outer_half_extent_mm");
    analysis->CreateNtupleIColumn("fanout_row");
      analysis->CreateNtupleIColumn("fanout_column");
      analysis->CreateNtupleIColumn("physical_fibre_present");
      analysis->FinishNtuple();
    }
    if(fDetector->GetRecordGeometryControlNtuples()) {
      const auto triggerCrossingNtupleId=analysis->CreateNtuple(
        "primary_trigger_crossing",
        "Stage-2 primary truth inside each trigger volume");
      if(triggerCrossingNtupleId!=kPrimaryTriggerCrossingNtupleId) {
        G4ExceptionDescription message;
        message << "primary_trigger_crossing ntuple ID shifted: expected="
                << kPrimaryTriggerCrossingNtupleId
                << " actual=" << triggerCrossingNtupleId << ".";
        G4Exception("OpticalRunAction::BeginOfRunAction",
                    "HODO_GATE2_NTUPLE_ID",FatalException,message);
      }
      analysis->CreateNtupleIColumn("event");
      analysis->CreateNtupleIColumn("track");
      analysis->CreateNtupleIColumn("parent");
      analysis->CreateNtupleIColumn("pdg");
      analysis->CreateNtupleIColumn("trigger");
      analysis->CreateNtupleIColumn("reached");
      analysis->CreateNtupleDColumn("entry_x_mm");
      analysis->CreateNtupleDColumn("entry_y_mm");
      analysis->CreateNtupleDColumn("entry_z_mm");
      analysis->CreateNtupleDColumn("exit_x_mm");
      analysis->CreateNtupleDColumn("exit_y_mm");
      analysis->CreateNtupleDColumn("exit_z_mm");
      analysis->CreateNtupleDColumn("entry_local_x_mm");
      analysis->CreateNtupleDColumn("entry_local_y_mm");
      analysis->CreateNtupleDColumn("entry_local_z_mm");
      analysis->CreateNtupleDColumn("exit_local_x_mm");
      analysis->CreateNtupleDColumn("exit_local_y_mm");
      analysis->CreateNtupleDColumn("exit_local_z_mm");
      analysis->CreateNtupleDColumn("entry_dir_x");
      analysis->CreateNtupleDColumn("entry_dir_y");
      analysis->CreateNtupleDColumn("entry_dir_z");
      analysis->CreateNtupleDColumn("exit_dir_x");
      analysis->CreateNtupleDColumn("exit_dir_y");
      analysis->CreateNtupleDColumn("exit_dir_z");
      analysis->CreateNtupleDColumn("entry_time_ns");
      analysis->CreateNtupleDColumn("exit_time_ns");
      analysis->CreateNtupleDColumn("path_mm");
      analysis->CreateNtupleDColumn("edep_mev");
      analysis->CreateNtupleIColumn("steps");
      analysis->FinishNtuple();

      const auto triggerMapNtupleId=analysis->CreateNtuple(
        "validation_trigger_map",
        "Stage-2 configured trigger centre, axes, and half lengths");
      if(triggerMapNtupleId!=kValidationTriggerMapNtupleId) {
        G4ExceptionDescription message;
        message << "validation_trigger_map ntuple ID shifted: expected="
                << kValidationTriggerMapNtupleId
                << " actual=" << triggerMapNtupleId << ".";
        G4Exception("OpticalRunAction::BeginOfRunAction",
                    "HODO_GATE2_NTUPLE_ID",FatalException,message);
      }
      analysis->CreateNtupleIColumn("trigger");
      analysis->CreateNtupleDColumn("centre_x_mm");
      analysis->CreateNtupleDColumn("centre_y_mm");
      analysis->CreateNtupleDColumn("centre_z_mm");
      analysis->CreateNtupleDColumn("axis_x_x");
      analysis->CreateNtupleDColumn("axis_x_y");
      analysis->CreateNtupleDColumn("axis_x_z");
      analysis->CreateNtupleDColumn("axis_y_x");
      analysis->CreateNtupleDColumn("axis_y_y");
      analysis->CreateNtupleDColumn("axis_y_z");
      analysis->CreateNtupleDColumn("axis_z_x");
      analysis->CreateNtupleDColumn("axis_z_y");
      analysis->CreateNtupleDColumn("axis_z_z");
      analysis->CreateNtupleDColumn("half_x_mm");
      analysis->CreateNtupleDColumn("half_y_mm");
      analysis->CreateNtupleDColumn("half_z_mm");
      analysis->CreateNtupleDColumn("rotation_deg");
      analysis->CreateNtupleDColumn("yaw_command_deg");
      analysis->FinishNtuple();

      const auto provenanceNtupleId=analysis->CreateNtuple(
        "validation_geometry_provenance",
        "Stage-2 source-table and seed provenance");
      if(provenanceNtupleId!=kValidationGeometryProvenanceNtupleId) {
        G4ExceptionDescription message;
        message << "validation_geometry_provenance ntuple ID shifted: expected="
                << kValidationGeometryProvenanceNtupleId
                << " actual=" << provenanceNtupleId << ".";
        G4Exception("OpticalRunAction::BeginOfRunAction",
                    "HODO_GATE2_NTUPLE_ID",FatalException,message);
      }
      analysis->CreateNtupleIColumn("schema_version");
      analysis->CreateNtupleIColumn("geometry_probe_mode");
      analysis->CreateNtupleDColumn("cosmic_seed");
      analysis->CreateNtupleIColumn("row_offset");
      analysis->CreateNtupleIColumn("expected_rows");
      analysis->CreateNtupleSColumn("probe_file");
      analysis->FinishNtuple();
    }
    if(fDetector->GetRecordGenerationNtuples()) {
      const auto stepId=analysis->CreateNtuple(
          "validation_generation_step",
          "Stage-3 charged fibre steps and direct optical secondaries");
      if(stepId!=kValidationGenerationStepNtupleId) {
        G4ExceptionDescription message;
        message << "validation_generation_step ntuple ID shifted: expected="
                << kValidationGenerationStepNtupleId
                << " actual=" << stepId << ".";
        G4Exception("OpticalRunAction::BeginOfRunAction",
                    "HODO_GATE3_NTUPLE_ID",FatalException,message);
      }
      const std::array<G4String,20> stepIntegerColumns{{
          "event","source_row","cohort","target_plane","target_channel",
          "expected_volume_class","expected_region","chord_bin",
          "repeat_index","case_index","track","parent","pdg",
          "step_number","plane","physical_copy","local_channel",
          "global_channel","volume_class","material_code"}};
      for(const auto& name:stepIntegerColumns)
        analysis->CreateNtupleIColumn(name);
      const std::array<G4String,14> stepDoubleColumns{{
          "step_length_mm","active_path_mm","total_edep_mev",
          "nonionizing_edep_mev","pre_x_mm","pre_y_mm","pre_z_mm",
          "post_x_mm","post_y_mm","post_z_mm","actual_yield_per_mev",
          "actual_birks_mm_per_mev",
          "active_optical_birth_max_abs_u_mm",
          "fanout_optical_birth_min_abs_u_mm"}};
      for(const auto& name:stepDoubleColumns)
        analysis->CreateNtupleDColumn(name);
      const std::array<G4String,9> stepPhotonColumns{{
          "scintillation_count","cerenkov_count","other_optical_count",
          "scint_active_count","scint_fanout_count",
          "cerenkov_active_count","cerenkov_fanout_count",
          "other_active_count","other_fanout_count"}};
      for(const auto& name:stepPhotonColumns)
        analysis->CreateNtupleIColumn(name);
      analysis->FinishNtuple();

      const auto sourceId=analysis->CreateNtuple(
          "validation_generation_source",
          "Stage-3 event-to-frozen-TSV source bridge");
      if(sourceId!=kValidationGenerationSourceNtupleId) {
        G4ExceptionDescription message;
        message << "validation_generation_source ntuple ID shifted: expected="
                << kValidationGenerationSourceNtupleId
                << " actual=" << sourceId << ".";
        G4Exception("OpticalRunAction::BeginOfRunAction",
                    "HODO_GATE3_NTUPLE_ID",FatalException,message);
      }
      const std::array<G4String,9> sourceIntegerPrefix{{
          "event","source_row","cohort","target_plane","target_channel",
          "expected_volume_class","expected_region","chord_bin",
          "repeat_index"}};
      for(const auto& name:sourceIntegerPrefix)
        analysis->CreateNtupleIColumn(name);
      const std::array<G4String,9> sourceDoubleColumns{{
          "x_mm","y_mm","z_mm","dx","dy","dz","energy_gev",
          "analytic_core_chord_mm","analytic_target_path_mm"}};
      for(const auto& name:sourceDoubleColumns)
        analysis->CreateNtupleDColumn(name);
      analysis->CreateNtupleIColumn("case_index");
      analysis->FinishNtuple();

      const auto eventId=analysis->CreateNtuple(
          "stage3_generation_event",
          "Stage-3 target-plane step and creator ledger");
      if(eventId!=kStage3GenerationEventNtupleId) {
        G4ExceptionDescription message;
        message << "stage3_generation_event ntuple ID shifted: expected="
                << kStage3GenerationEventNtupleId
                << " actual=" << eventId << ".";
        G4Exception("OpticalRunAction::BeginOfRunAction",
                    "HODO_GATE3_NTUPLE_ID",FatalException,message);
      }
      const std::array<G4String,10> eventIntegerPrefix{{
          "event","source_row","cohort","target_plane","target_channel",
          "expected_volume_class","expected_region","chord_bin",
          "repeat_index","case_index"}};
      for(const auto& name:eventIntegerPrefix)
        analysis->CreateNtupleIColumn(name);
      const std::array<G4String,17> eventDoubleColumns{{
          "analytic_core_chord_mm","analytic_target_path_mm",
          "primary_core_path_mm","primary_active_core_path_mm",
          "primary_fanout_core_path_mm","primary_cladding_path_mm",
          "primary_active_cladding_path_mm",
          "primary_fanout_cladding_path_mm",
          "primary_core_total_edep_mev",
          "primary_core_nonionizing_edep_mev",
          "primary_core_birks_visible_edep_mev",
          "all_core_total_edep_mev","all_core_nonionizing_edep_mev",
          "primary_cladding_total_edep_mev",
          "primary_cladding_nonionizing_edep_mev",
          "all_cladding_total_edep_mev",
          "all_cladding_nonionizing_edep_mev"}};
      for(const auto& name:eventDoubleColumns)
        analysis->CreateNtupleDColumn(name);
      const std::array<G4String,19> eventIntegerSuffix{{
          "primary_core_nonzero_niel_steps",
          "primary_scintillation_count","secondary_scintillation_count",
          "primary_cerenkov_count","secondary_cerenkov_count",
          "primary_other_optical_count","secondary_other_optical_count",
          "scint_active_count","scint_fanout_count",
          "cerenkov_active_count","cerenkov_fanout_count",
          "other_active_count","other_fanout_count",
          "cladding_scintillation_count","cladding_cerenkov_count",
          "cladding_other_optical_count","direct_optical_total",
          "charged_fibre_steps","primary_core_steps"}};
      for(const auto& name:eventIntegerSuffix)
        analysis->CreateNtupleIColumn(name);
      const std::array<G4String,15> targetChannelDoubleColumns{{
          "target_channel_primary_core_path_mm",
          "target_channel_primary_active_core_path_mm",
          "target_channel_primary_fanout_core_path_mm",
          "target_channel_primary_cladding_path_mm",
          "target_channel_primary_active_cladding_path_mm",
          "target_channel_primary_fanout_cladding_path_mm",
          "target_channel_primary_core_total_edep_mev",
          "target_channel_primary_core_nonionizing_edep_mev",
          "target_channel_primary_core_birks_visible_edep_mev",
          "target_channel_all_core_total_edep_mev",
          "target_channel_all_core_nonionizing_edep_mev",
          "target_channel_primary_cladding_total_edep_mev",
          "target_channel_primary_cladding_nonionizing_edep_mev",
          "target_channel_all_cladding_total_edep_mev",
          "target_channel_all_cladding_nonionizing_edep_mev"}};
      for(const auto& name:targetChannelDoubleColumns)
        analysis->CreateNtupleDColumn(name);
      const std::array<G4String,19> targetChannelIntegerColumns{{
          "target_channel_primary_core_nonzero_niel_steps",
          "target_channel_primary_scintillation_count",
          "target_channel_secondary_scintillation_count",
          "target_channel_primary_cerenkov_count",
          "target_channel_secondary_cerenkov_count",
          "target_channel_primary_other_optical_count",
          "target_channel_secondary_other_optical_count",
          "target_channel_scint_active_count",
          "target_channel_scint_fanout_count",
          "target_channel_cerenkov_active_count",
          "target_channel_cerenkov_fanout_count",
          "target_channel_other_active_count",
          "target_channel_other_fanout_count",
          "target_channel_cladding_scintillation_count",
          "target_channel_cladding_cerenkov_count",
          "target_channel_cladding_other_optical_count",
          "target_channel_direct_optical_total",
          "target_channel_charged_fibre_steps",
          "target_channel_primary_core_steps"}};
      for(const auto& name:targetChannelIntegerColumns)
        analysis->CreateNtupleIColumn(name);
      analysis->FinishNtuple();

      const auto provenanceId=analysis->CreateNtuple(
          "validation_generation_provenance",
          "Stage-3 actual material, process, source, and run state");
      if(provenanceId!=kValidationGenerationProvenanceNtupleId) {
        G4ExceptionDescription message;
        message << "validation_generation_provenance ntuple ID shifted: "
                << "expected=" << kValidationGenerationProvenanceNtupleId
                << " actual=" << provenanceId << ".";
        G4Exception("OpticalRunAction::BeginOfRunAction",
                    "HODO_GATE3_NTUPLE_ID",FatalException,message);
      }
      analysis->CreateNtupleIColumn("schema_version");
      analysis->CreateNtupleIColumn("generation_probe_mode");
      analysis->CreateNtupleIColumn("record_generation_ntuples");
      analysis->CreateNtupleSColumn("source_file");
      analysis->CreateNtupleSColumn("source_sha256");
      analysis->CreateNtupleIColumn("source_rows");
      analysis->CreateNtupleIColumn("source_offset");
      analysis->CreateNtupleDColumn("cosmic_seed");
      analysis->CreateNtupleDColumn("clhep_seed_0");
      analysis->CreateNtupleDColumn("clhep_seed_1");
      analysis->CreateNtupleDColumn("requested_max_step_mm");
      analysis->CreateNtupleDColumn("observed_max_core_step_mm");
      analysis->CreateNtupleIColumn("generation_only");
      analysis->CreateNtupleIColumn("scintillation_active");
      analysis->CreateNtupleIColumn("cerenkov_active");
      analysis->CreateNtupleIColumn("em_saturation_attached");
      analysis->CreateNtupleIColumn("scint_by_particle_type");
      analysis->CreateNtupleDColumn("h1_yield_per_mev");
      analysis->CreateNtupleDColumn("h2_yield_per_mev");
      analysis->CreateNtupleDColumn("h1_birks_mm_per_mev");
      analysis->CreateNtupleDColumn("h2_birks_mm_per_mev");
      analysis->CreateNtupleDColumn("h1_resolution_scale");
      analysis->CreateNtupleDColumn("h2_resolution_scale");
      analysis->CreateNtupleDColumn("production_cut_mm");
      analysis->CreateNtupleIColumn("completed_events");
      analysis->CreateNtupleIColumn(
          "optical_photon_scintillation_absent");
      analysis->CreateNtupleIColumn("muon_msc_active");
      analysis->CreateNtupleIColumn(
          "muon_coulomb_scattering_active");
      analysis->CreateNtupleIColumn("muon_ionisation_active");
      analysis->CreateNtupleIColumn("muon_bremsstrahlung_active");
      analysis->CreateNtupleIColumn("muon_pair_production_active");
      analysis->CreateNtupleIColumn("muon_decay_active");
      analysis->FinishNtuple();
    }
    fAlignmentConfigurationNtupleId=-1;
    if(fDetector->GetRecordAlignmentStudyNtuple()) {
      fAlignmentConfigurationNtupleId=analysis->CreateNtuple(
          "alignment_configuration",
          "Actual alignment-study centres, rotations, and relative errors");
      analysis->CreateNtupleIColumn("schema_version");
      analysis->CreateNtupleIColumn("alignment_reference_code");
      analysis->CreateNtupleDColumn("t1_centre_x_mm");
      analysis->CreateNtupleDColumn("t1_centre_y_mm");
      analysis->CreateNtupleDColumn("t1_centre_z_mm");
      analysis->CreateNtupleDColumn("t2_centre_x_mm");
      analysis->CreateNtupleDColumn("t2_centre_y_mm");
      analysis->CreateNtupleDColumn("t2_centre_z_mm");
      analysis->CreateNtupleDColumn("t1_yaw_command_deg");
      analysis->CreateNtupleDColumn("t2_yaw_command_deg");
      analysis->CreateNtupleDColumn("t1_placement_angle_deg");
      analysis->CreateNtupleDColumn("t2_placement_angle_deg");
      analysis->CreateNtupleDColumn("h1_centre_x_mm");
      analysis->CreateNtupleDColumn("h1_centre_y_mm");
      analysis->CreateNtupleDColumn("h2_centre_x_mm");
      analysis->CreateNtupleDColumn("h2_centre_y_mm");
      analysis->CreateNtupleDColumn("h1x_yaw_deg");
      analysis->CreateNtupleDColumn("h1y_yaw_deg");
      analysis->CreateNtupleDColumn("h2x_yaw_deg");
      analysis->CreateNtupleDColumn("h2y_yaw_deg");
      analysis->CreateNtupleDColumn("trigger_relative_dx_mm");
      analysis->CreateNtupleDColumn("trigger_relative_dy_mm");
      analysis->CreateNtupleDColumn("trigger_relative_r_mm");
      analysis->CreateNtupleDColumn("trigger_relative_yaw_deg");
      analysis->CreateNtupleDColumn("hodoscope_relative_dx_mm");
      analysis->CreateNtupleDColumn("hodoscope_relative_dy_mm");
      analysis->CreateNtupleDColumn("hodoscope_relative_r_mm");
      analysis->CreateNtupleDColumn("hodoscope_relative_yaw_deg");
      analysis->FinishNtuple();
    }
    G4String rootName = fDetector->GetDiagnosticRootFile();
    if (rootName.size() < 5 || rootName.substr(rootName.size()-5) != ".root")
      rootName += ".root";
    if (!analysis->OpenFile(rootName)) {
      G4ExceptionDescription message;
      message << "Cannot open V21 analysis output: " << rootName
              << ". Verify that its parent directory exists and is writable.";
      G4Exception("OpticalRunAction::BeginOfRunAction",
                  "HODO_V21_ANALYSIS_OPEN", FatalException, message);
    }
    // With ROOT ntuple merging the master owns no fillable ntuple instance.
    // One identical provenance row per worker is intentional and merge-safe.
    const auto* currentRunManager=G4RunManager::GetRunManager();
    const G4bool fillGeometryProvenance =
        !IsMaster() ||
        (currentRunManager &&
         currentRunManager->GetRunManagerType()==G4RunManager::sequentialRM);
    if(fillGeometryProvenance) {
      analysis->FillNtupleIColumn(5,0,fDetector->GetH1RoundSides());
      analysis->FillNtupleIColumn(
          5,1,fDetector->GetH1FanoutHierarchyCode());
      for (G4int plane=0; plane<4; ++plane)
        analysis->FillNtupleIColumn(
            5,2+plane,fDetector->GetPlaneLayoutCode(plane));
      for (G4int plane=0; plane<4; ++plane)
        analysis->FillNtupleIColumn(
            5,6+plane,fDetector->GetPlaneChannelOrderCode(plane));
      for (G4int plane=0; plane<4; ++plane)
        analysis->FillNtupleDColumn(
            5,10+plane,fDetector->GetPlaneRowZSpacing(plane)/mm);
      analysis->FillNtupleIColumn(
          5,14,fDetector->GetMixedLayoutPattern());
      analysis->FillNtupleIColumn(
          5,15,fDetector->GetRecordActiveCoreChordBranches());
      analysis->FillNtupleIColumn(5,16,1);
      analysis->FillNtupleIColumn(5,17,0);
      analysis->FillNtupleIColumn(5,18,0);
      analysis->FillNtupleIColumn(5,19,fDetector->GetFanoutPreview());
      analysis->FillNtupleDColumn(
          5,20,fDetector->GetHodoscopeOffsetX(0)/mm);
      analysis->FillNtupleDColumn(
          5,21,fDetector->GetHodoscopeOffsetY(0)/mm);
      analysis->FillNtupleDColumn(
          5,22,fDetector->GetHodoscopeOffsetX(1)/mm);
      analysis->FillNtupleDColumn(
          5,23,fDetector->GetHodoscopeOffsetY(1)/mm);
      const auto& geometryTolerance =
          hodoscope::v21::GetGeometryToleranceSnapshot();
      analysis->FillNtupleDColumn(
          5,24,geometryTolerance.requestedWorldMaximumExtent/mm);
      analysis->FillNtupleDColumn(
          5,25,geometryTolerance.surfaceTolerance/mm);
      analysis->FillNtupleDColumn(
          5,26,geometryTolerance.radialTolerance/mm);
      analysis->FillNtupleIColumn(5,27,0);  // 0 = SerialOnly
      analysis->FillNtupleDColumn(
          5,28,fDetector->GetH1RoundPolygonPhase());
      analysis->FillNtupleIColumn(
          5,29,fDetector->GetH1RoundSweepRepresentationCode());
      analysis->FillNtupleDColumn(
          5,30,fDetector->GetH1RoundMiterJointOverlap()/mm);
      analysis->FillNtupleIColumn(
          5,31,fDetector->GetH1ActiveFacetRefinement());
      analysis->FillNtupleIColumn(
          5,32,fDetector->GetH2FanoutHierarchyCode());
      analysis->FillNtupleDColumn(5,33,fDetector->GetFibreFarEndExtension()/mm);
      analysis->FillNtupleIColumn(5,34,fDetector->GetActivePlaneOverlap());
      analysis->FillNtupleIColumn(5,35,fDetector->GetH2SquareFibreCount());
      analysis->FillNtupleIColumn(5,36,fDetector->GetRecordSelectedDiagnosticsOnly());
      analysis->AddNtupleRow(5);
      if(fAlignmentConfigurationNtupleId>=0) {
        const auto t1=fDetector->GetTriggerCentre(0);
        const auto t2=fDetector->GetTriggerCentre(1);
        const auto h1x=fDetector->GetHodoscopeOffsetX(0);
        const auto h1y=fDetector->GetHodoscopeOffsetY(0);
        const auto h2x=fDetector->GetHodoscopeOffsetX(1);
        const auto h2y=fDetector->GetHodoscopeOffsetY(1);
        const auto tdx=t2.x()-t1.x();
        const auto tdy=t2.y()-t1.y();
        const auto hdx=h2x-h1x;
        const auto hdy=h2y-h1y;
        const auto id=fAlignmentConfigurationNtupleId;
        analysis->FillNtupleIColumn(id,0,1);
        analysis->FillNtupleIColumn(
            id,1,fDetector->GetAlignmentReferenceCode());
        analysis->FillNtupleDColumn(id,2,t1.x()/mm);
        analysis->FillNtupleDColumn(id,3,t1.y()/mm);
        analysis->FillNtupleDColumn(id,4,t1.z()/mm);
        analysis->FillNtupleDColumn(id,5,t2.x()/mm);
        analysis->FillNtupleDColumn(id,6,t2.y()/mm);
        analysis->FillNtupleDColumn(id,7,t2.z()/mm);
        analysis->FillNtupleDColumn(id,8,fDetector->GetTriggerYaw(0)/deg);
        analysis->FillNtupleDColumn(id,9,fDetector->GetTriggerYaw(1)/deg);
        analysis->FillNtupleDColumn(
            id,10,fDetector->GetTriggerRotationAngle(0)/deg);
        analysis->FillNtupleDColumn(
            id,11,fDetector->GetTriggerRotationAngle(1)/deg);
        analysis->FillNtupleDColumn(id,12,h1x/mm);
        analysis->FillNtupleDColumn(id,13,h1y/mm);
        analysis->FillNtupleDColumn(id,14,h2x/mm);
        analysis->FillNtupleDColumn(id,15,h2y/mm);
        for(G4int plane=0;plane<4;++plane)
          analysis->FillNtupleDColumn(
              id,16+plane,fDetector->GetPlaneYaw(plane)/deg);
        analysis->FillNtupleDColumn(id,20,tdx/mm);
        analysis->FillNtupleDColumn(id,21,tdy/mm);
        analysis->FillNtupleDColumn(id,22,std::hypot(tdx,tdy)/mm);
        analysis->FillNtupleDColumn(
            id,23,std::abs(fDetector->GetTriggerYaw(1)-
                           fDetector->GetTriggerYaw(0))/deg);
        analysis->FillNtupleDColumn(id,24,hdx/mm);
        analysis->FillNtupleDColumn(id,25,hdy/mm);
        analysis->FillNtupleDColumn(id,26,std::hypot(hdx,hdy)/mm);
        analysis->FillNtupleDColumn(
            id,27,std::abs(fDetector->GetPlaneYaw(2)-
                           fDetector->GetPlaneYaw(0))/deg);
        analysis->AddNtupleRow(id);
      }
      for (G4int plane=0; plane<4; ++plane) {
        for (G4int channel=0; channel<32; ++channel) {
          analysis->FillNtupleIColumn(6,0,plane);
          analysis->FillNtupleIColumn(6,1,channel);
          analysis->FillNtupleDColumn(
              6,2,fDetector->GetPlaneSourceTransverse(plane,channel)/mm);
          analysis->FillNtupleDColumn(
              6,3,fDetector->GetPlaneSourceZ(plane,channel)/mm);
          analysis->FillNtupleIColumn(
              6,4,fDetector->GetPlaneFanoutRow(plane,channel));
          analysis->FillNtupleIColumn(
              6,5,fDetector->GetPlaneFanoutColumn(plane,channel));
          const auto activeBounds=fDetector->GetActiveLongitudinalBounds(plane);
          analysis->FillNtupleDColumn(6,6,activeBounds[0]/mm);
          analysis->FillNtupleDColumn(6,7,activeBounds[1]/mm);
          analysis->FillNtupleDColumn(6,8,fDetector->GetFibreFarEndExtension()>0.
              ? (activeBounds[0]-fDetector->GetFibreFarEndExtension())/mm : -16.);
          analysis->FillNtupleDColumn(6,9,fDetector->GetFanoutOuterSize(plane)/mm);
          analysis->FillNtupleIColumn(6,10,fDetector->HasPhysicalFibre(plane,channel));
          analysis->AddNtupleRow(6);
        }
      }
      if(fDetector->GetRecordSpatialClosureNtuples()) {
        constexpr std::array<G4int,4> mid{{1,9,15,13}};
        for(G4int plane=0;plane<4;++plane) {
          const G4bool xAxis=plane==0||plane==2;
          const auto axis=fDetector->PlaneNominalToWorldDirection(
              plane,xAxis?G4ThreeVector(1.,0.,0.)
                         :G4ThreeVector(0.,1.,0.)).unit();
          const auto transverse=fDetector->PlaneNominalToWorldDirection(
              plane,xAxis?G4ThreeVector(0.,1.,0.)
                         :G4ThreeVector(1.,0.,0.)).unit();
          const auto normal=fDetector->PlaneNominalToWorldDirection(
              plane,G4ThreeVector(0.,0.,1.)).unit();
          const G4double coreHalf=plane<2
              ? (fDetector->GetH1SingleClad()?.4795*mm:.4695*mm)
              : .4695*mm;
          for(G4int channel=0;channel<32;++channel) {
            const G4double sourceTransverse=
                fDetector->GetPlaneSourceTransverse(plane,channel);
            const G4double sourceZ=
                fDetector->GetPlaneSourceZ(plane,channel);
            const auto nominalCentre=xAxis
                ? G4ThreeVector(0.,sourceTransverse,sourceZ)
                : G4ThreeVector(sourceTransverse,0.,sourceZ);
            const auto centre=fDetector->PlaneNominalToWorldPosition(
                plane,nominalCentre);
            const G4int physicalCopy=(plane%2)*32+channel;
            const G4int globalChannel=plane*32+channel;
            analysis->FillNtupleIColumn(11,0,plane);
            analysis->FillNtupleIColumn(11,1,mid[plane]);
            analysis->FillNtupleIColumn(11,2,physicalCopy);
            analysis->FillNtupleIColumn(11,3,channel);
            analysis->FillNtupleIColumn(11,4,globalChannel);
            analysis->FillNtupleDColumn(11,5,centre.x()/mm);
            analysis->FillNtupleDColumn(11,6,centre.y()/mm);
            analysis->FillNtupleDColumn(11,7,centre.z()/mm);
            analysis->FillNtupleDColumn(11,8,axis.x());
            analysis->FillNtupleDColumn(11,9,axis.y());
            analysis->FillNtupleDColumn(11,10,axis.z());
            analysis->FillNtupleDColumn(11,11,transverse.x());
            analysis->FillNtupleDColumn(11,12,transverse.y());
            analysis->FillNtupleDColumn(11,13,transverse.z());
            analysis->FillNtupleDColumn(11,14,normal.x());
            analysis->FillNtupleDColumn(11,15,normal.y());
            analysis->FillNtupleDColumn(11,16,normal.z());
            analysis->FillNtupleDColumn(11,17,coreHalf/mm);
            analysis->FillNtupleDColumn(11,18,.4995);
            analysis->FillNtupleIColumn(
                11,19,fDetector->GetPlaneFanoutRow(plane,channel));
            analysis->FillNtupleIColumn(
                11,20,fDetector->GetPlaneFanoutColumn(plane,channel));
            analysis->FillNtupleIColumn(11,21,fDetector->HasPhysicalFibre(plane,channel));
            analysis->AddNtupleRow(11);
          }
        }
      }
      if(fDetector->GetRecordGeometryControlNtuples()) {
        for(G4int trigger=0;trigger<4;++trigger) {
          if(!fDetector->GetTriggerEnabled(trigger))continue;
          const auto centre=fDetector->GetTriggerCentre(trigger);
          const auto angle=fDetector->GetTriggerRotationAngle(trigger);
          const G4ThreeVector axisX(std::cos(angle),std::sin(angle),0.);
          const G4ThreeVector axisY(-std::sin(angle),std::cos(angle),0.);
          const G4ThreeVector axisZ(0.,0.,1.);
          analysis->FillNtupleIColumn(
              kValidationTriggerMapNtupleId,0,trigger);
          analysis->FillNtupleDColumn(
              kValidationTriggerMapNtupleId,1,centre.x()/mm);
          analysis->FillNtupleDColumn(
              kValidationTriggerMapNtupleId,2,centre.y()/mm);
          analysis->FillNtupleDColumn(
              kValidationTriggerMapNtupleId,3,centre.z()/mm);
          analysis->FillNtupleDColumn(
              kValidationTriggerMapNtupleId,4,axisX.x());
          analysis->FillNtupleDColumn(
              kValidationTriggerMapNtupleId,5,axisX.y());
          analysis->FillNtupleDColumn(
              kValidationTriggerMapNtupleId,6,axisX.z());
          analysis->FillNtupleDColumn(
              kValidationTriggerMapNtupleId,7,axisY.x());
          analysis->FillNtupleDColumn(
              kValidationTriggerMapNtupleId,8,axisY.y());
          analysis->FillNtupleDColumn(
              kValidationTriggerMapNtupleId,9,axisY.z());
          analysis->FillNtupleDColumn(
              kValidationTriggerMapNtupleId,10,axisZ.x());
          analysis->FillNtupleDColumn(
              kValidationTriggerMapNtupleId,11,axisZ.y());
          analysis->FillNtupleDColumn(
              kValidationTriggerMapNtupleId,12,axisZ.z());
          analysis->FillNtupleDColumn(
              kValidationTriggerMapNtupleId,13,
              fDetector->GetTriggerHalfLength(trigger)/mm);
          analysis->FillNtupleDColumn(
              kValidationTriggerMapNtupleId,14,
              fDetector->GetTriggerHalfWidth(trigger)/mm);
          analysis->FillNtupleDColumn(
              kValidationTriggerMapNtupleId,15,
              fDetector->GetTriggerHalfThickness(trigger)/mm);
          analysis->FillNtupleDColumn(
              kValidationTriggerMapNtupleId,16,angle/deg);
          analysis->FillNtupleDColumn(
              kValidationTriggerMapNtupleId,17,
              fDetector->GetTriggerYaw(trigger)/deg);
          analysis->AddNtupleRow(kValidationTriggerMapNtupleId);
        }
        analysis->FillNtupleIColumn(
            kValidationGeometryProvenanceNtupleId,0,1);
        analysis->FillNtupleIColumn(
            kValidationGeometryProvenanceNtupleId,1,
            fDetector->GetGeometryProbeMode());
        analysis->FillNtupleDColumn(
            kValidationGeometryProvenanceNtupleId,2,
            static_cast<G4double>(fDetector->GetCosmicSeed()));
        analysis->FillNtupleIColumn(
            kValidationGeometryProvenanceNtupleId,3,
            fDetector->GetGeometryProbeRowOffset());
        analysis->FillNtupleIColumn(
            kValidationGeometryProvenanceNtupleId,4,
            DetectorConstruction::kGeometryProbeExpectedRows);
        analysis->FillNtupleSColumn(
            kValidationGeometryProvenanceNtupleId,5,
            fDetector->GetGeometryProbeFile());
        analysis->AddNtupleRow(kValidationGeometryProvenanceNtupleId);
      }
    }
  }
  void EndOfRunAction(const G4Run* base) override {
    if (IsMaster()) {
      const auto* run = static_cast<const OpticalRun*>(base);
      const auto denominator = run->H1();
      const G4double efficiency = denominator > 0
          ? static_cast<G4double>(run->Both())/denominator : 0.;
      const G4double collection = run->Produced() > 0
          ? static_cast<G4double>(run->Detected())/run->Produced() : 0.;
      const auto triggered = run->Triggered();
      const auto triggeredEfficiency = [triggered](G4long numerator) {
        return triggered > 0 ? static_cast<G4double>(numerator)/triggered : 0.;
      };
      G4cout << "\n========== OPTICAL HODOSCOPE V21-MT ==========\n"
             << "Generated events                 = " << base->GetNumberOfEvent() << "\n"
             << "N(T1 && T2)                      = " << run->Triggered() << "\n"
             << "N(T1 && T2 && T3)                = " << run->TriggeredT123() << "\n"
             << "N(T1 && T2 && T3 && T4)          = " << run->TriggeredT1234() << "\n"
             << "triggerCoincidenceMode           = "
             << fDetector->GetTriggerCoincidenceMode() << "\n"
             << "----- STRUCTURAL PATH TRUTH: primary core chord, T1&&T2 denominator -----\n"
             << "active core T-only               = "
             << run->StructuralActiveCore(0) << "\n"
             << "active core H1-only              = "
             << run->StructuralActiveCore(1) << "\n"
             << "active core H2-only              = "
             << run->StructuralActiveCore(2) << "\n"
             << "active core H1&H2                = "
             << run->StructuralActiveCore(3) << "\n"
             << "active plane MID1/MID9/MID15/MID13 = "
             << run->StructuralPlaneActiveCore(0) << "/"
             << run->StructuralPlaneActiveCore(1) << "/"
             << run->StructuralPlaneActiveCore(2) << "/"
             << run->StructuralPlaneActiveCore(3) << "\n"
             << "total core T-only                = "
             << run->StructuralTotalCore(0) << "\n"
             << "total core H1-only               = "
             << run->StructuralTotalCore(1) << "\n"
             << "total core H2-only               = "
             << run->StructuralTotalCore(2) << "\n"
             << "total core H1&H2                 = "
             << run->StructuralTotalCore(3) << "\n"
             << "total plane MID1/MID9/MID15/MID13 = "
             << run->StructuralPlaneTotalCore(0) << "/"
             << run->StructuralPlaneTotalCore(1) << "/"
             << run->StructuralPlaneTotalCore(2) << "/"
             << run->StructuralPlaneTotalCore(3) << "\n"
             << "----- STRUCTURAL PATH TRUTH: active core chord, T1&&T2&&T3 denominator -----\n"
             << "T123 active core T-only          = "
             << run->StructuralT123ActiveCore(0) << "\n"
             << "T123 active core H1-only         = "
             << run->StructuralT123ActiveCore(1) << "\n"
             << "T123 active core H2-only         = "
             << run->StructuralT123ActiveCore(2) << "\n"
             << "T123 active core H1&H2           = "
             << run->StructuralT123ActiveCore(3) << "\n"
             << "T123 active plane MID1/MID9/MID15/MID13 = "
             << run->StructuralT123PlaneActiveCore(0) << "/"
             << run->StructuralT123PlaneActiveCore(1) << "/"
             << run->StructuralT123PlaneActiveCore(2) << "/"
             << run->StructuralT123PlaneActiveCore(3) << "\n"
             << "----- STRUCTURAL PATH TRUTH: active core chord, T1&&T2&&T3&&T4 denominator -----\n"
             << "T1234 active core T-only         = "
             << run->StructuralT1234ActiveCore(0) << "\n"
             << "T1234 active core H1-only        = "
             << run->StructuralT1234ActiveCore(1) << "\n"
             << "T1234 active core H2-only        = "
             << run->StructuralT1234ActiveCore(2) << "\n"
             << "T1234 active core H1&H2          = "
             << run->StructuralT1234ActiveCore(3) << "\n"
             << "T1234 active plane MID1/MID9/MID15/MID13 = "
             << run->StructuralT1234PlaneActiveCore(0) << "/"
             << run->StructuralT1234PlaneActiveCore(1) << "/"
             << run->StructuralT1234PlaneActiveCore(2) << "/"
             << run->StructuralT1234PlaneActiveCore(3) << "\n"
             << "----- VIRTUAL ALL-ROUND HODOSCOPE PATH ORACLE -----\n"
             << "virtualRoundModelEnabled         = "
             << fDetector->GetVirtualRoundModelEnabled() << "\n"
             << "virtualRoundH3Enabled            = "
             << fDetector->GetVirtualRoundH3Enabled() << "\n"
             << "virtualRoundH3Z [mm]             = "
             << fDetector->GetVirtualRoundH3Z()/mm << "\n"
             << "T1T2 virtual H1R/H2R/H3R         = "
             << run->VirtualT12H1R() << "/"
             << run->VirtualT12H2R() << "/"
             << run->VirtualT12H3R() << "\n"
             << "T1T2 virtual H1R&&H2R            = "
             << run->VirtualT12H1H2() << "\n"
             << "T1T2 virtual H1R&&H2R&&H3R       = "
             << run->VirtualT12H1H2H3() << "\n"
             << "T123 virtual H1R/H2R/H3R         = "
             << run->VirtualT123H1R() << "/"
             << run->VirtualT123H2R() << "/"
             << run->VirtualT123H3R() << "\n"
             << "T123 virtual H1R&&H2R            = "
             << run->VirtualT123H1H2() << "\n"
             << "T123 virtual H1R&&H2R&&H3R       = "
             << run->VirtualT123H1H2H3() << "\n"
             << "----- PRIMARY: Park/Yewon [200,600), channel mean+5RMS, plane-AND -----\n"
             << "N(Park MID1 round X | T1T2)      = " << run->ParkRoundX() << "\n"
             << "N(Park MID9 round Y | T1T2)      = " << run->ParkRoundY() << "\n"
             << "N(Park MID15 square X | T1T2)    = " << run->ParkSquareX() << "\n"
             << "N(Park MID13 square Y | T1T2)    = " << run->ParkSquareY() << "\n"
             << "N(Park MID1 && MID9 | T1T2)      = " << run->ParkRoundPair() << "\n"
             << "N(Park MID15 && MID13 | T1T2)    = " << run->ParkSquarePair() << "\n"
             << "efficiency(Park round XY | T1T2) = "
             << triggeredEfficiency(run->ParkRoundPair()) << "\n"
             << "efficiency(Park square XY | T1T2)= "
             << triggeredEfficiency(run->ParkSquarePair()) << "\n"
             << "----- DIAGNOSTIC: timing/noise -----\n"
             << "efficiency(Park round timed)      = "
             << triggeredEfficiency(run->ParkRoundTimed()) << "\n"
             << "efficiency(Park square timed)     = "
             << triggeredEfficiency(run->ParkSquareTimed()) << "\n"
             << "false-hit(Park round noise-only)  = "
             << triggeredEfficiency(run->ParkRoundNoise()) << "\n"
             << "false-hit(Park square noise-only) = "
             << triggeredEfficiency(run->ParkSquareNoise()) << "\n"
             << "----- DIAGNOSTIC LEGACY: fixed 120 ADC OR/pair -----\n"
             << "N(T1 && T2 && H1)                = " << run->H1() << "\n"
             << "N(T1 && T2 && H2)                = " << run->H2() << "\n"
             << "N(T1 && T2 && H1 && H2)          = " << run->Both() << "\n"
             << "efficiency(H2 | H1, triggered)   = " << efficiency << "\n"
             << "N(T1 && T2 && MID1 round X)      = " << run->RoundX() << "\n"
             << "N(T1 && T2 && MID9 round Y)      = " << run->RoundY() << "\n"
             << "N(T1 && T2 && MID15 "
             << (fDetector->GetH2PhysicalRound() ? "H2 round X)  = "
                                                : "square X)    = ")
             << run->SquareX() << "\n"
             << "N(T1 && T2 && MID13 "
             << (fDetector->GetH2PhysicalRound() ? "H2 round Y)  = "
                                                : "square Y)    = ")
             << run->SquareY() << "\n"
             << "N(T1 && T2 && MID1 && MID9)      = " << run->RoundPair() << "\n"
             << "N(T1 && T2 && MID13 && MID15)    = " << run->SquarePair() << "\n"
             << "efficiency(round XY | T1T2)      = "
             << triggeredEfficiency(run->RoundPair()) << "\n"
             << "efficiency("
             << (fDetector->GetH2PhysicalRound() ? "H2 round XY" : "square XY")
             << " | T1T2)     = "
             << triggeredEfficiency(run->SquarePair()) << "\n"
             << "efficiency(all four | T1T2)      = "
             << triggeredEfficiency(run->AllPlanes()) << "\n"
             << "optical photons produced         = " << run->Produced() << "\n"
             << "photons detected at SiPM surface = " << run->Detected() << "\n"
             << "global optical detection fraction= " << collection << "\n"
             << "H1 round edep [MeV]             = " << run->EdepH1()/MeV << "\n"
             << "H1 round photons produced       = " << run->ProducedH1() << "\n"
             << "H1 round photons detected       = " << run->DetectedH1() << "\n"
             << "H1 stage active/fanout-end/window/SiPM boundary = "
             << run->StageH1(0) << "/" << run->StageH1(1) << "/"
             << run->StageH1(2) << "/" << run->StageH1(3) << "\n"
             << "H1 PDE rejected                 = " << run->RejectedH1() << "\n"
             << "H1 round detected/produced      = "
             << (run->ProducedH1() ? G4double(run->DetectedH1())/run->ProducedH1() : 0.)
             << "\n"
             << "H2 "
             << (fDetector->GetH2PhysicalRound() ? "round" : "square")
             << " edep [MeV]            = " << run->EdepH2()/MeV << "\n"
             << "H2 "
             << (fDetector->GetH2PhysicalRound() ? "round" : "square")
             << " photons produced      = " << run->ProducedH2() << "\n"
             << "H2 "
             << (fDetector->GetH2PhysicalRound() ? "round" : "square")
             << " photons detected      = " << run->DetectedH2() << "\n"
             << "H2 stage active/fanout-end/window/SiPM boundary = "
             << run->StageH2(0) << "/" << run->StageH2(1) << "/"
             << run->StageH2(2) << "/" << run->StageH2(3) << "\n"
             << "H2 PDE rejected                 = " << run->RejectedH2() << "\n"
             << "H2 "
             << (fDetector->GetH2PhysicalRound() ? "round" : "square")
             << " detected/produced     = "
             << (run->ProducedH2() ? G4double(run->DetectedH2())/run->ProducedH2() : 0.)
             << "\n"
             << "============================================\n" << G4endl;
    }
    auto* analysis = G4AnalysisManager::Instance();
    if(fDetector->GetRecordGenerationNtuples()) {
      if(gGenerationProbeLoadedRows.load(std::memory_order_relaxed)!=
             DetectorConstruction::kGenerationProbeExpectedRows)
        G4Exception("OpticalRunAction::EndOfRunAction",
                    "HODO_GATE3_SOURCE_NOT_LOADED",FatalException,
                    "Stage-3 did not parse exactly 196960 frozen source rows "
                    "before finalization.");
      const auto readMaterial=[](const G4Material* material,
                                 const char* label) {
        struct Values { G4double yield;G4double birks;G4double resolution; };
        if(!material||!material->GetMaterialPropertiesTable()||
           !material->GetIonisation()) {
          G4ExceptionDescription message;
          message << "Missing initialized Stage-3 " << label
                  << " material/MPT/ionisation state.";
          G4Exception("OpticalRunAction::EndOfRunAction",
                      "HODO_GATE3_MATERIAL_PROVENANCE",FatalException,
                      message);
        }
        const auto* table=material->GetMaterialPropertiesTable();
        if(!table->ConstPropertyExists("SCINTILLATIONYIELD")||
           !table->ConstPropertyExists("RESOLUTIONSCALE")) {
          G4ExceptionDescription message;
          message << "Stage-3 " << label
                  << " lacks SCINTILLATIONYIELD or RESOLUTIONSCALE.";
          G4Exception("OpticalRunAction::EndOfRunAction",
                      "HODO_GATE3_MATERIAL_PROPERTIES",FatalException,
                      message);
        }
        return Values{
            table->GetConstProperty("SCINTILLATIONYIELD")*MeV,
            material->GetIonisation()->GetBirksConstant()/(mm/MeV),
            table->GetConstProperty("RESOLUTIONSCALE")};
      };
      const auto h1=readMaterial(fDetector->GetH1CoreMaterial(),"H1 core");
      const auto h2=readMaterial(fDetector->GetH2CoreMaterial(),"H2 core");

      G4int scintillationActive=-1,cerenkovActive=-1;
      G4int muonMscActive=-1,muonCoulombScatteringActive=-1;
      G4int muonIonisationActive=-1,muonBremsstrahlungActive=-1;
      G4int muonPairProductionActive=-1,muonDecayActive=-1;
      G4int emSaturationAttached=-1,scintByParticleType=-1;
      auto* muon=G4ParticleTable::GetParticleTable()->FindParticle("mu-");
      auto* processManager=muon?muon->GetProcessManager():nullptr;
      if(processManager) {
        auto* processes=processManager->GetProcessList();
        for(G4int i=0;i<processManager->GetProcessListLength();++i) {
          auto* process=(*processes)[i];
          if(!process)continue;
          if(process->GetProcessName()=="Scintillation") {
            scintillationActive=processManager->GetProcessActivation(process);
            const auto* scintillation=dynamic_cast<G4Scintillation*>(process);
            if(scintillation) {
              emSaturationAttached=scintillation->GetSaturation()!=nullptr;
              scintByParticleType=
                  scintillation->GetScintillationByParticleType();
            }
          } else if(process->GetProcessName()=="Cerenkov") {
            cerenkovActive=processManager->GetProcessActivation(process);
          } else if(process->GetProcessName()=="msc") {
            muonMscActive=processManager->GetProcessActivation(process);
          } else if(process->GetProcessName()=="CoulombScat") {
            muonCoulombScatteringActive=
                processManager->GetProcessActivation(process);
          } else if(process->GetProcessName()=="muIoni") {
            muonIonisationActive=
                processManager->GetProcessActivation(process);
          } else if(process->GetProcessName()=="muBrems") {
            muonBremsstrahlungActive=
                processManager->GetProcessActivation(process);
          } else if(process->GetProcessName()=="muPairProd") {
            muonPairProductionActive=
                processManager->GetProcessActivation(process);
          } else if(process->GetProcessName()=="Decay") {
            muonDecayActive=processManager->GetProcessActivation(process);
          }
        }
      }
      // Geant4 may omit a process that was disabled in PreInit from the
      // particle's process vector entirely. Absence plus an inactive
      // initialized optical-parameter state is an observed inactive process,
      // not missing provenance. An active parameter without a process remains
      // fatal below.
      const auto* opticalParameters=G4OpticalParameters::Instance();
      if(cerenkovActive<0&&opticalParameters&&
         !opticalParameters->GetProcessActivation("Cerenkov"))
        cerenkovActive=0;
      if(scintillationActive<0&&opticalParameters&&
         !opticalParameters->GetProcessActivation("Scintillation"))
        scintillationActive=0;
      const auto* runManager=G4RunManager::GetRunManager();
      const auto* physics=runManager?runManager->GetUserPhysicsList():nullptr;
      if(scintillationActive<0||cerenkovActive<0||
         emSaturationAttached<0||scintByParticleType<0||
         muonMscActive<0||muonCoulombScatteringActive<0||
         muonIonisationActive<0||muonBremsstrahlungActive<0||
         muonPairProductionActive<0||muonDecayActive<0||!physics) {
        G4ExceptionDescription message;
        message << "Cannot recover initialized mu- optical process state: "
                << "scintillation_active=" << scintillationActive
                << " cerenkov_active=" << cerenkovActive
                << " em_saturation_attached=" << emSaturationAttached
                << " scint_by_particle_type=" << scintByParticleType
                << " muon_msc_active=" << muonMscActive
                << " muon_coulomb_scattering_active="
                << muonCoulombScatteringActive
                << " muon_ionisation_active=" << muonIonisationActive
                << " muon_bremsstrahlung_active="
                << muonBremsstrahlungActive
                << " muon_pair_production_active="
                << muonPairProductionActive
                << " muon_decay_active=" << muonDecayActive
                << " physics_list=" << (physics?1:0) << ".";
        G4Exception("OpticalRunAction::EndOfRunAction",
                    "HODO_GATE3_PROCESS_PROVENANCE",FatalException,message);
      }
      // This is deliberately distinct from mu- Scintillation activation.
      // The predecessor physics contract removes the Scintillation process
      // from the optical-photon process manager itself; retaining it changed
      // the approved default-off replay.  Read the initialized process list,
      // and refuse to write a Stage-3 artifact if that cleanup regresses.
      auto* opticalPhoton=G4OpticalPhoton::Definition();
      auto* opticalPhotonProcessManager=
          opticalPhoton?opticalPhoton->GetProcessManager():nullptr;
      G4int opticalPhotonScintillationAbsent=-1;
      if(opticalPhotonProcessManager) {
        opticalPhotonScintillationAbsent=1;
        auto* processes=opticalPhotonProcessManager->GetProcessList();
        for(G4int i=0;
            i<opticalPhotonProcessManager->GetProcessListLength();++i) {
          const auto* process=(*processes)[i];
          if(process&&process->GetProcessName()=="Scintillation") {
            opticalPhotonScintillationAbsent=0;
            break;
          }
        }
      }
      if(opticalPhotonScintillationAbsent!=1)
        G4Exception("OpticalRunAction::EndOfRunAction",
                    "HODO_GATE3_OPTICAL_PHOTON_SCINTILLATION",FatalException,
                    "Initialized optical-photon ProcessManager is missing or "
                    "still contains Scintillation; predecessor cleanup was "
                    "not preserved.");

      constexpr auto id=kValidationGenerationProvenanceNtupleId;
      analysis->FillNtupleIColumn(id,0,2);
      analysis->FillNtupleIColumn(id,1,fDetector->GetGenerationProbeMode());
      analysis->FillNtupleIColumn(
          id,2,fDetector->GetRecordGenerationNtuples());
      analysis->FillNtupleSColumn(id,3,fDetector->GetGenerationProbeFile());
      analysis->FillNtupleSColumn(
          id,4,fDetector->GetGenerationProbeSha256());
      analysis->FillNtupleIColumn(
          id,5,gGenerationProbeLoadedRows.load(std::memory_order_relaxed));
      analysis->FillNtupleIColumn(
          id,6,fDetector->GetGenerationProbeRowOffset());
      analysis->FillNtupleDColumn(
          id,7,static_cast<G4double>(fDetector->GetCosmicSeed()));
      analysis->FillNtupleDColumn(
          id,8,static_cast<G4double>(fInitialRandomSeed0));
      analysis->FillNtupleDColumn(
          id,9,static_cast<G4double>(fInitialRandomSeed1));
      analysis->FillNtupleDColumn(id,10,fDetector->GetFibreMaxStep()/mm);
      analysis->FillNtupleDColumn(
          id,11,gGenerationObservedMaxCoreStep.load(
              std::memory_order_relaxed)/mm);
      analysis->FillNtupleIColumn(id,12,fDetector->GetGenerationOnlyKill());
      analysis->FillNtupleIColumn(id,13,scintillationActive);
      analysis->FillNtupleIColumn(id,14,cerenkovActive);
      analysis->FillNtupleIColumn(id,15,emSaturationAttached);
      analysis->FillNtupleIColumn(id,16,scintByParticleType);
      analysis->FillNtupleDColumn(id,17,h1.yield);
      analysis->FillNtupleDColumn(id,18,h2.yield);
      analysis->FillNtupleDColumn(id,19,h1.birks);
      analysis->FillNtupleDColumn(id,20,h2.birks);
      analysis->FillNtupleDColumn(id,21,h1.resolution);
      analysis->FillNtupleDColumn(id,22,h2.resolution);
      analysis->FillNtupleDColumn(id,23,physics->GetDefaultCutValue()/mm);
      analysis->FillNtupleIColumn(id,24,base->GetNumberOfEvent());
      analysis->FillNtupleIColumn(
          id,25,opticalPhotonScintillationAbsent);
      analysis->FillNtupleIColumn(id,26,muonMscActive);
      analysis->FillNtupleIColumn(
          id,27,muonCoulombScatteringActive);
      analysis->FillNtupleIColumn(id,28,muonIonisationActive);
      analysis->FillNtupleIColumn(id,29,muonBremsstrahlungActive);
      analysis->FillNtupleIColumn(id,30,muonPairProductionActive);
      analysis->FillNtupleIColumn(id,31,muonDecayActive);
      analysis->AddNtupleRow(id);
    }
    const G4bool writeOk = analysis->Write();
    const G4bool closeOk = analysis->CloseFile();
    if (!writeOk || !closeOk) {
      G4ExceptionDescription message;
      message << "V15 analysis finalization failed: write=" << writeOk
              << ", close=" << closeOk << ".";
      G4Exception("OpticalRunAction::EndOfRunAction",
                  "HODO_V21_ANALYSIS_FINALIZE", FatalException,
                  message);
    }
  }
 private:
  const DetectorConstruction* fDetector;
  G4int fAlignmentConfigurationNtupleId=-1;
  G4long fInitialRandomSeed0=0;
  G4long fInitialRandomSeed1=0;
};

class OpticalEventAction final : public G4UserEventAction {
 public:
  explicit OpticalEventAction(const DetectorConstruction* detector) : fDetector(detector) {}
  const DetectorConstruction* Detector() const { return fDetector; }
  void BeginOfEventAction(const G4Event* event) override {
    // G4RunManager builds user actions before a macro can set the Stage-2
    // controls.  Defer external electronics input loading until the first
    // event, when the detector configuration is final.  This also guarantees
    // that a geometry-only process never depends on DAQ input files.
    if(!fDetector->GetGate4OpticalAudit())PrepareElectronicsInputs();

    // Reset the transport RNG per event so optical branching in one event
    // cannot alter the primary/trigger physics of later events in a
    // same-cosmic-seed sensitivity comparison.
    std::uint64_t value=static_cast<std::uint64_t>(fDetector->GetCosmicSeed());
    const auto replayEventId =
        static_cast<std::uint64_t>(event->GetEventID())+
        static_cast<std::uint64_t>(
            std::max(0,fDetector->GetTransportSeedEventOffset()));
    value += 0x9e3779b97f4a7c15ULL*
             (replayEventId+1);
    value=(value^(value>>30))*0xbf58476d1ce4e5b9ULL;
    value=(value^(value>>27))*0x94d049bb133111ebULL;
    value^=value>>31;
    G4Random::setTheSeed(static_cast<long>((value & 0x7ffffffeULL)+1));
    fT1=fT2=fT3=fT4=0.;
    fDetected.fill(0);
    fFirstTime.fill(std::numeric_limits<G4double>::infinity());
    for(auto& times:fArrivalTimes)times.clear();
    fProduced = fProducedFiber = fProducedTrigger = fEscaped = fPdeRejected = 0;
    fProducedH1=fProducedH2=0;
    fProducedActivePlane.fill(0);
    fProducedFanoutPlane.fill(0);
    fProducedFarEndPlane.fill(0);
    fStageH1.fill(0);fStageH2.fill(0);fPdeRejectedPlane.fill(0);
    fPdeRejectedByPlane.fill(0);
    for(auto& plane:fStagePlane)plane.fill(0);
    for(auto& plane:fStageSeen)for(auto& stage:plane)stage.clear();
    fBulkAbsCore = fBulkAbsClad = fBulkAbsCoupling = fBulkAbsOther = 0;
    fLossStress=fLossCrack=fLossPacking=fLossAdhesive=0;
    fLossH1Cladding=fLossH2Cladding=0;
    fLossFanoutExit=fLossEndFace=fLossGreaseCoupling=0;
    fInitialEnergy = fAtTrigger1 = -1.;
    fTrigger1Time = std::numeric_limits<G4double>::infinity();
    fEdepConcrete = fEdepBlackout = fEdepTrigger = fEdepH1 = fEdepH2 = 0.;
    fEdepH1Core = fEdepH2Core = 0.;
    fPrimaryCoreTrack.fill(0.);
    fPrimaryActiveCoreTrack.fill(0.);
    fCoreEdepPlane.fill(0.);
    fProducedPlane.fill(0);
    fPrimaryCoreEdepPlane.fill(0.);
    fSecondaryCoreEdepPlane.fill(0.);
    fPrimaryProducedPlane.fill(0);
    fSecondaryProducedPlane.fill(0);
    for(auto& cores:fPrimaryCoreCopies)cores.clear();
    fPhotonAudit.clear();
    fFibreCrossings.clear();
    fLastCrossingByTrack.clear();
    fNextCrossingIndex.clear();
    fPendingPhotonBirth.clear();
    for(auto& surface:fPrimarySurface)surface=PrimarySurfaceState{};
    for(auto& crossing:fPrimaryTriggerCrossing)
      crossing=PrimaryTriggerCrossing{};
    for(auto& ledger:fGenerationPlaneLedger)
      ledger=GenerationPlaneLedger{};
    for(auto& ledger:fGenerationChannelLedger)
      ledger=GenerationPlaneLedger{};
    fGenerationSource=GenerationProbeRow{};
    fGenerationSourceValid=false;
    if(fDetector->GetGenerationProbeMode()) {
      const auto* vertex=event->GetPrimaryVertex();
      const auto* primary=vertex?vertex->GetPrimary():nullptr;
      const auto* info=primary?dynamic_cast<const GenerationProbePrimaryInfo*>(
          primary->GetUserInformation()):nullptr;
      if(!info)
        G4Exception("OpticalEventAction::BeginOfEventAction",
                    "HODO_GATE3_SOURCE_INFO",FatalException,
                    "Stage-3 primary lacks its frozen source-row metadata.");
      fGenerationSource=info->value;
      fGenerationSourceValid=true;
    }
    if(fDetector->GetGate6ElectronicsProbeMode()) {
      const auto channel=std::clamp(
          fDetector->GetGate6InjectedChannelForEvent(event->GetEventID()),0,127);
      const auto count=fDetector->GetGate6InjectedPe();
      const auto start=fDetector->GetGate6InjectedStartTime();
      const auto spacing=fDetector->GetGate6InjectedSpacing();
      fDetected[channel]=count;
      if(count>0)fFirstTime[channel]=start;
      fArrivalTimes[channel].reserve(static_cast<std::size_t>(count));
      for(G4int index=0;index<count;++index)
        fArrivalTimes[channel].push_back(start+index*spacing);
    }
  }
  void PrepareElectronicsInputs() {
    if(fElectronicsInputsPrepared)return;
    if(fDetector->GetGeometryProbeMode()||
       fDetector->GetGenerationProbeMode()) {
      // Geometry probes never create optical avalanches and do not judge ADC.
      // Keep the legacy event/channel schemas writable with deterministic zero
      // electronics state, without depending on external DAQ input files.
      fGain.fill(0.);
      fNoiseRms.fill(0.);
      fParkThreshold.fill(120.);
      G4cout << "[validation probe] electronics input loading bypassed"
             << G4endl;
    } else {
      // Preserve the predecessor input-loading order exactly when Stage-2 is
      // off; loading occurs later, but consumes no random numbers.
      LoadGains();
      LoadNoise();
      LoadParkThresholds();
      LoadSinglePeTemplate();
    }
    fElectronicsInputsPrepared=true;
  }
  void AddEnergy(const G4String& volume, G4double edep) {
    if (volume == "Trigger1") fT1 += edep;
    else if (volume == "Trigger2") fT2 += edep;
    else if (volume == "Trigger3") fT3 += edep;
    else if (volume == "Trigger4") fT4 += edep;
  }
  void AddProducedPhoton(const G4String& material, G4int plane=-1,
                         G4int parentTrackId=-1,
                         const G4ThreeVector& position={}) {
    ++fProduced;
    if (material == "Trigger_vinyltoluene_scintillator") {
      ++fProducedTrigger;
    } else if (material == "Kuraray_SCSF78_core") {
      ++fProducedFiber;
      ++fProducedH1;
      if(plane>=0&&plane<2) {
        ++fProducedPlane[plane];
        const auto nominalPosition=
            fDetector->PlaneWorldToNominalPosition(plane,position);
        const auto longitudinal=
            plane==0?nominalPosition.x():nominalPosition.y();
        const auto bounds=fDetector->GetActiveLongitudinalBounds(plane);
        if(longitudinal>=bounds[0]&&longitudinal<=bounds[1])
          ++fProducedActivePlane[plane];
        else if(longitudinal<bounds[0])++fProducedFarEndPlane[plane];
        else ++fProducedFanoutPlane[plane];
        if(parentTrackId==1)++fPrimaryProducedPlane[plane];
        else ++fSecondaryProducedPlane[plane];
      }
    } else if (material == "BCF20_core") {
      ++fProducedFiber;
      ++fProducedH2;
      if(plane>=2&&plane<4) {
        ++fProducedPlane[plane];
        const auto nominalPosition=
            fDetector->PlaneWorldToNominalPosition(plane,position);
        const auto longitudinal=
            plane==2?nominalPosition.x():nominalPosition.y();
        const auto bounds=fDetector->GetActiveLongitudinalBounds(plane);
        if(longitudinal>=bounds[0]&&longitudinal<=bounds[1])
          ++fProducedActivePlane[plane];
        else if(longitudinal<bounds[0])++fProducedFarEndPlane[plane];
        else ++fProducedFanoutPlane[plane];
        if(parentTrackId==1)++fPrimaryProducedPlane[plane];
        else ++fSecondaryProducedPlane[plane];
      }
    }
  }
  void AddEscapedPhoton(const G4Track* track) {
    ++fEscaped;
    MarkPhotonTerminal(track,TerminalReason::world_escape,
                       AuditRegion::world);
  }
  void SetInitialEnergy(G4double value) { if (fInitialEnergy < 0.) fInitialEnergy = value; }
  void SetAtTrigger1(G4double value, G4double time) {
    if (fAtTrigger1 < 0.) { fAtTrigger1 = value; fTrigger1Time = time; }
  }
  void AddMaterialEnergy(const G4String& volume, G4double edep) {
    if (volume == "ConcreteFloor") fEdepConcrete += edep;
    else if (volume == "BlackoutSheet") fEdepBlackout += edep;
    else if (volume == "Trigger1" || volume == "Trigger2" ||
             volume == "Trigger3" || volume == "Trigger4")
      fEdepTrigger += edep;
    else if (volume.rfind("H1", 0) == 0) fEdepH1 += edep;
    else if (volume.rfind("H2", 0) == 0) fEdepH2 += edep;
  }
  void AddCoreEnergy(const G4String& material, G4double edep) {
    if (material == "Kuraray_SCSF78_core") fEdepH1Core += edep;
    else if (material == "BCF20_core") fEdepH2Core += edep;
  }
  void AddCoreTraversal(G4int plane,G4int copy,G4double length,
                        G4double activeLength,G4double edep,G4int trackId) {
    if(plane<0||plane>=4)return;
    fCoreEdepPlane[plane]+=edep;
    if(trackId==1) {
      fPrimaryCoreEdepPlane[plane]+=edep;
      fPrimaryCoreTrack[plane]+=length;
      fPrimaryActiveCoreTrack[plane]+=std::clamp(activeLength,0.,length);
      fPrimaryCoreCopies[plane].insert(copy);
    } else fSecondaryCoreEdepPlane[plane]+=edep;
  }
  void ObserveChargedStep(const G4Step* step) {
    if((!fDetector->GetRecordSpatialClosureNtuples()&&
        !fDetector->GetRecordGeometryControlNtuples())||!step)return;
    const auto* track=step->GetTrack();
    if(!track||track->GetDefinition()==G4OpticalPhoton::Definition())return;

    if(fDetector->GetRecordGeometryControlNtuples()&&
       track->GetTrackID()==1) {
      const auto* preVolume=step->GetPreStepPoint()->GetPhysicalVolume();
      const auto& name=preVolume?preVolume->GetName():G4String{};
      const G4int trigger=name=="Trigger1"?0:name=="Trigger2"?1:
                          name=="Trigger3"?2:name=="Trigger4"?3:-1;
      if(trigger>=0) {
        auto& crossing=fPrimaryTriggerCrossing[trigger];
        const auto* pre=step->GetPreStepPoint();
        const auto* post=step->GetPostStepPoint();
        if(!crossing.reached) {
          crossing.reached=true;
          crossing.track=track->GetTrackID();
          crossing.parent=track->GetParentID();
          crossing.pdg=track->GetDefinition()->GetPDGEncoding();
          crossing.entryWorld=pre->GetPosition();
          crossing.entryDirection=pre->GetMomentumDirection();
          crossing.entryTime=pre->GetGlobalTime();
        }
        crossing.exitWorld=post->GetPosition();
        crossing.exitDirection=post->GetMomentumDirection();
        crossing.exitTime=post->GetGlobalTime();
        crossing.path+=step->GetStepLength();
        crossing.edep+=step->GetTotalEnergyDeposit();
        ++crossing.steps;
      }
    }
    if(!fDetector->GetRecordSpatialClosureNtuples())return;

    if(track->GetTrackID()==1) {
      const auto* pre=step->GetPreStepPoint();
      const auto* post=step->GetPostStepPoint();
      const G4double a=pre->GetPosition().z();
      const G4double b=post->GetPosition().z();
      const G4double delta=b-a;
      for(G4int surface=0;surface<6;++surface) {
        auto& state=fPrimarySurface[surface];
        if(state.reached)continue;
        const G4double z=fDetector->GetValidationSurfaceZ(surface);
        const G4bool bracketed=(a-z)*(b-z)<=0.&&
            (std::abs(delta)>1.e-15*mm||std::abs(a-z)<=1.e-12*mm);
        if(!bracketed)continue;
        const G4double fraction=std::abs(delta)>1.e-15*mm
            ? std::clamp((z-a)/delta,0.,1.) : 0.;
        state.reached=true;
        state.track=track->GetTrackID();
        state.position=pre->GetPosition()+
            fraction*(post->GetPosition()-pre->GetPosition());
        state.direction=pre->GetMomentumDirection();
        state.kineticEnergy=pre->GetKineticEnergy()+
            fraction*(post->GetKineticEnergy()-pre->GetKineticEnergy());
        state.time=pre->GetGlobalTime()+
            fraction*(post->GetGlobalTime()-pre->GetGlobalTime());
      }
    }

    const auto* material=step->GetPreStepPoint()->GetMaterial();
    const auto materialName=material?material->GetName():G4String{};
    FibreVolumeClass volumeClass=FibreVolumeClass::unknown;
    if(materialName=="Kuraray_SCSF78_core"||materialName=="BCF20_core")
      volumeClass=FibreVolumeClass::core;
    else if(materialName=="PMMA_inner_cladding")
      volumeClass=FibreVolumeClass::inner_cladding;
    else if(materialName=="Fluorinated_outer_cladding")
      volumeClass=FibreVolumeClass::outer_cladding;
    if(volumeClass==FibreVolumeClass::unknown)return;

    G4int plane=-1,globalChannel=-1;
    ResolvePlaneChannel(step->GetPreStepPoint()->GetTouchableHandle(),
                        plane,globalChannel);
    const G4int localChannel=plane>=0&&plane<4&&globalChannel>=0
        ? globalChannel-plane*32 : -1;
    const G4int physicalCopy=plane>=0&&plane<4&&localChannel>=0
        ? (plane%2)*32+localChannel : -1;
    const auto local=[this,plane,localChannel](const G4ThreeVector& world) {
      if(plane<0||plane>=4||localChannel<0||localChannel>=32)return world;
      const auto nominal=fDetector->PlaneWorldToNominalPosition(plane,world);
      const G4bool xAxis=plane==0||plane==2;
      const G4double centreTransverse=
          fDetector->GetPlaneSourceTransverse(plane,localChannel);
      const G4double centreZ=fDetector->GetPlaneSourceZ(plane,localChannel);
      return xAxis
          ? G4ThreeVector(nominal.x(),nominal.y()-centreTransverse,
                          nominal.z()-centreZ)
          : G4ThreeVector(nominal.y(),nominal.x()-centreTransverse,
                          nominal.z()-centreZ);
    };
    const auto* pre=step->GetPreStepPoint();
    const auto* post=step->GetPostStepPoint();
    const G4int trackId=track->GetTrackID();
    const G4int stepNumber=track->GetCurrentStepNumber();
    const G4double length=std::max(0.,step->GetStepLength());
    G4double activeLength=0.,farEndLength=0.;
    if(volumeClass==FibreVolumeClass::core&&plane>=0&&plane<4) {
      const auto before=fDetector->PlaneWorldToNominalPosition(
          plane,pre->GetPosition());
      const auto after=fDetector->PlaneWorldToNominalPosition(
          plane,post->GetPosition());
      const G4double first=plane==0||plane==2?before.x():before.y();
      const G4double last=plane==0||plane==2?after.x():after.y();
      const auto bounds=fDetector->GetActiveLongitudinalBounds(plane);
      const G4double low=bounds[0],high=bounds[1];
      const G4double deltaLong=last-first;
      G4double fraction=0.;
      if(std::abs(deltaLong)<=1.e-12*mm)
        fraction=first>=low&&first<=high?1.:0.;
      else {
        G4double enter=(low-first)/deltaLong;
        G4double leave=(high-first)/deltaLong;
        if(enter>leave)std::swap(enter,leave);
        fraction=std::max(0.,std::min(1.,leave)-std::max(0.,enter));
      }
      activeLength=length*fraction;
      if(fDetector->GetActivePlaneOverlap()||fDetector->GetFibreFarEndExtension()>0.) {
        G4double farFraction=0.;
        if(std::abs(deltaLong)<=1.e-12*mm)farFraction=first<low?1.:0.;
        else {
          const auto split=std::clamp((low-first)/deltaLong,0.,1.);
          farFraction=deltaLong>0.?split:1.-split;
        }
        farEndLength=length*farFraction;
      }
    }

    std::size_t index=std::numeric_limits<std::size_t>::max();
    const auto found=fLastCrossingByTrack.find(trackId);
    if(found!=fLastCrossingByTrack.end()) {
      auto& previous=fFibreCrossings[found->second];
      const G4double continuity=(previous.exitWorld-pre->GetPosition()).mag();
      if(previous.lastStep+1==stepNumber&&previous.plane==plane&&
         previous.physicalCopy==physicalCopy&&
         previous.volumeClass==volumeClass&&continuity<=1.e-9*mm)
        index=found->second;
    }
    if(index==std::numeric_limits<std::size_t>::max()) {
      FibreCrossing crossing;
      crossing.track=trackId;
      crossing.parent=track->GetParentID();
      crossing.pdg=track->GetDefinition()->GetPDGEncoding();
      crossing.crossingIndex=fNextCrossingIndex[trackId]++;
      crossing.plane=plane;
      crossing.physicalCopy=physicalCopy;
      crossing.localChannel=localChannel;
      crossing.globalChannel=globalChannel;
      crossing.volumeClass=volumeClass;
      crossing.firstStep=stepNumber;
      crossing.lastStep=stepNumber;
      crossing.entryWorld=pre->GetPosition();
      crossing.exitWorld=post->GetPosition();
      crossing.entryLocal=local(pre->GetPosition());
      crossing.exitLocal=local(post->GetPosition());
      crossing.entryDirection=pre->GetMomentumDirection();
      crossing.exitDirection=post->GetMomentumDirection();
      crossing.entryEnergy=pre->GetKineticEnergy();
      crossing.exitEnergy=post->GetKineticEnergy();
      crossing.entryTime=pre->GetGlobalTime();
      crossing.exitTime=post->GetGlobalTime();
      crossing.path=length;
      crossing.corePath=volumeClass==FibreVolumeClass::core?length:0.;
      crossing.activeCorePath=activeLength;
      crossing.farEndCorePath=farEndLength;
      crossing.edep=std::max(0.,step->GetTotalEnergyDeposit());
      fFibreCrossings.push_back(crossing);
      index=fFibreCrossings.size()-1;
      fLastCrossingByTrack[trackId]=index;
    } else {
      auto& crossing=fFibreCrossings[index];
      crossing.lastStep=stepNumber;
      crossing.exitWorld=post->GetPosition();
      crossing.exitLocal=local(post->GetPosition());
      crossing.exitDirection=post->GetMomentumDirection();
      crossing.exitEnergy=post->GetKineticEnergy();
      crossing.exitTime=post->GetGlobalTime();
      crossing.path+=length;
      if(volumeClass==FibreVolumeClass::core)crossing.corePath+=length;
      crossing.activeCorePath+=activeLength;
      crossing.farEndCorePath+=farEndLength;
      crossing.edep+=std::max(0.,step->GetTotalEnergyDeposit());
    }
  }
  void RecordGenerationStep(const G4Step* step,G4int plane,
                            G4int globalChannel,
                            const GenerationSecondaryCounts& photons) {
    if(!fDetector->GetRecordGenerationNtuples()||!step)return;
    if(!fGenerationSourceValid)
      G4Exception("OpticalEventAction::RecordGenerationStep",
                  "HODO_GATE3_SOURCE_INFO",FatalException,
                  "A Stage-3 step was observed without source metadata.");
    const auto* track=step->GetTrack();
    const auto* material=step->GetPreStepPoint()->GetMaterial();
    if(!track||!material||
       track->GetDefinition()==G4OpticalPhoton::Definition()||
       track->GetDefinition()->GetPDGCharge()==0.)return;
    const auto& materialName=material->GetName();
    FibreVolumeClass volumeClass=FibreVolumeClass::unknown;
    G4int materialCode=0;
    if(materialName=="Kuraray_SCSF78_core") {
      volumeClass=FibreVolumeClass::core;
      materialCode=1;
    } else if(materialName=="BCF20_core") {
      volumeClass=FibreVolumeClass::core;
      materialCode=2;
    } else if(materialName=="PMMA_inner_cladding") {
      volumeClass=FibreVolumeClass::inner_cladding;
      materialCode=3;
    } else if(materialName=="Fluorinated_outer_cladding") {
      volumeClass=FibreVolumeClass::outer_cladding;
      materialCode=4;
    } else return;

    if(plane<0||plane>=4||globalChannel<plane*32||
       globalChannel>=(plane+1)*32) {
      G4ExceptionDescription message;
      message << "Charged fibre step cannot be mapped to a Stage-3 plane/"
              << "channel: material=" << materialName
              << " plane=" << plane
              << " global_channel=" << globalChannel << ".";
      G4Exception("OpticalEventAction::RecordGenerationStep",
                  "HODO_GATE3_STEP_ORPHAN",FatalException,message);
    }
    const G4int localChannel=globalChannel-plane*32;
    const G4int physicalCopy=(plane%2)*32+localChannel;
    const auto* pre=step->GetPreStepPoint();
    const auto* post=step->GetPostStepPoint();
    const G4double length=std::max(0.,step->GetStepLength());
    const auto nominalPre=fDetector->PlaneWorldToNominalPosition(
        plane,pre->GetPosition());
    const auto nominalPost=fDetector->PlaneWorldToNominalPosition(
        plane,post->GetPosition());
    const G4bool xAxis=plane==0||plane==2;
    const G4double first=xAxis?nominalPre.x():nominalPre.y();
    const G4double last=xAxis?nominalPost.x():nominalPost.y();
    const auto bounds=fDetector->GetActiveLongitudinalBounds(plane);
    const G4double low=bounds[0],high=bounds[1];
    G4double activeFraction=0.;
    const G4double delta=last-first;
    if(std::abs(delta)<=1.e-12*mm)
      activeFraction=first>=low&&first<=high?1.:0.;
    else {
      G4double enter=(low-first)/delta;
      G4double leave=(high-first)/delta;
      if(enter>leave)std::swap(enter,leave);
      activeFraction=std::max(0.,std::min(1.,leave)-std::max(0.,enter));
    }
    const G4double activePath=length*activeFraction;
    const G4double edep=std::max(0.,step->GetTotalEnergyDeposit());
    const G4double niel=std::max(0.,step->GetNonIonizingEnergyDeposit());
    const auto* mpt=material->GetMaterialPropertiesTable();
    const G4double yield=mpt&&mpt->ConstPropertyExists("SCINTILLATIONYIELD")
        ? mpt->GetConstProperty("SCINTILLATIONYIELD")*MeV : 0.;
    const auto* ionisation=material->GetIonisation();
    const G4double birks=ionisation
        ? ionisation->GetBirksConstant()/(mm/MeV) : 0.;

    if(volumeClass==FibreVolumeClass::core)
      ObserveGenerationCoreStepLength(length);
    const auto accumulate=[&](GenerationPlaneLedger& ledger) {
      ++ledger.chargedFibreSteps;
      ledger.directOpticalTotal+=photons.Total();
      if(track->GetTrackID()==1) {
        ledger.primaryScintillation+=photons.scintillation;
        ledger.primaryCerenkov+=photons.cerenkov;
        ledger.primaryOther+=photons.other;
      } else {
        ledger.secondaryScintillation+=photons.scintillation;
        ledger.secondaryCerenkov+=photons.cerenkov;
        ledger.secondaryOther+=photons.other;
      }
      if(volumeClass==FibreVolumeClass::core) {
        ledger.allCoreEdep+=edep;
        ledger.allCoreNiel+=niel;
        ledger.scintillationActive+=photons.scintillationActive;
        ledger.scintillationFanout+=photons.scintillationFanout;
        ledger.cerenkovActive+=photons.cerenkovActive;
        ledger.cerenkovFanout+=photons.cerenkovFanout;
        ledger.otherActive+=photons.otherActive;
        ledger.otherFanout+=photons.otherFanout;
        if(track->GetTrackID()==1) {
          ++ledger.primaryCoreSteps;
          ledger.primaryCorePath+=length;
          ledger.primaryActiveCorePath+=activePath;
          ledger.primaryCoreEdep+=edep;
          ledger.primaryCoreNiel+=niel;
          if(niel==0.) {
            const G4double birksInternal=ionisation
                ? ionisation->GetBirksConstant() : 0.;
            ledger.primaryCoreVisibleEdep+=length>0.
                ? edep/(1.+birksInternal*edep/length) : 0.;
          } else ++ledger.primaryCoreNonzeroNielSteps;
        }
      } else {
        ledger.allCladdingEdep+=edep;
        ledger.allCladdingNiel+=niel;
        ledger.claddingScintillation+=photons.scintillation;
        ledger.claddingCerenkov+=photons.cerenkov;
        ledger.claddingOther+=photons.other;
        if(track->GetTrackID()==1) {
          ledger.primaryCladdingPath+=length;
          ledger.primaryActiveCladdingPath+=activePath;
          ledger.primaryCladdingEdep+=edep;
          ledger.primaryCladdingNiel+=niel;
        }
      }
    };
    accumulate(fGenerationPlaneLedger[plane]);
    accumulate(fGenerationChannelLedger[globalChannel]);

    constexpr auto id=kValidationGenerationStepNtupleId;
    const std::array<G4int,20> integers{{
        G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID(),
        fGenerationSource.event,fGenerationSource.cohort,
        fGenerationSource.targetPlane,fGenerationSource.targetChannel,
        fGenerationSource.expectedVolumeClass,
        fGenerationSource.expectedRegion,fGenerationSource.chordBin,
        fGenerationSource.repeatIndex,fGenerationSource.caseIndex,
        track->GetTrackID(),track->GetParentID(),
        track->GetDefinition()->GetPDGEncoding(),
        track->GetCurrentStepNumber(),plane,physicalCopy,localChannel,
        globalChannel,static_cast<G4int>(volumeClass),materialCode}};
    auto* analysis=G4AnalysisManager::Instance();
    for(std::size_t column=0;column<integers.size();++column)
      analysis->FillNtupleIColumn(id,column,integers[column]);
    const std::array<G4double,14> doubles{{
        length/mm,activePath/mm,edep/MeV,niel/MeV,
        pre->GetPosition().x()/mm,pre->GetPosition().y()/mm,
        pre->GetPosition().z()/mm,post->GetPosition().x()/mm,
        post->GetPosition().y()/mm,post->GetPosition().z()/mm,yield,birks,
        photons.activeMaxAbsLongitudinal<0.
            ? -1. : photons.activeMaxAbsLongitudinal/mm,
        photons.fanoutMinAbsLongitudinal<0.
            ? -1. : photons.fanoutMinAbsLongitudinal/mm}};
    for(std::size_t column=0;column<doubles.size();++column)
      analysis->FillNtupleDColumn(id,20+column,doubles[column]);
    const std::array<G4int,9> secondaryCounts{{
        photons.scintillation,photons.cerenkov,photons.other,
        photons.scintillationActive,photons.scintillationFanout,
        photons.cerenkovActive,photons.cerenkovFanout,
        photons.otherActive,photons.otherFanout}};
    for(std::size_t column=0;column<secondaryCounts.size();++column)
      analysis->FillNtupleIColumn(id,34+column,secondaryCounts[column]);
    analysis->AddNtupleRow(id);
  }
  void ApplyBirthContext(PhotonAudit& audit,
                         const PhotonBirthContext& context) {
    audit.parentStep=context.parentStep;
    const G4double scale=std::max(1.*mm,context.position.mag());
    const G4double tolerance=
        32.*std::numeric_limits<G4double>::epsilon()*scale;
    audit.birthContextMatched=
        audit.parent==context.parentTrack&&audit.plane==context.plane&&
        audit.originChannel==context.globalChannel&&
        (audit.origin-context.position).mag()<=tolerance;
  }
  void RegisterPhotonBirth(const G4Track* secondary,const G4Step* parentStep,
                           G4int plane,G4int globalChannel) {
    if(!fDetector->GetRecordSpatialClosureNtuples()||!secondary||!parentStep||
       secondary->GetDefinition()!=G4OpticalPhoton::Definition())return;
    const auto* parent=parentStep->GetTrack();
    const auto* material=parentStep->GetPreStepPoint()->GetMaterial();
    if(!parent||!material||
       (material->GetName()!="Kuraray_SCSF78_core"&&
        material->GetName()!="BCF20_core")||
       plane<0||plane>=4||globalChannel<plane*32||
       globalChannel>=(plane+1)*32)return;
    PhotonBirthContext context;
    context.parentTrack=parent->GetTrackID();
    context.parentStep=parent->GetCurrentStepNumber();
    context.plane=plane;
    context.globalChannel=globalChannel;
    context.position=secondary->GetPosition();
    const auto audit=fPhotonAudit.find(secondary->GetTrackID());
    if(audit!=fPhotonAudit.end()) {
      ApplyBirthContext(audit->second,context);
      return;
    }
    fPendingPhotonBirth[secondary]=context;
  }
  void BeginOpticalTrack(const G4Track* track) {
    if(!fDetector->GetRecordPhotonAuditNtuple() || !track ||
       track->GetDefinition()!=G4OpticalPhoton::Definition() ||
       fPhotonAudit.count(track->GetTrackID())) return;
    const auto* volume=track->GetVolume();
    const auto* vertexLogical=track->GetLogicalVolumeAtVertex();
    const auto* material=vertexLogical
        ? vertexLogical->GetMaterial()
        : volume&&volume->GetLogicalVolume()
            ? volume->GetLogicalVolume()->GetMaterial() : nullptr;
    if(!material ||
       (material->GetName()!="Kuraray_SCSF78_core" &&
        material->GetName()!="BCF20_core")) return;
    PhotonAudit audit;
    audit.track=track->GetTrackID();
    audit.parent=track->GetParentID();
    LocatePlaneChannel(track->GetTouchableHandle(),
                       audit.plane,audit.originChannel);
    if(audit.plane>=0&&audit.plane<4&&audit.originChannel>=0) {
      audit.originLocalChannel=audit.originChannel-audit.plane*32;
      audit.originPhysicalCopy=(audit.plane%2)*32+
                               audit.originLocalChannel;
    }
    const auto* creator=track->GetCreatorProcess();
    if(!creator) audit.creatorProcess=0;
    else if(creator->GetProcessName()=="Scintillation")audit.creatorProcess=1;
    else if(creator->GetProcessName()=="Cerenkov")audit.creatorProcess=2;
    else if(creator->GetProcessName()=="OpWLS")audit.creatorProcess=3;
    audit.origin=track->GetPosition();
    audit.initialDirection=track->GetMomentumDirection();
    audit.originTime=track->GetGlobalTime();
    audit.terminalPosition=audit.origin;
    audit.terminalTime=audit.originTime;
    if(track->GetTotalEnergy()>0.)
      audit.wavelength=h_Planck*c_light/track->GetTotalEnergy();
    audit.originRegion=RegionAt(
        volume,material,audit.origin,audit.plane);
    audit.lastRegion=audit.originRegion;
    audit.wasInCore=IsCore(audit.originRegion);
    auto inserted=fPhotonAudit.emplace(audit.track,std::move(audit));
    const auto pending=fPendingPhotonBirth.find(track);
    if(pending!=fPendingPhotonBirth.end()) {
      ApplyBirthContext(inserted.first->second,pending->second);
      fPendingPhotonBirth.erase(pending);
    }
  }
  void ObserveOpticalStep(const G4Step* step,
                          const G4OpBoundaryProcess* boundary) {
    if(!fDetector->GetRecordPhotonAuditNtuple() || !step)return;
    const auto* track=step->GetTrack();
    BeginOpticalTrack(track);
    auto found=fPhotonAudit.find(track->GetTrackID());
    if(found==fPhotonAudit.end())return;
    auto& audit=found->second;
    const auto* prePoint=step->GetPreStepPoint();
    const auto* postPoint=step->GetPostStepPoint();
    const auto* prePV=prePoint->GetPhysicalVolume();
    const auto* postPV=postPoint->GetPhysicalVolume();
    const auto preRegion=RegionAt(
        prePV,prePoint->GetMaterial(),prePoint->GetPosition(),audit.plane);
    const auto candidatePostRegion=RegionAt(
        postPV,postPoint->GetMaterial(),postPoint->GetPosition(),audit.plane);
    const G4bool atBoundary=postPoint->GetStepStatus()==fGeomBoundary;
    const auto status=atBoundary&&boundary
        ? boundary->GetStatus() : Undefined;
    const G4bool reflected=atBoundary&&IsReflection(status);
    const auto actualPostRegion=reflected ? preRegion : candidatePostRegion;
    const auto length=std::max(0.,step->GetStepLength());
    AddAuditPath(audit,preRegion,length);
    audit.totalTrackLength+=length;
    ++audit.steps;
    if(fDetector->GetGate4OpticalAudit()) {
      const G4bool enteredCladding=
          IsCore(preRegion)&&IsCladding(actualPostRegion);
      const G4bool traversedCladding=IsCladding(preRegion)&&length>0.;
      // At a reflected daughter/parent boundary the navigator may briefly
      // report the parent cladding as the pre-volume of a zero-length recovery
      // step.  Gate 4 needs physical cladding travel, not that ownership-only
      // transition, for its analytic core-only trapping denominator.
      if(enteredCladding||traversedCladding)audit.everCladding=true;
      if(enteredCladding) {
        ++audit.coreToClad;
        audit.wasInCladding=true;
      }
      if(IsCladding(preRegion)&&IsCore(actualPostRegion)&&
         (audit.wasInCladding||traversedCladding)) {
        ++audit.cladToCore;
        audit.reenteredCore=true;
      }
      audit.wasInCore=audit.wasInCore||IsCore(actualPostRegion);
      audit.wasInCladding=audit.wasInCladding||enteredCladding||
                           traversedCladding;
    } else {
      // Preserve the approved predecessor bookkeeping bit-for-bit whenever
      // the Gate-4 optical audit is disabled.  These legacy counters include
      // navigator ownership transitions and are retained for lineage only.
      if(IsCladding(preRegion)||IsCladding(actualPostRegion))
        audit.everCladding=true;
      if(IsCore(preRegion)&&IsCladding(actualPostRegion)) {
        ++audit.coreToClad;
        audit.wasInCladding=true;
      }
      if(IsCladding(preRegion)&&IsCore(actualPostRegion)) {
        ++audit.cladToCore;
        audit.reenteredCore=true;
      }
      audit.wasInCore=audit.wasInCore||IsCore(actualPostRegion);
      audit.wasInCladding=audit.wasInCladding||IsCladding(actualPostRegion);
    }
    audit.lastRegion=actualPostRegion;
    audit.terminalPosition=postPoint->GetPosition();
    audit.terminalTime=track->GetGlobalTime();

    if(atBoundary) {
      if(status==TotalInternalReflection)++audit.tirCount;
      else if(status==FresnelReflection)++audit.fresnelReflectionCount;
      else if(reflected)++audit.otherReflectionCount;
    }

    const G4bool xAxis=audit.plane==0||audit.plane==2;
    const auto nominalBefore=fDetector->PlaneWorldToNominalPosition(
        audit.plane,prePoint->GetPosition());
    const auto nominalAfter=fDetector->PlaneWorldToNominalPosition(
        audit.plane,postPoint->GetPosition());
    const auto before=xAxis?nominalBefore.x():nominalBefore.y();
    const auto after=xAxis?nominalAfter.x():nominalAfter.y();
    const auto activeExit=fDetector->GetActiveLongitudinalBounds(audit.plane)[1];
    if(audit.plane>=0&&before<activeExit&&after>=activeExit)
      audit.stageMask|=stage_active_exit;

    const auto preName=prePV?prePV->GetName():G4String{};
    const auto postName=postPV?postPV->GetName():G4String{};
    const G4bool greaseTarget=postName=="H1_Grease"||
                              postName=="H2_Grease";
    // Gate 4 needs the actual result of exactly one side-wall encounter per
    // photon.  This observer is deliberately stdout-only: it draws no random
    // numbers, changes no track, and leaves the established ROOT schema
    // untouched.  Requiring the candidate post region to be cladding excludes
    // the readout/far ends and makes the trace specific to the core/cladding
    // optical-law fixture.
    if(fDetector->GetGate4BoundaryTrace()&&!audit.gate4BoundaryTraceRecorded&&
       atBoundary&&IsCore(preRegion)&&IsCladding(candidatePostRegion)) {
      const auto preDirection=prePoint->GetMomentumDirection();
      const auto postDirection=postPoint->GetMomentumDirection();
      const auto* currentEvent=
          G4RunManager::GetRunManager()->GetCurrentEvent();
      std::ostringstream trace;
      trace << std::setprecision(17)
            << "GATE4_BOUNDARY"
            << " event=" << (currentEvent?currentEvent->GetEventID():-1)
            << " track=" << track->GetTrackID()
            << " plane=" << audit.plane
            << " channel=" << audit.originChannel
            << " status=" << static_cast<G4int>(status)
            << " tir=" << (status==TotalInternalReflection)
            << " fresnel_reflection=" << (status==FresnelReflection)
            << " reflected=" << reflected
            << " transmitted="
            << (!reflected&&status!=Absorption&&status!=Detection&&
                status!=NoRINDEX)
            << " origin_x_mm=" << audit.origin.x()/mm
            << " origin_y_mm=" << audit.origin.y()/mm
            << " origin_z_mm=" << audit.origin.z()/mm
            << " boundary_x_mm=" << postPoint->GetPosition().x()/mm
            << " boundary_y_mm=" << postPoint->GetPosition().y()/mm
            << " boundary_z_mm=" << postPoint->GetPosition().z()/mm
            << " pre_dx=" << preDirection.x()
            << " pre_dy=" << preDirection.y()
            << " pre_dz=" << preDirection.z()
            << " post_dx=" << postDirection.x()
            << " post_dy=" << postDirection.y()
            << " post_dz=" << postDirection.z();
      G4cout << trace.str() << G4endl;
      audit.gate4BoundaryTraceRecorded=true;
    }
    const G4bool crossed=atBoundary&&!reflected&&status!=Absorption&&
                         status!=Detection&&status!=NoRINDEX;
    if(atBoundary&&greaseTarget&&IsCore(preRegion)) {
      ++audit.coreGreaseAttempts;
      if(crossed)audit.stageMask|=stage_core_to_grease;
    }
    if(atBoundary&&greaseTarget&&IsCladding(preRegion)) {
      ++audit.cladGreaseAttempts;
      if(crossed)audit.stageMask|=stage_clad_to_grease;
    }
    if(atBoundary&&(preName=="H1_Grease"||preName=="H2_Grease")&&
       postName=="SiPMWindow") {
      ++audit.greaseWindowAttempts;
      if(crossed)audit.stageMask|=stage_window_entry;
    }
    if(atBoundary&&postName=="SiPM") {
      ++audit.sipmAttempts;
      audit.stageMask|=stage_sipm_boundary;
      if(reflected)++audit.sipmReflections;
    }
  }
  void MarkPhotonTerminal(const G4Track* track,TerminalReason reason,
                          AuditRegion region=AuditRegion::unknown) {
    if(!fDetector->GetRecordPhotonAuditNtuple()||!track)return;
    auto found=fPhotonAudit.find(track->GetTrackID());
    if(found==fPhotonAudit.end())return;
    auto& audit=found->second;
    if(audit.terminal!=TerminalReason::unclassified) {
      ++audit.duplicateTerminalMarks;
      return;
    }
    audit.terminal=reason;
    audit.terminalRegion=region==AuditRegion::unknown
        ? audit.lastRegion : region;
    audit.terminalPosition=track->GetPosition();
    audit.terminalTime=track->GetGlobalTime();
  }
  void FinalizeOpticalTrack(const G4Track* track) {
    if(!fDetector->GetRecordPhotonAuditNtuple()||!track)return;
    auto found=fPhotonAudit.find(track->GetTrackID());
    if(found==fPhotonAudit.end())return;
    auto& audit=found->second;
    if(audit.terminal==TerminalReason::unclassified) {
      audit.terminal=audit.lastRegion==AuditRegion::world
          ? TerminalReason::world_escape
          : TerminalReason::other_terminal;
      audit.terminalRegion=audit.lastRegion;
      audit.terminalPosition=track->GetPosition();
      audit.terminalTime=track->GetGlobalTime();
    }
  }
  void AddDetectedPhoton(G4int channel, G4double time, G4double wavelength,
                         const G4Track* track) {
    if (channel < 0 || channel >= 128) return;
    ++fDetected[channel];
    fFirstTime[channel] = std::min(fFirstTime[channel], time);
    fArrivalTimes[channel].push_back(time);
    const auto found=fPhotonAudit.find(track?track->GetTrackID():-1);
    if(found!=fPhotonAudit.end())found->second.terminalGlobalChannel=channel;
    RecordSipmPhoton(channel, 1, wavelength, time);
    MarkPhotonTerminal(track,TerminalReason::sipm_detected,
                       AuditRegion::sipm);
  }
  void AddPdeRejectedPhoton(G4int channel, G4double time, G4double wavelength,
                            const G4Track* track) {
    ++fPdeRejected;
    if(channel>=0&&channel<128) {
      ++fPdeRejectedPlane[channel<64?0:1];
      ++fPdeRejectedByPlane[channel/32];
    }
    const auto found=fPhotonAudit.find(track?track->GetTrackID():-1);
    if(found!=fPhotonAudit.end())found->second.terminalGlobalChannel=channel;
    RecordSipmPhoton(channel, 0, wavelength, time);
    MarkPhotonTerminal(track,TerminalReason::sipm_pde_rejected,
                       AuditRegion::sipm);
  }
  void TrackOpticalStage(const G4Step* step) {
    const auto* track=step->GetTrack();
    const auto* prePV=step->GetPreStepPoint()->GetPhysicalVolume();
    const auto* postPV=step->GetPostStepPoint()->GetPhysicalVolume();
    if(!track||!prePV)return;
    const auto& preName=prePV->GetName();
    const auto planeFromName=[](const G4String& name) {
      if(name=="H1_X")return 0;
      if(name=="H1_Y")return 1;
      if(name=="H2_X")return 2;
      if(name=="H2_Y")return 3;
      if(name.rfind("H1X",0)==0)return 0;
      if(name.rfind("H1Y",0)==0)return 1;
      if(name.rfind("H2X",0)==0)return 2;
      if(name.rfind("H2Y",0)==0)return 3;
      return -1;
    };
    G4int plane=planeFromName(preName);
    // Count the active-area exit through the complete fibre assembly. A
    // photon may leave the core into a cladding before x/y=16 mm and remain
    // guided; restricting this counter to Core made the stages non-nested,
    // especially for square fibres.
    if(plane>=0) {
      const G4bool xAxis=plane==0||plane==2;
      const auto a=fDetector->PlaneWorldToNominalPosition(
          plane,step->GetPreStepPoint()->GetPosition());
      const auto b=fDetector->PlaneWorldToNominalPosition(
          plane,step->GetPostStepPoint()->GetPosition());
      const G4double before=xAxis?a.x():a.y(),after=xAxis?b.x():b.y();
      const auto activeExit=fDetector->GetActiveLongitudinalBounds(plane)[1];
      if(before<activeExit&&after>=activeExit)AddStage(plane,0,track->GetTrackID());
    }
    if(!postPV)return;
    const auto& postName=postPV->GetName();
    if(plane>=0&&preName.find("FanoutCore")!=std::string::npos&&
       (postName=="H1_Grease"||postName=="H2_Grease"))
      AddStage(plane,1,track->GetTrackID());
    if((preName=="H1_Grease"||preName=="H2_Grease")&&postName=="SiPMWindow") {
      const auto touch = step->GetPreStepPoint()->GetTouchableHandle();
      if (touch && touch->GetVolume()) {
        const G4int localCopy = touch->GetCopyNumber();
        if(localCopy>=0&&localCopy<64) {
          plane=(preName=="H1_Grease"?0:2)+(localCopy/32);
          AddStage(plane,2,track->GetTrackID());
        }
      }
    }
    if(postName=="SiPM") {
      const auto touch = step->GetPostStepPoint()->GetTouchableHandle();
      if (touch && touch->GetVolume()) {
        const G4int channel = touch->GetCopyNumber();
        if(channel>=0&&channel<128)AddStage(channel/32,3,track->GetTrackID());
      }
    }
  }
  void AddBulkAbsorption(const G4String& volume,const G4Track* track,
                         const G4StepPoint* prePoint) {
    const auto region=prePoint
        ? RegionAt(prePoint->GetPhysicalVolume(),prePoint->GetMaterial(),
                   prePoint->GetPosition(),AuditPlane(track))
        : AuditRegion::unknown;

    // Keep every legacy V13 event counter bit-for-bit invariant.  The more
    // precise material/position classification belongs only to the new V15
    // per-photon lineage below.
    if (volume.find("Core") != std::string::npos) {
      ++fBulkAbsCore;
    } else if (volume.find("Clad") != std::string::npos) {
      ++fBulkAbsClad;
    } else if (volume.find("Grease") != std::string::npos ||
               volume == "SiPMWindow") {
      ++fBulkAbsCoupling;
    } else {
      ++fBulkAbsOther;
    }

    const auto reason = IsCore(region) ||
                        volume.find("Core") != std::string::npos
        ? TerminalReason::bulk_core
        : IsCladding(region)
            ? (region==AuditRegion::active_outer_clad ||
               region==AuditRegion::fanout_outer_clad ||
               region==AuditRegion::far_end_outer_clad
                  ? TerminalReason::bulk_outer_clad
                  : TerminalReason::bulk_inner_clad)
            : region==AuditRegion::grease ||
              region==AuditRegion::window ||
              volume.find("Grease") != std::string::npos ||
              volume=="SiPMWindow"
                ? TerminalReason::bulk_coupling
                : TerminalReason::bulk_other;
    MarkPhotonTerminal(track,reason,region);
  }
  void AddDegradationLoss(G4int category,const G4Track* track) {
    if (category==0) ++fLossStress;
    else if (category==1) ++fLossCrack;
    else if (category==2) ++fLossPacking;
    else if (category==3) ++fLossAdhesive;
    else if (category==4) ++fLossH1Cladding;
    else if (category==5) ++fLossH2Cladding;
    const auto reason=category==0 ? TerminalReason::degradation_stress
        : category==1 ? TerminalReason::degradation_crack
        : category==2 ? TerminalReason::degradation_packing
        : category==3 ? TerminalReason::degradation_adhesive
        : category==4 ? TerminalReason::degradation_h1_cladding
                      : TerminalReason::degradation_h2_cladding;
    MarkPhotonTerminal(track,reason);
  }
  G4bool ApplyFanoutDegradation(const G4Step* step) {
    const auto* prePV=step->GetPreStepPoint()->GetPhysicalVolume();
    if (!prePV) return false;
    const auto& name=prePV->GetName();
    if (name.find("FanoutCore")==std::string::npos) return false;
    const auto& position=step->GetPreStepPoint()->GetPosition();
    const G4bool xAxis=name.find("XFanoutCore")!=std::string::npos;
    const G4bool isH1=name.rfind("H1",0)==0;
    const G4int plane=(isH1 ? 0 : 2)+(xAxis ? 0 : 1);
    const auto nominalPosition=
        fDetector->PlaneWorldToNominalPosition(plane,position);
    const G4double longitudinal=
        xAxis ? nominalPosition.x() : nominalPosition.y();
    if (longitudinal<=16.*mm) return false;
    std::array<G4double,6> alpha{
      std::max(0.,fDetector->GetStressLossPerM()),
      std::max(0.,fDetector->GetCrackLossPerM()),
      std::max(0.,fDetector->GetPackingLossPerM()),
      0., 0., 0.
    };
    if (longitudinal>=86.*mm-fDetector->GetAdhesiveIntrusionLength())
      alpha[3]=std::max(0.,fDetector->GetAdhesiveLossPerM());
    alpha[isH1 ? 4 : 5]=std::max(
        0.,isH1 ? fDetector->GetH1CladdingDamagePerM()
                : fDetector->GetH2CladdingDamagePerM());
    G4double total=0.;
    for (const auto value : alpha) total+=value;
    if (total<=0. ||
        G4UniformRand()>=1.-std::exp(-total*step->GetStepLength()))
      return false;
    G4double selector=G4UniformRand()*total;
    G4int category=0;
    while (category<5 && selector>=alpha[category]) {
      selector-=alpha[category];
      ++category;
    }
    AddDegradationLoss(category,step->GetTrack());
    const_cast<G4Track*>(step->GetTrack())->SetTrackStatus(fStopAndKill);
    return true;
  }
  G4bool ApplyV15InterfaceLoss(const G4Step* step) {
    if (step->GetPostStepPoint()->GetStepStatus() != fGeomBoundary) return false;
    const auto* pre=step->GetPreStepPoint()->GetPhysicalVolume();
    const auto* post=step->GetPostStepPoint()->GetPhysicalVolume();
    if(!pre||!post)return false;
    const auto& a=pre->GetName();
    const auto& b=post->GetName();
    const auto reject=[step](G4double efficiency) {
      const auto keep=std::clamp(efficiency,0.,1.);
      if(G4UniformRand()<keep)return false;
      const_cast<G4Track*>(step->GetTrack())->SetTrackStatus(fStopAndKill);
      return true;
    };
    if(a.find("FanoutCore")!=std::string::npos&&
       (b=="H1_Grease"||b=="H2_Grease")) {
      G4int plane=-1,globalChannel=-1;
      ResolvePlaneChannel(
          step->GetPreStepPoint()->GetTouchableHandle(),
          plane,globalChannel);
      const auto channelTransmission=
          fDetector->GetFanoutChannelTransmission(globalChannel);
      if(reject(fDetector->GetBendTransmission()*channelTransmission)) {
        ++fLossFanoutExit;
        MarkPhotonTerminal(step->GetTrack(),
                           TerminalReason::fanout_rejected);
        return true;
      }
      if(reject(fDetector->GetEndFaceEfficiency())) {
        ++fLossEndFace;
        MarkPhotonTerminal(step->GetTrack(),
                           TerminalReason::endface_rejected);
        return true;
      }
    }
    if((a=="H1_Grease"||a=="H2_Grease")&&b=="SiPMWindow"&&
       reject(fDetector->GetCouplingEfficiency())) {
      ++fLossGreaseCoupling;
      MarkPhotonTerminal(step->GetTrack(),
                         TerminalReason::coupling_rejected);
      return true;
    }
    return false;
  }
  void WriteSpatialClosure(
      G4int eventId,const std::array<G4double,128>& adc,
      const std::array<G4double,128>& noiselessAdc,
      const std::array<G4int,4>& parkBestChannel) {
    if(!fDetector->GetRecordSpatialClosureNtuples())return;
    if(fDetector->GetRecordSelectedDiagnosticsOnly()&&!HasTriggerCoincidence())return;
    auto* analysis=G4AnalysisManager::Instance();
    const G4double missing=std::numeric_limits<G4double>::quiet_NaN();
    for(G4int surface=0;surface<6;++surface) {
      const auto& state=fPrimarySurface[surface];
      analysis->FillNtupleIColumn(7,0,eventId);
      analysis->FillNtupleIColumn(7,1,state.track);
      analysis->FillNtupleIColumn(7,2,surface);
      analysis->FillNtupleIColumn(7,3,state.reached);
      analysis->FillNtupleDColumn(7,4,state.reached?state.position.x()/mm:missing);
      analysis->FillNtupleDColumn(7,5,state.reached?state.position.y()/mm:missing);
      analysis->FillNtupleDColumn(7,6,state.reached?state.position.z()/mm:missing);
      analysis->FillNtupleDColumn(7,7,state.reached?state.direction.x():missing);
      analysis->FillNtupleDColumn(7,8,state.reached?state.direction.y():missing);
      analysis->FillNtupleDColumn(7,9,state.reached?state.direction.z():missing);
      analysis->FillNtupleDColumn(
          7,10,state.reached?state.kineticEnergy/GeV:missing);
      analysis->FillNtupleDColumn(7,11,state.reached?state.time/ns:missing);
      analysis->FillNtupleDColumn(
          7,12,fDetector->GetValidationSurfaceZ(surface)/mm);
      analysis->AddNtupleRow(7);
    }

    std::vector<std::size_t> order(fFibreCrossings.size());
    std::iota(order.begin(),order.end(),0);
    std::sort(order.begin(),order.end(),[this](std::size_t left,
                                               std::size_t right) {
      const auto& a=fFibreCrossings[left];
      const auto& b=fFibreCrossings[right];
      return std::tie(a.track,a.crossingIndex,a.plane,a.physicalCopy) <
             std::tie(b.track,b.crossingIndex,b.plane,b.physicalCopy);
    });
    for(const auto index:order) {
      const auto& value=fFibreCrossings[index];
      analysis->FillNtupleIColumn(8,0,eventId);
      analysis->FillNtupleIColumn(8,1,value.track);
      analysis->FillNtupleIColumn(8,2,value.parent);
      analysis->FillNtupleIColumn(8,3,value.pdg);
      analysis->FillNtupleIColumn(8,4,value.crossingIndex);
      analysis->FillNtupleIColumn(8,5,value.plane);
      analysis->FillNtupleIColumn(8,6,value.physicalCopy);
      analysis->FillNtupleIColumn(8,7,value.localChannel);
      analysis->FillNtupleIColumn(8,8,value.globalChannel);
      analysis->FillNtupleIColumn(
          8,9,static_cast<G4int>(value.volumeClass));
      analysis->FillNtupleIColumn(8,10,value.firstStep);
      analysis->FillNtupleIColumn(8,11,value.lastStep);
      analysis->FillNtupleDColumn(8,12,value.entryWorld.x()/mm);
      analysis->FillNtupleDColumn(8,13,value.entryWorld.y()/mm);
      analysis->FillNtupleDColumn(8,14,value.entryWorld.z()/mm);
      analysis->FillNtupleDColumn(8,15,value.exitWorld.x()/mm);
      analysis->FillNtupleDColumn(8,16,value.exitWorld.y()/mm);
      analysis->FillNtupleDColumn(8,17,value.exitWorld.z()/mm);
      analysis->FillNtupleDColumn(8,18,value.entryLocal.x()/mm);
      analysis->FillNtupleDColumn(8,19,value.entryLocal.y()/mm);
      analysis->FillNtupleDColumn(8,20,value.entryLocal.z()/mm);
      analysis->FillNtupleDColumn(8,21,value.exitLocal.x()/mm);
      analysis->FillNtupleDColumn(8,22,value.exitLocal.y()/mm);
      analysis->FillNtupleDColumn(8,23,value.exitLocal.z()/mm);
      analysis->FillNtupleDColumn(8,24,value.entryDirection.x());
      analysis->FillNtupleDColumn(8,25,value.entryDirection.y());
      analysis->FillNtupleDColumn(8,26,value.entryDirection.z());
      analysis->FillNtupleDColumn(8,27,value.exitDirection.x());
      analysis->FillNtupleDColumn(8,28,value.exitDirection.y());
      analysis->FillNtupleDColumn(8,29,value.exitDirection.z());
      analysis->FillNtupleDColumn(8,30,value.entryEnergy/GeV);
      analysis->FillNtupleDColumn(8,31,value.exitEnergy/GeV);
      analysis->FillNtupleDColumn(8,32,value.entryTime/ns);
      analysis->FillNtupleDColumn(8,33,value.exitTime/ns);
      analysis->FillNtupleDColumn(8,34,value.path/mm);
      analysis->FillNtupleDColumn(8,35,value.corePath/mm);
      analysis->FillNtupleDColumn(8,36,value.activeCorePath/mm);
      analysis->FillNtupleDColumn(8,37,value.edep/MeV);
      analysis->FillNtupleDColumn(8,38,value.farEndCorePath/mm);
      analysis->FillNtupleDColumn(8,39,
          std::max(0.,value.corePath-value.activeCorePath-value.farEndCorePath)/mm);
      analysis->AddNtupleRow(8);
    }

    constexpr std::array<G4int,4> mid{{1,9,15,13}};
    for(G4int plane=0;plane<4;++plane) {
      const G4int bestGlobal=parkBestChannel[plane];
      const G4long accepted=std::accumulate(
          fDetected.begin()+plane*32,fDetected.begin()+(plane+1)*32,0L);
      G4bool hit=false;
      for(G4int channel=plane*32;channel<(plane+1)*32;++channel)
        hit=hit||adc[channel]>fParkThreshold[channel];
      analysis->FillNtupleIColumn(10,0,eventId);
      analysis->FillNtupleIColumn(10,1,plane);
      analysis->FillNtupleIColumn(10,2,mid[plane]);
      analysis->FillNtupleIColumn(
          10,3,static_cast<G4int>(fPrimaryCoreCopies[plane].size()));
      analysis->FillNtupleDColumn(10,4,fPrimaryCoreTrack[plane]/mm);
      analysis->FillNtupleDColumn(10,5,fPrimaryActiveCoreTrack[plane]/mm);
      analysis->FillNtupleDColumn(10,6,fCoreEdepPlane[plane]/MeV);
      analysis->FillNtupleIColumn(10,7,fProducedPlane[plane]);
      analysis->FillNtupleIColumn(10,8,fStagePlane[plane][3]);
      analysis->FillNtupleIColumn(10,9,accepted);
      analysis->FillNtupleIColumn(10,10,fPdeRejectedByPlane[plane]);
      analysis->FillNtupleIColumn(10,11,bestGlobal);
      analysis->FillNtupleIColumn(10,12,bestGlobal-plane*32);
      analysis->FillNtupleIColumn(10,13,fDetected[bestGlobal]);
      analysis->FillNtupleDColumn(10,14,noiselessAdc[bestGlobal]);
      analysis->FillNtupleDColumn(10,15,adc[bestGlobal]);
      analysis->FillNtupleDColumn(10,16,fParkThreshold[bestGlobal]);
      analysis->FillNtupleIColumn(10,17,hit);
      analysis->AddNtupleRow(10);
    }
    if(fDetector->GetRecordGeometryControlNtuples()) {
      constexpr auto ntupleId=kPrimaryTriggerCrossingNtupleId;
      for(G4int trigger=0;trigger<4;++trigger) {
        if(!fDetector->GetTriggerEnabled(trigger))continue;
        const auto& crossing=fPrimaryTriggerCrossing[trigger];
        const auto entryLocal=crossing.reached
            ? fDetector->TriggerWorldToLocalPosition(
                  trigger,crossing.entryWorld)
            : G4ThreeVector();
        const auto exitLocal=crossing.reached
            ? fDetector->TriggerWorldToLocalPosition(
                  trigger,crossing.exitWorld)
            : G4ThreeVector();
        const auto entryDirection=crossing.reached
            ? fDetector->TriggerWorldToLocalDirection(
                  trigger,crossing.entryDirection)
            : G4ThreeVector();
        const auto exitDirection=crossing.reached
            ? fDetector->TriggerWorldToLocalDirection(
                  trigger,crossing.exitDirection)
            : G4ThreeVector();
        analysis->FillNtupleIColumn(ntupleId,0,eventId);
        analysis->FillNtupleIColumn(ntupleId,1,crossing.track);
        analysis->FillNtupleIColumn(ntupleId,2,crossing.parent);
        analysis->FillNtupleIColumn(ntupleId,3,crossing.pdg);
        analysis->FillNtupleIColumn(ntupleId,4,trigger);
        analysis->FillNtupleIColumn(ntupleId,5,crossing.reached);
        const auto coordinate=[&](G4double value) {
          return crossing.reached?value/mm:missing;
        };
        analysis->FillNtupleDColumn(
            ntupleId,6,coordinate(crossing.entryWorld.x()));
        analysis->FillNtupleDColumn(
            ntupleId,7,coordinate(crossing.entryWorld.y()));
        analysis->FillNtupleDColumn(
            ntupleId,8,coordinate(crossing.entryWorld.z()));
        analysis->FillNtupleDColumn(
            ntupleId,9,coordinate(crossing.exitWorld.x()));
        analysis->FillNtupleDColumn(
            ntupleId,10,coordinate(crossing.exitWorld.y()));
        analysis->FillNtupleDColumn(
            ntupleId,11,coordinate(crossing.exitWorld.z()));
        analysis->FillNtupleDColumn(ntupleId,12,coordinate(entryLocal.x()));
        analysis->FillNtupleDColumn(ntupleId,13,coordinate(entryLocal.y()));
        analysis->FillNtupleDColumn(ntupleId,14,coordinate(entryLocal.z()));
        analysis->FillNtupleDColumn(ntupleId,15,coordinate(exitLocal.x()));
        analysis->FillNtupleDColumn(ntupleId,16,coordinate(exitLocal.y()));
        analysis->FillNtupleDColumn(ntupleId,17,coordinate(exitLocal.z()));
        analysis->FillNtupleDColumn(
            ntupleId,18,crossing.reached?entryDirection.x():missing);
        analysis->FillNtupleDColumn(
            ntupleId,19,crossing.reached?entryDirection.y():missing);
        analysis->FillNtupleDColumn(
            ntupleId,20,crossing.reached?entryDirection.z():missing);
        analysis->FillNtupleDColumn(
            ntupleId,21,crossing.reached?exitDirection.x():missing);
        analysis->FillNtupleDColumn(
            ntupleId,22,crossing.reached?exitDirection.y():missing);
        analysis->FillNtupleDColumn(
            ntupleId,23,crossing.reached?exitDirection.z():missing);
        analysis->FillNtupleDColumn(
            ntupleId,24,crossing.reached?crossing.entryTime/ns:missing);
        analysis->FillNtupleDColumn(
            ntupleId,25,crossing.reached?crossing.exitTime/ns:missing);
        analysis->FillNtupleDColumn(ntupleId,26,crossing.path/mm);
        analysis->FillNtupleDColumn(ntupleId,27,crossing.edep/MeV);
        analysis->FillNtupleIColumn(ntupleId,28,crossing.steps);
        analysis->AddNtupleRow(ntupleId);
      }
    }
  }
  void WriteGenerationLedger(const G4Event* event) {
    if(!fDetector->GetRecordGenerationNtuples())return;
    if(!event||!fGenerationSourceValid||fGenerationSource.targetPlane<0||
       fGenerationSource.targetPlane>=4)
      G4Exception("OpticalEventAction::WriteGenerationLedger",
                  "HODO_GATE3_EVENT_METADATA",FatalException,
                  "Stage-3 event/source/target-plane metadata is invalid.");
    const auto* vertex=event->GetPrimaryVertex();
    const auto* primary=vertex?vertex->GetPrimary():nullptr;
    if(!vertex||!primary)
      G4Exception("OpticalEventAction::WriteGenerationLedger",
                  "HODO_GATE3_PRIMARY",FatalException,
                  "Stage-3 event has no primary vertex/particle.");
    const G4ThreeVector actualDirection(
        primary->GetPx(),primary->GetPy(),primary->GetPz());
    const auto direction=actualDirection.mag2()>0.
        ? actualDirection.unit() : G4ThreeVector();
    const G4ThreeVector position(vertex->GetX0(),vertex->GetY0(),
                                 vertex->GetZ0());
    if((position-fGenerationSource.position).mag()>1.e-12*mm||
       (direction-fGenerationSource.direction).mag()>1.e-12||
       std::abs(primary->GetKineticEnergy()-fGenerationSource.energy)>
           1.e-9*GeV)
      G4Exception("OpticalEventAction::WriteGenerationLedger",
                  "HODO_GATE3_PRIMARY_SOURCE_MISMATCH",FatalException,
                  "Actual Stage-3 primary does not match its frozen TSV row.");

    auto* analysis=G4AnalysisManager::Instance();
    constexpr auto sourceId=kValidationGenerationSourceNtupleId;
    const std::array<G4int,9> sourceIntegers{{
        event->GetEventID(),fGenerationSource.event,
        fGenerationSource.cohort,fGenerationSource.targetPlane,
        fGenerationSource.targetChannel,
        fGenerationSource.expectedVolumeClass,
        fGenerationSource.expectedRegion,fGenerationSource.chordBin,
        fGenerationSource.repeatIndex}};
    for(std::size_t column=0;column<sourceIntegers.size();++column)
      analysis->FillNtupleIColumn(sourceId,column,sourceIntegers[column]);
    const std::array<G4double,9> sourceDoubles{{
        position.x()/mm,position.y()/mm,position.z()/mm,
        direction.x(),direction.y(),direction.z(),
        primary->GetKineticEnergy()/GeV,
        fGenerationSource.analyticCoreChord/mm,
        fGenerationSource.analyticTargetPath/mm}};
    for(std::size_t column=0;column<sourceDoubles.size();++column)
      analysis->FillNtupleDColumn(sourceId,9+column,sourceDoubles[column]);
    analysis->FillNtupleIColumn(sourceId,18,fGenerationSource.caseIndex);
    analysis->AddNtupleRow(sourceId);

    const G4int plane=fGenerationSource.targetPlane;
    const auto& ledger=fGenerationPlaneLedger[plane];
    constexpr auto eventId=kStage3GenerationEventNtupleId;
    const std::array<G4int,10> eventIntegers{{
        event->GetEventID(),fGenerationSource.event,
        fGenerationSource.cohort,fGenerationSource.targetPlane,
        fGenerationSource.targetChannel,
        fGenerationSource.expectedVolumeClass,
        fGenerationSource.expectedRegion,fGenerationSource.chordBin,
        fGenerationSource.repeatIndex,fGenerationSource.caseIndex}};
    for(std::size_t column=0;column<eventIntegers.size();++column)
      analysis->FillNtupleIColumn(eventId,column,eventIntegers[column]);
    const std::array<G4double,17> eventDoubles{{
        fGenerationSource.analyticCoreChord/mm,
        fGenerationSource.analyticTargetPath/mm,
        ledger.primaryCorePath/mm,ledger.primaryActiveCorePath/mm,
        (ledger.primaryCorePath-ledger.primaryActiveCorePath)/mm,
        ledger.primaryCladdingPath/mm,
        ledger.primaryActiveCladdingPath/mm,
        (ledger.primaryCladdingPath-
             ledger.primaryActiveCladdingPath)/mm,
        ledger.primaryCoreEdep/MeV,ledger.primaryCoreNiel/MeV,
        ledger.primaryCoreVisibleEdep/MeV,ledger.allCoreEdep/MeV,
        ledger.allCoreNiel/MeV,ledger.primaryCladdingEdep/MeV,
        ledger.primaryCladdingNiel/MeV,ledger.allCladdingEdep/MeV,
        ledger.allCladdingNiel/MeV}};
    for(std::size_t column=0;column<eventDoubles.size();++column)
      analysis->FillNtupleDColumn(eventId,10+column,eventDoubles[column]);
    const std::array<G4int,19> eventSuffix{{
        ledger.primaryCoreNonzeroNielSteps,
        ledger.primaryScintillation,ledger.secondaryScintillation,
        ledger.primaryCerenkov,ledger.secondaryCerenkov,
        ledger.primaryOther,ledger.secondaryOther,
        ledger.scintillationActive,ledger.scintillationFanout,
        ledger.cerenkovActive,ledger.cerenkovFanout,
        ledger.otherActive,ledger.otherFanout,
        ledger.claddingScintillation,ledger.claddingCerenkov,
        ledger.claddingOther,ledger.directOpticalTotal,
        ledger.chargedFibreSteps,ledger.primaryCoreSteps}};
    for(std::size_t column=0;column<eventSuffix.size();++column)
      analysis->FillNtupleIColumn(eventId,27+column,eventSuffix[column]);
    if(fGenerationSource.targetChannel<-1||
       fGenerationSource.targetChannel>=32)
      G4Exception("OpticalEventAction::WriteGenerationLedger",
                  "HODO_GATE3_TARGET_CHANNEL",FatalException,
                  "Stage-3 target channel must be -1 or [0,31].");
    const GenerationPlaneLedger emptyTargetLedger{};
    const auto& targetLedger=fGenerationSource.targetChannel>=0
        ? fGenerationChannelLedger[
              plane*32+fGenerationSource.targetChannel]
        : emptyTargetLedger;
    const std::array<G4double,15> targetChannelDoubles{{
        targetLedger.primaryCorePath/mm,
        targetLedger.primaryActiveCorePath/mm,
        (targetLedger.primaryCorePath-
             targetLedger.primaryActiveCorePath)/mm,
        targetLedger.primaryCladdingPath/mm,
        targetLedger.primaryActiveCladdingPath/mm,
        (targetLedger.primaryCladdingPath-
             targetLedger.primaryActiveCladdingPath)/mm,
        targetLedger.primaryCoreEdep/MeV,
        targetLedger.primaryCoreNiel/MeV,
        targetLedger.primaryCoreVisibleEdep/MeV,
        targetLedger.allCoreEdep/MeV,targetLedger.allCoreNiel/MeV,
        targetLedger.primaryCladdingEdep/MeV,
        targetLedger.primaryCladdingNiel/MeV,
        targetLedger.allCladdingEdep/MeV,
        targetLedger.allCladdingNiel/MeV}};
    for(std::size_t column=0;column<targetChannelDoubles.size();++column)
      analysis->FillNtupleDColumn(eventId,46+column,
                                  targetChannelDoubles[column]);
    const std::array<G4int,19> targetChannelIntegers{{
        targetLedger.primaryCoreNonzeroNielSteps,
        targetLedger.primaryScintillation,
        targetLedger.secondaryScintillation,
        targetLedger.primaryCerenkov,targetLedger.secondaryCerenkov,
        targetLedger.primaryOther,targetLedger.secondaryOther,
        targetLedger.scintillationActive,
        targetLedger.scintillationFanout,targetLedger.cerenkovActive,
        targetLedger.cerenkovFanout,targetLedger.otherActive,
        targetLedger.otherFanout,targetLedger.claddingScintillation,
        targetLedger.claddingCerenkov,targetLedger.claddingOther,
        targetLedger.directOpticalTotal,targetLedger.chargedFibreSteps,
        targetLedger.primaryCoreSteps}};
    for(std::size_t column=0;column<targetChannelIntegers.size();++column)
      analysis->FillNtupleIColumn(eventId,61+column,
                                  targetChannelIntegers[column]);
    analysis->AddNtupleRow(eventId);
  }
  void EndOfEventAction(const G4Event* event) override {
    // Gate 4 ends at optical transport.  Writing the established photon and
    // closure ledgers here avoids downstream ADC/noise work and the wide
    // event row without changing a photon or consuming an electronics random
    // number.  The trees remain booked, so the legacy schema is untouched.
    if(fDetector->GetGate4OpticalAudit()) {
      WritePhotonAudit(event->GetEventID());
      return;
    }
    std::array<G4double,128> adc{};
    std::array<G4double,128> noiselessAdc{};
    std::array<G4double,128> promptAdc{};
    std::array<G4double,128> waveformAdc{};
    std::array<G4double,128> noiseOnlyAdc{};
    std::array<G4long,128> effectiveDetected{};
    std::array<G4double,128> effectiveFirstTime{};
    effectiveFirstTime.fill(std::numeric_limits<G4double>::infinity());
    G4long totalDetected = 0;
    G4long detectedH1=0,detectedH2=0;
    if(fDetector->GetGeometryProbeMode()||
       fDetector->GetGenerationProbeMode()||
       fDetector->GetGate4OpticalAudit() ||
       (fDetector->GetSkipUntriggeredDigitization() && !HasTriggerCoincidence())) {
      // Probe lanes do not judge ADC. Do not draw electronics random numbers;
      // Stage-3/Gate-4 full optical transport still contributes detected
      // counts.  Gate 4 must remain isolated from downstream electronics.
      for(G4int channel=0;channel<128;++channel) {
        effectiveDetected[channel]=fDetected[channel];
        effectiveFirstTime[channel]=fFirstTime[channel];
        totalDetected+=fDetected[channel];
        if(channel<64)detectedH1+=fDetected[channel];
        else detectedH2+=fDetected[channel];
      }
    } else for (G4int channel=0; channel<128; ++channel) {
      auto avalancheTimes=fArrivalTimes[channel];
      const G4double effectiveScale=std::clamp(
          channel<64 ? fDetector->GetEffectiveLightScaleH1()
                     : fDetector->GetEffectiveLightScaleH2(),
          0.,1.);
      if(effectiveScale<=0.) {
        avalancheTimes.clear();
      } else if(effectiveScale<1.) {
        std::vector<G4double> keptTimes;
        keptTimes.reserve(avalancheTimes.size());
        for(std::size_t photonIndex=0; photonIndex<avalancheTimes.size();
            ++photonIndex) {
          const auto time=avalancheTimes[photonIndex];
          if(KeepEffectiveLightPe(
                 event->GetEventID(),channel,photonIndex,effectiveScale))
            keptTimes.push_back(time);
        }
        avalancheTimes.swap(keptTimes);
      }
      effectiveDetected[channel]=static_cast<G4long>(avalancheTimes.size());
      for(const auto time:avalancheTimes)
        if(time<effectiveFirstTime[channel])
          effectiveFirstTime[channel]=time;
      totalDetected += effectiveDetected[channel];
      if (channel<64) detectedH1+=effectiveDetected[channel];
      else detectedH2+=effectiveDetected[channel];
      const G4long dark=G4Poisson(
          fDetector->GetDarkCountRate()*fDetector->GetReadoutGate());
      for(G4long i=0;i<dark;++i)
        avalancheTimes.push_back(G4UniformRand()*fDetector->GetReadoutGate());
      G4long avalanches=static_cast<G4long>(avalancheTimes.size());

      // Galton-Watson prompt optical-crosstalk cascade, capped for malformed input.
      const auto pCross = std::clamp(fDetector->GetCrossTalkProbability(), 0., .95);
      G4long parent=0;
      while (parent<avalanches && avalanches < 100000) {
        if (G4UniformRand() < pCross) {
          avalancheTimes.push_back(avalancheTimes[parent]);
          ++avalanches;
        }
        ++parent;
      }
      const auto pAfter = std::clamp(fDetector->GetAfterPulseProbability(), 0., 1.);
      const G4long beforeAfterpulse = avalanches;
      for (G4long i=0; i<beforeAfterpulse; ++i) {
        if (G4UniformRand() < pAfter) {
          const auto u=std::clamp(G4UniformRand(),1.e-12,1.-1.e-12);
          avalancheTimes.push_back(avalancheTimes[i]-20.*ns*std::log(u));
          ++avalanches;
        }
      }

      const G4double cells = std::max(1, fDetector->GetMicrocellCount());
      const G4double fired = cells*(1.-std::exp(-static_cast<G4double>(avalanches)/cells));
      promptAdc[channel]=fired*fGain[channel];
      const G4double saturationScale=avalanches>0 ? fired/avalanches : 0.;
      waveformAdc[channel]=WaveformPeak(avalancheTimes,channel)*saturationScale*fGain[channel];
      noiselessAdc[channel]=fDetector->GetUseWaveformDigitization()
          ? waveformAdc[channel] : promptAdc[channel];
      if (fDetector->GetUseMeasuredNoise()) {
        const G4double sigma=std::max(
            0.,fNoiseRms[channel]*fDetector->GetAdcNoiseScale());
        if(fDetector->GetUseWaveformDigitization()) {
          const G4double pedestalError=G4RandGauss::shoot(
              0.,sigma/std::sqrt(std::max(1,fDetector->GetPedestalWindowSamples())));
          adc[channel]=WaveformPeakWithNoise(avalancheTimes,channel,
              saturationScale,sigma,pedestalError,noiseOnlyAdc[channel]);
          continue;
        }
        const G4double n=std::max(2,fDetector->GetPeakWindowSamples());
        const G4double root=std::sqrt(2.*std::log(n));
        const G4double b=root-(std::log(std::log(n))+
            std::log(4.*pi))/(2.*root);
        const G4double a=1./root;
        const G4double u=std::clamp(G4UniformRand(),1.e-12,1.-1.e-12);
        noiseOnlyAdc[channel]=sigma*(b-a*std::log(-std::log(u)));
        const G4double pedestalError=G4RandGauss::shoot(
            0.,sigma/std::sqrt(std::max(1,fDetector->GetPedestalWindowSamples())));
        const G4double signalSample=noiselessAdc[channel]+
            G4RandGauss::shoot(0.,sigma);
        adc[channel]=std::max(0.,std::max(noiseOnlyAdc[channel],signalSample)
                                  -pedestalError);
      } else {
        const G4double noise=G4RandGauss::shoot(
            0.,std::max(0.,fDetector->GetElectronicsNoise()));
        adc[channel]=std::max(0.,noiselessAdc[channel]+noise*fGain[channel]);
        noiseOnlyAdc[channel]=std::max(0.,noise)*fGain[channel];
      }
    }
    const auto roundXIt=std::max_element(adc.begin(),adc.begin()+32);
    const auto roundYIt=std::max_element(adc.begin()+32,adc.begin()+64);
    const auto squareXIt=std::max_element(adc.begin()+64,adc.begin()+96);
    const auto squareYIt=std::max_element(adc.begin()+96,adc.end());
    const G4double roundXMax=*roundXIt;
    const G4double roundYMax=*roundYIt;
    const G4double squareXMax=*squareXIt;
    const G4double squareYMax=*squareYIt;
    const G4int roundXChannel=std::distance(adc.begin(),roundXIt);
    const G4int roundYChannel=std::distance(adc.begin(),roundYIt);
    const G4int squareXChannel=std::distance(adc.begin(),squareXIt);
    const G4int squareYChannel=std::distance(adc.begin(),squareYIt);
    const auto h1It = std::max_element(adc.begin(), adc.begin()+64);
    const auto h2It = std::max_element(adc.begin()+64, adc.end());
    const G4double h1Max = *h1It;
    const G4double h2Max = *h2It;
    const G4int h1Channel = std::distance(adc.begin(), h1It);
    const G4int h2Channel = std::distance(adc.begin(), h2It);
    const G4double h1Noiseless=*std::max_element(
        noiselessAdc.begin(),noiselessAdc.begin()+64);
    const G4double h2Noiseless=*std::max_element(
        noiselessAdc.begin()+64,noiselessAdc.end());
    const G4double h1Prompt=*std::max_element(promptAdc.begin(),promptAdc.begin()+64);
    const G4double h2Prompt=*std::max_element(promptAdc.begin()+64,promptAdc.end());
    const G4double h1Waveform=*std::max_element(
        waveformAdc.begin(),waveformAdc.begin()+64);
    const G4double h2Waveform=*std::max_element(
        waveformAdc.begin()+64,waveformAdc.end());
    const G4double h1NoiseOnly=*std::max_element(
        noiseOnlyAdc.begin(),noiseOnlyAdc.begin()+64);
    const G4double h2NoiseOnly=*std::max_element(
        noiseOnlyAdc.begin()+64,noiseOnlyAdc.end());
    const G4long roundXPe=std::accumulate(
        effectiveDetected.begin(),effectiveDetected.begin()+32,0L);
    const G4long roundYPe=std::accumulate(
        effectiveDetected.begin()+32,effectiveDetected.begin()+64,0L);
    const G4long squareXPe=std::accumulate(
        effectiveDetected.begin()+64,effectiveDetected.begin()+96,0L);
    const G4long squareYPe=std::accumulate(
        effectiveDetected.begin()+96,effectiveDetected.end(),0L);
    const G4double roundXWaveform=*std::max_element(
        waveformAdc.begin(),waveformAdc.begin()+32);
    const G4double roundYWaveform=*std::max_element(
        waveformAdc.begin()+32,waveformAdc.begin()+64);
    const G4double squareXWaveform=*std::max_element(
        waveformAdc.begin()+64,waveformAdc.begin()+96);
    const G4double squareYWaveform=*std::max_element(
        waveformAdc.begin()+96,waveformAdc.end());
    const G4bool h1 = h1Max >= fDetector->GetH1AdcThreshold();
    const G4bool h2 = h2Max >= fDetector->GetH2AdcThreshold();
    const G4bool roundX = roundXMax >= fDetector->GetH1AdcThreshold();
    const G4bool roundY = roundYMax >= fDetector->GetH1AdcThreshold();
    const G4bool squareX = squareXMax >= fDetector->GetH2AdcThreshold();
    const G4bool squareY = squareYMax >= fDetector->GetH2AdcThreshold();
    const auto parkDecision = hodoscope::hit::Classify(
        fT1, fT2, fDetector->GetHitThreshold(), adc, fParkThreshold);
    const G4bool parkRoundX=
        parkDecision.plane[hodoscope::hit::kRoundX];
    const G4bool parkRoundY=
        parkDecision.plane[hodoscope::hit::kRoundY];
    const G4bool parkSquareX=
        parkDecision.plane[hodoscope::hit::kSquareX];
    const G4bool parkSquareY=
        parkDecision.plane[hodoscope::hit::kSquareY];
    std::array<G4int,4> parkBestChannel{};
    for(G4int plane=0;plane<4;++plane) {
      const G4int first=plane*32;
      parkBestChannel[plane]=first;
      for(G4int channel=first+1;channel<first+32;++channel)
        if(adc[channel]-fParkThreshold[channel] >
           adc[parkBestChannel[plane]]-fParkThreshold[parkBestChannel[plane]])
          parkBestChannel[plane]=channel;
    }
    const auto parkNoiseDecision =
        hodoscope::hit::ClassifyPlanes(noiseOnlyAdc, fParkThreshold);
    const G4bool parkRoundNoise=parkNoiseDecision.round;
    const G4bool parkSquareNoise=parkNoiseDecision.square;
    const auto pairDt = [this,&effectiveFirstTime](G4int xChannel,G4int yChannel) {
      const auto x=effectiveFirstTime[xChannel],y=effectiveFirstTime[yChannel];
      return std::isfinite(x)&&std::isfinite(y) ? std::abs(x-y) : -1.;
    };
    const G4double parkRoundDt=pairDt(roundXChannel,roundYChannel);
    const G4double parkSquareDt=pairDt(squareXChannel,squareYChannel);
    const G4bool parkRoundTimed=parkRoundX&&parkRoundY&&parkRoundDt>=0.&&
        parkRoundDt<=fDetector->GetXYCoincidenceWindow();
    const G4bool parkSquareTimed=parkSquareX&&parkSquareY&&parkSquareDt>=0.&&
        parkSquareDt<=fDetector->GetXYCoincidenceWindow();
    const G4bool t1 = parkDecision.trigger1;
    const G4bool t2 = parkDecision.trigger2;
    const G4bool t3 = !fDetector->GetTriggerEnabled(2) ? false
        : fT3 >= fDetector->GetHitThreshold();
    const G4bool t4 = !fDetector->GetTriggerEnabled(3) ? false
        : fT4 >= fDetector->GetHitThreshold();
    const auto virtualRoundPairs = VirtualRoundPairs(event);
    const G4bool selectedTrigger = HasSelectedTriggerCoincidence();
    if (event->GetEventID() < 10) {
      const auto printTime = [this](G4int channel) {
        return std::isfinite(fFirstTime[channel]) ? fFirstTime[channel]/ns : -1.;
      };
      G4cout << "[optical event] id=" << event->GetEventID()
             << " produced=" << fProduced << "(fiber=" << fProducedFiber
             << ",trigger=" << fProducedTrigger << ") detected=" << totalDetected
             << " pdeRejected=" << fPdeRejected << " escaped=" << fEscaped
             << " bulkAbs(core=" << fBulkAbsCore << ",clad=" << fBulkAbsClad
             << ",coupling=" << fBulkAbsCoupling << ",other=" << fBulkAbsOther << ")"
             << " degradation(stress=" << fLossStress << ",crack=" << fLossCrack
             << ",packing=" << fLossPacking << ",adhesive=" << fLossAdhesive
             << ",h1clad=" << fLossH1Cladding
             << ",h2clad=" << fLossH2Cladding
             << ",fanoutExit=" << fLossFanoutExit
             << ",endFace=" << fLossEndFace
             << ",greaseWindow=" << fLossGreaseCoupling << ")"
             << " H1(ch=" << h1Channel << ",n=" << fDetected[h1Channel]
             << ",first_ns=" << printTime(h1Channel) << ",ADC=" << h1Max << ")"
             << " H2(ch=" << h2Channel-64 << ",n=" << fDetected[h2Channel]
             << ",first_ns=" << printTime(h2Channel) << ",ADC=" << h2Max << ")"
             << " T1=" << t1 << " T2=" << t2
             << " T3=" << t3 << " T4=" << t4
             << " selectedTrigger=" << selectedTrigger
             << " H1=" << h1 << " H2=" << h2
             << G4endl;
      G4cout << "[muon energy] initial_MeV=" << fInitialEnergy/MeV
             << " at_T1_MeV=" << fAtTrigger1/MeV
             << " edep_concrete_MeV=" << fEdepConcrete/MeV
             << " edep_blackout_MeV=" << fEdepBlackout/MeV
             << " edep_triggers_MeV=" << fEdepTrigger/MeV
             << " edep_H1_MeV=" << fEdepH1/MeV
             << " edep_H2_MeV=" << fEdepH2/MeV << G4endl;
    }
    auto* run = static_cast<OpticalRun*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    run->Count(t1, t2, h1, h2, t3, t4, roundX, roundY, squareX, squareY,
               parkRoundX,parkRoundY,parkSquareX,parkSquareY,
               parkRoundTimed,parkSquareTimed,
               parkRoundNoise,parkSquareNoise,
               fProduced, totalDetected,
               fProducedH1,fProducedH2,detectedH1,detectedH2,fEdepH1,fEdepH2,
               fStageH1,fStageH2,fPdeRejectedPlane[0],fPdeRejectedPlane[1],
               fPrimaryCoreTrack,fPrimaryActiveCoreTrack);
    run->CountVirtualRound(t1, t2, t3, virtualRoundPairs[0],
                           virtualRoundPairs[1], virtualRoundPairs[2]);
    auto* analysis = G4AnalysisManager::Instance();
    analysis->FillNtupleIColumn(0, 0, event->GetEventID());
    analysis->FillNtupleDColumn(0, 1, fDetector->GetH1AttenuationLength()/m);
    analysis->FillNtupleIColumn(0, 2, fDetector->GetH1SingleClad());
    analysis->FillNtupleIColumn(0, 3, fDetector->GetH1Uses22PointEmission());
    analysis->FillNtupleIColumn(0, 4, fProduced);
    analysis->FillNtupleIColumn(0, 5, totalDetected);
    analysis->FillNtupleIColumn(0, 6, fPdeRejected);
    analysis->FillNtupleIColumn(0, 7, fEscaped);
    analysis->FillNtupleIColumn(0, 8, fBulkAbsCore);
    analysis->FillNtupleIColumn(0, 9, fBulkAbsClad);
    analysis->FillNtupleIColumn(0, 10, fBulkAbsCoupling);
    analysis->FillNtupleIColumn(0, 11, fBulkAbsOther);
    analysis->FillNtupleIColumn(0, 12, fDetector->GetMicrocellCount());
    analysis->FillNtupleDColumn(0, 13, fDetector->GetDarkCountRate()/(1.e3/second));
    analysis->FillNtupleDColumn(0, 14, fDetector->GetPDEScale());
    analysis->FillNtupleDColumn(0, 15, fDetector->GetCouplingEfficiency());
    analysis->FillNtupleDColumn(0, 16, fDetector->GetBendTransmission());
    analysis->FillNtupleDColumn(0, 17, fDetector->GetScintillationYieldScale());
    const G4int acceptedAfter = t1 && t2
        ? gAcceptedT1T2.fetch_add(1, std::memory_order_relaxed)+1
        : gAcceptedT1T2.load(std::memory_order_relaxed);
    const G4int selectedAcceptedAfter = selectedTrigger
        ? gAcceptedTriggerCoincidence.fetch_add(1, std::memory_order_relaxed)+1
        : gAcceptedTriggerCoincidence.load(std::memory_order_relaxed);
    const G4int progressBucket = selectedAcceptedAfter/100;
    G4int previousBucket = gLastProgressBucket.load(std::memory_order_relaxed);
    const G4bool reportAccepted =
        progressBucket > previousBucket &&
        gLastProgressBucket.compare_exchange_strong(
            previousBucket, progressBucket, std::memory_order_relaxed);
    if ((event->GetEventID()+1)%1000 == 0 || reportAccepted) {
      G4cout << "[optical-v21 progress] event=" << event->GetEventID()
             << " processed=" << event->GetEventID()+1
             << " accepted_t1t2=" << acceptedAfter
             << " accepted_selected=" << selectedAcceptedAfter
             << " mode=" << fDetector->GetTriggerCoincidenceMode()
             << " target=" << fDetector->GetTargetCoincidences()
             << G4endl;
    }
    analysis->FillNtupleIColumn(0, 18, t1);
    analysis->FillNtupleIColumn(0, 19, t2);
    analysis->FillNtupleIColumn(0, 20, h1);
    analysis->FillNtupleIColumn(0, 21, h2);
    analysis->FillNtupleDColumn(0, 22, fT1/MeV);
    analysis->FillNtupleDColumn(0, 23, fT2/MeV);
    analysis->FillNtupleIColumn(0, 24, acceptedAfter);
    analysis->FillNtupleIColumn(0, 25, fDetector->GetTargetCoincidences());
    const auto* vertex = event->GetPrimaryVertex();
    const auto* primary = vertex ? vertex->GetPrimary() : nullptr;
    const G4ThreeVector momentum = primary
        ? G4ThreeVector(primary->GetPx(), primary->GetPy(), primary->GetPz()).unit()
        : G4ThreeVector();
    analysis->FillNtupleDColumn(0, 26, vertex ? vertex->GetX0()/mm : 0.);
    analysis->FillNtupleDColumn(0, 27, vertex ? vertex->GetY0()/mm : 0.);
    analysis->FillNtupleDColumn(0, 28, vertex ? vertex->GetZ0()/mm : 0.);
    analysis->FillNtupleDColumn(0, 29, momentum.x());
    analysis->FillNtupleDColumn(0, 30, momentum.y());
    analysis->FillNtupleDColumn(0, 31, momentum.z());
    analysis->FillNtupleDColumn(0, 32,
        primary ? std::acos(std::clamp(-momentum.z(), -1., 1.))/deg : 0.);
    analysis->FillNtupleDColumn(0, 33, h1Max);
    analysis->FillNtupleDColumn(0, 34, h2Max);
    analysis->FillNtupleIColumn(0, 35, fLossStress);
    analysis->FillNtupleIColumn(0, 36, fLossCrack);
    analysis->FillNtupleIColumn(0, 37, fLossPacking);
    analysis->FillNtupleIColumn(0, 38, fLossAdhesive);
    analysis->FillNtupleDColumn(0, 39, fInitialEnergy/GeV);
    analysis->FillNtupleDColumn(0, 40, fAtTrigger1/GeV);
    analysis->FillNtupleDColumn(0, 41, fEdepConcrete/MeV);
    analysis->FillNtupleDColumn(0, 42, fEdepBlackout/MeV);
    analysis->FillNtupleIColumn(0, 43, fLossH1Cladding);
    analysis->FillNtupleIColumn(0, 44, fLossH2Cladding);
    analysis->FillNtupleIColumn(0, 45, fDetector->GetUseMeasuredNoise());
    analysis->FillNtupleDColumn(0, 46, fDetector->GetAdcNoiseScale());
    analysis->FillNtupleDColumn(0, 47, h1Noiseless);
    analysis->FillNtupleDColumn(0, 48, h2Noiseless);
    analysis->FillNtupleDColumn(0, 49, h1NoiseOnly);
    analysis->FillNtupleDColumn(0, 50, h2NoiseOnly);
    analysis->FillNtupleIColumn(0, 51,
        h1 && h1Noiseless < fDetector->GetH1AdcThreshold());
    analysis->FillNtupleIColumn(0, 52,
        h2 && h2Noiseless < fDetector->GetH2AdcThreshold());
    analysis->FillNtupleIColumn(0, 53,
        !h1 && h1Noiseless >= fDetector->GetH1AdcThreshold());
    analysis->FillNtupleIColumn(0, 54,
        !h2 && h2Noiseless >= fDetector->GetH2AdcThreshold());
    analysis->FillNtupleIColumn(0,55,fDetector->GetUseWaveformDigitization());
    analysis->FillNtupleDColumn(0,56,fDetector->GetAdcSamplePeriod()/ns);
    analysis->FillNtupleDColumn(0,57,fDetector->GetSinglePeFwhmSamples());
    analysis->FillNtupleDColumn(0,58,h1Prompt);
    analysis->FillNtupleDColumn(0,59,h2Prompt);
    analysis->FillNtupleDColumn(0,60,h1Waveform);
    analysis->FillNtupleDColumn(0,61,h2Waveform);
    analysis->FillNtupleIColumn(0,62,fDetector->GetH2SingleClad());
    analysis->FillNtupleDColumn(0,63,fEdepH1/MeV);
    analysis->FillNtupleDColumn(0,64,fEdepH2/MeV);
    analysis->FillNtupleDColumn(0,65,fEdepH1Core/MeV);
    analysis->FillNtupleDColumn(0,66,fEdepH2Core/MeV);
    analysis->FillNtupleIColumn(0,67,fProducedH1);
    analysis->FillNtupleIColumn(0,68,fProducedH2);
    analysis->FillNtupleIColumn(0,69,detectedH1);
    analysis->FillNtupleIColumn(0,70,detectedH2);
    const G4double dz = momentum.z();
    const G4double projectScale = vertex && std::abs(dz)>1.e-12
        ? -vertex->GetZ0()/dz : 0.;
    analysis->FillNtupleDColumn(0,71,
        vertex ? (vertex->GetX0()+projectScale*momentum.x())/mm : 0.);
    analysis->FillNtupleDColumn(0,72,
        vertex ? (vertex->GetY0()+projectScale*momentum.y())/mm : 0.);
    const auto h1First=*std::min_element(
        effectiveFirstTime.begin(),effectiveFirstTime.begin()+64);
    const auto h2First=*std::min_element(
        effectiveFirstTime.begin()+64,effectiveFirstTime.end());
    const auto relative=[this](G4double first) {
      return std::isfinite(first) && std::isfinite(fTrigger1Time)
          ? (first-fTrigger1Time)/ns : -1.;
    };
    const G4double h1Relative=relative(h1First),h2Relative=relative(h2First);
    const G4double delay=fDetector->GetAdcPeakWindowDelay()/ns;
    const G4double samplePeriod=fDetector->GetAdcSamplePeriod()/ns;
    analysis->FillNtupleDColumn(0,73,
        std::isfinite(fTrigger1Time)?fTrigger1Time/ns:-1.);
    analysis->FillNtupleDColumn(0,74,delay);
    analysis->FillNtupleDColumn(0,75,h1Relative);
    analysis->FillNtupleDColumn(0,76,h2Relative);
    analysis->FillNtupleDColumn(0,77,h1Relative>=0.?200.+(h1Relative-delay)/samplePeriod:-1.);
    analysis->FillNtupleDColumn(0,78,h2Relative>=0.?200.+(h2Relative-delay)/samplePeriod:-1.);
    analysis->FillNtupleDColumn(0,79,roundXMax);
    analysis->FillNtupleDColumn(0,80,roundYMax);
    analysis->FillNtupleDColumn(0,81,squareXMax);
    analysis->FillNtupleDColumn(0,82,squareYMax);
    analysis->FillNtupleIColumn(0,83,roundX);
    analysis->FillNtupleIColumn(0,84,roundY);
    analysis->FillNtupleIColumn(0,85,squareX);
    analysis->FillNtupleIColumn(0,86,squareY);
    analysis->FillNtupleIColumn(0,87,roundX&&roundY);
    analysis->FillNtupleIColumn(0,88,squareX&&squareY);
    analysis->FillNtupleIColumn(0,89,roundX&&roundY&&squareX&&squareY);
    analysis->FillNtupleIColumn(0,90,parkRoundX);
    analysis->FillNtupleIColumn(0,91,parkRoundY);
    analysis->FillNtupleIColumn(0,92,parkSquareX);
    analysis->FillNtupleIColumn(0,93,parkSquareY);
    analysis->FillNtupleIColumn(0,94,parkRoundX&&parkRoundY);
    analysis->FillNtupleIColumn(0,95,parkSquareX&&parkSquareY);
    analysis->FillNtupleIColumn(0,96,parkRoundNoise);
    analysis->FillNtupleIColumn(0,97,parkSquareNoise);
    analysis->FillNtupleDColumn(0,98,parkRoundDt>=0.?parkRoundDt/ns:-1.);
    analysis->FillNtupleDColumn(0,99,parkSquareDt>=0.?parkSquareDt/ns:-1.);
    analysis->FillNtupleIColumn(0,100,parkRoundTimed);
    analysis->FillNtupleIColumn(0,101,parkSquareTimed);
    analysis->FillNtupleIColumn(0,102,fStageH1[3]);
    analysis->FillNtupleIColumn(0,103,fStageH2[3]);
    analysis->FillNtupleIColumn(0,104,fPdeRejectedPlane[0]);
    analysis->FillNtupleIColumn(0,105,fPdeRejectedPlane[1]);
    analysis->FillNtupleIColumn(0,106,roundXPe);
    analysis->FillNtupleIColumn(0,107,roundYPe);
    analysis->FillNtupleIColumn(0,108,squareXPe);
    analysis->FillNtupleIColumn(0,109,squareYPe);
    analysis->FillNtupleDColumn(0,110,roundXWaveform);
    analysis->FillNtupleDColumn(0,111,roundYWaveform);
    analysis->FillNtupleDColumn(0,112,squareXWaveform);
    analysis->FillNtupleDColumn(0,113,squareYWaveform);
    for(G4int plane=0;plane<4;++plane) {
      analysis->FillNtupleDColumn(0,114+plane,fPrimaryCoreTrack[plane]/mm);
      analysis->FillNtupleIColumn(0,118+plane,fPrimaryCoreCopies[plane].size());
      analysis->FillNtupleDColumn(0,122+plane,fCoreEdepPlane[plane]/MeV);
      analysis->FillNtupleIColumn(0,126+plane,fProducedPlane[plane]);
      analysis->FillNtupleDColumn(0,130+plane,fPrimaryCoreEdepPlane[plane]/MeV);
      analysis->FillNtupleDColumn(0,134+plane,fSecondaryCoreEdepPlane[plane]/MeV);
      analysis->FillNtupleIColumn(0,138+plane,fPrimaryProducedPlane[plane]);
      analysis->FillNtupleIColumn(0,142+plane,fSecondaryProducedPlane[plane]);
    }
    analysis->FillNtupleIColumn(0,146,fLossFanoutExit);
    analysis->FillNtupleIColumn(0,147,fLossEndFace);
    analysis->FillNtupleIColumn(0,148,fLossGreaseCoupling);
    analysis->FillNtupleIColumn(0,149,fStageH1[0]);
    analysis->FillNtupleIColumn(0,150,fStageH1[1]);
    analysis->FillNtupleIColumn(0,151,fStageH1[2]);
    analysis->FillNtupleIColumn(0,152,fStageH2[0]);
    analysis->FillNtupleIColumn(0,153,fStageH2[1]);
    analysis->FillNtupleIColumn(0,154,fStageH2[2]);
    analysis->FillNtupleIColumn(0,155,fDetector->GetFixedPhotonTest());
    const G4int fixedPhotonChannel=
        fDetector->GetFixedPhotonChannelForEvent(event->GetEventID());
    analysis->FillNtupleIColumn(0,156,fixedPhotonChannel);
    const G4double actualFixedWavelength =
        primary && primary->GetG4code()==G4OpticalPhoton::Definition() &&
        primary->GetTotalEnergy()>0.
        ? h_Planck*c_light/primary->GetTotalEnergy()/nm
        : fDetector->GetFixedPhotonWavelength()/nm;
    analysis->FillNtupleDColumn(0,157,actualFixedWavelength);
    const auto fixedPlane=fDetector->GetFixedPhotonPlane();
    const G4int fixedOffset=
        fixedPlane=="h1y"?32:fixedPlane=="h2x"?64:fixedPlane=="h2y"?96:0;
    const G4int fixedGlobal=fixedOffset+fixedPhotonChannel;
    analysis->FillNtupleIColumn(0,158,fixedGlobal);
    analysis->FillNtupleDColumn(0,159,
        fDetector->GetFanoutMinRadius(fixedGlobal)/mm);
    analysis->FillNtupleDColumn(0,160,
        fDetector->GetFanoutChannelTransmission(fixedGlobal));
    for(G4int plane=0;plane<4;++plane) {
      const G4int channel=parkBestChannel[plane];
      const G4int first=161+plane*9;
      analysis->FillNtupleIColumn(0,first,channel-plane*32);
      analysis->FillNtupleIColumn(0,first+1,effectiveDetected[channel]);
      analysis->FillNtupleDColumn(0,first+2,noiselessAdc[channel]);
      analysis->FillNtupleDColumn(0,first+3,adc[channel]);
      analysis->FillNtupleDColumn(0,first+4,fParkThreshold[channel]);
      analysis->FillNtupleDColumn(0,first+5,
          adc[channel]-fParkThreshold[channel]);
      analysis->FillNtupleDColumn(0,first+6,fGain[channel]);
      analysis->FillNtupleDColumn(0,first+7,fNoiseRms[channel]);
      analysis->FillNtupleDColumn(0,first+8,noiseOnlyAdc[channel]);
    }
    for(G4int plane=0;plane<4;++plane) {
      analysis->FillNtupleIColumn(0,197+plane*2,fProducedActivePlane[plane]);
      analysis->FillNtupleIColumn(0,198+plane*2,fProducedFanoutPlane[plane]);
    }
    analysis->FillNtupleDColumn(0,205,
        fDetector->GetH2AttenuationLength()/m);
    analysis->FillNtupleDColumn(0,206,
        fDetector->GetH1ScintillationYield());
    analysis->FillNtupleDColumn(0,207,
        fDetector->GetH2ScintillationYield());
    analysis->FillNtupleDColumn(0,208,fDetector->GetH1BirksConstant());
    analysis->FillNtupleDColumn(0,209,fDetector->GetH2BirksConstant());
    analysis->FillNtupleDColumn(0,210,fDetector->GetTriggerBirksConstant());
    for(G4int plane=0;plane<4;++plane)
      for(G4int stage=0;stage<4;++stage)
        analysis->FillNtupleIColumn(
            0,211+plane*4+stage,fStagePlane[plane][stage]);
    if (fDetector->GetRecordActiveCoreChordBranches())
      for (G4int plane=0; plane<4; ++plane)
        analysis->FillNtupleDColumn(
            0,227+plane,fPrimaryActiveCoreTrack[plane]/mm);
    const G4int triggerColumnBase =
        fDetector->GetRecordActiveCoreChordBranches() ? 231 : 227;
    analysis->FillNtupleIColumn(0,triggerColumnBase,t3);
    analysis->FillNtupleIColumn(0,triggerColumnBase+1,t4);
    analysis->FillNtupleDColumn(0,triggerColumnBase+2,fT3/MeV);
    analysis->FillNtupleDColumn(0,triggerColumnBase+3,fT4/MeV);
    analysis->FillNtupleIColumn(0,triggerColumnBase+4,selectedTrigger);
    analysis->FillNtupleIColumn(0,triggerColumnBase+5,selectedAcceptedAfter);
    if(fDetector->GetActivePlaneOverlap()||fDetector->GetFibreFarEndExtension()>0.)
      for(G4int plane=0;plane<4;++plane)
        analysis->FillNtupleIColumn(0,triggerColumnBase+6+plane,fProducedFarEndPlane[plane]);
    if(selectedTrigger || event->GetEventID()%fDetector->GetUnselectedEventStride()==0)
      analysis->AddNtupleRow(0);
    if(fDetector->GetRecordChannelNtuple() &&
       (!fDetector->GetRecordSelectedDiagnosticsOnly()||HasTriggerCoincidence())) {
      for(G4int channel=0;channel<128;++channel) {
        analysis->FillNtupleIColumn(2,0,event->GetEventID());
        analysis->FillNtupleIColumn(2,1,channel);
        analysis->FillNtupleIColumn(2,2,channel/32);
        analysis->FillNtupleIColumn(2,3,channel%32);
        analysis->FillNtupleDColumn(2,4,adc[channel]);
        analysis->FillNtupleDColumn(2,5,noiseOnlyAdc[channel]);
        analysis->FillNtupleDColumn(2,6,noiselessAdc[channel]);
        analysis->FillNtupleDColumn(2,7,fParkThreshold[channel]);
        analysis->FillNtupleDColumn(2,8,fGain[channel]);
        analysis->FillNtupleDColumn(2,9,fNoiseRms[channel]);
        analysis->AddNtupleRow(2);
      }
    }
    WriteGenerationLedger(event);
    WriteSpatialClosure(
        event->GetEventID(),adc,noiselessAdc,parkBestChannel);
    WritePhotonAudit(event->GetEventID());
    if (fDetector->GetTargetCoincidences() > 0 &&
        selectedTrigger &&
        selectedAcceptedAfter >= fDetector->GetTargetCoincidences()) {
      G4cout << "[optical-v21 stop] event=" << event->GetEventID()
             << " accepted_t1t2=" << acceptedAfter
             << " accepted_selected=" << selectedAcceptedAfter
             << " mode=" << fDetector->GetTriggerCoincidenceMode()
             << " target=" << fDetector->GetTargetCoincidences() << G4endl;
      G4RunManager::GetRunManager()->AbortRun(true);
    }
  }
 private:
  std::array<G4bool,3> VirtualRoundPairs(const G4Event* event) const {
    std::array<G4bool,3> pairs{{false,false,false}};
    if(!fDetector->GetVirtualRoundModelEnabled() || !event) return pairs;
    const auto* vertex = event->GetPrimaryVertex();
    const auto* primary = vertex ? vertex->GetPrimary() : nullptr;
    if(!primary) return pairs;
    const G4ThreeVector p0(vertex->GetX0(),vertex->GetY0(),vertex->GetZ0());
    const G4ThreeVector direction(primary->GetPx(),primary->GetPy(),
                                  primary->GetPz());
    if(direction.mag2()==0.) return pairs;
    const G4ThreeVector unit = direction.unit();
    pairs[0]=VirtualRoundModulePair(p0,unit,fDetector->GetHodoscopeFanoutZ(0),
                                    0.*mm,0.*mm,0.*deg);
    pairs[1]=VirtualRoundModulePair(p0,unit,fDetector->GetHodoscopeFanoutZ(1),
                                    0.*mm,0.*mm,0.*deg);
    if(fDetector->GetVirtualRoundH3Enabled()) {
      pairs[2]=VirtualRoundModulePair(
          p0,unit,fDetector->GetVirtualRoundH3Z(),
          fDetector->GetVirtualRoundH3OffsetX(),
          fDetector->GetVirtualRoundH3OffsetY(),
          fDetector->GetVirtualRoundH3Yaw());
    }
    return pairs;
  }

  G4bool VirtualRoundModulePair(const G4ThreeVector& p0,
                                const G4ThreeVector& unit,
                                G4double moduleBaseZ,
                                G4double offsetX,
                                G4double offsetY,
                                G4double yaw) const {
    return VirtualRoundPlaneHit(p0,unit,0,moduleBaseZ,offsetX,offsetY,yaw) &&
           VirtualRoundPlaneHit(p0,unit,1,moduleBaseZ,offsetX,offsetY,yaw);
  }

  G4bool VirtualRoundPlaneHit(const G4ThreeVector& p0,
                              const G4ThreeVector& unit,
                              G4int templatePlane,
                              G4double moduleBaseZ,
                              G4double offsetX,
                              G4double offsetY,
                              G4double yaw) const {
    constexpr G4double coreRadius = 0.4795*mm;
    constexpr G4double activeHalfLength = 16.*mm;
    const G4bool xAxis = templatePlane == 0;
    const G4double c = std::cos(yaw);
    const G4double s = std::sin(yaw);
    const G4ThreeVector axisLocal = xAxis ? G4ThreeVector(1.,0.,0.)
                                          : G4ThreeVector(0.,1.,0.);
    const G4ThreeVector axis(c*axisLocal.x()-s*axisLocal.y(),
                             s*axisLocal.x()+c*axisLocal.y(),0.);
    const auto toWorld = [&](const G4ThreeVector& local) {
      return G4ThreeVector(offsetX + c*local.x()-s*local.y(),
                           offsetY + s*local.x()+c*local.y(),
                           local.z());
    };
    for(G4int channel=0;channel<32;++channel) {
      const G4double transverse =
          fDetector->GetPlaneSourceTransverse(templatePlane,channel);
      const G4double relativeZ =
          fDetector->GetPlaneSourceZ(templatePlane,channel) -
          fDetector->GetHodoscopeFanoutZ(0);
      const G4ThreeVector localCentre =
          xAxis ? G4ThreeVector(0.,transverse,moduleBaseZ+relativeZ)
                : G4ThreeVector(transverse,0.,moduleBaseZ+relativeZ);
      const G4ThreeVector centre = toWorld(localCentre);
      G4double axisPosition = 0.;
      G4double trackPosition = 0.;
      if(!ClosestLineParameters(centre,axis,p0,unit,
                                axisPosition,trackPosition)) {
        continue;
      }
      if(std::abs(axisPosition) > activeHalfLength) continue;
      const auto axisPoint = centre + axisPosition*axis;
      const auto trackPoint = p0 + trackPosition*unit;
      if((axisPoint-trackPoint).mag() <= coreRadius) return true;
    }
    return false;
  }

  static G4bool ClosestLineParameters(const G4ThreeVector& p1,
                                      const G4ThreeVector& d1,
                                      const G4ThreeVector& p2,
                                      const G4ThreeVector& d2,
                                      G4double& s1,
                                      G4double& s2) {
    const G4double a=d1.dot(d1);
    const G4double b=d1.dot(d2);
    const G4double c=d2.dot(d2);
    const G4ThreeVector w=p1-p2;
    const G4double d=d1.dot(w);
    const G4double e=d2.dot(w);
    const G4double denom=a*c-b*b;
    if(std::abs(denom)<1.e-12) return false;
    s1=(b*e-c*d)/denom;
    s2=(a*e-b*d)/denom;
    return std::isfinite(s1)&&std::isfinite(s2);
  }

  static G4bool IsCore(AuditRegion region) {
    return region==AuditRegion::active_core||
           region==AuditRegion::fanout_core||region==AuditRegion::far_end_core;
  }
  static G4bool IsCladding(AuditRegion region) {
    return region==AuditRegion::active_inner_clad||
           region==AuditRegion::fanout_inner_clad||
           region==AuditRegion::active_outer_clad||
           region==AuditRegion::fanout_outer_clad||
           region==AuditRegion::far_end_inner_clad||region==AuditRegion::far_end_outer_clad;
  }
  static G4bool IsReflection(G4OpBoundaryProcessStatus status) {
    return status==FresnelReflection||
           status==TotalInternalReflection||
           status==LambertianReflection||
           status==LobeReflection||
           status==SpikeReflection||
           status==BackScattering;
  }
  static void LocatePlaneChannel(const G4TouchableHandle& touch,
                                 G4int& plane,G4int& channel) {
    ResolvePlaneChannel(touch,plane,channel);
  }
  AuditRegion RegionAt(const G4VPhysicalVolume* volume,
                       const G4Material* material,
                       const G4ThreeVector& position,
                       G4int plane) const {
    if(!volume)return AuditRegion::world;
    const auto& name=volume->GetName();
    const auto materialName=material?material->GetName():G4String{};
    const G4bool xAxis=plane==0||plane==2;
    const auto nominalPosition=
        fDetector->PlaneWorldToNominalPosition(plane,position);
    const auto longitudinal=
        xAxis?nominalPosition.x():nominalPosition.y();
    const auto activeBounds=fDetector->GetActiveLongitudinalBounds(plane);
    const G4bool active=longitudinal>=activeBounds[0]&&longitudinal<=activeBounds[1];
    if(longitudinal < activeBounds[0]) {
      if(materialName=="Kuraray_SCSF78_core"||materialName=="BCF20_core")
        return AuditRegion::far_end_core;
      if(materialName=="PMMA_inner_cladding")return AuditRegion::far_end_inner_clad;
      if(materialName=="Fluorinated_outer_cladding")return AuditRegion::far_end_outer_clad;
    }
    if(materialName=="Kuraray_SCSF78_core"||
       materialName=="BCF20_core")
      return active?AuditRegion::active_core:AuditRegion::fanout_core;
    if(materialName=="PMMA_inner_cladding")
      return active?AuditRegion::active_inner_clad
                   :AuditRegion::fanout_inner_clad;
    if(materialName=="Fluorinated_outer_cladding")
      return active?AuditRegion::active_outer_clad
                   :AuditRegion::fanout_outer_clad;
    if(name=="H1_Grease"||name=="H2_Grease")
      return AuditRegion::grease;
    if(name=="SiPMWindow")return AuditRegion::window;
    if(name=="SiPM")return AuditRegion::sipm;
    if(name=="World")return AuditRegion::world;
    return AuditRegion::other;
  }
  static void AddAuditPath(PhotonAudit& audit,AuditRegion region,
                           G4double length) {
    switch(region) {
      case AuditRegion::active_core:audit.activeCorePath+=length;break;
      case AuditRegion::fanout_core:audit.fanoutCorePath+=length;break;
      case AuditRegion::active_inner_clad:
        audit.activeInnerCladPath+=length;break;
      case AuditRegion::fanout_inner_clad:
        audit.fanoutInnerCladPath+=length;break;
      case AuditRegion::active_outer_clad:
        audit.activeOuterCladPath+=length;break;
      case AuditRegion::fanout_outer_clad:
        audit.fanoutOuterCladPath+=length;break;
      case AuditRegion::grease:audit.greasePath+=length;break;
      case AuditRegion::window:audit.windowPath+=length;break;
      default:audit.otherPath+=length;break;
    }
  }
  G4int AuditPlane(const G4Track* track) const {
    if(!track)return -1;
    const auto found=fPhotonAudit.find(track->GetTrackID());
    return found==fPhotonAudit.end()?-1:found->second.plane;
  }
  void WritePhotonAudit(G4int eventId) {
    if(!fDetector->GetRecordPhotonAuditNtuple())return;
    if(fDetector->GetRecordSelectedDiagnosticsOnly()&&!HasTriggerCoincidence())return;
    struct Closure {
      G4long registered=0,terminal=0,detected=0,rejected=0;
      G4long bulkCore=0,bulkInner=0,bulkOuter=0;
      G4long bulkCoupling=0,bulkOther=0,world=0,degradation=0;
      G4long fanoutRejected=0,endfaceRejected=0,couplingRejected=0;
      G4long boundaryAbsorption=0,otherTerminal=0,unclassified=0;
      G4long activeExit=0,coreGrease=0,cladGrease=0,window=0,sipm=0;
      G4long coreToCladTracks=0,cladToCoreTracks=0,duplicates=0;
    };
    std::array<Closure,4> closure{};
    std::vector<G4int> trackIds;
    trackIds.reserve(fPhotonAudit.size());
    for(const auto& item:fPhotonAudit)trackIds.push_back(item.first);
    std::sort(trackIds.begin(),trackIds.end());
    auto* analysis=G4AnalysisManager::Instance();
    for(const auto trackId:trackIds) {
      auto& audit=fPhotonAudit.at(trackId);
      // A secondary can reach the stack already stopped at a geometry
      // boundary. Such a track is a real generated photon but never enters
      // PreUserTrackingAction/PostUserTrackingAction. It is registered by the
      // stacking action and closed here without altering transport.
      if(audit.terminal==TerminalReason::unclassified) {
        audit.terminal=TerminalReason::other_terminal;
        audit.terminalRegion=audit.lastRegion;
        audit.terminalPosition=audit.origin;
        audit.terminalTime=audit.originTime;
      }
      if(fDetector->GetRecordPhotonAuditRows()) {
        analysis->FillNtupleIColumn(3,0,eventId);
        analysis->FillNtupleIColumn(3,1,audit.track);
        analysis->FillNtupleIColumn(3,2,audit.parent);
        analysis->FillNtupleIColumn(3,3,audit.plane);
        analysis->FillNtupleIColumn(3,4,audit.originChannel);
        analysis->FillNtupleIColumn(3,5,audit.creatorProcess);
        analysis->FillNtupleIColumn(
            3,6,static_cast<G4int>(audit.originRegion));
        analysis->FillNtupleIColumn(
            3,7,static_cast<G4int>(audit.terminal));
        analysis->FillNtupleIColumn(
            3,8,static_cast<G4int>(audit.terminalRegion));
        analysis->FillNtupleIColumn(
            3,9,static_cast<G4int>(audit.stageMask));
        analysis->FillNtupleIColumn(3,10,audit.everCladding);
        analysis->FillNtupleIColumn(3,11,audit.reenteredCore);
        analysis->FillNtupleDColumn(3,12,audit.wavelength/nm);
        analysis->FillNtupleDColumn(3,13,audit.origin.x()/mm);
        analysis->FillNtupleDColumn(3,14,audit.origin.y()/mm);
        analysis->FillNtupleDColumn(3,15,audit.origin.z()/mm);
        analysis->FillNtupleDColumn(3,16,audit.initialDirection.x());
        analysis->FillNtupleDColumn(3,17,audit.initialDirection.y());
        analysis->FillNtupleDColumn(3,18,audit.initialDirection.z());
        analysis->FillNtupleDColumn(3,19,audit.activeCorePath/mm);
        analysis->FillNtupleDColumn(3,20,audit.fanoutCorePath/mm);
        analysis->FillNtupleDColumn(3,21,audit.activeInnerCladPath/mm);
        analysis->FillNtupleDColumn(3,22,audit.fanoutInnerCladPath/mm);
        analysis->FillNtupleDColumn(3,23,audit.activeOuterCladPath/mm);
        analysis->FillNtupleDColumn(3,24,audit.fanoutOuterCladPath/mm);
        analysis->FillNtupleDColumn(3,25,audit.greasePath/mm);
        analysis->FillNtupleDColumn(3,26,audit.windowPath/mm);
        analysis->FillNtupleDColumn(3,27,audit.otherPath/mm);
        analysis->FillNtupleIColumn(3,28,audit.coreToClad);
        analysis->FillNtupleIColumn(3,29,audit.cladToCore);
        analysis->FillNtupleIColumn(3,30,audit.coreGreaseAttempts);
        analysis->FillNtupleIColumn(3,31,audit.cladGreaseAttempts);
        analysis->FillNtupleIColumn(3,32,audit.greaseWindowAttempts);
        analysis->FillNtupleIColumn(3,33,audit.sipmAttempts);
        analysis->FillNtupleIColumn(3,34,audit.sipmReflections);
        analysis->FillNtupleIColumn(3,35,audit.tirCount);
        analysis->FillNtupleIColumn(
            3,36,audit.fresnelReflectionCount);
        analysis->FillNtupleIColumn(3,37,audit.otherReflectionCount);
        analysis->FillNtupleIColumn(3,38,audit.duplicateTerminalMarks);
        analysis->FillNtupleDColumn(
            3,39,audit.terminalPosition.x()/mm);
        analysis->FillNtupleDColumn(
            3,40,audit.terminalPosition.y()/mm);
        analysis->FillNtupleDColumn(
            3,41,audit.terminalPosition.z()/mm);
        analysis->FillNtupleDColumn(3,42,audit.totalTrackLength/mm);
        analysis->FillNtupleIColumn(3,43,audit.steps);
        analysis->FillNtupleDColumn(3,44,audit.originTime/ns);
        analysis->FillNtupleDColumn(3,45,audit.terminalTime/ns);
        analysis->AddNtupleRow(3);
      }
      if(fDetector->GetRecordSpatialClosureNtuples()) {
        analysis->FillNtupleIColumn(9,0,eventId);
        analysis->FillNtupleIColumn(9,1,audit.track);
        analysis->FillNtupleIColumn(9,2,audit.parent);
        analysis->FillNtupleIColumn(9,3,audit.plane);
        analysis->FillNtupleIColumn(9,4,audit.originPhysicalCopy);
        analysis->FillNtupleIColumn(9,5,audit.originLocalChannel);
        analysis->FillNtupleIColumn(9,6,audit.originChannel);
        analysis->FillNtupleIColumn(9,7,audit.terminalGlobalChannel);
        analysis->FillNtupleIColumn(9,8,audit.creatorProcess);
        analysis->FillNtupleIColumn(
            9,9,static_cast<G4int>(audit.terminal));
        analysis->FillNtupleIColumn(
            9,10,static_cast<G4int>(audit.terminalRegion));
        analysis->FillNtupleDColumn(9,11,audit.originTime/ns);
        analysis->FillNtupleDColumn(9,12,audit.terminalTime/ns);
        analysis->FillNtupleIColumn(9,13,audit.birthContextMatched);
        analysis->FillNtupleIColumn(9,14,audit.parentStep);
        analysis->AddNtupleRow(9);
      }

      if(audit.plane<0||audit.plane>=4)continue;
      auto& value=closure[audit.plane];
      ++value.registered;
      value.duplicates+=audit.duplicateTerminalMarks;
      if(audit.terminal==TerminalReason::unclassified)++value.unclassified;
      else ++value.terminal;
      switch(audit.terminal) {
        case TerminalReason::sipm_detected:++value.detected;break;
        case TerminalReason::sipm_pde_rejected:++value.rejected;break;
        case TerminalReason::bulk_core:++value.bulkCore;break;
        case TerminalReason::bulk_inner_clad:++value.bulkInner;break;
        case TerminalReason::bulk_outer_clad:++value.bulkOuter;break;
        case TerminalReason::bulk_coupling:++value.bulkCoupling;break;
        case TerminalReason::bulk_other:++value.bulkOther;break;
        case TerminalReason::world_escape:++value.world;break;
        case TerminalReason::degradation_stress:
        case TerminalReason::degradation_crack:
        case TerminalReason::degradation_packing:
        case TerminalReason::degradation_adhesive:
        case TerminalReason::degradation_h1_cladding:
        case TerminalReason::degradation_h2_cladding:
          ++value.degradation;break;
        case TerminalReason::fanout_rejected:
          ++value.fanoutRejected;break;
        case TerminalReason::endface_rejected:
          ++value.endfaceRejected;break;
        case TerminalReason::coupling_rejected:
          ++value.couplingRejected;break;
        case TerminalReason::boundary_absorption:
          ++value.boundaryAbsorption;break;
        case TerminalReason::other_terminal:
          ++value.otherTerminal;break;
        case TerminalReason::unclassified:break;
      }
      value.activeExit+=(audit.stageMask&stage_active_exit)!=0;
      value.coreGrease+=(audit.stageMask&stage_core_to_grease)!=0;
      value.cladGrease+=(audit.stageMask&stage_clad_to_grease)!=0;
      value.window+=(audit.stageMask&stage_window_entry)!=0;
      value.sipm+=(audit.stageMask&stage_sipm_boundary)!=0;
      value.coreToCladTracks+=audit.coreToClad>0;
      value.cladToCoreTracks+=audit.cladToCore>0;
    }
    for(G4int plane=0;plane<4;++plane) {
      const auto& value=closure[plane];
      analysis->FillNtupleIColumn(4,0,eventId);
      analysis->FillNtupleIColumn(4,1,plane);
      analysis->FillNtupleIColumn(4,2,value.registered);
      analysis->FillNtupleIColumn(4,3,value.terminal);
      analysis->FillNtupleIColumn(4,4,value.detected);
      analysis->FillNtupleIColumn(4,5,value.rejected);
      analysis->FillNtupleIColumn(4,6,value.bulkCore);
      analysis->FillNtupleIColumn(4,7,value.bulkInner);
      analysis->FillNtupleIColumn(4,8,value.bulkOuter);
      analysis->FillNtupleIColumn(4,9,value.bulkCoupling);
      analysis->FillNtupleIColumn(4,10,value.bulkOther);
      analysis->FillNtupleIColumn(4,11,value.world);
      analysis->FillNtupleIColumn(4,12,value.degradation);
      analysis->FillNtupleIColumn(4,13,value.fanoutRejected);
      analysis->FillNtupleIColumn(4,14,value.endfaceRejected);
      analysis->FillNtupleIColumn(4,15,value.couplingRejected);
      analysis->FillNtupleIColumn(4,16,value.boundaryAbsorption);
      analysis->FillNtupleIColumn(4,17,value.otherTerminal);
      analysis->FillNtupleIColumn(4,18,value.unclassified);
      analysis->FillNtupleIColumn(4,19,value.activeExit);
      analysis->FillNtupleIColumn(4,20,value.coreGrease);
      analysis->FillNtupleIColumn(4,21,value.cladGrease);
      analysis->FillNtupleIColumn(4,22,value.window);
      analysis->FillNtupleIColumn(4,23,value.sipm);
      analysis->FillNtupleIColumn(4,24,value.coreToCladTracks);
      analysis->FillNtupleIColumn(4,25,value.cladToCoreTracks);
      analysis->FillNtupleIColumn(4,26,value.duplicates);
      analysis->FillNtupleIColumn(
          4,27,value.terminal-value.registered);
      analysis->FillNtupleIColumn(4,28,fProducedPlane[plane]);
      analysis->FillNtupleIColumn(
          4,29,value.registered-fProducedPlane[plane]);
      analysis->AddNtupleRow(4);
    }
  }
  G4double WaveformPeak(const std::vector<G4double>& times,G4int channel) const {
    if(times.empty())return 0.;
    const G4double dt=std::max(1.e-6*ns,fDetector->GetAdcSamplePeriod());
    const G4int plane=std::clamp(channel/32,0,3);
    const G4int group=std::clamp(fDetector->GetGainBoard(plane)-1,0,3);
    if(!fTemplateTime[group].empty()) {
      const G4double start=ReadoutWindowStart();
      const G4int samples=std::max(1,fDetector->GetPeakWindowSamples());
      G4double peak=0.;
      for(G4int sample=0;sample<samples;++sample) {
        const G4double t=start+sample*dt;
        G4double amplitude=0.;
        for(const auto arrival:times)amplitude+=TemplateAmplitude(group,t-arrival);
        peak=std::max(peak,amplitude);
      }
      return peak;
    }
    const G4double sigma=std::max(1.e-6*ns,
        fDetector->GetSinglePeFwhmSamples()*dt/2.354820045);
    const G4double start=ReadoutWindowStart();
    const G4int samples=std::max(1,fDetector->GetPeakWindowSamples());
    G4double peak=0.;
    for(G4int sample=0;sample<samples;++sample) {
      const G4double t=start+sample*dt;
      G4double amplitude=0.;
      for(const auto arrival:times) {
        const G4double z=(t-arrival)/sigma;
        if(std::abs(z)<5.)amplitude+=std::exp(-.5*z*z);
      }
      peak=std::max(peak,amplitude);
    }
    return peak;
  }
  G4double WaveformPeakWithNoise(const std::vector<G4double>& times,G4int channel,
      G4double saturationScale,G4double sigma,G4double pedestalError,
      G4double& noiseOnlyPeak) const {
    const G4double dt=std::max(1.e-6*ns,fDetector->GetAdcSamplePeriod());
    const G4int plane=std::clamp(channel/32,0,3);
    const G4int group=std::clamp(fDetector->GetGainBoard(plane)-1,0,3);
    const G4double start=ReadoutWindowStart();
    const G4int samples=std::max(1,fDetector->GetPeakWindowSamples());
    G4double peak=0.;noiseOnlyPeak=0.;
    for(G4int sample=0;sample<samples;++sample) {
      const G4double t=start+sample*dt;
      G4double photoelectrons=0.;
      for(const auto arrival:times)photoelectrons+=TemplateAmplitude(group,t-arrival);
      const G4double noise=G4RandGauss::shoot(0.,sigma);
      noiseOnlyPeak=std::max(noiseOnlyPeak,noise-pedestalError);
      peak=std::max(peak,photoelectrons*saturationScale*fGain[channel]+
                          noise-pedestalError);
    }
    return std::max(0.,peak);
  }
  G4double TemplateAmplitude(G4int group,G4double relative) const {
    const auto& time=fTemplateTime[group];const auto& amplitude=fTemplateAmplitude[group];
    if(time.empty()||relative<time.front()||relative>time.back())return 0.;
    const auto upper=std::upper_bound(time.begin(),time.end(),relative);
    if(upper==time.begin())return amplitude.front();
    if(upper==time.end())return amplitude.back();
    const auto hi=std::distance(time.begin(),upper),lo=hi-1;
    const auto fraction=(relative-time[lo])/(time[hi]-time[lo]);
    return amplitude[lo]+fraction*(amplitude[hi]-amplitude[lo]);
  }
  void LoadSinglePeTemplate() {
    for(G4int group=0;group<4;++group) {
      const auto& filename=fDetector->GetSinglePeTemplateFile(group);
      std::ifstream input(filename);
      if(!input)input.open(("../"+filename).c_str());
      std::string line;G4double time=0.,amplitude=0.;
      while(std::getline(input,line)) {
        if(line.empty()||line.rfind("time_ns",0)==0)continue;
        std::replace(line.begin(),line.end(),',',' ');
        std::istringstream row(line);
        if(row>>time>>amplitude) {
          fTemplateTime[group].push_back(time*ns);
          fTemplateAmplitude[group].push_back(amplitude);
        }
      }
      G4cout<<"[optical-v21] measured single-p.e. template group="<<group+1
            <<" points="<<fTemplateTime[group].size()<<" source="<<filename<<G4endl;
      if(fDetector->GetUseWaveformDigitization() &&
         (fTemplateTime[group].empty() ||
          fTemplateTime[group].size()!=fTemplateAmplitude[group].size())) {
        G4Exception("OpticalEventAction::LoadSinglePeTemplate",
                    "HODO_V21_TEMPLATE",FatalException,
                    ("Missing or malformed single-p.e. template: "+filename).c_str());
      }
    }
  }
  void AddStage(G4int plane,G4int stage,G4int trackId) {
    if(plane<0||plane>3||stage<0||stage>3)return;
    if(!fStageSeen[plane][stage].insert(trackId).second)return;
    ++fStagePlane[plane][stage];
    if(plane<2)++fStageH1[stage];
    else ++fStageH2[stage];
  }
  void RecordSipmPhoton(G4int channel, G4int outcome, G4double wavelength,
                        G4double time) {
    if (!fDetector->GetRecordPhotonNtuple()) return;
    if(fDetector->GetRecordSelectedDiagnosticsOnly()&&!HasTriggerCoincidence())return;
    auto* analysis = G4AnalysisManager::Instance();
    auto* event = G4RunManager::GetRunManager()->GetCurrentEvent();
    analysis->FillNtupleIColumn(1, 0, event ? event->GetEventID() : -1);
    analysis->FillNtupleIColumn(1, 1, channel);
    analysis->FillNtupleIColumn(1, 2, outcome);
    analysis->FillNtupleDColumn(1, 3, wavelength/nm);
    analysis->FillNtupleDColumn(1, 4, time/ns);
    analysis->AddNtupleRow(1);
  }
  void LoadNoise() {
    // Plane order is H1X(MID1), H1Y(MID9), H2X(MID15), H2Y(MID13).
    constexpr std::array<G4double,4> fallback{4.3572,3.7715,3.7912,9.4214};
    for (G4int plane=0; plane<4; ++plane)
      std::fill(fNoiseRms.begin()+plane*32,fNoiseRms.begin()+(plane+1)*32,
                fallback[plane]);
    std::ifstream input("channel_noise_summary.csv");
    if (!input)
      input.open("../reference/experimental/run15105/channel_noise_summary.csv");
    std::string line;
    G4int loaded=0;
    while (std::getline(input,line)) {
      if (line.empty() || line.rfind("mid,",0)==0) continue;
      std::replace(line.begin(),line.end(),',',' ');
      std::istringstream row(line);
      G4int mid=0,channel=0;
      G4double rms=0.;
      if (!(row>>mid>>channel>>rms)) continue;
      const G4int plane=mid==1 ? 0 : mid==9 ? 1 : mid==15 ? 2 :
                        mid==13 ? 3 : -1;
      if (plane>=0 && channel>=1 && channel<=32 && rms>0.) {
        fNoiseRms[plane*32+channel-1]=rms;
        ++loaded;
      }
    }
    G4cout<<"[optical-v21] measured noise RMS channels loaded="<<loaded
          <<" (CONFIRMED MID1,9,15,13 -> H1X,H1Y,H2X,H2Y)"<<G4endl;
    if(fDetector->GetUseMeasuredNoise()&&loaded!=128)
      G4Exception("OpticalEventAction::LoadNoise","HODO_V21_NOISE",
                  FatalException,"Park/Yewon primary requires all 128 measured-noise channels.");
  }
  void LoadParkThresholds() {
    fParkThreshold.fill(120.);
    const std::array<std::string,2> filenames{
      "park_yewon_thresholds_mid1_mid9.csv",
      "park_yewon_thresholds_mid13_mid15.csv"
    };
    G4int loaded=0;
    G4double minimum=std::numeric_limits<G4double>::infinity();
    G4double maximum=0.,sum=0.;
    for(const auto& filename:filenames) {
      std::ifstream input(filename);
      if(!input)input.open(("../reference/experimental/run15105/"+filename).c_str());
      std::string line;
      while(std::getline(input,line)) {
        if(line.empty()||line.rfind("mid,",0)==0)continue;
        std::replace(line.begin(),line.end(),',',' ');
        std::istringstream row(line);
        G4int mid=0,channel=0;
        G4double mean=0.,rms=0.,threshold=0.;
        if(!(row>>mid>>channel>>mean>>rms>>threshold))continue;
        const G4int plane=mid==1?0:mid==9?1:mid==15?2:mid==13?3:-1;
        if(plane<0||channel<1||channel>32||threshold<=0.)continue;
        fParkThreshold[plane*32+channel-1]=threshold;
        minimum=std::min(minimum,threshold);
        maximum=std::max(maximum,threshold);
        sum+=threshold;++loaded;
      }
    }
    G4cout<<"[optical-v21] Park/Yewon per-channel mean+5RMS thresholds loaded="
          <<loaded<<" mean_ADC="<<(loaded?sum/loaded:0.)
          <<" range_ADC="<<(loaded?minimum:0.)<<"--"<<maximum<<G4endl;
    if(loaded!=128)
      G4Exception("OpticalEventAction::LoadParkThresholds",
                  "HODO_V21_PARK_THRESHOLD",FatalException,
                  "Park/Yewon primary requires all 128 channel thresholds.");
  }
  void LoadGains() {
    const auto voltage = fDetector->GetBiasVoltage();
    fGain.fill(std::abs(voltage-43.) < .25 ? 18.4
               : std::abs(voltage-44.) < .25 ? 23.4 : 29.3556);
    std::ifstream input("SiPM_channel_gain_vs_bias.dat");
    if (!input) input.open("../SiPM_channel_gain_vs_bias.dat");
    std::array<std::array<G4double,32>,4> measured{};
    for (auto& boardValues : measured) boardValues.fill(-1.);
    G4int board=0, channel=0;
    G4double gain45=0., gain44=0., gain43=0.;
    std::string line;
    G4int loaded=0;
    while (std::getline(input, line)) {
      if (line.empty() || line.front() == '#') continue;
      std::istringstream row(line);
      if (!(row >> board >> channel >> gain45 >> gain44 >> gain43)) continue;
      if (board >= 1 && board <= 4 && channel >= 1 && channel <= 32) {
        const auto gain = std::abs(voltage-43.) < .25 ? gain43
                        : std::abs(voltage-44.) < .25 ? gain44 : gain45;
        measured[board-1][channel-1] = gain;
        ++loaded;
      }
    }
    for (G4int plane=0; plane<4; ++plane) {
      const G4int mappedBoard = fDetector->GetGainBoard(plane);
      if (mappedBoard < 1 || mappedBoard > 4) {
        G4cerr << "[optical-v21] Invalid gain board " << mappedBoard
               << " for plane " << plane << "; retaining mean gains." << G4endl;
        continue;
      }
      for (G4int channel=0; channel<32; ++channel) {
        const auto gain = measured[mappedBoard-1][channel];
        if (gain > 0.) fGain[plane*32 + channel] = gain;
      }
    }
    G4cout << "[optical-v21] photo/layout gain-board mapping "
           << "H1X=" << fDetector->GetGainBoard(0)
           << " H1Y=" << fDetector->GetGainBoard(1)
           << " H2X=" << fDetector->GetGainBoard(2)
           << " H2Y=" << fDetector->GetGainBoard(3)
           << " loaded_channels=" << loaded
           << " (not established by supplied files)" << G4endl;
    if(loaded!=128)
      G4Exception("OpticalEventAction::LoadGains","HODO_V21_GAIN",
                  FatalException,"Waveform digitization requires all 128 measured gain channels.");
  }
  const DetectorConstruction* fDetector;
  G4bool fElectronicsInputsPrepared=false;
  G4double fT1=0., fT2=0., fT3=0., fT4=0.;
  G4long fProduced=0, fProducedFiber=0, fProducedTrigger=0, fEscaped=0, fPdeRejected=0;
  G4long fProducedH1=0, fProducedH2=0;
  std::array<G4long,4> fStageH1{},fStageH2{};
  std::array<std::array<G4long,4>,4> fStagePlane{};
  std::array<G4long,2> fPdeRejectedPlane{};
  std::array<G4long,4> fPdeRejectedByPlane{};
  std::array<G4double,4> fPrimaryCoreTrack{},fPrimaryActiveCoreTrack{},
      fCoreEdepPlane{};
  std::array<G4long,4> fProducedPlane{};
  std::array<G4long,4> fProducedActivePlane{},fProducedFanoutPlane{};
  std::array<G4long,4> fProducedFarEndPlane{};
  std::array<G4double,4> fPrimaryCoreEdepPlane{},fSecondaryCoreEdepPlane{};
  std::array<G4long,4> fPrimaryProducedPlane{},fSecondaryProducedPlane{};
  std::array<std::unordered_set<G4int>,4> fPrimaryCoreCopies;
  std::array<std::array<std::unordered_set<G4int>,4>,4> fStageSeen;
  G4long fBulkAbsCore=0, fBulkAbsClad=0, fBulkAbsCoupling=0, fBulkAbsOther=0;
  G4long fLossStress=0, fLossCrack=0, fLossPacking=0, fLossAdhesive=0;
  G4long fLossH1Cladding=0, fLossH2Cladding=0;
  G4long fLossFanoutExit=0, fLossEndFace=0, fLossGreaseCoupling=0;
  G4double ReadoutWindowStart() const {
    const G4double trigger=std::isfinite(fTrigger1Time)?fTrigger1Time:0.;
    return trigger+fDetector->GetAdcPeakWindowDelay();
  }
  G4double fInitialEnergy=-1., fAtTrigger1=-1.;
  G4double fTrigger1Time=std::numeric_limits<G4double>::infinity();
  G4double fEdepConcrete=0., fEdepBlackout=0., fEdepTrigger=0., fEdepH1=0., fEdepH2=0.;
  G4double fEdepH1Core=0., fEdepH2Core=0.;
  std::array<G4long,128> fDetected{};
  std::array<G4double,128> fFirstTime{};
  std::array<std::vector<G4double>,128> fArrivalTimes{};
  std::array<G4double,128> fGain{};
  std::array<G4double,128> fNoiseRms{};
  std::array<G4double,128> fParkThreshold{};
  std::array<std::vector<G4double>,4> fTemplateTime{};
  std::array<std::vector<G4double>,4> fTemplateAmplitude{};
  std::unordered_map<G4int,PhotonAudit> fPhotonAudit;
  std::vector<FibreCrossing> fFibreCrossings;
  std::unordered_map<G4int,std::size_t> fLastCrossingByTrack;
  std::unordered_map<G4int,G4int> fNextCrossingIndex;
  std::unordered_map<const G4Track*,PhotonBirthContext> fPendingPhotonBirth;
  std::array<PrimarySurfaceState,6> fPrimarySurface{};
  std::array<PrimaryTriggerCrossing,4> fPrimaryTriggerCrossing{};
  GenerationProbeRow fGenerationSource{};
  G4bool fGenerationSourceValid=false;
  std::array<GenerationPlaneLedger,4> fGenerationPlaneLedger{};
  std::array<GenerationPlaneLedger,128> fGenerationChannelLedger{};

 public:
  G4bool HasTriggerCoincidence() const {
    return fT1>=fDetector->GetHitThreshold() &&
           fT2>=fDetector->GetHitThreshold();
  }
  G4bool HasSelectedTriggerCoincidence() const {
    if(!HasTriggerCoincidence())return false;
    const auto& mode=fDetector->GetTriggerCoincidenceMode();
    if(mode=="t1t2")return true;
    if(mode=="t1t2t3")
      return fDetector->GetTriggerEnabled(2) &&
             fT3>=fDetector->GetHitThreshold();
    if(mode=="t1t2t3t4")
      return fDetector->GetTriggerEnabled(2) &&
             fDetector->GetTriggerEnabled(3) &&
             fT3>=fDetector->GetHitThreshold() &&
             fT4>=fDetector->GetHitThreshold();
    return true;
  }
};

G4ClassificationOfNewTrack OpticalLastStackingAction::ClassifyNewTrack(
    const G4Track* track) {
  if(track->GetDefinition()==G4OpticalPhoton::Definition()) {
    // G4EventManager assigns the secondary track ID and origin touchable
    // before G4StackManager invokes this callback. Register fibre photons here
    // so even already-stopped boundary secondaries remain in lineage closure.
    fEvent->BeginOpticalTrack(track);
    if(fEvent->Detector()->GetGenerationOnlyKill())
      return fKill;
    if(fEvent->Detector()->GetGate4SpectrumOnlyKill())return fKill;
    return fWaiting;
  }
  return fUrgent;
}

class OpticalTrackingAction final : public G4UserTrackingAction {
 public:
  explicit OpticalTrackingAction(OpticalEventAction* event) : fEvent(event) {}
  void PreUserTrackingAction(const G4Track* track) override {
    fEvent->BeginOpticalTrack(track);
    if(track->GetDefinition()==G4OpticalPhoton::Definition() &&
       fEvent->Detector()->GetSkipOpticalUntilT1T2() &&
       !fEvent->HasTriggerCoincidence()) {
      fEvent->MarkPhotonTerminal(track,TerminalReason::other_terminal);
      const_cast<G4Track*>(track)->SetTrackStatus(fStopAndKill);
    }
  }
  void PostUserTrackingAction(const G4Track* track) override {
    fEvent->FinalizeOpticalTrack(track);
  }
 private:
  OpticalEventAction* fEvent;
};

class OpticalSteppingAction final : public G4UserSteppingAction {
 public:
  explicit OpticalSteppingAction(OpticalEventAction* event) : fEvent(event) {}
  void UserSteppingAction(const G4Step* step) override {
    const auto* track = step->GetTrack();
    if (track->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) {
      // Gate-1 instrumentation is passive and executes before any early return
      // so zero-edep core/cladding crossings remain observable.
      fEvent->ObserveChargedStep(step);
      const auto* parentMaterial = step->GetPreStepPoint()->GetMaterial();
      const auto touch = step->GetPreStepPoint()->GetTouchableHandle();
      G4int corePlane=-1,coreCopy=-1;
      G4int globalChannel=-1;
      ResolvePlaneChannel(touch,corePlane,globalChannel);
      coreCopy=corePlane<2?globalChannel:globalChannel-64;
      GenerationSecondaryCounts generationPhotons;
      const auto parentMaterialName=parentMaterial
          ? parentMaterial->GetName() : G4String{};
      const G4bool generationFibreStep=
          fEvent->Detector()->GetGenerationProbeMode()&&
          track->GetDefinition()->GetPDGCharge()!=0.&&
          (parentMaterialName=="Kuraray_SCSF78_core"||
           parentMaterialName=="BCF20_core"||
           parentMaterialName=="PMMA_inner_cladding"||
           parentMaterialName=="Fluorinated_outer_cladding")&&
          corePlane>=0&&corePlane<4&&globalChannel>=corePlane*32&&
          globalChannel<(corePlane+1)*32;
      if (const auto* secondaries = step->GetSecondaryInCurrentStep()) {
        for (const auto* secondary : *secondaries) {
          if (secondary &&
              secondary->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
            if(generationFibreStep) {
              const auto* creator=secondary->GetCreatorProcess();
              const auto creatorName=
                  creator?creator->GetProcessName():G4String{};
              const auto nominal=
                  fEvent->Detector()->PlaneWorldToNominalPosition(
                      corePlane,secondary->GetPosition());
              const G4double longitudinal=corePlane==0||corePlane==2
                  ? nominal.x() : nominal.y();
              const G4double absoluteLongitudinal=std::abs(longitudinal);
              const G4bool active=longitudinal>=-16.*mm&&
                                  longitudinal<=16.*mm;
              if(active)
                generationPhotons.activeMaxAbsLongitudinal=std::max(
                    generationPhotons.activeMaxAbsLongitudinal,
                    absoluteLongitudinal);
              else if(generationPhotons.fanoutMinAbsLongitudinal<0.||
                      absoluteLongitudinal<
                          generationPhotons.fanoutMinAbsLongitudinal)
                generationPhotons.fanoutMinAbsLongitudinal=
                    absoluteLongitudinal;
              if(creatorName=="Scintillation") {
                ++generationPhotons.scintillation;
                if(active)++generationPhotons.scintillationActive;
                else ++generationPhotons.scintillationFanout;
              } else if(creatorName=="Cerenkov") {
                ++generationPhotons.cerenkov;
                if(active)++generationPhotons.cerenkovActive;
                else ++generationPhotons.cerenkovFanout;
              } else {
                ++generationPhotons.other;
                if(active)++generationPhotons.otherActive;
                else ++generationPhotons.otherFanout;
              }
            }
            fEvent->RegisterPhotonBirth(
                secondary,step,corePlane,globalChannel);
            fEvent->AddProducedPhoton(
                parentMaterial ? parentMaterial->GetName() : "unknown",
                corePlane,track->GetTrackID(),
                fEvent->Detector()->GetGenerationProbeMode()
                    ? secondary->GetPosition()
                    : step->GetPreStepPoint()->GetPosition());
          }
        }
      }
      fEvent->RecordGenerationStep(
          step,corePlane,globalChannel,generationPhotons);
      if (track->GetTrackID() == 1 && track->GetCurrentStepNumber() == 1)
        fEvent->SetInitialEnergy(step->GetPreStepPoint()->GetKineticEnergy());
      if (track->GetTrackID() == 1 && step->GetPostStepPoint()->GetPhysicalVolume() &&
          step->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "Trigger1")
        fEvent->SetAtTrigger1(step->GetPostStepPoint()->GetKineticEnergy(),
                             step->GetPostStepPoint()->GetGlobalTime());
      const auto edep = step->GetTotalEnergyDeposit();
      if(parentMaterial &&
         (parentMaterial->GetName()=="Kuraray_SCSF78_core"||
          parentMaterial->GetName()=="BCF20_core")) {
        G4double activeLength=0.;
        if(corePlane>=0&&corePlane<4) {
          const auto before=fEvent->Detector()->PlaneWorldToNominalPosition(
              corePlane,step->GetPreStepPoint()->GetPosition());
          const auto after=fEvent->Detector()->PlaneWorldToNominalPosition(
              corePlane,step->GetPostStepPoint()->GetPosition());
          const auto longitudinal=[](G4int plane,const G4ThreeVector& value) {
            return plane==0||plane==2 ? value.x() : value.y();
          };
          const G4double a=longitudinal(corePlane,before);
          const G4double b=longitudinal(corePlane,after);
          const auto bounds=fEvent->Detector()->GetActiveLongitudinalBounds(corePlane);
          const G4double low=bounds[0],high=bounds[1];
          G4double fraction=0.;
          const G4double delta=b-a;
          if(std::abs(delta)<=1.e-12*mm) {
            fraction=(a>=low&&a<=high)?1.:0.;
          } else {
            G4double first=(low-a)/delta;
            G4double last=(high-a)/delta;
            if(first>last)std::swap(first,last);
            const G4double enter=std::max(0.,first);
            const G4double leave=std::min(1.,last);
            fraction=std::max(0.,leave-enter);
          }
          activeLength=step->GetStepLength()*fraction;
        }
        fEvent->AddCoreTraversal(corePlane,coreCopy,step->GetStepLength(),
                                 activeLength,std::max(0.,edep),
                                 track->GetTrackID());
      }
      if (edep <= 0.) return;
      if (touch && touch->GetVolume())
        fEvent->AddMaterialEnergy(touch->GetVolume()->GetName(), edep);
      if (parentMaterial) fEvent->AddCoreEnergy(parentMaterial->GetName(), edep);
      if (touch) {
        for (G4int depth=0; depth<=touch->GetHistoryDepth(); ++depth) {
          const auto* volume = touch->GetVolume(depth);
          if (!volume) continue;
          const auto& name = volume->GetName();
          if (name == "Trigger1" || name == "Trigger2" ||
              name == "Trigger3" || name == "Trigger4") {
            fEvent->AddEnergy(name, edep);
            break;
          }
        }
      }
      return;
    }

    auto* boundary=BoundaryProcess();
    const G4double opticalMaxTime=fEvent->Detector()->GetOpticalPhotonMaxTime();
    if (opticalMaxTime>0. && track->GetGlobalTime()>opticalMaxTime) {
      fEvent->MarkPhotonTerminal(track,TerminalReason::other_terminal);
      const_cast<G4Track*>(track)->SetTrackStatus(fStopAndKill);
      return;
    }
    fEvent->ObserveOpticalStep(step,boundary);
    fEvent->TrackOpticalStage(step);
    if (fEvent->ApplyFanoutDegradation(step)) return;
    if (fEvent->ApplyV15InterfaceLoss(step)) return;
    const auto* process = step->GetPostStepPoint()->GetProcessDefinedStep();
    if (process && process->GetProcessName() == "OpAbsorption") {
      const auto* prePV = step->GetPreStepPoint()->GetPhysicalVolume();
      fEvent->AddBulkAbsorption(
          prePV ? prePV->GetName() : "unknown",track,
          step->GetPreStepPoint());
      // OpAbsorption is terminal. At a geometry endpoint Geant4 may still
      // expose World as the post-step volume; do not also classify that same
      // terminal step as a world escape.
      return;
    }
    if (step->GetPostStepPoint()->GetStepStatus() != fGeomBoundary) return;
    const auto* prePV = step->GetPreStepPoint()->GetPhysicalVolume();
    const auto* postPV = step->GetPostStepPoint()->GetPhysicalVolume();
    if (prePV && postPV && postPV->GetName() == "World" &&
        (prePV->GetName().rfind("H1",0) == 0 || prePV->GetName().rfind("H2",0) == 0)) {
      fEvent->AddEscapedPhoton(track);
      const_cast<G4Track*>(track)->SetTrackStatus(fStopAndKill);
      return;
    }
    if (!boundary) return;
    const auto post = step->GetPostStepPoint()->GetTouchableHandle();
    const G4bool atSipm=post&&post->GetVolume()&&
                        post->GetVolume()->GetName()=="SiPM";
    if(boundary->GetStatus()==Absorption&&!atSipm) {
      fEvent->MarkPhotonTerminal(
          track,TerminalReason::boundary_absorption);
      return;
    }
    if (!atSipm) return;
    const G4double wavelength = h_Planck*c_light/track->GetTotalEnergy();
    if (boundary->GetStatus() == Detection)
      fEvent->AddDetectedPhoton(
          post->GetCopyNumber(),track->GetGlobalTime(),wavelength,track);
    else if (boundary->GetStatus() == Absorption)
      fEvent->AddPdeRejectedPhoton(
          post->GetCopyNumber(),track->GetGlobalTime(),wavelength,track);
  }
 private:
  static G4OpBoundaryProcess* BoundaryProcess() {
    static G4ThreadLocal G4OpBoundaryProcess* process = nullptr;
    if (process) return process;
    auto* manager = G4OpticalPhoton::OpticalPhotonDefinition()->GetProcessManager();
    auto* vector = manager->GetPostStepProcessVector(typeDoIt);
    for (G4int i=0; i<vector->entries(); ++i) {
      process = dynamic_cast<G4OpBoundaryProcess*>((*vector)[i]);
      if (process) break;
    }
    return process;
  }
  OpticalEventAction* fEvent;
};
}

void ActionInitialization::BuildForMaster() const { SetUserAction(new OpticalRunAction(fDetector)); }

void ActionInitialization::Build() const {
  SetUserAction(new PrimaryGenerator(fDetector));
  SetUserAction(new OpticalRunAction(fDetector));
  auto* event = new OpticalEventAction(fDetector);
  SetUserAction(event);
  SetUserAction(new OpticalLastStackingAction(event));
  SetUserAction(new OpticalTrackingAction(event));
  SetUserAction(new OpticalSteppingAction(event));
}
