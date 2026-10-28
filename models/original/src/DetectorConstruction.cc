#include "DetectorConstruction.hh"
#include "H2CoreBoundaryTessellatedSolid.hh"

#include "G4Box.hh"
#include "G4ApplicationState.hh"
#include "G4Colour.hh"
#include "G4CutTubs.hh"
#include "G4DisplacedSolid.hh"
#include "G4Exception.hh"
#include "G4GenericMessenger.hh"
#include "G4GeometryTolerance.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4MultiUnion.hh"
#include "G4NistManager.hh"
#include "G4Navigator.hh"
#include "G4OpticalSurface.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SolidStore.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4TessellatedSolid.hh"
#include "G4Transform3D.hh"
#include "G4TriangularFacet.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

#include <algorithm>
#include <array>
#include <cfloat>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {
const std::vector<G4double> kPhotonEnergy = {
    2.38*eV, 2.53*eV, 2.70*eV, 2.88*eV, 3.10*eV}; // 521 ... 400 nm

void AddConstantOptics(G4Material* material, G4double rindex, G4double absorption) {
  auto* table = new G4MaterialPropertiesTable;
  table->AddProperty("RINDEX", kPhotonEnergy,
                     std::vector<G4double>(kPhotonEnergy.size(), rindex));
  table->AddProperty("ABSLENGTH", kPhotonEnergy,
                     std::vector<G4double>(kPhotonEnergy.size(), absorption));
  material->SetMaterialPropertiesTable(table);
}

void AddScintillatorOptics(G4Material* material, G4double rindex,
                           G4double attenuation, G4double yield,
                           G4double decay,
                           const std::vector<G4double>& spectrum,
                           G4double birks) {
  auto* table = new G4MaterialPropertiesTable;
  table->AddProperty("RINDEX", kPhotonEnergy,
                     std::vector<G4double>(kPhotonEnergy.size(), rindex));
  table->AddProperty("ABSLENGTH", kPhotonEnergy,
                     std::vector<G4double>(kPhotonEnergy.size(), attenuation));
  table->AddProperty("SCINTILLATIONCOMPONENT1", kPhotonEnergy, spectrum);
  table->AddConstProperty("SCINTILLATIONYIELD", yield/MeV);
  table->AddConstProperty("SCINTILLATIONTIMECONSTANT1", decay);
  table->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
  table->AddConstProperty("RESOLUTIONSCALE", 1.0);
  material->SetMaterialPropertiesTable(table);
  material->GetIonisation()->SetBirksConstant(std::max(0., birks));
}

void AddSCSF78Optics(G4Material* material, G4double attenuation, G4double yield,
                     G4double birks) {
  const std::vector<G4double> energy = {
      2.066*eV, 2.103*eV, 2.137*eV, 2.175*eV, 2.214*eV, 2.255*eV,
      2.296*eV, 2.340*eV, 2.384*eV, 2.431*eV, 2.480*eV, 2.530*eV,
      2.583*eV, 2.637*eV, 2.695*eV, 2.755*eV, 2.787*eV, 2.818*eV,
      2.850*eV, 2.883*eV, 2.952*eV, 3.024*eV};
  const std::vector<G4double> spectrum = {
      0.00, 0.01, 0.02, 0.04, 0.06, 0.10, 0.15, 0.22, 0.32, 0.42, 0.48,
      0.50, 0.52, 0.48, 0.58, 0.82, 0.95, 1.00, 0.75, 0.35, 0.05, 0.00};
  auto* table = new G4MaterialPropertiesTable;
  table->AddProperty("RINDEX", energy, std::vector<G4double>(energy.size(), 1.59));
  table->AddProperty("ABSLENGTH", energy,
                     std::vector<G4double>(energy.size(), attenuation));
  table->AddProperty("SCINTILLATIONCOMPONENT1", energy, spectrum, true);
  table->AddConstProperty("SCINTILLATIONYIELD", yield/MeV);
  table->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.8*ns);
  table->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
  table->AddConstProperty("RESOLUTIONSCALE", 1.0);
  material->SetMaterialPropertiesTable(table);
  material->GetIonisation()->SetBirksConstant(std::max(0., birks));
}

void PrintGate4MptProperty(const G4String& family,
                           const G4Material* material,
                           const G4String& property) {
  const auto* table=material ? material->GetMaterialPropertiesTable() : nullptr;
  const auto* values=table ? table->GetProperty(property) : nullptr;
  if(!values) {
    G4cout << "GATE4_MPT_MISSING family=" << family
           << " property=" << property << G4endl;
    return;
  }
  for(std::size_t index=0;index<values->GetVectorLength();++index) {
    std::ostringstream line;
    line << std::setprecision(17)
         << "GATE4_MPT family=" << family
         << " property=" << property
         << " index=" << index
         << " nodes=" << values->GetVectorLength()
         << " energy_eV=" << values->Energy(index)/eV
         << " value=" << (*values)[index];
    if(property=="ABSLENGTH")line << " unit=mm";
    else line << " unit=dimensionless";
    G4cout << line.str() << G4endl;
  }
}

struct LayoutSlot {
  G4int row = 0;
  G4int column = 0;
  G4int columns = 32;
};

LayoutSlot ResolveLayoutSlot(const G4String& layout, const G4String& order,
                             G4int requestedChannel) {
  G4int channel = std::clamp(requestedChannel, 0, 31);
  if (order == "descending") channel = 31-channel;
  const G4bool twoRows = layout != "single1x32";
  LayoutSlot slot;
  slot.columns = twoRows ? 16 : 32;
  slot.row = twoRows ? channel/16 : 0;
  slot.column = twoRows ? channel%16 : channel;
  if (twoRows && order == "serpentine" && slot.row == 1)
    slot.column = 15-slot.column;
  return slot;
}

// Gate-3 H2 outer-boundary correction, restricted to that boundary only.
// The inherited tessellation, voxelization, Inside(), DistanceToOut(), safety,
// extent, visualization, and geometry-type behavior stay untouched.  Only a
// zero directional DistanceToIn at a point already classified as kSurface is
// disambiguated by probing a bounded distance along the supplied unit ray.
class H2OuterBoundaryTessellatedSolid final
    : public G4TessellatedSolid {
 public:
  explicit H2OuterBoundaryTessellatedSolid(const G4String& name)
      : G4TessellatedSolid(name) {}

  H2OuterBoundaryTessellatedSolid(
      const H2OuterBoundaryTessellatedSolid& right)
      : G4TessellatedSolid(right) {}

  H2OuterBoundaryTessellatedSolid& operator=(
      const H2OuterBoundaryTessellatedSolid& right) {
    G4TessellatedSolid::operator=(right);
    return *this;
  }

  using G4TessellatedSolid::DistanceToIn;

  G4double DistanceToIn(const G4ThreeVector& p,
                        const G4ThreeVector& v) const override {
    const G4double baseDistance=
        G4TessellatedSolid::DistanceToIn(p,v);
    if(baseDistance!=0.||G4TessellatedSolid::Inside(p)!=kSurface)
      return baseDistance;

    G4double epsilon=std::max(
        128.*G4GeometryTolerance::GetInstance()->GetSurfaceTolerance(),
        1.*nm);
    constexpr G4int maximumAttempts=8;
    for(G4int attempt=0;attempt<maximumAttempts;++attempt) {
      const auto q=p+epsilon*v;
      const auto qState=G4TessellatedSolid::Inside(q);
      if(qState==kInside)
        return 0.;
      if(qState==kOutside) {
        const G4double reentry=
            G4TessellatedSolid::DistanceToIn(q,v);
        if(reentry>=kInfinity)
          return kInfinity;
        if(reentry>=0.&&std::isfinite(reentry))
          return epsilon+reentry;
        return baseDistance;
      }
      epsilon*=2.;
    }
    epsilon=std::max(
        128.*G4GeometryTolerance::GetInstance()->GetSurfaceTolerance(),
        1.*nm);
    for(G4int attempt=0;attempt<maximumAttempts;++attempt) {
      const auto q=p-epsilon*v;
      const auto qState=G4TessellatedSolid::Inside(q);
      if(qState==kOutside)
        return 0.;
      if(qState==kInside)
        return kInfinity;
      epsilon*=2.;
    }
    return kInfinity;
  }

  G4double DistanceToOut(const G4ThreeVector& p,
                         const G4ThreeVector& v,
                         const G4bool calcNorm=false,
                         G4bool* validNorm=nullptr,
                         G4ThreeVector* n=nullptr) const override {
    const G4double baseDistance=
        G4TessellatedSolid::DistanceToOut(p,v,calcNorm,validNorm,n);
    if(baseDistance!=0.||G4TessellatedSolid::Inside(p)!=kSurface)
      return baseDistance;

    G4double epsilon=std::max(
        128.*G4GeometryTolerance::GetInstance()->GetSurfaceTolerance(),
        1.*nm);
    constexpr G4int maximumAttempts=8;
    for(G4int attempt=0;attempt<maximumAttempts;++attempt) {
      const auto q=p+epsilon*v;
      const auto qState=G4TessellatedSolid::Inside(q);
      if(qState==kOutside)
        return baseDistance;
      if(qState==kInside) {
        const G4double interiorDistance=
            G4TessellatedSolid::DistanceToOut(
                q,v,false,nullptr,nullptr);
        if(interiorDistance>0.&&std::isfinite(interiorDistance))
          return epsilon+interiorDistance;
        if(interiorDistance>=kInfinity)
          return kInfinity;
        return epsilon;
      }
      epsilon*=2.;
    }
    return std::max(
        128.*G4GeometryTolerance::GetInstance()->GetSurfaceTolerance(),
        1.*nm);
  }

  G4VSolid* Clone() const override {
    return new H2OuterBoundaryTessellatedSolid(*this);
  }
};

}

DetectorConstruction::DetectorConstruction() {
  fFanoutChannelTransmission.fill(1.);
  fFanoutMinRadius.fill(0.);
  fMessenger = new G4GenericMessenger(this, "/hodo/", "Optical hodoscope controls");
  fMessenger->DeclarePropertyWithUnit("hitThreshold", "MeV", fHitThreshold);
  fMessenger->DeclareProperty("useBuilding", fUseBuilding,
      "Add three assumed concrete floors and a blackout sheet.");
  fMessenger->DeclarePropertyWithUnit("concreteThickness", "cm", fConcreteThickness,
      "ASSUMED thickness of each of three concrete floors.");
  fMessenger->DeclarePropertyWithUnit("blackoutThickness", "mm", fBlackoutThickness,
      "ASSUMED polyethylene-equivalent blackout-sheet thickness.");
  fMessenger->DeclareProperty("cosmicBeam", fCosmicBeam);
  fMessenger->DeclareProperty("cosmicSeed", fCosmicSeed,
      "Seed for event-indexed cosmic positions, directions, and energies.");
  fMessenger->DeclareProperty("replayEventIdOffset", fReplayEventIdOffset,
      "Nonnegative global event ID offset for disjoint serial replay chunks; "
      "shifts primary, transport RNG and truth IDs together. Default zero.")
      .SetParameterName("replayEventIdOffset",false)
      .SetRange("replayEventIdOffset>=0 && replayEventIdOffset<=1000000000")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty(
      "transportSeedEventOffset", fTransportSeedEventOffset,
      "Diagnostic-only nonnegative offset for event-indexed transport RNG "
      "replay; production must use 0.");
  fMessenger->DeclareProperty("fixedPhotonTest", fFixedPhotonTest,
      "Inject one deterministic optical photon per event inside a fibre core.");
  fMessenger->DeclareProperty("fixedPhotonPlane", fFixedPhotonPlane,
      "Fixed-photon plane: h1x, h1y, h2x, or h2y.");
  fMessenger->DeclareProperty("fixedPhotonAngularMode", fFixedPhotonAngularMode,
      "Fixed-photon direction: isotropic, forward, axis, or polar.");
  fMessenger->DeclareProperty("fixedPhotonSpectrum", fFixedPhotonSpectrum,
      "Fixed-photon replay spectrum: mono or bcf20_5point.");
  fMessenger->DeclareProperty("fixedPhotonSource", fFixedPhotonSource,
      "Gate-1 fixed-photon source: active or terminal.")
      .SetCandidates("active terminal");
  fMessenger->DeclareProperty("fixedPhotonChannel", fFixedPhotonChannel,
      "Local channel 0..31, or -1 to sweep channels 0..31 by event ID.");
  fMessenger->DeclarePropertyWithUnit("fixedPhotonWavelength", "nm",
      fFixedPhotonWavelength, "Fixed-photon wavelength.");
  fMessenger->DeclarePropertyWithUnit("fixedPhotonPolarAngle", "deg",
      fFixedPhotonPolarAngle,
      "Gate-4 polar mode angle from the positive fibre axis.");
  fMessenger->DeclarePropertyWithUnit("fixedPhotonAzimuth", "deg",
      fFixedPhotonAzimuth,
      "Gate-4 polar mode azimuth around the fibre axis.");
  fMessenger->DeclarePropertyWithUnit("fixedPhotonLongitudinalOffset", "mm",
      fFixedPhotonLongitudinalOffset,
      "Gate-4 source offset along the fibre axis from its centre.");
  fMessenger->DeclareProperty("gate4OpticalAudit", fGate4OpticalAudit,
      "Print runtime optical MPT nodes.");
  fMessenger->DeclareProperty("gate4BoundaryTrace", fGate4BoundaryTrace,
      "Print the first core-side boundary direction/status per photon.");
  fMessenger->DeclareProperty("gate4SpectrumOnlyKill", fGate4SpectrumOnlyKill,
      "Gate-4 emission fixture: register scintillation photons then kill "
      "them before transport.");
  fMessenger->DeclareProperty("cosmicEnergySpectrum", fCosmicEnergySpectrum,
      "Cosmic energy model: modified_gaisser or fixed.");
  fMessenger->DeclareProperty("cosmicAngularDistribution",
      fCosmicAngularDistribution,
      "Cosmic zenith model: cos2, flat_plane_cos2, cos_power, or "
      "uniform_solid_angle. modified_gaisser_joint samples E and theta "
      "from the horizontal-plane modified-Gaisser flux.")
      .SetCandidates("cos2 flat_plane_cos2 cos_power uniform_solid_angle modified_gaisser_joint");
  fMessenger->DeclareProperty("cosmicZenithPower", fCosmicZenithPower,
      "Exponent n for cos_power zenith sampling: "
      "p(cos(theta)) is proportional to cos^n(theta). "
      "cos2 uses n=2 regardless of this value.");
  fMessenger->DeclarePropertyWithUnit("cosmicMinEnergy", "GeV", fCosmicMinEnergy,
      "Lower kinetic-energy bound for modified_gaisser sampling.");
  fMessenger->DeclarePropertyWithUnit("cosmicMaxEnergy", "GeV", fCosmicMaxEnergy,
      "Upper kinetic-energy bound for modified_gaisser sampling.");
  fMessenger->DeclareProperty("fanoutPreview", fFanoutPreview,
      "Use the photo-derived H1/H2 4x8 fan-out as the detector geometry.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("fibreFarEndExtension", "mm",
      fFibreFarEndExtension,
      "Extra continuous fibre length beyond the active boundary, opposite the fan-out, "
      "for all four planes. Default 0 preserves historical runs.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("activePlaneOverlap",fActivePlaneOverlap,
      "Active is the X/Y packed-plane overlap, for zero individual plane yaw.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("fanoutLinearFanout", fFanoutLinearFanout,
      "Use the original full fan-out geometry with linear interpolation "
      "from the active packed fibres to the SiPM-side grid. Default false "
      "keeps the smootherstep fan-out.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("fanoutPreviewLite", fFanoutPreviewLite,
      "Visualization-only fan-out preview: representative simple fibres, "
      "not physics-production geometry.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("fanoutPreviewLiteChannelStride",
      fFanoutPreviewLiteChannelStride,
      "Visualization-only channel stride for fanoutPreviewLite. "
      "Channels divisible by this value and channel 31 are drawn.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("fanoutPreviewLiteSegments",
      fFanoutPreviewLiteSegments,
      "Visualization-only number of curved fan-out segments per fibre.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("fibreOverlapClearance", "mm",
      fFibreOverlapClearance,
      "Optional fan-out construction clearance subtracted from the modelled "
      "outer fibre half-width/radius. Default 0 preserves the predecessor; "
      "structural overlap diagnostics use at most 0.015 mm.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("maxZenith", "deg", fMaxZenith);
  fMessenger->DeclarePropertyWithUnit("cosmicTargetHalfX", "mm", fCosmicTargetHalfX,
      "Half-width of uniform cosmic target plane at Trigger1 z.");
  fMessenger->DeclarePropertyWithUnit("cosmicTargetHalfY", "mm", fCosmicTargetHalfY,
      "Half-height of uniform cosmic target plane at Trigger1 z.");
  fMessenger->DeclareProperty("cosmicPositionModel", fCosmicPositionModel,
      "target_backproject (legacy) or uniform_source_plane (position independent of direction).")
      .SetCandidates("target_backproject uniform_source_plane").SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("cosmicSourceHalfX", "mm", fCosmicSourceHalfX)
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("cosmicSourceHalfY", "mm", fCosmicSourceHalfY)
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("unselectedEventStride", fUnselectedEventStride,
      "Keep all trigger-selected events plus every Nth unselected event; default 1 keeps all.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("skipUntriggeredDigitization", fSkipUntriggeredDigitization,
      "Skip ADC/noise for events failing the energy-deposition trigger. Changes downstream RNG history, not the trigger definition.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger1Yaw", "deg", fTriggerYaw[0],
      "Legacy sign: positive command subtracts from Trigger1's nominal zero "
      "degree global-Z placement angle.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger2Yaw", "deg", fTriggerYaw[1],
      "Legacy sign: positive command subtracts from Trigger2's nominal 90 "
      "degree global-Z placement angle.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger3Yaw", "deg", fTriggerYaw[2],
      "Optional Trigger3 global-Z yaw for trigger-acceptance studies.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger4Yaw", "deg", fTriggerYaw[3],
      "Optional Trigger4 global-Z yaw for trigger-acceptance studies.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("h1xYaw", "deg", fPlaneYaw[0],
      "V15 global-Z yaw delta about the H1X active-plane centre.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("h1yYaw", "deg", fPlaneYaw[1],
      "V15 global-Z yaw delta about the H1Y active-plane centre.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("h2xYaw", "deg", fPlaneYaw[2],
      "V15 global-Z yaw delta about the H2X active-plane centre.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("h2yYaw", "deg", fPlaneYaw[3],
      "V15 global-Z yaw delta about the H2Y active-plane centre.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("h2ModuleRotate180", fH2ModuleRotate180,
      "Structural layout option: rotate the whole H2 module by 180 degrees "
      "about global Z while keeping H2 yaw-delta commands at zero.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("h1OffsetX", "mm",
      fHodoscopeOffsetX[0],
      "V21/V16-compatible rigid H1 module global-X translation.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("h1OffsetY", "mm",
      fHodoscopeOffsetY[0],
      "V21/V16-compatible rigid H1 module global-Y translation.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("h2OffsetX", "mm",
      fHodoscopeOffsetX[1],
      "V21/V16-compatible rigid H2 module global-X translation.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("h2OffsetY", "mm",
      fHodoscopeOffsetY[1],
      "V21/V16-compatible rigid H2 module global-Y translation.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("targetCoincidences", fTargetCoincidences);
  fMessenger->DeclareProperty("triggerCoincidenceMode",
      fTriggerCoincidenceMode,
      "Run-stop/progress denominator: t1t2, t1t2t3, or t1t2t3t4.")
      .SetCandidates("t1t2 t1t2t3 t1t2t3t4");
  fMessenger->DeclarePropertyWithUnit("triggerLength", "mm",
      fTriggerLength,
      "Trigger scintillator full local-X length; converted to a G4Box "
      "half-length internally.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("triggerWidth", "mm",
      fTriggerWidth,
      "Trigger scintillator full local-Y width; converted to a G4Box "
      "half-width internally.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("triggerThickness", "mm",
      fTriggerThickness,
      "Trigger scintillator full local-Z thickness; converted to a G4Box "
      "half-thickness internally.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger1Z", "mm", fTriggerZ[0],
      "Structural-acceptance Trigger 1 centre Z position.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger2Z", "mm", fTriggerZ[1],
      "Structural-acceptance Trigger 2 centre Z position.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("trigger3Enabled", fTriggerEnabled[2],
      "Enable optional Trigger3 volume for trigger-acceptance studies.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("trigger4Enabled", fTriggerEnabled[3],
      "Enable optional Trigger4 volume for trigger-acceptance studies.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger3Z", "mm", fTriggerZ[2],
      "Optional Trigger3 centre Z position.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger4Z", "mm", fTriggerZ[3],
      "Optional Trigger4 centre Z position.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger3Length", "mm",
      fTrigger3Length, "Optional Trigger3 full local-X length.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger3Width", "mm",
      fTrigger3Width, "Optional Trigger3 full local-Y width.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger3Thickness", "mm",
      fTrigger3Thickness, "Optional Trigger3 full local-Z thickness.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger4Length", "mm",
      fTrigger4Length, "Optional Trigger4 full local-X length.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger4Width", "mm",
      fTrigger4Width, "Optional Trigger4 full local-Y width.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger4Thickness", "mm",
      fTrigger4Thickness, "Optional Trigger4 full local-Z thickness.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("alignmentReference", fAlignmentReference,
      "Alignment reference frame: legacy or ideal_axis.")
      .SetCandidates("legacy ideal_axis")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger1OffsetX", "mm",
      fTriggerOffsetX[0],
      "Trigger 1 global-X translation relative to alignmentReference.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger1OffsetY", "mm",
      fTriggerOffsetY[0],
      "Trigger 1 global-Y translation relative to alignmentReference.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger2OffsetX", "mm",
      fTriggerOffsetX[1],
      "Trigger 2 global-X translation relative to alignmentReference.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger2OffsetY", "mm",
      fTriggerOffsetY[1],
      "Trigger 2 global-Y translation relative to alignmentReference.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger3OffsetX", "mm",
      fTriggerOffsetX[2],
      "Optional Trigger3 global-X translation.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger3OffsetY", "mm",
      fTriggerOffsetY[2],
      "Optional Trigger3 global-Y translation.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger4OffsetX", "mm",
      fTriggerOffsetX[3],
      "Optional Trigger4 global-X translation.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("trigger4OffsetY", "mm",
      fTriggerOffsetY[3],
      "Optional Trigger4 global-Y translation.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("h1FanoutZ", "mm",
      fHodoscopeFanoutZ[0],
      "Structural-acceptance H1 fan-out preview module centre Z position.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("h2FanoutZ", "mm",
      fHodoscopeFanoutZ[1],
      "Structural-acceptance H2 fan-out preview module centre Z position.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("virtualRoundModelEnabled",
      fVirtualRoundModelEnabled,
      "Enable analysis-only all-round hodoscope path oracle in run summary.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("virtualRoundH3Enabled",
      fVirtualRoundH3Enabled,
      "Enable virtual third round hodoscope station in the path oracle.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("virtualRoundH3Z", "mm",
      fVirtualRoundH3Z,
      "Virtual H3 round module reference Z position for the path oracle.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("virtualRoundH3OffsetX", "mm",
      fVirtualRoundH3OffsetX,
      "Virtual H3 global-X offset for the path oracle.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("virtualRoundH3OffsetY", "mm",
      fVirtualRoundH3OffsetY,
      "Virtual H3 global-Y offset for the path oracle.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("virtualRoundH3Yaw", "deg",
      fVirtualRoundH3Yaw,
      "Virtual H3 yaw angle for the path oracle.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("h2PhysicalRound", fH2PhysicalRound,
      "Construct H2 as physical round fibres instead of square fibres. "
      "Default false preserves the validated mixed detector geometry.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("h2SquareFibreCount", fH2SquareFibreCount,
      "Physical fibres per H2 square plane: 29 or legacy 32. "
      "29 uses centred single-row fibres on provisional local channels 0..28.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("h1AdcThreshold", fH1AdcThreshold);
  fMessenger->DeclareProperty("h2AdcThreshold", fH2AdcThreshold);
  fMessenger->DeclareProperty("fiberSurfaceRoughness", fFiberSurfaceRoughness,
      "UNIFIED sigma_alpha in radians; 0 is ideally polished.");
  fMessenger->DeclareProperty("endFaceRoughness", fEndFaceRoughness,
      "V15 UNIFIED sigma_alpha at the sanded fibre-end/grease boundary.");
  fMessenger->DeclareProperty("farEndReflectivity", fFarEndReflectivity);
  fMessenger->DeclareProperty("couplingEfficiency", fCouplingEfficiency,
      "Residual contact-loss multiplier after explicit V15 optical interfaces.");
  fMessenger->DeclarePropertyWithUnit("greaseGap", "mm", fGreaseGap,
      "V15 axial fibre-end to SiPM-window gap, completely filled with grease.");
  fMessenger->DeclarePropertyWithUnit("sipmActiveSize", "mm", fSiPMActiveSize,
      "V15 square SiPM optical active-side length.");
  fMessenger->DeclarePropertyWithUnit("sipmOffsetU", "mm", fSiPMOffsetU,
      "V15 systematic SiPM lateral offset in fan-out transverse direction.");
  fMessenger->DeclarePropertyWithUnit("sipmOffsetV", "mm", fSiPMOffsetV,
      "V15 systematic SiPM lateral offset in fan-out Z direction.");
  fMessenger->DeclareProperty("bendTransmission", fBendTransmission,
      "Transmission through the photographed fibre fan-out and bends.");
  fMessenger->DeclareProperty("fanoutReferenceTransmissionH1",
      fFanoutReferenceTransmissionH1,
      "V15 sensitivity: H1 survival at the reference bend radius.");
  fMessenger->DeclareProperty("fanoutReferenceTransmissionH2",
      fFanoutReferenceTransmissionH2,
      "V15 sensitivity: H2 survival at the reference bend radius.");
  fMessenger->DeclarePropertyWithUnit("fanoutReferenceRadius", "mm",
      fFanoutReferenceRadius,
      "Reference radius for curvature-dependent fan-out survival.");
  fMessenger->DeclareProperty("fanoutCurvaturePower", fFanoutCurvaturePower,
      "Exponent controlling curvature dependence of channel transmission.");
  fMessenger->DeclareProperty("fanoutChannelRms", fFanoutChannelRms,
      "Log-normal channel-to-channel RMS sensitivity for packing/micro-bends.");
  fMessenger->DeclareProperty("pdeScale", fPDEScale,
      "Scale factor applied to the provisional wavelength-dependent PDE.");
  fMessenger->DeclareProperty("scintillationYieldScale", fScintillationYieldScale,
      "Sensitivity scale on manufacturer light yield; 1 keeps the nominal value.");
  fMessenger->DeclareProperty("effectiveLightScaleH1",
      fEffectiveLightScaleH1,
      "Post-transport/PDE unresolved effective light scale for H1 readout; 1 applies no extra reduction.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("effectiveLightScaleH2",
      fEffectiveLightScaleH2,
      "Post-transport/PDE unresolved effective light scale for H2 readout; 1 applies no extra reduction.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("h1BirksConstant", fH1BirksConstant,
      "SENSITIVITY: H1 SCSF-78 Birks coefficient in mm/MeV; generic PS default.");
  fMessenger->DeclareProperty("h2BirksConstant", fH2BirksConstant,
      "SENSITIVITY: H2 BCF-20 Birks coefficient in mm/MeV; generic PS default.");
  fMessenger->DeclareProperty("triggerBirksConstant", fTriggerBirksConstant,
      "SENSITIVITY: trigger Birks coefficient in mm/MeV; material unresolved.");
  fMessenger->DeclarePropertyWithUnit("microbendAmplitude", "mm",
      fMicrobendAmplitude, "ASSUMED sinusoidal centreline perturbation.");
  fMessenger->DeclarePropertyWithUnit("microbendPeriod", "mm",
      fMicrobendPeriod, "ASSUMED micro-bend wavelength.");
  fMessenger->DeclareProperty("bendOvality", fBendOvality,
      "ASSUMED maximum fractional bend-region cross-section ovality.");
  fMessenger->DeclareProperty("claddingAsymmetry", fCladdingAsymmetry,
      "ASSUMED 0..1 core/cladding offset fraction in bends.");
  fMessenger->DeclareProperty("stressLossPerM",
      fStressLossPerM, "ASSUMED bend-stress photon-loss coefficient.");
  fMessenger->DeclareProperty("crackLossPerM",
      fCrackLossPerM, "ASSUMED micro-crack photon-loss coefficient.");
  fMessenger->DeclareProperty("packingLossPerM",
      fPackingLossPerM, "ASSUMED packing-pressure photon-loss coefficient.");
  fMessenger->DeclarePropertyWithUnit("adhesiveIntrusionLength", "mm",
      fAdhesiveIntrusionLength, "ASSUMED EJ-550 intrusion length from fibre end.");
  fMessenger->DeclareProperty("adhesiveLossPerM",
      fAdhesiveLossPerM, "ASSUMED loss coefficient in adhesive-intruded end.");
  fMessenger->DeclareProperty("h1CladdingDamagePerM", fH1CladdingDamagePerM,
      "SENSITIVITY effective H1 cladding breach/leakage coefficient per metre.");
  fMessenger->DeclareProperty("h2CladdingDamagePerM", fH2CladdingDamagePerM,
      "SENSITIVITY effective H2 cladding breach/leakage coefficient per metre.");
  fMessenger->DeclareProperty("endFaceEfficiency", fEndFaceEfficiency,
      "ASSUMED fibre-end to grease survival probability, applied before SiPM PDE.");
  fMessenger->DeclareProperty("crossTalkProbability", fCrossTalkProbability);
  fMessenger->DeclareProperty("afterPulseProbability", fAfterPulseProbability);
  fMessenger->DeclarePropertyWithUnit("darkCountRate", "Hz", fDarkCountRate);
  fMessenger->DeclarePropertyWithUnit("readoutGate", "ns", fReadoutGate);
  fMessenger->DeclareProperty("electronicsNoise", fElectronicsNoise,
      "Gaussian electronics noise in photoelectron-equivalent units RMS.");
  fMessenger->DeclareProperty("useMeasuredNoise", fUseMeasuredNoise,
      "Use run15105 channel pedestal RMS and peak-window maximum model.");
  fMessenger->DeclareProperty("adcNoiseScale", fAdcNoiseScale,
      "Sensitivity multiplier on measured per-channel ADC noise RMS.");
  fMessenger->DeclareProperty("pedestalWindowSamples", fPedestalWindowSamples);
  fMessenger->DeclareProperty("peakWindowSamples", fPeakWindowSamples);
  fMessenger->DeclareProperty("useWaveformDigitization", fUseWaveformDigitization,
      "Build peak ADC by summing time-shifted single-p.e. pulses.");
  fMessenger->DeclarePropertyWithUnit("adcSamplePeriod", "ns", fAdcSamplePeriod,
      "Nominal fixed DAQ sample period; recovered DRS source uses about 0.2 ns/bin.");
  fMessenger->DeclareProperty("singlePeFwhmSamples", fSinglePeFwhmSamples,
      "PROVISIONAL normalized pulse FWHM in samples; raw cosmic average is about 80.");
  fMessenger->DeclarePropertyWithUnit("adcPeakWindowDelay", "ns", fAdcPeakWindowDelay,
      "Fixed delay from Trigger1 crossing to the first sample of the peak window.");
  fMessenger->DeclarePropertyWithUnit("xyCoincidenceWindow", "ns",
      fXYCoincidenceWindow,
      "V11 diagnostic maximum X/Y hit-channel first-arrival time difference.");
  fMessenger->DeclareProperty("singlePeTemplateFile1", fSinglePeTemplateFile[0]);
  fMessenger->DeclareProperty("singlePeTemplateFile2", fSinglePeTemplateFile[1]);
  fMessenger->DeclareProperty("singlePeTemplateFile3", fSinglePeTemplateFile[2]);
  fMessenger->DeclareProperty("singlePeTemplateFile4", fSinglePeTemplateFile[3]);
  fMessenger->DeclareProperty("recordPhotonNtuple", fRecordPhotonNtuple,
      "Write per-photon ROOT rows; disable for compact sensitivity scans.");
  fMessenger->DeclareProperty("recordChannelNtuple", fRecordChannelNtuple,
      "Write per-event, per-channel ADC/noise diagnostics (128 rows/event).");
  fMessenger->DeclareProperty("recordSelectedDiagnosticsOnly", fRecordSelectedDiagnosticsOnly,
      "Write long diagnostics only for T1&&T2 events; retain all event rows and physics.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("recordPhotonAuditNtuple",
      fRecordPhotonAuditNtuple,
      "Track V21 per-fibre-photon lineage and write per-plane loss_summary.");
  fMessenger->DeclareProperty("recordPhotonAuditRows",
      fRecordPhotonAuditRows,
      "Write photon_audit rows when lineage tracking is enabled; disable only "
      "these rows for compact production statistics.");
  fMessenger->DeclareProperty("recordActiveCoreChordBranches",
      fRecordActiveCoreChordBranches,
      "Append the exact V16 four-plane primary active-straight chord fields. "
      "Legacy mid*_primary_core_track_mm retains active+fan-out total.");
  fMessenger->DeclareProperty("recordSpatialClosureNtuples",
      fRecordSpatialClosureNtuples,
      "Write passive Gate-1 surface, fibre-crossing, photon-endpoint, "
      "event-plane, and extended geometry-map relational diagnostics.");
  fMessenger->DeclareProperty("pathOnlyMuonTruth", fPathOnlyMuonTruth,
      "Allow compact structural-acceptance runs to write primary muon "
      "trigger/fibre truth while suppressing photon/channel/audit rows.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("recordGeometryControlNtuples",
      fRecordGeometryControlNtuples,
      "Append passive Stage-2 trigger-volume and trigger-map geometry truth; "
      "requires recordSpatialClosureNtuples=true.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("geometryProbeMode", fGeometryProbeMode,
      "Use a deterministic event-indexed geantino geometry-probe table.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("geometryProbeFile", fGeometryProbeFile,
      "Stage-2 deterministic geometry-probe TSV input.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("geometryProbeRowOffset", fGeometryProbeRowOffset,
      "First TSV row used as geometry-probe event zero.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("generationProbeMode", fGenerationProbeMode,
      "Stage-3 deterministic event-indexed 4 GeV muon source table.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("generationProbeFile", fGenerationProbeFile,
      "Stage-3 deterministic generation-probe TSV input.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("generationProbeSha256",
      fGenerationProbeSha256,
      "Expected SHA-256 of the exact Stage-3 generation-probe TSV.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("generationProbeRowOffset",
      fGenerationProbeRowOffset,
      "First TSV row used as Stage-3 generation-probe event zero.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("recordGenerationNtuples",
      fRecordGenerationNtuples,
      "Append Stage-3 charged-step, photon-birth, and provenance ntuples.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("fibreMaxStep", "mm", fFibreMaxStep,
      "Stage-3 maximum charged-particle step in every fibre core and "
      "cladding logical volume; 0 disables the limit.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("generationOnlyKill", fGenerationOnlyKill,
      "Stage-3 optimization: classify generated optical secondaries, then "
      "kill their transport without changing charged-particle transport.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("opticalPhotonMaxTime", "ns",
      fOpticalPhotonMaxTime,
      "Production safety valve: kill optical photons after this global time; "
      "0 disables the limit and preserves validation defaults.");
  fMessenger->DeclareProperty("skipOpticalUnlessT1T2", fSkipOpticalUntilT1T2,
      "Production optimization for T1&&T2 conditional studies: skip optical "
      "transport in events that do not pass trigger energy threshold.");
  fMessenger->DeclareProperty("gate6ElectronicsProbeMode",
      fGate6ElectronicsProbeMode,
      "Gate-6 default-OFF exact-p.e. injection fixture for digitizer tests.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("gate6InjectedChannel", fGate6InjectedChannel,
      "Gate-6 global channel; -1 cycles event modulo 128.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("gate6InjectedPe", fGate6InjectedPe,
      "Gate-6 exact primary photoelectrons injected into one channel.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("gate6InjectedStartTime", "ns",
      fGate6InjectedStartTime,
      "Gate-6 first injected primary-p.e. arrival time.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("gate6InjectedSpacing", "ns",
      fGate6InjectedSpacing,
      "Gate-6 constant spacing between injected primary p.e.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("microcellCount", fMicrocellCount);
  fMessenger->DeclareProperty("biasVoltage", fBiasVoltage,
      "Measured gain-table voltage: 43, 44, or 45 V.");
  fMessenger->DeclareProperty("sipmModel", fSiPMModel,
      "CONFIRMED Hamamatsu S14160-1310PS.");
  fMessenger->DeclarePropertyWithUnit("h1AttenuationLength", "m", fH1AttenuationLength,
      "SCSF-78 flat attenuation sensitivity: compare 4.0 m and 4.67 m.");
  fMessenger->DeclarePropertyWithUnit("h2AttenuationLength", "m", fH2AttenuationLength,
      "BCF-20 nominal attenuation length; independently configurable.");
  fMessenger->DeclareProperty("h1ScintillationYield", fH1ScintillationYield,
      "SCSF-78 nominal photons/MeV before the legacy global yield scale.");
  fMessenger->DeclareProperty("h2ScintillationYield", fH2ScintillationYield,
      "BCF-20 nominal photons/MeV before the legacy global yield scale.");
  fMessenger->DeclareProperty("triggerScintillationYield", fTriggerScintillationYield,
      "Trigger paddle photons/MeV. Default preserves validation behaviour; "
      "set to 0 for production speed when T1/T2 use edep only.");
  fMessenger->DeclareProperty("h1CladdingModel", fH1CladdingModel,
      "H1 construction sensitivity: single or double.");
  fMessenger->DeclareProperty("h2CladdingModel", fH2CladdingModel,
      "H2 BCF-20 square construction sensitivity: single or double.");
  fMessenger->DeclareProperty("h1EmissionModel", fH1EmissionModel,
      "H1 emission sensitivity: legacy5 or scsf78_22.");
  fMessenger->DeclareProperty("roundTessellatedSides", fH1RoundSides,
      "V21 round-fibre tessellation sensitivity: validated 16 or 64; "
      "production candidate is 64.");
  fMessenger->DeclareProperty(
      "roundActiveFacetRefinement", fH1ActiveFacetRefinement,
      "V21 mesh-only OAT: conforming barycentric subdivision factor "
      "1, 2, 4, or 8 for each existing H1 active-span triangle.");
  fMessenger->DeclareProperty("roundCladdingHierarchy", fH1FanoutHierarchy,
      "V21 round-cladding hierarchy: nested or sibling_shells; production "
      "candidate is sibling_shells.");
  fMessenger->DeclareProperty("squareCladdingHierarchy", fH2FanoutHierarchy,
      "V21 square-cladding hierarchy: nested or sibling_shells. Production "
      "must explicitly select sibling_shells after the H2 navigation gate.")
      .SetCandidates("nested sibling_shells");
  fMessenger->DeclareProperty(
      "roundSweepRepresentation", fH1RoundSweepRepresentation,
      "V21 geometry-only OAT: tessellated, native_miter_union, or "
      "native_capsule_union; native_miter_subtraction_union is a Boolean "
      "control, native_miter_public_union is a local-normal control, and "
      "native_balanced_union is the public G4UnionSolid candidate. Native "
      "modes apply only to H1 round fibres.")
      .SetCandidates(
          "tessellated native_miter_union native_capsule_union "
          "native_miter_subtraction_union native_miter_public_union "
          "native_balanced_union");
  fMessenger->DeclarePropertyWithUnit(
      "roundMiterJointOverlap", "nm", fH1RoundMiterJointOverlap,
      "NUMERICAL OAT only: one-sided axial extension at each internal "
      "G4CutTubs miter end; primary default 0 nm, not an efficiency fit.");
  fMessenger->DeclareProperty("h1RoundSides", fH1RoundSides,
      "Compatibility alias for roundTessellatedSides.");
  fMessenger->DeclareProperty("h1FanoutHierarchy", fH1FanoutHierarchy,
      "Compatibility alias for roundCladdingHierarchy.");
  fMessenger->DeclareProperty("h2FanoutHierarchy", fH2FanoutHierarchy,
      "Compatibility alias for squareCladdingHierarchy.");
  fMessenger->DeclareMethod(
      "runNavigatorSolidProbe",
      &DetectorConstruction::RunNavigatorSolidProbe,
      "V21 GEOMETRY OAT ONLY: after /run/initialize, replay a frozen TSV "
      "directly against named G4VSolid DistanceToIn/Inside methods.")
      .SetParameterName("path",false);
  fMessenger->DeclareProperty(
      "navigatorSolidProbeMode", fNavigatorSolidProbeMode,
      "Fail-closed direct-solid mode: frozen_baseline preserves the raw "
      "captured ray; candidate_gate constructs an Outside-to-Inside query.")
      .SetCandidates("frozen_baseline candidate_gate");
  fMessenger->DeclareProperty(
      "navigatorSolidProbeExpectedCondition",
      fNavigatorSolidProbeExpectedCondition,
      "Exact condition label required on every TSV data row.");
  fMessenger->DeclareProperty(
      "navigatorSolidProbeExpectedRows", fNavigatorSolidProbeExpectedRows,
      "Exact non-zero TSV data-row count required by the probe.");
  fMessenger->DeclareProperty("h1xLayout", fPlaneLayout[0],
      "ASSUMED layout sensitivity: single1x32, 2x16_aligned, or "
      "2x16_staggered.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("h1yLayout", fPlaneLayout[1],
      "ASSUMED layout sensitivity: single1x32, 2x16_aligned, or "
      "2x16_staggered.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("h2xLayout", fPlaneLayout[2],
      "ASSUMED layout sensitivity: single1x32, 2x16_aligned, or "
      "2x16_staggered.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("h2yLayout", fPlaneLayout[3],
      "ASSUMED layout sensitivity: single1x32, 2x16_aligned, or "
      "2x16_staggered.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("h1xChannelOrder", fPlaneChannelOrder[0],
      "ASSUMED channel-order sensitivity: ascending, serpentine, descending.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("h1yChannelOrder", fPlaneChannelOrder[1],
      "ASSUMED channel-order sensitivity: ascending, serpentine, descending.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("h2xChannelOrder", fPlaneChannelOrder[2],
      "ASSUMED channel-order sensitivity: ascending, serpentine, descending.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("h2yChannelOrder", fPlaneChannelOrder[3],
      "ASSUMED channel-order sensitivity: ascending, serpentine, descending.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("h1xRowZSpacing", "mm",
      fPlaneRowZSpacing[0], "ASSUMED row-Z spacing; unused by single1x32.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("h1yRowZSpacing", "mm",
      fPlaneRowZSpacing[1], "ASSUMED row-Z spacing; unused by single1x32.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("h2xRowZSpacing", "mm",
      fPlaneRowZSpacing[2], "ASSUMED row-Z spacing; unused by single1x32.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclarePropertyWithUnit("h2yRowZSpacing", "mm",
      fPlaneRowZSpacing[3], "ASSUMED row-Z spacing; unused by single1x32.")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareMethod("layoutPreset",
      &DetectorConstruction::SetLayoutPreset,
      "Set an exact source-level layout preset: mixed_assumed or "
      "v18d_compat. Process default is v18d_compat; plane-specific commands "
      "may override either preset afterwards.")
      .SetCandidates("mixed_assumed v18d_compat")
      .SetStates(G4State_PreInit);
  fMessenger->DeclareProperty("diagnosticRootFile", fDiagnosticRootFile,
      "ROOT diagnostic output basename; set before /run/initialize.");
  fMessenger->DeclareProperty("gainBoardH1X", fGainBoard[0],
      "PROVISIONAL photo/layout mapping: round H1 X uses hodoscope/group 3.");
  fMessenger->DeclareProperty("gainBoardH1Y", fGainBoard[1],
      "PROVISIONAL photo/layout mapping: round H1 Y uses hodoscope/group 1.");
  fMessenger->DeclareProperty("gainBoardH2X", fGainBoard[2],
      "PROVISIONAL photo/layout mapping: square H2 X uses hodoscope/group 2.");
  fMessenger->DeclareProperty("gainBoardH2Y", fGainBoard[3],
      "PHOTO-CONFIRMED mapping: square H2 Y uses hodoscope/group 4.");

  // Keep the predecessor messenger allocation sequence unchanged. Alignment
  // scans are independent processes configured by one explicit environment
  // record: reference,t1x,t1y,t1yaw,h1x,h1y,h1yaw,h2x,h2y,h2yaw,
  // t2x,t2y,t2yaw (mm,deg). A missing variable takes this zero-allocation
  // fast path and therefore preserves the approved default run.
  const char* alignmentSpec=std::getenv("HODO_ALIGNMENT_CONFIG");
  if(alignmentSpec&&*alignmentSpec) {
    std::vector<G4String> token;
    std::istringstream parser(alignmentSpec);
    G4String value;
    while(std::getline(parser,value,','))token.push_back(value);
    if(token.size()!=13)
      G4Exception("DetectorConstruction::DetectorConstruction",
                  "HODO_ALIGNMENT_ENV",FatalException,
                  "HODO_ALIGNMENT_CONFIG needs 13 comma-separated fields.");
    try {
      fAlignmentReference=token[0];
      const auto number=[&token](std::size_t index) {
        return std::stod(std::string(token[index]));
      };
      fTriggerOffsetX[0]=number(1)*mm;
      fTriggerOffsetY[0]=number(2)*mm;
      fTriggerYaw[0]=number(3)*deg;
      fHodoscopeOffsetX[0]=number(4)*mm;
      fHodoscopeOffsetY[0]=number(5)*mm;
      fPlaneYaw[0]=fPlaneYaw[1]=number(6)*deg;
      fHodoscopeOffsetX[1]=number(7)*mm;
      fHodoscopeOffsetY[1]=number(8)*mm;
      fPlaneYaw[2]=fPlaneYaw[3]=number(9)*deg;
      fTriggerOffsetX[1]=number(10)*mm;
      fTriggerOffsetY[1]=number(11)*mm;
      fTriggerYaw[1]=number(12)*deg;
      fRecordAlignmentStudyNtuple=true;
    } catch(const std::exception&) {
      G4Exception("DetectorConstruction::DetectorConstruction",
                  "HODO_ALIGNMENT_ENV_VALUE",FatalException,
                  "HODO_ALIGNMENT_CONFIG contains a non-numeric field.");
    }
  }
}

DetectorConstruction::~DetectorConstruction() {
  delete fFibreUserLimits;
  delete fMessenger;
}

void DetectorConstruction::ApplyFibreUserLimits(
    G4LogicalVolume* logical) const {
  if (logical && fFibreUserLimits) logical->SetUserLimits(fFibreUserLimits);
}

G4double DetectorConstruction::GetH1RoundPolygonPhase() const {
  // V21 polygon-phase OAT: only the production 64-gon is rotated so that
  // frame axes cross facet centres instead of polygon vertices.  Preserve
  // the historical 16-gon sensitivity and every square-fibre coordinate.
  return fH1RoundSides == 64 ? CLHEP::pi/64. : 0.;
}

void DetectorConstruction::SetLayoutPreset(G4String preset) {
  if(preset=="mixed_assumed") {
    fPlaneLayout={{
        "single1x32","2x16_staggered",
        "2x16_aligned","single1x32"}};
    fPlaneChannelOrder={{
        "ascending","ascending","ascending","ascending"}};
    fPlaneRowZSpacing={{1.*mm,.8660254*mm,1.*mm,1.*mm}};
    return;
  }
  if(preset=="v18d_compat") {
    fPlaneLayout={{
        "2x16_staggered","2x16_staggered",
        "single1x32","single1x32"}};
    fPlaneChannelOrder={{
        "ascending","ascending","ascending","ascending"}};
    fPlaneRowZSpacing={{.8660254*mm,.8660254*mm,1.*mm,1.*mm}};
    return;
  }
  G4ExceptionDescription message;
  message << "layoutPreset must be mixed_assumed or v18d_compat, got: "
          << preset;
  G4Exception("DetectorConstruction::SetLayoutPreset",
              "HODO_V21_LAYOUT_PRESET",FatalException,message);
}

void DetectorConstruction::RunNavigatorSolidProbe(G4String path) {
  const auto fail=[&](const char* code,const G4String& text) {
    G4ExceptionDescription message;
    message << text;
    G4Exception("DetectorConstruction::RunNavigatorSolidProbe",
                code,FatalException,message);
  };
  if(fNavigatorSolidProbeExpectedCondition.empty() ||
     fNavigatorSolidProbeExpectedRows<=0) {
    fail("HODO_V21_SOLID_PROBE_CONTRACT",
         "The probe requires a non-empty expected condition and a positive "
         "expected row count.");
    return;
  }
  const G4bool frozenMode=
      fNavigatorSolidProbeMode=="frozen_baseline";
  const G4bool candidateMode=
      fNavigatorSolidProbeMode=="candidate_gate";
  if(!frozenMode&&!candidateMode) {
    fail("HODO_V21_SOLID_PROBE_MODE",
         "navigatorSolidProbeMode must be frozen_baseline or candidate_gate.");
    return;
  }
  std::ifstream input(path);
  if(!input) {
    fail("HODO_V21_SOLID_PROBE_OPEN",
         "Cannot open navigator solid-probe TSV: "+path);
    return;
  }

  struct ProbeRow {
    G4String condition;
    G4int ordinal=0;
    G4String solidName;
    G4double expectedDistance=0.;
    G4String expectedState;
    G4ThreeVector original;
    G4ThreeVector direction;
    G4ThreeVector capturedIntersection;
  };
  const G4String expectedHeader =
      "condition\tordinal\tsolid\tdistance_to_in_mm\treturned_state\t"
      "original_x_mm\toriginal_y_mm\toriginal_z_mm\tdirection_x\t"
      "direction_y\tdirection_z\tintersection_x_mm\tintersection_y_mm\t"
      "intersection_z_mm";
  G4String header;
  if(!std::getline(input,header)) {
    fail("HODO_V21_SOLID_PROBE_SCHEMA",
         "The probe TSV is empty: "+path);
    return;
  }
  if(!header.empty()&&header.back()=='\r')header.pop_back();
  if(header!=expectedHeader) {
    fail("HODO_V21_SOLID_PROBE_SCHEMA",
         "Unexpected TSV header in "+path+": "+header);
    return;
  }
  std::vector<ProbeRow> rows;
  G4String line;
  G4int lineNumber=1;
  while(std::getline(input,line)) {
    ++lineNumber;
    if(!line.empty()&&line.back()=='\r')line.pop_back();
    if(line.empty()) {
      fail("HODO_V21_SOLID_PROBE_SCHEMA",
           "Blank data row at line "+std::to_string(lineNumber)+
           " in "+path);
      return;
    }
    ProbeRow row;
    G4double ox=0.,oy=0.,oz=0.;
    G4double dx=0.,dy=0.,dz=0.;
    G4double ix=0.,iy=0.,iz=0.;
    std::istringstream fields(line);
    if(!(fields>>row.condition>>row.ordinal>>row.solidName>>
         row.expectedDistance>>row.expectedState>>
         ox>>oy>>oz>>dx>>dy>>dz>>ix>>iy>>iz)) {
      fail("HODO_V21_SOLID_PROBE_SCHEMA",
           "Malformed TSV row at line "+std::to_string(lineNumber)+
           " in "+path);
      return;
    }
    fields>>std::ws;
    if(!fields.eof()) {
      fail("HODO_V21_SOLID_PROBE_SCHEMA",
           "Extra TSV field at line "+std::to_string(lineNumber)+
           " in "+path);
      return;
    }
    row.original=G4ThreeVector(ox*mm,oy*mm,oz*mm);
    row.direction=G4ThreeVector(dx,dy,dz);
    row.capturedIntersection=G4ThreeVector(ix*mm,iy*mm,iz*mm);
    const auto finiteVector=[](const G4ThreeVector& value) {
      return std::isfinite(value.x())&&std::isfinite(value.y())&&
             std::isfinite(value.z());
    };
    if(row.condition!=fNavigatorSolidProbeExpectedCondition) {
      fail("HODO_V21_SOLID_PROBE_CONDITION",
           "Condition mismatch at line "+std::to_string(lineNumber)+
           ": expected "+fNavigatorSolidProbeExpectedCondition+
           ", got "+row.condition);
      return;
    }
    if(row.ordinal!=static_cast<G4int>(rows.size())+1) {
      fail("HODO_V21_SOLID_PROBE_ORDINAL",
           "Ordinals must be contiguous from one; got "+
           std::to_string(row.ordinal)+" at line "+
           std::to_string(lineNumber));
      return;
    }
    if(!std::isfinite(row.expectedDistance)||row.expectedDistance<0.||
       !finiteVector(row.original)||!finiteVector(row.direction)||
       !finiteVector(row.capturedIntersection)||
       (row.expectedState!="kOutside"&&
        row.expectedState!="kSurface"&&row.expectedState!="kInside")) {
      fail("HODO_V21_SOLID_PROBE_VALUES",
           "Non-finite or invalid value at line "+
           std::to_string(lineNumber));
      return;
    }
    const G4double directionMagnitude2=row.direction.mag2();
    if(!std::isfinite(directionMagnitude2)||
       std::abs(directionMagnitude2-1.)>2.e-12) {
      fail("HODO_V21_SOLID_PROBE_DIRECTION",
           "The frozen direction is not unit length at line "+
           std::to_string(lineNumber)+": mag2="+
           std::to_string(directionMagnitude2));
      return;
    }
    rows.push_back(row);
  }
  if(static_cast<G4int>(rows.size())!=
     fNavigatorSolidProbeExpectedRows) {
    fail("HODO_V21_SOLID_PROBE_ROW_COUNT",
         "Expected "+std::to_string(fNavigatorSolidProbeExpectedRows)+
         " rows for "+fNavigatorSolidProbeExpectedCondition+", got "+
         std::to_string(rows.size()));
    return;
  }

  auto* store=G4SolidStore::GetInstance();
  store->UpdateMap();
  const auto stateName=[](EInside state) {
    return state==kOutside ? "kOutside"
         : state==kSurface ? "kSurface" : "kInside";
  };
  const std::array<G4double,3> perturbations{{1.*nm,10.*nm,100.*nm}};
  G4int total=static_cast<G4int>(rows.size());
  G4int missing=0;
  G4int nonUnique=0;
  G4int finiteDistance=0;
  G4int inconsistent=0;
  G4int expectedMatch=0;
  G4int capturedStateMatch=0;
  G4int capturedNonSurface=0;
  G4int deterministicFailure=0;
  G4int directQueries=0;
  G4int adaptiveQueries=0;
  G4int adaptiveMinusQueries=0;
  G4int adaptivePlusQueries=0;
  G4int unmappable=0;
  G4int surfaceFailure=0;
  G4int bracket10Failure=0;
  G4int bracket100Failure=0;
  G4int localityFailure=0;
  G4int dtoEvaluated=0;
  G4int dtoDiagnosticNonlocal=0;
  G4int dtoRequired=0;
  G4int dtoPass=0;
  G4int dtoFiniteFailure=0;
  G4int dtoDeterministicFailure=0;
  G4int dtoSurfaceFailure=0;
  G4int dtoBracket10Failure=0;
  G4int dtoBracket100Failure=0;
  G4int dtoLocalityFailure=0;
  G4int candidatePass=0;
  // The largest intentional radial/envelope change is about 1.2 um.  The
  // frozen minimum |ray dot facet-normal| is 0.034454, so the corresponding
  // along-ray displacement is about 35 um.  A 100-um bound is conservative
  // but still three orders of magnitude shorter than the 102-mm fibre and
  // therefore cannot silently remap a query to a remote segment.
  const G4double candidateLocalityLimit=100.*um;

  for(const auto& row:rows) {
    auto bucket=store->GetMap().find(row.solidName);
    if(bucket==store->GetMap().end()||bucket->second.empty()) {
      ++missing;
      G4cout << "V21_SOLID_PROBE id=" << row.condition << ":"
             << row.ordinal << " solid=" << row.solidName
             << " status=MISSING" << G4endl;
      continue;
    }
    if(bucket->second.size()!=1||bucket->second.front()==nullptr) {
      ++nonUnique;
      G4cout << "V21_SOLID_PROBE id=" << row.condition << ":"
             << row.ordinal << " solid=" << row.solidName
             << " status=NON_UNIQUE count=" << bucket->second.size()
             << G4endl;
      continue;
    }
    auto* solid=bucket->second.front();
    if(!solid) {
      ++missing;
      G4cout << "V21_SOLID_PROBE id=" << row.condition << ":"
             << row.ordinal << " solid=" << row.solidName
             << " status=NULL" << G4endl;
      continue;
    }

    const auto transformEntry=
        fNavigatorSolidLocalToWorld.find(row.solidName);
    const G4bool transformed=
        transformEntry!=fNavigatorSolidLocalToWorld.end();
    const G4Transform3D identity;
    const auto& localToWorld=
        transformed ? transformEntry->second : identity;
    const auto toLocalPoint=[&](const G4ThreeVector& worldPoint) {
      if(!transformed)return worldPoint;
      return localToWorld.getRotation().inverse()*
          (worldPoint-localToWorld.getTranslation());
    };
    const auto toLocalDirection=[&](const G4ThreeVector& worldDirection) {
      if(!transformed)return worldDirection;
      return localToWorld.getRotation().inverse()*worldDirection;
    };
    const auto toWorldPoint=[&](const G4ThreeVector& localPoint) {
      if(!transformed)return localPoint;
      return localToWorld.getRotation()*localPoint+
          localToWorld.getTranslation();
    };
    const auto insideWorld=[&](const G4ThreeVector& worldPoint) {
      return solid->Inside(toLocalPoint(worldPoint));
    };
    const auto localDirection=toLocalDirection(row.direction);
    const auto originalState=insideWorld(row.original);

    // A candidate may classify a tolerance-singular frozen Outside point as
    // Surface or Inside.  Diagnose the nearest forward exit in either case
    // and require it for every adaptive-plus route.  This DTO query does not
    // replace the DTI contract below: the latter is still run from a
    // certified Outside origin with the original, unnormalised frozen
    // direction.
    const G4bool dtoEvaluatedThis=candidateMode&&
        originalState!=kOutside;
    G4bool dtoRequiredThis=false;
    G4bool dtoFinite=false;
    G4bool dtoDeterministic=false;
    G4bool dtoSurface=false;
    G4bool dtoBracket10=false;
    G4bool dtoBracket100=false;
    G4bool dtoLocal=false;
    G4bool dtoPassedThis=false;
    G4double dtoDistance1=kInfinity;
    G4ThreeVector dtoIntersection1=row.original;
    if(dtoEvaluatedThis) {
      ++dtoEvaluated;
      const auto localOriginal=toLocalPoint(row.original);
      dtoDistance1=solid->DistanceToOut(localOriginal,localDirection);
      const G4double dtoDistance2=
          solid->DistanceToOut(localOriginal,localDirection);
      const G4bool dtoFinite1=std::isfinite(dtoDistance1)&&
                              dtoDistance1>=0.&&
                              dtoDistance1<.5*kInfinity;
      const G4bool dtoFinite2=std::isfinite(dtoDistance2)&&
                              dtoDistance2>=0.&&
                              dtoDistance2<.5*kInfinity;
      dtoFinite=dtoFinite1&&dtoFinite2;
      const auto localDtoIntersection1=
          dtoFinite1 ? localOriginal+dtoDistance1*localDirection
                     : localOriginal;
      const auto localDtoIntersection2=
          dtoFinite2 ? localOriginal+dtoDistance2*localDirection
                     : localOriginal;
      dtoIntersection1=toWorldPoint(localDtoIntersection1);
      const auto dtoIntersection2=toWorldPoint(localDtoIntersection2);
      const auto dtoState1=solid->Inside(localDtoIntersection1);
      const auto dtoState2=solid->Inside(localDtoIntersection2);
      dtoDeterministic=dtoFinite&&(dtoDistance1==dtoDistance2)&&
          (dtoIntersection1==dtoIntersection2)&&dtoState1==dtoState2;
      dtoSurface=dtoFinite&&dtoState1==kSurface;
      if(dtoFinite) {
        dtoBracket10=
            insideWorld(dtoIntersection1-perturbations[1]*row.direction)
                ==kInside&&
            insideWorld(dtoIntersection1+perturbations[1]*row.direction)
                ==kOutside;
        dtoBracket100=
            insideWorld(dtoIntersection1-perturbations[2]*row.direction)
                ==kInside&&
            insideWorld(dtoIntersection1+perturbations[2]*row.direction)
                ==kOutside;
        dtoLocal=(dtoIntersection1-row.capturedIntersection).mag()<=
            candidateLocalityLimit;
      }
      dtoPassedThis=dtoFinite&&dtoDeterministic&&dtoSurface&&
          dtoBracket10&&dtoBracket100&&dtoLocal;
      if(!dtoLocal)++dtoDiagnosticNonlocal;
    }

    G4ThreeVector queryOrigin=row.original;
    G4double adaptiveOffset=0.;
    G4int adaptiveSign=0;
    G4bool mapped=true;
    if(candidateMode&&originalState!=kOutside) {
      mapped=false;
      G4double trialOffset=1.*nm;
      while(trialOffset<=candidateLocalityLimit) {
        const auto minusOrigin=
            row.capturedIntersection-trialOffset*row.direction;
        const auto plusOrigin=
            row.capturedIntersection+trialOffset*row.direction;
        const G4bool minusOutside=insideWorld(minusOrigin)==kOutside;
        const G4bool plusOutside=insideWorld(plusOrigin)==kOutside;
        if(minusOutside||plusOutside) {
          // Both signs are tested at every distance.  The first successful
          // distance is therefore the nearest sampled Outside origin; a
          // same-distance tie is resolved toward minus deterministically.
          adaptiveOffset=trialOffset;
          adaptiveSign=minusOutside ? -1 : 1;
          queryOrigin=adaptiveSign<0 ? minusOrigin : plusOrigin;
          mapped=true;
          break;
        }
        if(trialOffset==candidateLocalityLimit)break;
        trialOffset=std::min(2.*trialOffset,candidateLocalityLimit);
      }
      if(mapped) {
        ++adaptiveQueries;
        if(adaptiveSign<0)++adaptiveMinusQueries;
        else ++adaptivePlusQueries;
      }
    } else {
      ++directQueries;
    }
    // The local DTO is a required supplemental gate only when the remapper
    // had to cross the nearby forward exit and select an Outside point on the
    // plus side.  For adaptive-minus, the original +direction DTI below
    // already certifies the local Outside -> Surface -> Inside transition;
    // its +direction DTO can legitimately be the remote far-end exit.
    dtoRequiredThis=mapped&&adaptiveSign>0;
    if(dtoRequiredThis) {
      ++dtoRequired;
      if(dtoPassedThis)++dtoPass;
      if(!dtoFinite)++dtoFiniteFailure;
      if(!dtoDeterministic)++dtoDeterministicFailure;
      if(!dtoSurface)++dtoSurfaceFailure;
      if(!dtoBracket10)++dtoBracket10Failure;
      if(!dtoBracket100)++dtoBracket100Failure;
      if(!dtoLocal)++dtoLocalityFailure;
    }
    if(!mapped) {
      ++unmappable;
      G4cout << "V21_SOLID_PROBE id=" << row.condition << ":"
             << row.ordinal << " solid=" << row.solidName
             << " mode=" << fNavigatorSolidProbeMode
             << " route=unmappable"
             << " start_kind=captured_intersection_bilateral_search"
             << " status=UNMAPPABLE original_state="
             << stateName(originalState)
             << " dto_evaluated=" << dtoEvaluatedThis
             << " dto_required=" << dtoRequiredThis
             << " dto_diagnostic_pass=" << dtoPassedThis << G4endl;
      continue;
    }

    // Preserve the frozen raw direction exactly.  Every input direction is
    // already unit length to <2e-12 in mag2; calling unit() here changes the
    // original floating-point query and can hide or create a singular result.
    const auto localQueryOrigin=toLocalPoint(queryOrigin);
    const G4double distance1=
        solid->DistanceToIn(localQueryOrigin,localDirection);
    const G4double distance2=
        solid->DistanceToIn(localQueryOrigin,localDirection);
    const G4bool finite1=std::isfinite(distance1)&&distance1>=0.&&
                         distance1<.5*kInfinity;
    const G4bool finite2=std::isfinite(distance2)&&distance2>=0.&&
                         distance2<.5*kInfinity;
    const G4bool finite=finite1&&finite2;
    if(finite)++finiteDistance;
    const auto localIntersection1=
        finite1 ? localQueryOrigin+distance1*localDirection
                : localQueryOrigin;
    const auto localIntersection2=
        finite2 ? localQueryOrigin+distance2*localDirection
                : localQueryOrigin;
    const auto intersection1=toWorldPoint(localIntersection1);
    const auto intersection2=toWorldPoint(localIntersection2);
    const auto intersectionState1=
        finite1 ? solid->Inside(localIntersection1)
                : solid->Inside(localQueryOrigin);
    const auto intersectionState2=
        finite2 ? solid->Inside(localIntersection2)
                : solid->Inside(localQueryOrigin);
    const G4bool deterministic=
        finite&&(distance1==distance2)&&
        (intersection1==intersection2)&&
        intersectionState1==intersectionState2;
    if(!deterministic)++deterministicFailure;
    if(finite&&intersectionState1!=kSurface)++inconsistent;

    const auto capturedState=insideWorld(row.capturedIntersection);
    if(capturedState!=kSurface)++capturedNonSurface;
    if(row.expectedState==stateName(capturedState))++capturedStateMatch;
    const G4double distanceTolerance=std::max(
        1.e-10*mm,1.e-10*std::abs(row.expectedDistance*mm));
    const G4bool matches=frozenMode&&finite &&
        std::abs(distance1-row.expectedDistance*mm)<=distanceTolerance &&
        (intersection1-row.capturedIntersection).mag()<=1.e-8*mm &&
        row.expectedState==stateName(intersectionState1);
    if(matches)++expectedMatch;

    std::array<EInside,3> minusState{{kOutside,kOutside,kOutside}};
    std::array<EInside,3> plusState{{kOutside,kOutside,kOutside}};
    if(finite) {
      for(std::size_t i=0;i<perturbations.size();++i) {
        minusState[i]=insideWorld(
            intersection1-perturbations[i]*row.direction);
        plusState[i]=insideWorld(
            intersection1+perturbations[i]*row.direction);
      }
    }
    const G4bool bracket10=finite&&
        minusState[1]==kOutside&&plusState[1]==kInside;
    const G4bool bracket100=finite&&
        minusState[2]==kOutside&&plusState[2]==kInside;
    const G4bool localIntersection=finite&&
        (intersection1-row.capturedIntersection).mag()<=
            candidateLocalityLimit;
    if(candidateMode) {
      if(finite&&intersectionState1!=kSurface)++surfaceFailure;
      if(!bracket10)++bracket10Failure;
      if(!bracket100)++bracket100Failure;
      if(!localIntersection)++localityFailure;
      const G4bool passed=finite&&deterministic&&
          intersectionState1==kSurface&&bracket10&&bracket100&&
          localIntersection&&(!dtoRequiredThis||dtoPassedThis);
      if(passed)++candidatePass;
    }

    std::ostringstream perturbationReport;
    perturbationReport << std::setprecision(17);
    for(std::size_t i=0;i<perturbations.size();++i) {
      perturbationReport
          << " d" << perturbations[i]/nm << "nm=("
          << stateName(minusState[i]) << ","
          << stateName(plusState[i]) << ")";
    }
    G4cout << std::setprecision(17)
           << "V21_SOLID_PROBE id=" << row.condition << ":"
           << row.ordinal << " solid=" << row.solidName
           << " mode=" << fNavigatorSolidProbeMode
           << " solid_frame="
           << (transformed ? "placed_root_local" : "common_identity")
           << " route="
           << (candidateMode&&originalState!=kOutside
               ? (adaptiveSign<0 ? "adaptive_minus" : "adaptive_plus")
               : "direct")
           << " start_kind="
           << (candidateMode&&originalState!=kOutside
               ? "captured_intersection_bilateral_search" : "old_original")
           << " original_state=" << stateName(originalState)
           << " adaptive_offset_nm=" << adaptiveOffset/nm
           << " distance_mm=" << (finite ? distance1/mm : -1.)
           << " deterministic=" << deterministic
           << " intersection_state=" << stateName(intersectionState1)
           << " expected_distance_mm=" << row.expectedDistance
           << " captured_intersection_state=" << stateName(capturedState)
           << " expected_state=" << row.expectedState
           << " captured_delta_nm="
           << (finite
               ? (intersection1-row.capturedIntersection).mag()/nm : -1.)
           << " frozen_expected_match=" << matches
           << " locality_limit_um=" << candidateLocalityLimit/um
           << " local_intersection=" << localIntersection
           << " bracket10=" << bracket10
           << " bracket100=" << bracket100
           << " dto_required=" << dtoRequiredThis
           << " dto_evaluated=" << dtoEvaluatedThis
           << " dto_diagnostic_pass=" << dtoPassedThis
           << " dto_distance_nm="
           << (dtoFinite ? dtoDistance1/nm : -1.)
           << " dto_deterministic=" << dtoDeterministic
           << " dto_surface=" << dtoSurface
           << " dto_bracket10=" << dtoBracket10
           << " dto_bracket100=" << dtoBracket100
           << " dto_local=" << dtoLocal
           << " dto_captured_delta_nm="
           << (dtoFinite
               ? (dtoIntersection1-row.capturedIntersection).mag()/nm
               : -1.)
           << perturbationReport.str() << G4endl;
  }

  G4cout << "V21_SOLID_PROBE_COMPLETE path=" << path
         << " mode=" << fNavigatorSolidProbeMode
         << " condition=" << fNavigatorSolidProbeExpectedCondition
         << " total=" << total
         << " missing=" << missing
         << " non_unique=" << nonUnique
         << " finite_distance=" << finiteDistance
         << " deterministic_failure=" << deterministicFailure
         << " inconsistent=" << inconsistent
         << " expected_match=" << expectedMatch
         << " captured_non_surface=" << capturedNonSurface
         << " captured_state_match=" << capturedStateMatch
         << " direct_queries=" << directQueries
         << " adaptive_queries=" << adaptiveQueries
         << " adaptive_minus_queries=" << adaptiveMinusQueries
         << " adaptive_plus_queries=" << adaptivePlusQueries
         << " unmappable=" << unmappable
         << " surface_failure=" << surfaceFailure
         << " bracket10_failure=" << bracket10Failure
         << " bracket100_failure=" << bracket100Failure
         << " locality_failure=" << localityFailure
         << " dto_evaluated=" << dtoEvaluated
         << " dto_diagnostic_nonlocal=" << dtoDiagnosticNonlocal
         << " dto_required=" << dtoRequired
         << " dto_pass=" << dtoPass
         << " dto_finite_failure=" << dtoFiniteFailure
         << " dto_deterministic_failure=" << dtoDeterministicFailure
         << " dto_surface_failure=" << dtoSurfaceFailure
         << " dto_bracket10_failure=" << dtoBracket10Failure
         << " dto_bracket100_failure=" << dtoBracket100Failure
         << " dto_locality_failure=" << dtoLocalityFailure
         << " dto_route_valid="
         << (dtoRequired==adaptivePlusQueries)
         << " dto_gate_valid="
         << (dtoRequired==adaptivePlusQueries&&
             dtoPass==dtoRequired&&dtoFiniteFailure==0&&
             dtoDeterministicFailure==0&&dtoSurfaceFailure==0&&
             dtoBracket10Failure==0&&dtoBracket100Failure==0&&
             dtoLocalityFailure==0)
         << " candidate_pass=" << candidatePass
         << " route_count_valid="
         << (directQueries+adaptiveMinusQueries+adaptivePlusQueries+
                 unmappable==total&&
             adaptiveQueries==adaptiveMinusQueries+adaptivePlusQueries)
         << G4endl;

  const G4bool infrastructureFailure=
      missing!=0||nonUnique!=0||
      directQueries+adaptiveMinusQueries+adaptivePlusQueries+
              unmappable!=total||
      adaptiveQueries!=adaptiveMinusQueries+adaptivePlusQueries||
      dtoRequired!=adaptivePlusQueries||
      (candidateMode&&
       (finiteDistance!=total||deterministicFailure!=0));
  const G4bool transitionFailure=candidateMode&&
      (unmappable!=0||surfaceFailure!=0||bracket10Failure!=0||
       bracket100Failure!=0||localityFailure!=0||
       dtoPass!=dtoRequired||dtoFiniteFailure!=0||
       dtoDeterministicFailure!=0||dtoSurfaceFailure!=0||
       dtoBracket10Failure!=0||dtoBracket100Failure!=0||
       dtoLocalityFailure!=0||candidatePass!=total);
  if(infrastructureFailure||transitionFailure) {
    fail("HODO_V21_SOLID_PROBE_GATE",
         "Fail-closed solid-probe contract failed for "+
         fNavigatorSolidProbeExpectedCondition+" in "+
         fNavigatorSolidProbeMode+" mode.");
    return;
  }
}

G4double DetectorConstruction::GetPlaneSourceTransverse(
    G4int plane, G4int channel) const {
  if(GetPlaneFibreCount(plane)!=32)
    return (channel-.5*(GetPlaneFibreCount(plane)-1))*mm;
  const auto slot = ResolveLayoutSlot(
      GetPlaneLayout(plane), GetPlaneChannelOrder(plane), channel);
  const G4double centre = .5*(slot.columns-1);
  const G4double stagger =
      GetPlaneLayout(plane) == "2x16_staggered" ? .5*slot.row*mm : 0.;
  return (slot.column-centre)*mm+stagger;
}

std::array<G4double,2> DetectorConstruction::GetActiveLongitudinalBounds(G4int plane) const {
  if(!fActivePlaneOverlap||plane<0||plane>=4)return {{-16.*mm,16.*mm}};
  const G4int other=plane^1;
  G4double low=DBL_MAX,high=-DBL_MAX;
  const G4double radius=GetFanoutOuterSize(other);
  // Project the OTHER plane's packed transverse envelope into this plane's
  // longitudinal frame. This includes the existing module offsets/180-deg
  // rotation and plane pivots; it does not move either plane.
  for(G4int ch=0;ch<GetPlaneFibreCount(other);++ch) for(const auto sign:{-1.,1.}) {
    const auto v=GetPlaneSourceTransverse(other,ch)+sign*radius;
    const auto p=other%2==0?G4ThreeVector(0.,v,0.):G4ThreeVector(v,0.,0.);
    const auto local=PlaneWorldToNominalPosition(plane,
        PlaneNominalToWorldPosition(other,p));
    const auto coordinate=plane%2==0?local.x():local.y();
    low=std::min(low,coordinate);high=std::max(high,coordinate);
  }
  return {{low,high}};
}

G4double DetectorConstruction::GetPlaneSourceZ(
    G4int plane, G4int channel) const {
  // The photo-derived fan-out preview and the compact straight detector use
  // different plane centres.  Fixed-photon validation must start inside the
  // geometry that is actually constructed, rather than reusing the fan-out
  // source coordinates for both representations.
  if(!fFanoutPreview) {
    if(plane<0||plane>=4)return 0.;
    if(plane<2) {
      const auto slot=ResolveLayoutSlot(
          GetPlaneLayout(plane),GetPlaneChannelOrder(plane),channel);
      return 2.5*mm+slot.row*GetPlaneRowZSpacing(plane)+
             plane*2.2*mm;
    }
    return -2.5*mm-(plane-2)*1.*mm;
  }
  const std::array<G4double,4> baseZ{{
      fHodoscopeFanoutZ[0],
      fHodoscopeFanoutZ[0]+2.2*mm,
      fHodoscopeFanoutZ[1],
      fHodoscopeFanoutZ[1]-1.*mm}};
  if (plane < 0 || plane >= 4) return 0.;
  const auto slot = ResolveLayoutSlot(
      GetPlaneLayout(plane), GetPlaneChannelOrder(plane), channel);
  return baseZ[plane]+slot.row*GetPlaneRowZSpacing(plane);
}

G4double DetectorConstruction::GetValidationSurfaceZ(G4int surface) const {
  if (surface == 0) return fTriggerZ[0];  // Trigger 1 centre
  if (surface == 1) return fTriggerZ[1];  // Trigger 2 centre
  const G4int plane = surface-2;
  if (plane < 0 || plane >= 4) return 0.;
  G4double z = 0.;
  for (G4int channel=0; channel<GetPlaneFibreCount(plane); ++channel)
    z += GetPlaneSourceZ(plane,channel);
  return z/GetPlaneFibreCount(plane);
}

G4ThreeVector DetectorConstruction::GetTriggerCentre(G4int trigger) const {
  if(trigger<0||trigger>=4)return {};
  // Preserve the predecessor arithmetic path exactly.  This is deliberately
  // checked before constructing or adding an offset vector.
  if(trigger<2&&fAlignmentReference=="legacy"&&
     fTriggerOffsetX[trigger]==0.&&fTriggerOffsetY[trigger]==0.) {
    if(trigger==0)return {-7.*mm,0.,fTriggerZ[0]};
    return {0.,-7.*mm,fTriggerZ[1]};
  }
  G4ThreeVector centre;
  if(trigger<2&&fAlignmentReference=="legacy")
    centre=trigger==0 ? G4ThreeVector(-7.*mm,0.,fTriggerZ[0])
                      : G4ThreeVector(0.,-7.*mm,fTriggerZ[1]);
  else
    centre=G4ThreeVector(0.,0.,fTriggerZ[trigger]);
  centre.setX(centre.x()+fTriggerOffsetX[trigger]);
  centre.setY(centre.y()+fTriggerOffsetY[trigger]);
  return centre;
}

G4double DetectorConstruction::GetTriggerRotationAngle(G4int trigger) const {
  if(trigger==0)return -fTriggerYaw[0];
  if(trigger==1)return 90.*deg-fTriggerYaw[1];
  if(trigger==2||trigger==3)return -fTriggerYaw[trigger];
  return 0.;
}

G4ThreeVector DetectorConstruction::TriggerWorldToLocalPosition(
    G4int trigger,const G4ThreeVector& position) const {
  const auto relative=position-GetTriggerCentre(trigger);
  const auto angle=GetTriggerRotationAngle(trigger);
  const auto cosine=std::cos(angle);
  const auto sine=std::sin(angle);
  return {cosine*relative.x()+sine*relative.y(),
         -sine*relative.x()+cosine*relative.y(),relative.z()};
}

G4ThreeVector DetectorConstruction::TriggerWorldToLocalDirection(
    G4int trigger,const G4ThreeVector& direction) const {
  const auto angle=GetTriggerRotationAngle(trigger);
  const auto cosine=std::cos(angle);
  const auto sine=std::sin(angle);
  return {cosine*direction.x()+sine*direction.y(),
         -sine*direction.x()+cosine*direction.y(),direction.z()};
}

G4int DetectorConstruction::GetPlaneFanoutRow(
    G4int plane, G4int channel) const {
  const auto slot = ResolveLayoutSlot(
      GetPlaneLayout(plane), GetPlaneChannelOrder(plane), channel);
  return slot.columns == 32
      ? slot.column%4 : 2*slot.row+slot.column%2;
}

G4int DetectorConstruction::GetPlaneFanoutColumn(
    G4int plane, G4int channel) const {
  const auto slot = ResolveLayoutSlot(
      GetPlaneLayout(plane), GetPlaneChannelOrder(plane), channel);
  return slot.columns == 32 ? slot.column/4 : slot.column/2;
}

G4bool DetectorConstruction::GetMixedLayoutPattern() const {
  return fPlaneLayout[0] == "single1x32" &&
         fPlaneLayout[1] == "2x16_staggered" &&
         fPlaneLayout[2] == "2x16_aligned" &&
         fPlaneLayout[3] == "single1x32";
}

G4ThreeVector DetectorConstruction::PlaneNominalPivot(G4int plane) const {
  if (plane < 0 || plane >= 4) return {};
  G4double transverse = 0.;
  for (G4int channel=0; channel<GetPlaneFibreCount(plane); ++channel)
    transverse += GetPlaneSourceTransverse(plane,channel);
  transverse /= GetPlaneFibreCount(plane);
  return plane==0||plane==2
      ? G4ThreeVector(0.,transverse,0.)
      : G4ThreeVector(transverse,0.,0.);
}

G4ThreeVector DetectorConstruction::PlaneNominalToWorldPosition(
    G4int plane, const G4ThreeVector& position) const {
  if (plane < 0 || plane >= 4) return position;
  const auto module=plane/2;
  const auto offsetX=fHodoscopeOffsetX[module];
  const auto offsetY=fHodoscopeOffsetY[module];
  const G4double moduleYaw =
      fH2ModuleRotate180 && module == 1 ? CLHEP::pi : 0.;
  const G4double totalYaw = moduleYaw + fPlaneYaw[plane];
  if (totalYaw == 0. && offsetX == 0. && offsetY == 0.)
    return position;
  if (totalYaw == 0.)
    return position+G4ThreeVector(offsetX,offsetY,0.);
  const auto pivot = PlaneNominalPivot(plane);
  const auto relative = position - pivot;
  const auto cosine = std::cos(totalYaw);
  const auto sine = std::sin(totalYaw);
  return pivot+G4ThreeVector(offsetX,offsetY,0.)+
         G4ThreeVector(cosine*relative.x()-sine*relative.y(),
                       sine*relative.x()+cosine*relative.y(),
                       relative.z());
}

G4ThreeVector DetectorConstruction::PlaneWorldToNominalPosition(
    G4int plane, const G4ThreeVector& position) const {
  if (plane < 0 || plane >= 4) return position;
  const auto module=plane/2;
  const auto offsetX=fHodoscopeOffsetX[module];
  const auto offsetY=fHodoscopeOffsetY[module];
  const G4double moduleYaw =
      fH2ModuleRotate180 && module == 1 ? CLHEP::pi : 0.;
  const G4double totalYaw = moduleYaw + fPlaneYaw[plane];
  if (totalYaw == 0. && offsetX == 0. && offsetY == 0.)
    return position;
  if (totalYaw == 0.)
    return position-G4ThreeVector(offsetX,offsetY,0.);
  const auto pivot = PlaneNominalPivot(plane);
  const auto relative =
      position-pivot-G4ThreeVector(offsetX,offsetY,0.);
  const auto cosine = std::cos(totalYaw);
  const auto sine = std::sin(totalYaw);
  return pivot + G4ThreeVector(cosine*relative.x() + sine*relative.y(),
                             -sine*relative.x() + cosine*relative.y(),
                              relative.z());
}

G4ThreeVector DetectorConstruction::PlaneNominalToWorldDirection(
    G4int plane, const G4ThreeVector& direction) const {
  if (plane < 0 || plane >= 4) return direction;
  const auto module=plane/2;
  const G4double moduleYaw =
      fH2ModuleRotate180 && module == 1 ? CLHEP::pi : 0.;
  const G4double totalYaw = moduleYaw + fPlaneYaw[plane];
  if (totalYaw == 0.) return direction;
  const auto cosine = std::cos(totalYaw);
  const auto sine = std::sin(totalYaw);
  return {cosine*direction.x() - sine*direction.y(),
          sine*direction.x() + cosine*direction.y(), direction.z()};
}

G4ThreeVector DetectorConstruction::PlaneWorldToNominalDirection(
    G4int plane, const G4ThreeVector& direction) const {
  if (plane < 0 || plane >= 4) return direction;
  const auto module=plane/2;
  const G4double moduleYaw =
      fH2ModuleRotate180 && module == 1 ? CLHEP::pi : 0.;
  const G4double totalYaw = moduleYaw + fPlaneYaw[plane];
  if (totalYaw == 0.) return direction;
  const auto cosine = std::cos(totalYaw);
  const auto sine = std::sin(totalYaw);
  return {cosine*direction.x() + sine*direction.y(),
         -sine*direction.x() + cosine*direction.y(), direction.z()};
}

void DetectorConstruction::DefineMaterials() {
  auto* nist = G4NistManager::Instance();
  auto* polystyrene = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
  fSCSF78 = new G4Material("Kuraray_SCSF78_core", polystyrene->GetDensity(), 1);
  fSCSF78->AddMaterial(polystyrene, 1.0);
  if (GetH1Uses22PointEmission()) {
    AddSCSF78Optics(fSCSF78, fH1AttenuationLength,
                    fH1ScintillationYield*fScintillationYieldScale,
                    fH1BirksConstant*mm/MeV);
  } else {
    AddScintillatorOptics(fSCSF78, 1.59, fH1AttenuationLength,
                          fH1ScintillationYield*fScintillationYieldScale, 2.8*ns,
                          {0.01, 0.16, 0.88, 1.00, 0.18},
                          fH1BirksConstant*mm/MeV);
  }

  fBCF20 = new G4Material("BCF20_core", polystyrene->GetDensity(), 1);
  fBCF20->AddMaterial(polystyrene, 1.0);
  // Approximate BCF-20 emission: green, near 492 nm.
  AddScintillatorOptics(fBCF20, 1.60, fH2AttenuationLength,
                        fH2ScintillationYield*fScintillationYieldScale, 2.7*ns,
                        {0.30, 1.00, 0.42, 0.07, 0.01},
                        fH2BirksConstant*mm/MeV);

  auto* pmma = nist->FindOrBuildMaterial("G4_PLEXIGLASS");
  fInnerCladding = new G4Material("PMMA_inner_cladding", pmma->GetDensity(), 1);
  fInnerCladding->AddMaterial(pmma, 1.0);
  AddConstantOptics(fInnerCladding, 1.49, 10.*m);

  auto* polyethylene = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  fOuterCladding = new G4Material("Fluorinated_outer_cladding",
                                  polyethylene->GetDensity(), 1);
  fOuterCladding->AddMaterial(polyethylene, 1.0);
  AddConstantOptics(fOuterCladding, 1.42, 10.*m);

  fGrease = new G4Material("EJ550_coupling_grease", 1.0*g/cm3, 1);
  fGrease->AddMaterial(nist->FindOrBuildMaterial("G4_WATER"), 1.0);
  AddConstantOptics(fGrease, 1.46, 1.*m);

  fSiPMWindow = new G4Material("S14160_silicone_resin_window", 1.0*g/cm3, 1);
  fSiPMWindow->AddMaterial(nist->FindOrBuildMaterial("G4_WATER"), 1.0);
  AddConstantOptics(fSiPMWindow, 1.57, 10.*cm);

  fSilicon = nist->FindOrBuildMaterial("G4_Si");
  AddConstantOptics(fSilicon, 3.8, 1.*um);

  auto* vinyltoluene = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  fTriggerScint = new G4Material("Trigger_vinyltoluene_scintillator",
                                 vinyltoluene->GetDensity(), 1);
  fTriggerScint->AddMaterial(vinyltoluene, 1.0);
  // Representative plastic-paddle response; the actual trigger type is unresolved.
  AddScintillatorOptics(fTriggerScint, 1.58, 2.*m,
                        fTriggerScintillationYield, 2.1*ns,
                        {0.01, 0.10, 0.75, 1.00, 0.20},
                        fTriggerBirksConstant*mm/MeV);

  auto* air = nist->FindOrBuildMaterial("G4_AIR");
  AddConstantOptics(air, 1.0003, 100.*m);
}

void DetectorConstruction::DefineSurfaces() {
  fFiberSurface = new G4OpticalSurface("Fibre_core_surface");
  fFiberSurface->SetType(dielectric_dielectric);
  fFiberSurface->SetModel(unified);
  fFiberSurface->SetFinish(fFiberSurfaceRoughness > 0. ? ground : polished);
  fFiberSurface->SetSigmaAlpha(std::max(0., fFiberSurfaceRoughness));

  fEndFaceSurface = new G4OpticalSurface("Sanded_fibre_end_surface");
  fEndFaceSurface->SetType(dielectric_dielectric);
  fEndFaceSurface->SetModel(unified);
  fEndFaceSurface->SetFinish(fEndFaceRoughness > 0. ? ground : polished);
  fEndFaceSurface->SetSigmaAlpha(std::max(0.,fEndFaceRoughness));

  fSiPMSurface = new G4OpticalSurface("SiPM_photodetection_surface");
  fSiPMSurface->SetType(dielectric_metal);
  fSiPMSurface->SetModel(unified);
  fSiPMSurface->SetFinish(polished);
  auto* sipm = new G4MaterialPropertiesTable;
  // S14160-1310PS, 10 um-pixel family curve digitized from the supplied manual.
  // Active area is 1.3 x 1.3 mm; the optical column remains clipped to the
  // 1 mm fibre aperture until the measured fan-out is explicitly modelled.
  const std::vector<G4double> nominalPde{0.165, 0.178, 0.180, 0.170, 0.145};
  std::vector<G4double> pde;
  // V15 keeps the SiPM PDE independent of upstream optical losses. Residual
  // fan-out, end-face, and coupling efficiencies are applied at their own
  // physical boundaries by the stepping action and receive separate counters.
  for (auto value : nominalPde)
    pde.push_back(std::clamp(value*fPDEScale, 0., 1.));
  sipm->AddProperty("EFFICIENCY", kPhotonEnergy, pde);
  sipm->AddProperty("REFLECTIVITY", kPhotonEnergy,
                    std::vector<G4double>(kPhotonEnergy.size(), 0.05));
  fSiPMSurface->SetMaterialPropertiesTable(sipm);

  fMirrorSurface = new G4OpticalSurface("Far_end_reflector_surface");
  fMirrorSurface->SetType(dielectric_metal);
  fMirrorSurface->SetModel(unified);
  fMirrorSurface->SetFinish(fFiberSurfaceRoughness > 0. ? ground : polished);
  fMirrorSurface->SetSigmaAlpha(fFiberSurfaceRoughness);
  auto* mirror = new G4MaterialPropertiesTable;
  mirror->AddProperty("REFLECTIVITY", kPhotonEnergy,
                      std::vector<G4double>(kPhotonEnergy.size(),
                                            std::clamp(fFarEndReflectivity, 0., 1.)));
  mirror->AddProperty("EFFICIENCY", kPhotonEnergy,
                      std::vector<G4double>(kPhotonEnergy.size(), 0.));
  fMirrorSurface->SetMaterialPropertiesTable(mirror);

  fBlackWrapSurface = new G4OpticalSurface("Trigger_black_wrapping");
  fBlackWrapSurface->SetType(dielectric_metal);
  fBlackWrapSurface->SetModel(unified);
  fBlackWrapSurface->SetFinish(ground);
  auto* black = new G4MaterialPropertiesTable;
  black->AddProperty("REFLECTIVITY", kPhotonEnergy,
                     std::vector<G4double>(kPhotonEnergy.size(), 0.02));
  black->AddProperty("EFFICIENCY", kPhotonEnergy,
                     std::vector<G4double>(kPhotonEnergy.size(), 0.));
  fBlackWrapSurface->SetMaterialPropertiesTable(black);
}

G4VPhysicalVolume* DetectorConstruction::Construct() {
  fNavigatorSolidLocalToWorld.clear();
  fFanoutEndpoint.fill(G4ThreeVector());
  fFanoutTerminalDirection.fill(G4ThreeVector());
  fFanoutEndpointReady.fill(false);
  if(!std::isfinite(fFibreMaxStep) || fFibreMaxStep < 0.)
    G4Exception("DetectorConstruction::Construct",
                "HODO_GATE3_FIBRE_MAX_STEP",FatalException,
                "fibreMaxStep must be finite and nonnegative; 0 disables "
                "the Stage-3 user limit.");
  if(fGenerationProbeRowOffset < 0 ||
     (fGenerationProbeMode &&
      fGenerationProbeRowOffset >= kGenerationProbeExpectedRows)) {
    G4ExceptionDescription message;
    message << "generationProbeRowOffset must be nonnegative and, when "
            << "generationProbeMode is enabled, in [0,"
            << kGenerationProbeExpectedRows-1 << "].";
    G4Exception("DetectorConstruction::Construct",
                "HODO_GATE3_PROBE_OFFSET",FatalException,
                message);
  }
  if(fGeometryProbeMode && fGenerationProbeMode)
    G4Exception("DetectorConstruction::Construct",
                "HODO_GATE3_PROBE_MODE_CONFLICT",FatalException,
                "geometryProbeMode and generationProbeMode are mutually "
                "exclusive event-indexed primary sources.");
  if(fTriggerCoincidenceMode!="t1t2" &&
     fTriggerCoincidenceMode!="t1t2t3" &&
     fTriggerCoincidenceMode!="t1t2t3t4") {
    G4ExceptionDescription message;
    message << "triggerCoincidenceMode must be t1t2, t1t2t3, or "
            << "t1t2t3t4; got: " << fTriggerCoincidenceMode;
    G4Exception("DetectorConstruction::Construct",
                "HODO_TRIGGER_COINCIDENCE_MODE",FatalException,message);
  }
  if((fTriggerCoincidenceMode=="t1t2t3"&&!fTriggerEnabled[2]) ||
     (fTriggerCoincidenceMode=="t1t2t3t4"&&
      (!fTriggerEnabled[2]||!fTriggerEnabled[3]))) {
    G4Exception("DetectorConstruction::Construct",
                "HODO_TRIGGER_COINCIDENCE_DISABLED",FatalException,
                "The selected triggerCoincidenceMode requires the optional "
                "Trigger3/Trigger4 volumes to be enabled.");
  }
  if(fGate6ElectronicsProbeMode &&
     (fGeometryProbeMode || fGenerationProbeMode || fCosmicBeam ||
      fFixedPhotonTest || fGate4OpticalAudit || !fRecordChannelNtuple ||
      fRecordPhotonAuditNtuple || fRecordPhotonNtuple))
    G4Exception("DetectorConstruction::Construct",
                "HODO_GATE6_PROBE_CONTRACT",FatalException,
                "gate6ElectronicsProbeMode requires geometry/generation/"
                "cosmic/fixed-photon/Gate-4 modes OFF, recordChannelNtuple "
                "true, and photon/photon-audit ntuples false.");
  if(fGate6InjectedChannel < -1 || fGate6InjectedChannel >= 128 ||
     fGate6InjectedPe < 0 || fGate6InjectedPe > 100000 ||
     !std::isfinite(fGate6InjectedStartTime) ||
     !std::isfinite(fGate6InjectedSpacing) ||
     fGate6InjectedSpacing < 0. || fGate6InjectedSpacing > 1000.*ns ||
     std::abs(fGate6InjectedStartTime) > 1000.*ns)
    G4Exception("DetectorConstruction::Construct",
                "HODO_GATE6_PROBE_INPUT",FatalException,
                "Gate-6 injection requires channel -1..127, p.e. 0..100000, "
                "finite start time within +/-1000 ns, and spacing 0..1000 ns.");
  // generationOnlyKill was introduced for the frozen Gate-3 source-table
  // probe, but it is also useful for final-model cosmic diagnostics where the
  // charged-particle path and scintillation photon production must be kept
  // while expensive optical-photon transport is intentionally suppressed.  The
  // flag is false by default and therefore does not affect normal production.
  if(fGenerationProbeMode) {
    const G4bool digestIsHex=
        fGenerationProbeSha256.size()==64 &&
        std::all_of(fGenerationProbeSha256.begin(),
                    fGenerationProbeSha256.end(),[](char value) {
          return (value>='0'&&value<='9') ||
                 (value>='a'&&value<='f') ||
                 (value>='A'&&value<='F');
        });
    if(!fRecordGenerationNtuples || fGenerationProbeFile.empty() ||
       !digestIsHex || fCosmicBeam || fFixedPhotonTest)
      G4Exception("DetectorConstruction::Construct",
                  "HODO_GATE3_GENERATION_PROBE_CONTRACT",FatalException,
                  "generationProbeMode requires generation ntuples, a "
                  "non-empty TSV, an exact 64-hex SHA-256, cosmicBeam=false, "
                  "and fixedPhotonTest=false.");
  }
  if(fFibreMaxStep>0. && !fFibreUserLimits)
    fFibreUserLimits=new G4UserLimits(fFibreMaxStep);
  if(fFixedPhotonTest) {
    const G4bool validPlane=fFixedPhotonPlane=="h1x"||
        fFixedPhotonPlane=="h1y"||fFixedPhotonPlane=="h2x"||
        fFixedPhotonPlane=="h2y";
    if(!validPlane)
      G4Exception("DetectorConstruction::Construct",
                  "HODO_V21_FIXED_PHOTON_PLANE",FatalException,
                  "fixedPhotonPlane must be h1x, h1y, h2x, or h2y.");
    if(fFixedPhotonAngularMode!="isotropic"&&
       fFixedPhotonAngularMode!="forward"&&
       fFixedPhotonAngularMode!="axis"&&
       fFixedPhotonAngularMode!="polar")
      G4Exception("DetectorConstruction::Construct",
                  "HODO_V21_FIXED_PHOTON_ANGLE",FatalException,
                  "fixedPhotonAngularMode must be isotropic, forward, axis, "
                  "or polar.");
    if(fFixedPhotonSource!="active"&&fFixedPhotonSource!="terminal")
      G4Exception("DetectorConstruction::Construct",
                  "HODO_GATE1_FIXED_PHOTON_SOURCE",FatalException,
                  "fixedPhotonSource must be active or terminal.");
    if(fFixedPhotonChannel < -1 || fFixedPhotonChannel>=32)
      G4Exception("DetectorConstruction::Construct",
                  "HODO_V21_FIXED_PHOTON_CHANNEL",FatalException,
                  "fixedPhotonChannel must be -1 or in [0,31].");
    if(fFixedPhotonWavelength<=0.)
      G4Exception("DetectorConstruction::Construct",
                  "HODO_V21_FIXED_PHOTON_WAVELENGTH",FatalException,
                  "fixedPhotonWavelength must be positive.");
    if(fFixedPhotonAngularMode=="polar"&&
       (fFixedPhotonPolarAngle<0.||fFixedPhotonPolarAngle>=89.9*deg))
      G4Exception("DetectorConstruction::Construct",
                  "HODO_GATE4_FIXED_PHOTON_POLAR",FatalException,
                  "fixedPhotonPolarAngle must be in [0,89.9) deg.");
    const auto sourceLimit=fFanoutPreview ? 15.9*mm : 82.4*mm;
    if(std::abs(fFixedPhotonLongitudinalOffset)>sourceLimit)
      G4Exception("DetectorConstruction::Construct",
                  "HODO_GATE4_FIXED_PHOTON_OFFSET",FatalException,
                  "fixedPhotonLongitudinalOffset lies outside the active "
                  "fibre fixture.");
    if(!fFanoutPreview&&fFixedPhotonSource!="active")
      G4Exception("DetectorConstruction::Construct",
                  "HODO_V21_FIXED_PHOTON_GEOMETRY",FatalException,
                  "A terminal fixed-photon source requires "
                  "fanoutPreview=true; Gate-4 straight active sources do not.");
  }
  if(fGate4BoundaryTrace&&
     (!fGate4OpticalAudit||!fFixedPhotonTest||
      !fRecordPhotonAuditNtuple))
    G4Exception("DetectorConstruction::Construct",
                "HODO_GATE4_BOUNDARY_TRACE_CONTRACT",FatalException,
                "gate4BoundaryTrace requires gate4OpticalAudit=true, "
                "fixedPhotonTest=true, and recordPhotonAuditNtuple=true.");
  if(fGate4SpectrumOnlyKill&&
     (!fGate4OpticalAudit||fFixedPhotonTest||
      !fRecordPhotonAuditNtuple||!fRecordPhotonAuditRows))
    G4Exception("DetectorConstruction::Construct",
                "HODO_GATE4_SPECTRUM_KILL_CONTRACT",FatalException,
                "gate4SpectrumOnlyKill requires gate4OpticalAudit=true, "
                "fixedPhotonTest=false, recordPhotonAuditNtuple=true, and "
                "recordPhotonAuditRows=true.");
  if(fRecordSpatialClosureNtuples && !fFanoutPreview)
    G4Exception("DetectorConstruction::Construct",
                "HODO_GATE1_SPATIAL_FANOUT_CONTRACT",FatalException,
                "recordSpatialClosureNtuples requires fanoutPreview=true "
                "so active/fanout fibre truth is meaningful.");
  if(fRecordSpatialClosureNtuples && !fPathOnlyMuonTruth &&
     (!fRecordPhotonAuditNtuple || !fRecordPhotonAuditRows ||
      !fRecordPhotonNtuple || !fRecordChannelNtuple))
    G4Exception("DetectorConstruction::Construct",
                "HODO_GATE1_SPATIAL_AUDIT_CONTRACT",FatalException,
                "recordSpatialClosureNtuples requires fanoutPreview, "
                "recordPhotonNtuple, "
                "recordChannelNtuple, recordPhotonAuditNtuple, and "
                "recordPhotonAuditRows all true so the relational join is "
                "complete. Set pathOnlyMuonTruth=true only for compact "
                "muon-path structural-acceptance runs with photon production "
                "disabled.");
  if(fRecordGeometryControlNtuples && !fRecordSpatialClosureNtuples)
    G4Exception("DetectorConstruction::Construct",
                "HODO_GATE2_GEOMETRY_AUDIT_CONTRACT",FatalException,
                "recordGeometryControlNtuples requires "
                "recordSpatialClosureNtuples=true so trigger and fibre truth "
                "share one event key.");
  if(fGeometryProbeMode &&
     (!fRecordGeometryControlNtuples || fGeometryProbeFile.empty() ||
      fGeometryProbeRowOffset<0 ||
      fGeometryProbeRowOffset>=kGeometryProbeExpectedRows ||
      fCosmicBeam || fFixedPhotonTest))
    G4Exception("DetectorConstruction::Construct",
                "HODO_GATE2_GEOMETRY_PROBE_CONTRACT",FatalException,
                "geometryProbeMode requires geometry-control diagnostics, a "
                "non-empty canonical TSV, row offset in [0,532], "
                "cosmicBeam=false, and fixedPhotonTest=false.");
  if (fSiPMModel != "S14160-1310PS") {
    G4ExceptionDescription message;
    message << "V15 supports only the confirmed S14160-1310PS, got: " << fSiPMModel;
    G4Exception("DetectorConstruction::Construct", "HODO_V21_SIPM", FatalException, message);
  }
  if (fH1CladdingModel != "single" && fH1CladdingModel != "double") {
    G4ExceptionDescription message;
    message << "h1CladdingModel must be single or double, got: " << fH1CladdingModel;
    G4Exception("DetectorConstruction::Construct", "HODO_V21_CLADDING", FatalException, message);
  }
  if (fH2CladdingModel != "single" && fH2CladdingModel != "double") {
    G4ExceptionDescription message;
    message << "h2CladdingModel must be single or double, got: " << fH2CladdingModel;
    G4Exception("DetectorConstruction::Construct", "HODO_V21_H2_CLADDING",
                FatalException, message);
  }
  if (fH1EmissionModel != "legacy5" && fH1EmissionModel != "scsf78_22") {
    G4ExceptionDescription message;
    message << "h1EmissionModel must be legacy5 or scsf78_22, got: " << fH1EmissionModel;
    G4Exception("DetectorConstruction::Construct", "HODO_V21_EMISSION", FatalException, message);
  }
  if (fH1RoundSides != 16 && fH1RoundSides != 64) {
    G4ExceptionDescription message;
    message << "roundTessellatedSides must be 16 or 64, got: "
            << fH1RoundSides;
    G4Exception("DetectorConstruction::Construct", "HODO_V21_ROUND_SIDES",
                FatalException, message);
  }
  if (fH1ActiveFacetRefinement!=1 &&
      fH1ActiveFacetRefinement!=2 &&
      fH1ActiveFacetRefinement!=4 &&
      fH1ActiveFacetRefinement!=8) {
    G4ExceptionDescription message;
    message << "roundActiveFacetRefinement must be 1, 2, 4, or 8; got: "
            << fH1ActiveFacetRefinement;
    G4Exception("DetectorConstruction::Construct",
                "HODO_V21_ACTIVE_SUBDIVISION",FatalException,message);
  }
  if (fH1FanoutHierarchy != "nested" &&
      fH1FanoutHierarchy != "sibling_shells") {
    G4ExceptionDescription message;
    message << "roundCladdingHierarchy must be nested or sibling_shells, got: "
            << fH1FanoutHierarchy;
    G4Exception("DetectorConstruction::Construct", "HODO_V21_ROUND_HIERARCHY",
                FatalException, message);
  }
  if (fH2FanoutHierarchy != "nested" &&
      fH2FanoutHierarchy != "sibling_shells") {
    G4ExceptionDescription message;
    message << "squareCladdingHierarchy must be nested or sibling_shells, got: "
            << fH2FanoutHierarchy;
    G4Exception("DetectorConstruction::Construct",
                "HODO_V21_SQUARE_HIERARCHY", FatalException, message);
  }
  const G4bool nativeRoundSweep =
      fH1RoundSweepRepresentation == "native_miter_union" ||
      fH1RoundSweepRepresentation == "native_capsule_union" ||
      fH1RoundSweepRepresentation == "native_miter_subtraction_union" ||
      fH1RoundSweepRepresentation == "native_miter_public_union" ||
      fH1RoundSweepRepresentation == "native_balanced_union";
  if (fH1RoundSweepRepresentation != "tessellated" &&
      !nativeRoundSweep) {
    G4ExceptionDescription message;
    message << "roundSweepRepresentation must be tessellated, "
            << "native_miter_union, native_capsule_union, or "
            << "native_miter_subtraction_union/native_miter_public_union/"
            << "native_balanced_union; "
            << "got: "
            << fH1RoundSweepRepresentation;
    G4Exception("DetectorConstruction::Construct",
                "HODO_V21_ROUND_REPRESENTATION", FatalException, message);
  }
  if (!std::isfinite(fH1RoundMiterJointOverlap) ||
      fH1RoundMiterJointOverlap < 0. ||
      fH1RoundMiterJointOverlap > 1.e-3*mm) {
    G4ExceptionDescription message;
    message << "roundMiterJointOverlap must be finite in [0,0.001] mm; got: "
            << fH1RoundMiterJointOverlap/mm << " mm";
    G4Exception("DetectorConstruction::Construct",
                "HODO_V21_MITER_OVERLAP", FatalException, message);
  }
  if ((fH1RoundSweepRepresentation=="native_miter_union" ||
       fH1RoundSweepRepresentation=="native_miter_public_union" ||
       fH1RoundSweepRepresentation=="native_balanced_union") &&
      fH1RoundMiterJointOverlap!=0.) {
    G4Exception("DetectorConstruction::Construct",
                "HODO_V21_DIRECT_ANNULAR_OVERLAP", FatalException,
                "Direct-annular native_miter_union requires exact 0-nm "
                "joints: axial overlap can create cross-material overlap. "
                "Use native_miter_subtraction_union only as the explicit "
                "numerical-overlap control.");
  }
  if (nativeRoundSweep &&
      (fH1RoundSides != 64 ||
       fH1ActiveFacetRefinement != 1 ||
       fH1FanoutHierarchy != "sibling_shells" ||
       fBendOvality != 0. || fCladdingAsymmetry != 0.)) {
    G4ExceptionDescription message;
    message << "Native round-sweep OAT is fail-closed: it requires "
            << "roundTessellatedSides=64, "
            << "roundActiveFacetRefinement=1, "
            << "roundCladdingHierarchy=sibling_shells, bendOvality=0, and "
            << "claddingAsymmetry=0; got sides=" << fH1RoundSides
            << " active_facet_refinement=" << fH1ActiveFacetRefinement
            << " hierarchy=" << fH1FanoutHierarchy
            << " ovality=" << fBendOvality
            << " asymmetry=" << fCladdingAsymmetry;
    G4Exception("DetectorConstruction::Construct",
                "HODO_V21_NATIVE_OAT_DOMAIN", FatalException, message);
  }
  if(fH2SquareFibreCount!=29 && fH2SquareFibreCount!=32)
    G4Exception("DetectorConstruction::Construct","HODO_H2_FIBRE_COUNT",FatalException,
                "h2SquareFibreCount must be 29 or 32.");
  if(fH2SquareFibreCount==29 && (fH2PhysicalRound || !fFanoutPreview || fFanoutPreviewLite ||
      !fActivePlaneOverlap ||
      fPlaneLayout[2]!="single1x32" || fPlaneLayout[3]!="single1x32" ||
      fPlaneChannelOrder[2]!="ascending" || fPlaneChannelOrder[3]!="ascending"))
    G4Exception("DetectorConstruction::Construct","HODO_H2_29_DOMAIN",FatalException,
                "29-fibre validation supports full fanout, active-overlap, "
                "square single-row ascending mapping only. Other modes need separate validation.");
  if(fFixedPhotonTest && fFixedPhotonChannel>=0) {
    const G4int plane=fFixedPhotonPlane=="h1x"?0:fFixedPhotonPlane=="h1y"?1:
                      fFixedPhotonPlane=="h2x"?2:3;
    if(!HasPhysicalFibre(plane,fFixedPhotonChannel))
      G4Exception("DetectorConstruction::Construct","HODO_PHOTON_MISSING_FIBRE",FatalException,
                  "Fixed-photon source refers to a slot without a physical fibre.");
  }
  for (G4int plane=0; plane<4; ++plane) {
    const auto& layout=fPlaneLayout[plane];
    if(layout!="single1x32"&&layout!="2x16_aligned"&&
       layout!="2x16_staggered") {
      G4ExceptionDescription message;
      message << "plane " << plane << " layout must be single1x32, "
              << "2x16_aligned, or 2x16_staggered; got: " << layout;
      G4Exception("DetectorConstruction::Construct",
                  "HODO_V21_PLANE_LAYOUT",FatalException,message);
    }
    const auto& order=fPlaneChannelOrder[plane];
    if(order!="ascending"&&order!="serpentine"&&order!="descending") {
      G4ExceptionDescription message;
      message << "plane " << plane << " channel order must be ascending, "
              << "serpentine, or descending; got: " << order;
      G4Exception("DetectorConstruction::Construct",
                  "HODO_V21_CHANNEL_ORDER",FatalException,message);
    }
    if(!std::isfinite(fPlaneRowZSpacing[plane])||
       fPlaneRowZSpacing[plane]<=0.||
       fPlaneRowZSpacing[plane]>5.*mm) {
      G4ExceptionDescription message;
      message << "plane " << plane
              << " row-Z spacing must be finite in (0,5] mm; got: "
              << fPlaneRowZSpacing[plane]/mm << " mm";
      G4Exception("DetectorConstruction::Construct",
                  "HODO_V21_ROW_Z_SPACING",FatalException,message);
    }
  }
  if (fH1AttenuationLength <= 0.) {
    G4Exception("DetectorConstruction::Construct", "HODO_V21_ATTENUATION",
                FatalException, "h1AttenuationLength must be positive.");
  }
  if (fDiagnosticRootFile.empty()) {
    G4Exception("DetectorConstruction::Construct", "HODO_V21_OUTPUT",
                FatalException, "diagnosticRootFile must not be empty.");
  }
  if (fMicrobendPeriod<=0. || fFiberSurfaceRoughness<0. ||
      fEndFaceRoughness<0. || fBendOvality<0. || fBendOvality>=.25 ||
      fCladdingAsymmetry<0. || fCladdingAsymmetry>1. ||
      fEndFaceEfficiency<0. || fEndFaceEfficiency>1. ||
      fFanoutReferenceTransmissionH1<=0. ||
      fFanoutReferenceTransmissionH1>1. ||
      fFanoutReferenceTransmissionH2<=0. ||
      fFanoutReferenceTransmissionH2>1. ||
      fFanoutReferenceRadius<=0. || fFanoutCurvaturePower<0. ||
      fFanoutChannelRms<0.) {
    G4Exception("DetectorConstruction::Construct", "HODO_V21_DEGRADATION",
                FatalException, "Invalid V15 degradation parameter range.");
  }
  if (fGreaseGap <= 0. || fSiPMActiveSize <= 0. ||
      std::abs(fSiPMOffsetU)+0.999*mm/2. > fSiPMActiveSize/2. ||
      std::abs(fSiPMOffsetV)+0.999*mm/2. > fSiPMActiveSize/2.) {
    G4Exception("DetectorConstruction::Construct", "HODO_V21_COUPLING",
                FatalException,
                "Invalid V15 grease gap/active size/alignment: the 0.999 mm "
                "fibre face must remain inside the SiPM active face.");
  }
  if (fCosmicAngularDistribution=="cos_power" &&
      (!std::isfinite(fCosmicZenithPower) || fCosmicZenithPower <= -1.)) {
    G4Exception("DetectorConstruction::Construct",
                "HODO_COSMIC_ANGULAR_POWER", FatalException,
                "cosmicZenithPower must be finite and greater than -1 for "
                "cos_power sampling.");
  }
  if (fCosmicAngularDistribution=="modified_gaisser_joint" &&
      fCosmicEnergySpectrum!="modified_gaisser") {
    G4Exception("DetectorConstruction::Construct",
                "HODO_COSMIC_JOINT_REQUIRES_GAISSER", FatalException,
                "modified_gaisser_joint requires cosmicEnergySpectrum "
                "modified_gaisser.");
  }
  const auto validYaw = [](G4double yaw) {
    return std::isfinite(yaw) && std::abs(yaw) <= 60.*deg;
  };
  if (!std::all_of(fTriggerYaw.begin(), fTriggerYaw.end(), validYaw) ||
      !std::all_of(fPlaneYaw.begin(), fPlaneYaw.end(), validYaw)) {
    G4Exception("DetectorConstruction::Construct", "HODO_V21_YAW",
                FatalException,
                "All V15 yaw deltas must be finite and within +/-60 deg.");
  }
  const auto validOffset=[](G4double offset) {
    return std::isfinite(offset)&&std::abs(offset)<=20.*mm;
  };
  if(!std::all_of(fHodoscopeOffsetX.begin(),fHodoscopeOffsetX.end(),
                  validOffset)||
     !std::all_of(fHodoscopeOffsetY.begin(),fHodoscopeOffsetY.end(),
                  validOffset)) {
    G4Exception("DetectorConstruction::Construct","HODO_V21_MODULE_OFFSET",
                FatalException,
                "H1/H2 module offsets must be finite within +/-20 mm.");
  }
  if(fAlignmentReference!="legacy"&&fAlignmentReference!="ideal_axis") {
    G4Exception("DetectorConstruction::Construct","HODO_ALIGNMENT_REFERENCE",
                FatalException,
                "alignmentReference must be legacy or ideal_axis.");
  }
  const auto validTriggerOffset=[](G4double offset) {
    return std::isfinite(offset)&&std::abs(offset)<=60.*mm;
  };
  if(!std::all_of(fTriggerOffsetX.begin(),fTriggerOffsetX.end(),
                  validTriggerOffset)||
     !std::all_of(fTriggerOffsetY.begin(),fTriggerOffsetY.end(),
                  validTriggerOffset)) {
    G4Exception("DetectorConstruction::Construct","HODO_TRIGGER_OFFSET",
                FatalException,
                "Trigger offsets must be finite within +/-60 mm.");
  }
  if(fTriggerLength<=0.||fTriggerWidth<=0.||fTriggerThickness<=0.||
     fTriggerLength>200.*mm||fTriggerWidth>200.*mm||
     fTriggerThickness>50.*mm) {
    G4Exception("DetectorConstruction::Construct","HODO_TRIGGER_SIZE",
                FatalException,
                "Trigger full dimensions must be positive and within "
                "200 x 200 x 50 mm.");
  }
  const auto validStructuralZ=[](G4double value) {
    return std::isfinite(value)&&std::abs(value)<=1000.*mm;
  };
  if(!std::all_of(fTriggerZ.begin(),fTriggerZ.end(),validStructuralZ)||
     !std::all_of(fHodoscopeFanoutZ.begin(),fHodoscopeFanoutZ.end(),
                  validStructuralZ)) {
    G4Exception("DetectorConstruction::Construct","HODO_STRUCTURAL_Z",
                FatalException,
                "Trigger and fan-out preview centre Z positions must be "
                "finite within +/-1000 mm.");
  }
  if(fActivePlaneOverlap&&(!fFanoutPreview||
     std::any_of(fPlaneYaw.begin(),fPlaneYaw.end(),[](G4double yaw){return std::abs(yaw)>1.e-12;}))) {
    G4Exception("DetectorConstruction::Construct","HODO_ACTIVE_OVERLAP",FatalException,
                "activePlaneOverlap requires fanoutPreview with zero individual plane yaw.");
  }
  if((fActivePlaneOverlap||fFibreFarEndExtension>0.)&&
     (fFanoutPreviewLite||fGenerationProbeMode)) {
    G4Exception("DetectorConstruction::Construct","HODO_NEW_REGION_SCOPE",FatalException,
                "New overlap/tail geometry requires the full geometry, not the lite display. "
                "The legacy Gate3 generation-probe oracle has not been migrated to these regions.");
  }
  if(!std::isfinite(fFibreFarEndExtension)||fFibreFarEndExtension<0.||
     fFibreFarEndExtension>200.*mm||(!fFanoutPreview&&fFibreFarEndExtension>0.)) {
    G4Exception("DetectorConstruction::Construct","HODO_FAR_END_EXTENSION",
                FatalException,"fibreFarEndExtension must be finite in [0,200] mm and requires fanoutPreview.");
  }
  if(fFanoutPreviewLite &&
     (fFanoutPreviewLiteChannelStride<1 ||
      fFanoutPreviewLiteChannelStride>32 ||
      fFanoutPreviewLiteSegments<2 ||
      fFanoutPreviewLiteSegments>64)) {
    G4Exception("DetectorConstruction::Construct","HODO_FANOUT_LITE",
                FatalException,
                "fanoutPreviewLite requires channel stride in [1,32] and "
                "segments in [2,64].");
  }
  if(!std::isfinite(fFibreOverlapClearance) ||
     fFibreOverlapClearance < 0. ||
     fFibreOverlapClearance > 0.015*mm) {
    G4Exception("DetectorConstruction::Construct",
                "HODO_FIBRE_OVERLAP_CLEARANCE", FatalException,
                "fibreOverlapClearance must be finite in [0,0.015] mm.");
  }
  G4cout << "[optical-v21 config] SiPM=" << fSiPMModel
         << " cells=" << GetMicrocellCount()
         << " DCR_kHz=" << fDarkCountRate/(1.e3/second)
         << " gate6_probe=" << fGate6ElectronicsProbeMode
         << " gate6_channel=" << fGate6InjectedChannel
         << " gate6_pe=" << fGate6InjectedPe
         << " gate6_start_ns=" << fGate6InjectedStartTime/ns
         << " gate6_spacing_ns=" << fGate6InjectedSpacing/ns
         << " H1_att_m=" << fH1AttenuationLength/m
         << " H1_clad=" << fH1CladdingModel
         << " H1_emission=" << fH1EmissionModel
         << " round_tessellated_sides=" << fH1RoundSides
         << " round_active_facet_refinement="
         << fH1ActiveFacetRefinement
         << " round_polygon_phase_rad=" << GetH1RoundPolygonPhase()
         << " round_cladding_hierarchy=" << fH1FanoutHierarchy
         << " square_cladding_hierarchy=" << fH2FanoutHierarchy
         << " trigger_scint_yield=" << fTriggerScintillationYield
         << " round_sweep_representation=" << fH1RoundSweepRepresentation
         << " round_miter_joint_overlap_nm="
         << fH1RoundMiterJointOverlap/nm
         << " layout=(H1X:" << fPlaneLayout[0]
         << ",H1Y:" << fPlaneLayout[1]
         << ",H2X:" << fPlaneLayout[2]
         << ",H2Y:" << fPlaneLayout[3] << ")"
         << " h2_physical_round=" << fH2PhysicalRound
         << " channel_order=(H1X:" << fPlaneChannelOrder[0]
         << ",H1Y:" << fPlaneChannelOrder[1]
         << ",H2X:" << fPlaneChannelOrder[2]
         << ",H2Y:" << fPlaneChannelOrder[3] << ")"
         << " row_z_mm=(" << fPlaneRowZSpacing[0]/mm << ","
         << fPlaneRowZSpacing[1]/mm << ","
         << fPlaneRowZSpacing[2]/mm << ","
         << fPlaneRowZSpacing[3]/mm << ")"
         << " mixed_layout_pattern_ASSUMED=" << GetMixedLayoutPattern()
         << " active_core_chords=" << fRecordActiveCoreChordBranches
         << " grease_gap_mm=" << fGreaseGap/mm
         << " sipm_active_mm=" << fSiPMActiveSize/mm
         << " sipm_offset_uv_mm=(" << fSiPMOffsetU/mm << ","
         << fSiPMOffsetV/mm << ")"
         << " endface_sigma_alpha=" << fEndFaceRoughness
         << " yaw_deg=(T1:" << fTriggerYaw[0]/deg
         << ",T2:" << fTriggerYaw[1]/deg
         << ",H1X:" << fPlaneYaw[0]/deg
         << ",H1Y:" << fPlaneYaw[1]/deg
         << ",H2X:" << fPlaneYaw[2]/deg
         << ",H2Y:" << fPlaneYaw[3]/deg << ")"
         << " module_offset_mm=(H1:" << fHodoscopeOffsetX[0]/mm
         << "," << fHodoscopeOffsetY[0]/mm
         << ";H2:" << fHodoscopeOffsetX[1]/mm
         << "," << fHodoscopeOffsetY[1]/mm << ")"
         << " h2_module_rotate180=" << fH2ModuleRotate180
         << " alignment_reference=" << fAlignmentReference
         << " trigger_offset_mm=(T1:" << fTriggerOffsetX[0]/mm
         << "," << fTriggerOffsetY[0]/mm
         << ";T2:" << fTriggerOffsetX[1]/mm
         << "," << fTriggerOffsetY[1]/mm << ")"
         << " fanout_linear=" << fFanoutLinearFanout
         << " fibre_far_end_extension_mm=" << fFibreFarEndExtension/mm
         << " active_plane_overlap=" << fActivePlaneOverlap
         << " fanout_preview_lite=" << fFanoutPreviewLite
         << " fanout_preview_lite_stride="
         << fFanoutPreviewLiteChannelStride
         << " fanout_preview_lite_segments="
         << fFanoutPreviewLiteSegments
         << " fibre_overlap_clearance_mm="
         << fFibreOverlapClearance/mm
         << " root=" << fDiagnosticRootFile << G4endl;
  DefineMaterials();
  if(fGate4OpticalAudit) {
    for(const auto* property:{"RINDEX","ABSLENGTH",
                              "SCINTILLATIONCOMPONENT1"}) {
      PrintGate4MptProperty("H1",fSCSF78,property);
      PrintGate4MptProperty("H2",fBCF20,property);
    }
    for(const auto* property:{"RINDEX","ABSLENGTH"}) {
      PrintGate4MptProperty("INNER_CLAD",fInnerCladding,property);
      PrintGate4MptProperty("OUTER_CLAD",fOuterCladding,property);
    }
    G4cout << "GATE4_MPT_AUDIT_COMPLETE" << G4endl;
  }
  DefineSurfaces();
  auto* air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  const auto worldXY = fUseBuilding ? 2.5*m : 250.*mm;
  const auto worldZ = fUseBuilding ? 2.5*m : 250.*mm;
  if(fUnselectedEventStride<1 ||
     (fCosmicPositionModel=="uniform_source_plane" &&
      (fCosmicSourceHalfX<=0. || fCosmicSourceHalfY<=0. ||
       fCosmicSourceHalfX>=worldXY || fCosmicSourceHalfY>=worldXY)))
    G4Exception("DetectorConstruction::Construct", "HODO_SOURCE_PLANE_RANGE",
                FatalException,"Source half sizes must be positive and strictly inside world; event stride must be >=1.");
  G4cout << "COSMIC_POSITION model=" << fCosmicPositionModel
         << " source_half_x_mm=" << fCosmicSourceHalfX/mm
         << " source_half_y_mm=" << fCosmicSourceHalfY/mm
         << " unselected_event_stride=" << fUnselectedEventStride
         << " skip_untriggered_digitization=" << fSkipUntriggeredDigitization << G4endl;
  auto* worldSolid = new G4Box("World", worldXY, worldXY, worldZ);
  auto* world = new G4LogicalVolume(worldSolid, air, "WorldLV");
  auto* worldPV = new G4PVPlacement(nullptr, {}, world, "World", nullptr, false, 0, true);
  world->SetVisAttributes(G4VisAttributes::GetInvisible());

  BuildTrigger(world, 0);
  if (fTriggerEnabled[2]) BuildTrigger(world, 2);
  if (fFanoutPreview) {
    if (fFanoutPreviewLite) {
      BuildFanoutPreviewLite(world);
    } else {
      BuildFanoutPreview(world);
    }
  } else {
    BuildRoundHodoscope(world, +2.5*mm);
    if (fH2PhysicalRound) BuildRoundHodoscopeH2(world, -2.5*mm);
    else BuildSquareHodoscope(world, -2.5*mm);
  }
  if (fTriggerEnabled[3]) BuildTrigger(world, 3);
  BuildTrigger(world, 1);
  if (fUseBuilding) PlaceUpstreamMaterials(world);
  if(fH2SquareFibreCount==29) {
    std::array<G4int,4> counts{{0,0,0,0}};
    const std::array<G4String,4> names{{"H1_X","H1_Y","H2_X","H2_Y"}};
    for(G4int i=0;i<world->GetNoDaughters();++i)
      for(G4int plane=0;plane<4;++plane)
        if(world->GetDaughter(i)->GetName()==names[plane]) {
          ++counts[plane];
          const auto channel=world->GetDaughter(i)->GetCopyNo()%32;
          if(!HasPhysicalFibre(plane,channel))
            G4Exception("DetectorConstruction::Construct","HODO_H2_29_COPY",FatalException,
                        "A removed fibre was physically placed.");
        }
    if(counts!=std::array<G4int,4>{{32,32,29,29}})
      G4Exception("DetectorConstruction::Construct","HODO_H2_29_PLACEMENTS",FatalException,
                  "Physical daughter counts are not 32/32/29/29.");
    G4cout<<"H2_29_PHYSICAL_PLACEMENTS=PASS counts=32,32,29,29 "
          <<"mapping=PROVISIONAL_LOCAL_0_TO_28 readout_slots=32_per_mid"<<G4endl;
    // Independent navigator controls at known packed-core centres and beyond
    // the 29-fibre envelope. Only enabled for the new path-only configuration.
    G4Navigator auditNavigator;
    auditNavigator.SetWorldVolume(worldPV);
    G4int controls=0;
    for(G4int plane:{2,3}) {
      for(G4double u:{0.*mm,-30.*mm}) {
        for(G4double v:{-14.*mm,0.*mm,14.*mm,-15.*mm,15.*mm}) {
          const auto z=GetPlaneSourceZ(plane,0);
          const auto local=plane==2 ? G4ThreeVector(u,v,z) : G4ThreeVector(v,u,z);
          const auto p=PlaneNominalToWorldPosition(plane,local);
          auto* volume=auditNavigator.LocateGlobalPointAndSetup(p,nullptr,false);
          const auto core=volume && volume->GetLogicalVolume()->GetMaterial()==fBCF20;
          const auto expected=std::abs(v)<=14.*mm;
          if(core!=expected)
            G4Exception("DetectorConstruction::Construct","HODO_H2_29_NAV",FatalException,
                        "29-fibre packed/tail core versus outside navigator control failed.");
          ++controls;
        }
      }
    }
    G4cout<<"H2_29_NAVIGATOR_CONTROLS=PASS checks="<<controls<<G4endl;
  }
  return worldPV;
}

void DetectorConstruction::BuildFanoutPreview(G4LogicalVolume* world) {
  constexpr G4double straightExitLength = 20.*mm;
  constexpr G4double fanoutLength = 100.*mm;
  // Photo-derived preview assumption: the SiPM-side 4x8 grid is wider than
  // the tightly packed detector-side fibre array.
  constexpr G4double pitch = 4.*mm;
  constexpr G4double activeHalfLength = 16.*mm;
  constexpr G4int segments = 47;
  struct CouplingEnvelope {
    G4int plane;
    G4int channel;
    G4String prefix;
    G4ThreeVector centre;
    std::array<G4ThreeVector,3> axes;
    std::array<G4double,3> halfLengths;
  };
  std::vector<CouplingEnvelope> couplingEnvelopes;
  couplingEnvelopes.reserve(128);
  G4int terminalFrameAuditCount = 0;
  G4int h2OuterBoundaryAuditCount = 0;
  G4int h2OuterBoundaryCloneAuditCount = 0;
  G4int h2OuterFrozenRayAuditCount = 0;
  G4double minimumSamePlaneClearance=DBL_MAX;
  G4double minimumDifferentPlaneClearance=DBL_MAX;

  const auto buildAssembly = [&](G4bool round, G4int axis, G4double detectorZ,
                                 const G4String& prefix, const G4Colour& colour) {
    const G4int plane = (round ? 0 : 2) + axis;
    const G4bool roundFibre = round || (!round && fH2PhysicalRound);
    std::vector<std::vector<G4ThreeVector>> priorSquareCentrelines;
    priorSquareCentrelines.reserve(roundFibre ? 0 : 32);
    // Keep each endpoint near its real detector-plane height. The former
    // separated-holder display added an artificial vertical bend.
    const G4double planeZ = round ? detectorZ+axis*2.2*mm
                                  : (roundFibre ? detectorZ-axis*2.2*mm
                                                : detectorZ-axis*1.*mm);
    for (G4int channel=0; channel<GetPlaneFibreCount(plane); ++channel) {
      // V21 plane layouts and channel order are explicit sensitivities.
      // They determine both the detector-side slot and the route into the
      // frozen 4x8 SiPM grid; no efficiency-derived remapping is applied.
      const G4int row=GetPlaneFanoutRow(plane,channel);
      const G4int col=GetPlaneFanoutColumn(plane,channel);
      const G4double sourceTransverse=
          GetPlaneSourceTransverse(plane,channel);
      const G4double sourceZ=GetPlaneSourceZ(plane,channel);
      const G4double targetTransverse=(col-3.5)*pitch;
      const G4double targetZ=planeZ+(row-1.5)*pitch;

      std::vector<G4ThreeVector> centres;
      centres.reserve(segments+3);
      const auto activeBounds=GetActiveLongitudinalBounds(plane);
      const auto farEnd=fFibreFarEndExtension>0.
          ? activeBounds[0]-fFibreFarEndExtension : -activeHalfLength;
      centres.push_back(axis == 0
        ? G4ThreeVector(farEnd,sourceTransverse,sourceZ)
        : G4ThreeVector(sourceTransverse,farEnd,sourceZ));
      centres.push_back(axis == 0
        ? G4ThreeVector(activeHalfLength,sourceTransverse,sourceZ)
        : G4ThreeVector(sourceTransverse,activeHalfLength,sourceZ));
      // Keep the whole active detector region tightly packed and straight.
      // The fan-out starts only after the fibre has left that active region,
      // matching the intended "straight detector, gradual expansion to SiPM"
      // mechanical picture.
      centres.push_back(axis == 0
        ? G4ThreeVector(activeHalfLength+straightExitLength,
                        sourceTransverse,sourceZ)
        : G4ThreeVector(sourceTransverse,
                        activeHalfLength+straightExitLength,sourceZ));
      for (G4int step=1; step<=segments; ++step) {
        const G4double t=static_cast<G4double>(step)/segments;
        // Quintic smootherstep has zero slope and zero curvature at both ends,
        // avoiding visible kinks where the active fibre and SiPM face join.
        const auto smoother=[](G4double value) {
          const G4double u=std::clamp(value,0.,1.);
          return u*u*u*(10.-15.*u+6.*u*u);
        };
        // After the straight detector exit, round fibres use the same smooth
        // factor for z-row separation and transverse spreading.  Square H2
        // fibres are different: the single-row detector side is routed into a
        // 4x8 SiPM grid, so neighbouring channels that share a terminal column
        // can otherwise converge transversely before they have separated in z.
        // Delay only the square-fibre transverse convergence; this preserves
        // the active detector and terminal grid while avoiding an artificial
        // tessellated sweep overlap in the fan-out.
        const G4bool terminalStraight=step >= segments-1;
        const G4double smooth=
            terminalStraight ? 1. : (fFanoutLinearFanout ? t : smoother(t));
        const G4double transverseSmooth = roundFibre || terminalStraight
            ? smooth
            : (fFanoutLinearFanout
                ? std::pow(t,2.)
                : smoother((t-.25)/.75));
        const G4double longitudinal=
            activeHalfLength+straightExitLength+fanoutLength*t;
        const G4double transverse=sourceTransverse+
            (targetTransverse-sourceTransverse)*transverseSmooth;
        const G4double curveZ=sourceZ+(targetZ-sourceZ)*smooth;
        const G4double microWeight=std::sin(CLHEP::pi*t);
        const G4double microPhase=2.*CLHEP::pi*fanoutLength*t/
            std::max(fMicrobendPeriod,1.e-6*mm);
        const G4double microOffset=terminalStraight ? 0. :
            fMicrobendAmplitude*microWeight*std::sin(microPhase);
        centres.push_back(axis == 0
          ? G4ThreeVector(longitudinal,transverse+microOffset,curveZ)
          : G4ThreeVector(transverse+microOffset,longitudinal,curveZ));
      }

      // Geometry-derived minimum local curvature radius. The transmission
      // model below is a sensitivity model, not a manufacturer calibration.
      G4double minRadius=DBL_MAX;
      for(std::size_t station=1;station+1<centres.size();++station) {
        const auto ab=(centres[station]-centres[station-1]).mag();
        const auto bc=(centres[station+1]-centres[station]).mag();
        const auto ac=(centres[station+1]-centres[station-1]).mag();
        const auto twiceArea=
            ((centres[station]-centres[station-1]).cross(
              centres[station+1]-centres[station-1])).mag();
        if(twiceArea>1.e-12*mm*mm)
          minRadius=std::min(minRadius,ab*bc*ac/(2.*twiceArea));
      }
      const G4int transmissionChannel=
          (round ? 0 : 64)+axis*32+channel;
      fFanoutMinRadius[transmissionChannel]=
          std::isfinite(minRadius)?minRadius:0.;
      const auto referenceTransmission=roundFibre
          ? fFanoutReferenceTransmissionH1
          : fFanoutReferenceTransmissionH2;
      const auto severity=minRadius>0.&&std::isfinite(minRadius)
          ? std::pow(fFanoutReferenceRadius/minRadius,
                     fFanoutCurvaturePower) : 1.;
      // Fixed channel hash gives reproducible assembly variation without
      // consuming the Geant4 physics random stream.
      std::uint64_t hash=0x9e3779b97f4a7c15ULL*
          static_cast<std::uint64_t>(transmissionChannel+1);
      hash=(hash^(hash>>30))*0xbf58476d1ce4e5b9ULL;
      hash=(hash^(hash>>27))*0x94d049bb133111ebULL;
      hash^=hash>>31;
      const G4double u1=std::max(1.e-12,
          static_cast<G4double>((hash&0xffffffffULL)+1)/
          static_cast<G4double>(0x100000001ULL));
      const G4double u2=static_cast<G4double>(((hash>>32)&0xffffffffULL)+1)/
          static_cast<G4double>(0x100000001ULL);
      const G4double normal=std::sqrt(-2.*std::log(u1))*
          std::cos(2.*CLHEP::pi*u2);
      const G4double nominal=std::pow(referenceTransmission,severity);
      fFanoutChannelTransmission[transmissionChannel]=std::clamp(
          nominal*std::exp(fFanoutChannelRms*normal-
                           .5*fFanoutChannelRms*fFanoutChannelRms),0.,1.);
      // Curvature and transmission above are evaluated in the frozen nominal
      // frame. Rotate only after that calculation so yaw is a rigid alignment
      // change and cannot alter an optical coefficient through round-off.
      const auto module=plane/2;
      const G4bool moduleRotation =
          module==1 && GetH2ModuleRotate180();
      if (moduleRotation ||
          GetPlaneYaw(plane) != 0. ||
          GetHodoscopeOffsetX(module) != 0. ||
          GetHodoscopeOffsetY(module) != 0.)
        for (auto& centre : centres)
          centre = PlaneNominalToWorldPosition(plane, centre);

      const G4int sides = roundFibre ? fH1RoundSides : 4;
      const G4ThreeVector verticalAxis(0.,0.,1.);
      struct StationFrame {
        G4ThreeVector tangent;
        G4ThreeVector side;
        G4ThreeVector normal;
      };
      std::vector<StationFrame> stationFrames;
      stationFrames.reserve(centres.size());
      for (std::size_t station=0; station<centres.size(); ++station) {
        const G4ThreeVector tangent =
            (station == 0 ? centres[1]-centres[0]
             : station+1 == centres.size()
                 ? centres[station]-centres[station-1]
                 : centres[station+1]-centres[station-1]).unit();
        // One deterministic right-handed frame is shared by every material
        // layer and, at the terminal station, by the complete SiPM coupling
        // stack.  This prevents independently rounded rotations from opening
        // a picometre-scale wedge at the tessellated end cap.
        G4ThreeVector sideAxis=verticalAxis.cross(tangent);
        if (sideAxis.mag2() < 1.e-18)
          sideAxis=axis==0 ? G4ThreeVector(0.,1.,0.)
                           : G4ThreeVector(-1.,0.,0.);
        sideAxis=sideAxis.unit();
        const G4ThreeVector normalAxis=tangent.cross(sideAxis).unit();
        stationFrames.push_back({tangent,sideAxis,normalAxis});
      }
      const G4double frameTolerance =
          4096.*std::numeric_limits<G4double>::epsilon();
      for (std::size_t station=0; station<stationFrames.size(); ++station) {
        const auto& frame=stationFrames[station];
        const G4double defect=std::max({
            std::abs(frame.tangent.mag2()-1.),
            std::abs(frame.side.mag2()-1.),
            std::abs(frame.normal.mag2()-1.),
            std::abs(frame.tangent.dot(frame.side)),
            std::abs(frame.tangent.dot(frame.normal)),
            std::abs(frame.side.dot(frame.normal)),
            (frame.side.cross(frame.normal)-frame.tangent).mag()});
        if (defect > frameTolerance) {
          G4ExceptionDescription message;
          message << "Non-orthonormal fan-out station frame: prefix=" << prefix
                  << " channel=" << channel << " station=" << station
                  << " defect=" << defect
                  << " allowed=" << frameTolerance;
          G4Exception("DetectorConstruction::BuildFanoutPreview",
                      "HODO_V21_TERMINAL_FRAME",FatalException,message);
        }
      }
      const auto makeSwept = [&](const G4String& name, G4double halfWidth,
                                 G4double layerShift,
                                 G4bool h2OuterBoundaryFix=false) {
        std::vector<std::vector<G4ThreeVector>> rings(
            centres.size(), std::vector<G4ThreeVector>(sides));
        for (std::size_t station=0; station<centres.size(); ++station) {
          const auto& frame=stationFrames[station];
          const G4double bendT=station < 2 ? 0. :
              static_cast<G4double>(station-1)/segments;
          const G4double bendWeight=std::sin(CLHEP::pi*bendT);
          const G4double sideScale=1.;
          const G4double normalScale=1.-fBendOvality*bendWeight;
          const G4ThreeVector layerOffset=
              layerShift*fCladdingAsymmetry*bendWeight*frame.normal;
          for (G4int side=0; side<sides; ++side) {
            const G4double angle = roundFibre
                ? 2.*CLHEP::pi*side/sides+GetH1RoundPolygonPhase()
                : (45.+90.*side)*deg;
            const G4double scale =
                roundFibre ? halfWidth : halfWidth/std::sqrt(.5);
            rings[station][side] = centres[station]+layerOffset
                + scale*sideScale*std::cos(angle)*frame.side
                + scale*normalScale*std::sin(angle)*frame.normal;
          }
        }
        G4TessellatedSolid* swept=h2OuterBoundaryFix
            ? static_cast<G4TessellatedSolid*>(
                  new H2OuterBoundaryTessellatedSolid(name))
            : (!roundFibre && name.find("CoreSolid_")!=G4String::npos
                ? static_cast<G4TessellatedSolid*>(new H2CoreBoundaryTessellatedSolid(name))
                : new G4TessellatedSolid(name));
        G4double meshArea=0.;
        G4double signedMeshVolume=0.;
        G4double referenceMeshArea=0.;
        G4double signedReferenceMeshVolume=0.;
        std::uint64_t meshFacetCount=0;
        std::uint64_t referenceFacetCount=0;
        std::uint64_t uniformlyRefinedParentCount=0;
        std::uint64_t uniformChildFacetCount=0;
        std::uint64_t boundaryFanParentCount=0;
        std::uint64_t boundaryFanChildFacetCount=0;
        std::uint64_t interiorLatticePointCount=0;
        G4double maxChildPlaneResidual=0.;
        const auto addTriangle=[&](const G4ThreeVector& a,
                                   const G4ThreeVector& b,
                                   const G4ThreeVector& c) {
          swept->AddFacet(new G4TriangularFacet(a,b,c,ABSOLUTE));
          meshArea+=.5*((b-a).cross(c-a)).mag();
          signedMeshVolume+=a.dot(b.cross(c))/6.;
          ++meshFacetCount;
        };
        const auto addReferenceTriangle=[&](
            const G4ThreeVector& a,const G4ThreeVector& b,
            const G4ThreeVector& c) {
          referenceMeshArea+=.5*((b-a).cross(c-a)).mag();
          signedReferenceMeshVolume+=a.dot(b.cross(c))/6.;
          ++referenceFacetCount;
        };
        const G4int activeRefinement=
            roundFibre ? fH1ActiveFacetRefinement : 1;
        using MeshEdgeKey=std::pair<G4int,G4int>;
        std::map<MeshEdgeKey,std::vector<G4ThreeVector>>
            canonicalEdgePointCache;
        const auto originalVertexId=[&](std::size_t station,G4int side) {
          return static_cast<G4int>(station)*sides+side;
        };
        const auto originalVertex=[&](G4int id) -> const G4ThreeVector& {
          return rings[static_cast<std::size_t>(id/sides)]
                      [static_cast<std::size_t>(id%sides)];
        };
        const auto canonicalEdgePoint=[&](
            G4int fromId,G4int toId,G4int step) {
          const G4int lowId=std::min(fromId,toId);
          const G4int highId=std::max(fromId,toId);
          const MeshEdgeKey key(lowId,highId);
          auto found=canonicalEdgePointCache.find(key);
          if(found==canonicalEdgePointCache.end()) {
            std::vector<G4ThreeVector> points(
                static_cast<std::size_t>(activeRefinement+1));
            const auto& low=originalVertex(lowId);
            const auto& high=originalVertex(highId);
            points.front()=low;
            points.back()=high;
            for(G4int index=1;index<activeRefinement;++index) {
              const G4double fraction=
                  static_cast<G4double>(index)/activeRefinement;
              points[static_cast<std::size_t>(index)]=
                  low+(high-low)*fraction;
            }
            found=canonicalEdgePointCache.emplace(key,std::move(points)).first;
          }
          const G4int canonicalStep=
              fromId==lowId ? step : activeRefinement-step;
          return found->second[static_cast<std::size_t>(canonicalStep)];
        };
        const auto notePlanePoint=[&](
            const G4ThreeVector& point,const G4ThreeVector& origin,
            const G4ThreeVector& unnormalisedNormal) {
          const G4double normalMagnitude=unnormalisedNormal.mag();
          if(normalMagnitude>0.) {
            maxChildPlaneResidual=std::max(
                maxChildPlaneResidual,
                std::abs((point-origin).dot(unnormalisedNormal))/
                    normalMagnitude);
          }
        };
        const auto addUniformTriangleRefinement=[&](
            G4int aId,G4int bId,G4int cId,
            const G4ThreeVector& a,const G4ThreeVector& b,
            const G4ThreeVector& c) {
          const auto parentNormal=(b-a).cross(c-a);
          const auto point=[&](G4int i,G4int j) {
            if(j==0)
              return canonicalEdgePoint(aId,bId,i);
            if(i==0)
              return canonicalEdgePoint(aId,cId,j);
            if(i+j==activeRefinement)
              return canonicalEdgePoint(bId,cId,j);
            ++interiorLatticePointCount;
            return a+((b-a)*static_cast<G4double>(i)+
                      (c-a)*static_cast<G4double>(j))/
                         activeRefinement;
          };
          ++uniformlyRefinedParentCount;
          for(G4int i=0;i<activeRefinement;++i) {
            for(G4int j=0;j<activeRefinement-i;++j) {
              const auto p0=point(i,j);
              const auto p1=point(i+1,j);
              const auto p2=point(i,j+1);
              notePlanePoint(p0,a,parentNormal);
              notePlanePoint(p1,a,parentNormal);
              notePlanePoint(p2,a,parentNormal);
              addTriangle(p0,p1,p2);
              ++uniformChildFacetCount;
              if(i+j<activeRefinement-1) {
                const auto p3=point(i+1,j+1);
                notePlanePoint(p3,a,parentNormal);
                addTriangle(p1,p3,p2);
                ++uniformChildFacetCount;
              }
            }
          }
        };
        const auto addBoundaryEdgeFan=[&](
            const G4ThreeVector& parentA,
            const G4ThreeVector& parentB,
            const G4ThreeVector& parentC,
            G4int edgeFromId,G4int edgeToId,
            const G4ThreeVector& apex,G4bool apexFirst) {
          const auto parentNormal=(parentB-parentA).cross(parentC-parentA);
          ++boundaryFanParentCount;
          for(G4int step=0;step<activeRefinement;++step) {
            const auto edge0=
                canonicalEdgePoint(edgeFromId,edgeToId,step);
            const auto edge1=
                canonicalEdgePoint(edgeFromId,edgeToId,step+1);
            notePlanePoint(edge0,parentA,parentNormal);
            notePlanePoint(edge1,parentA,parentNormal);
            notePlanePoint(apex,parentA,parentNormal);
            if(apexFirst)
              addTriangle(apex,edge0,edge1);
            else
              addTriangle(edge0,edge1,apex);
            ++boundaryFanChildFacetCount;
          }
        };
        for (std::size_t station=0; station+1<rings.size(); ++station) {
          for (G4int side=0; side<sides; ++side) {
            const G4int nextSide=(side+1)%sides;
            const auto& a1=rings[station][side];
            const auto& b1=rings[station+1][nextSide];
            const auto& c1=rings[station+1][side];
            const auto& a2=rings[station][side];
            const auto& b2=rings[station][nextSide];
            const auto& c2=rings[station+1][nextSide];
            addReferenceTriangle(a1,b1,c1);
            addReferenceTriangle(a2,b2,c2);
            // Preserve the frozen refinement=1 vertices and AddFacet call
            // sequence bit-for-bit.  Even algebraically equivalent
            // barycentric endpoint reconstruction can move a vertex by a few
            // femtometres and flip a tolerance-singular negative control.
            if(activeRefinement==1) {
              addTriangle(a1,b1,c1);
              addTriangle(a2,b2,c2);
            } else if(station==0) {
              const auto a1Id=originalVertexId(station,side);
              const auto b1Id=originalVertexId(station+1,nextSide);
              const auto c1Id=originalVertexId(station+1,side);
              const auto a2Id=originalVertexId(station,side);
              const auto b2Id=originalVertexId(station,nextSide);
              const auto c2Id=originalVertexId(station+1,nextSide);
              addUniformTriangleRefinement(
                  a1Id,b1Id,c1Id,a1,b1,c1);
              addUniformTriangleRefinement(
                  a2Id,b2Id,c2Id,a2,b2,c2);
            } else if(station==1) {
              // The uniformly refined active interval subdivides the
              // station-1 ring edge.  Keep the next interval's first
              // triangle unchanged and fan only the adjacent second
              // triangle within its original plane.
              addTriangle(a1,b1,c1);
              addBoundaryEdgeFan(
                  a2,b2,c2,
                  originalVertexId(station,side),
                  originalVertexId(station,nextSide),
                  c2,false);
            } else {
              addTriangle(a1,b1,c1);
              addTriangle(a2,b2,c2);
            }
          }
        }
        const auto capCentreStart=centres.front();
        const auto capCentreEnd=centres.back();
        for (G4int side=0; side<sides; ++side) {
          const G4int nextSide=(side+1)%sides;
          addReferenceTriangle(
              capCentreStart,rings.front()[nextSide],rings.front()[side]);
          addReferenceTriangle(
              capCentreEnd,rings.back()[side],rings.back()[nextSide]);
          if(activeRefinement==1) {
            addTriangle(
                capCentreStart,rings.front()[nextSide],rings.front()[side]);
          } else {
            // The station-0 surface subdivision reaches the start-ring
            // boundary.  Reuse the same canonical edge points and fan the
            // cap triangle in its original plane, eliminating T-junctions.
            addBoundaryEdgeFan(
                capCentreStart,rings.front()[nextSide],rings.front()[side],
                originalVertexId(0,nextSide),originalVertexId(0,side),
                capCentreStart,true);
          }
          addTriangle(
              capCentreEnd,rings.back()[side],rings.back()[nextSide]);
        }
        const std::uint64_t baseFacetCount=
            2ULL*static_cast<std::uint64_t>(rings.size()-1)*
                static_cast<std::uint64_t>(sides)+
            2ULL*static_cast<std::uint64_t>(sides);
        const std::uint64_t expectedFacetCount=
            activeRefinement==1 ? baseFacetCount :
            baseFacetCount+
                2ULL*static_cast<std::uint64_t>(sides)*
                (static_cast<std::uint64_t>(activeRefinement)*
                     activeRefinement+
                 static_cast<std::uint64_t>(activeRefinement)-2ULL);
        const std::uint64_t expectedCanonicalEdges=
            activeRefinement==1 ? 0ULL :
            4ULL*static_cast<std::uint64_t>(sides);
        const std::uint64_t expectedUniformParents=
            activeRefinement==1 ? 0ULL :
            2ULL*static_cast<std::uint64_t>(sides);
        const std::uint64_t expectedUniformChildren=
            expectedUniformParents*
            static_cast<std::uint64_t>(activeRefinement)*
            static_cast<std::uint64_t>(activeRefinement);
        const std::uint64_t expectedBoundaryFanParents=
            activeRefinement==1 ? 0ULL :
            2ULL*static_cast<std::uint64_t>(sides);
        const std::uint64_t expectedBoundaryFanChildren=
            expectedBoundaryFanParents*
            static_cast<std::uint64_t>(activeRefinement);
        const G4bool topologyCountInvariant=
            referenceFacetCount==baseFacetCount &&
            meshFacetCount==expectedFacetCount &&
            canonicalEdgePointCache.size()==expectedCanonicalEdges &&
            uniformlyRefinedParentCount==expectedUniformParents &&
            uniformChildFacetCount==expectedUniformChildren &&
            boundaryFanParentCount==expectedBoundaryFanParents &&
            boundaryFanChildFacetCount==expectedBoundaryFanChildren;
        const G4double childPlaneTolerance=1.e-6*nm;
        const G4bool childPlaneInvariant=
            maxChildPlaneResidual<=childPlaneTolerance;
        swept->SetSolidClosed(true);
        if(h2OuterBoundaryFix) {
          auto* corrected=dynamic_cast<
              H2OuterBoundaryTessellatedSolid*>(swept);
          if(round||sides!=4||!corrected||
             swept->GetNumberOfFacets()!=400||meshFacetCount!=400||
             referenceFacetCount!=400||activeRefinement!=1||
             swept->GetEntityType()!="G4TessellatedSolid"||
             swept->GetVoxels().GetCountOfVoxels()<=1) {
            G4ExceptionDescription message;
            message << "Gate-3 H2 outer boundary type/geometry contract failed: "
                    << "name=" << name << " round=" << round
                    << " sides=" << sides
                    << " derived_type=" << (corrected!=nullptr)
                    << " facets=" << swept->GetNumberOfFacets()
                    << " mesh_facets=" << meshFacetCount
                    << " reference_facets=" << referenceFacetCount
                    << " active_refinement=" << activeRefinement
                    << " entity_type=" << swept->GetEntityType()
                    << " actual_voxels="
                    << swept->GetVoxels().GetCountOfVoxels();
            G4Exception("DetectorConstruction::BuildFanoutPreview",
                        "HODO_GATE3_H2_OUTER_BOUNDARY_CORRECTION_CONTRACT",
                        FatalException,message);
          }
          if(h2OuterBoundaryCloneAuditCount==0) {
            std::unique_ptr<G4VSolid> clone(swept->Clone());
            auto* correctedClone=dynamic_cast<
                H2OuterBoundaryTessellatedSolid*>(clone.get());
            auto* tessellatedClone=dynamic_cast<
                G4TessellatedSolid*>(clone.get());
            if(!correctedClone||!tessellatedClone||
               clone->GetEntityType()!="G4TessellatedSolid"||
               tessellatedClone->GetNumberOfFacets()!=400||
               !tessellatedClone->GetSolidClosed()||
               tessellatedClone->GetVoxels().GetCountOfVoxels()!=
                   swept->GetVoxels().GetCountOfVoxels()) {
              G4ExceptionDescription message;
              message << "Gate-3 H2 outer boundary Clone contract failed: "
                      << "derived=" << (correctedClone!=nullptr)
                      << " tessellated=" << (tessellatedClone!=nullptr)
                      << " entity_type=" << clone->GetEntityType()
                      << " facets="
                      << (tessellatedClone
                              ? tessellatedClone->GetNumberOfFacets() : -1)
                      << " closed="
                      << (tessellatedClone
                              ? tessellatedClone->GetSolidClosed() : false)
                      << " voxels="
                      << (tessellatedClone
                              ? tessellatedClone->GetVoxels()
                                    .GetCountOfVoxels() : -1)
                      << " source_voxels="
                      << swept->GetVoxels().GetCountOfVoxels();
              G4Exception("DetectorConstruction::BuildFanoutPreview",
                          "HODO_GATE3_H2_OUTER_BOUNDARY_CORRECTION_CLONE",
                          FatalException,message);
            }
            ++h2OuterBoundaryCloneAuditCount;
          }
          ++h2OuterBoundaryAuditCount;
        }
        const auto defects=swept->CheckStructure();
        const G4double areaDelta=
            std::abs(meshArea-referenceMeshArea);
        const G4double signedVolumeDelta=
            std::abs(signedMeshVolume-signedReferenceMeshVolume);
        const G4double refinementComplexity=
            std::pow(static_cast<G4double>(
                round?fH1ActiveFacetRefinement:1),2);
        const G4double metricFactor=
            1024.*std::numeric_limits<G4double>::epsilon()*
            refinementComplexity;
        const G4double areaTolerance=std::max(
            1.e-12*mm*mm,
            metricFactor*std::max(referenceMeshArea,1.*mm*mm));
        const G4double volumeTolerance=std::max(
            1.e-12*mm*mm*mm,
            metricFactor*std::max(
                std::abs(signedReferenceMeshVolume),1.*mm*mm*mm));
        const G4bool metricInvariant=
            areaDelta<=areaTolerance&&
            signedVolumeDelta<=volumeTolerance;
        G4cout << "[v21 mesh audit] name=" << name
               << " sides=" << sides << " defects=" << defects
               << " longitudinal_intervals=" << rings.size()-1
               << " active_facet_refinement="
               << activeRefinement
               << " refinement_method=conforming_uniform_barycentric"
               << " facets=" << meshFacetCount
               << " reference_facets=" << referenceFacetCount
               << " expected_facets=" << expectedFacetCount
               << " canonical_edges=" << canonicalEdgePointCache.size()
               << " canonical_edge_interior_points="
               << canonicalEdgePointCache.size()*
                  static_cast<std::size_t>(
                      std::max(0,activeRefinement-1))
               << " uniform_parents=" << uniformlyRefinedParentCount
               << " uniform_children=" << uniformChildFacetCount
               << " boundary_fan_parents=" << boundaryFanParentCount
               << " boundary_fan_children=" << boundaryFanChildFacetCount
               << " interior_lattice_point_evaluations="
               << interiorLatticePointCount
               << " max_child_plane_residual_nm="
               << maxChildPlaneResidual/nm
               << " child_plane_tolerance_nm="
               << childPlaneTolerance/nm
               << " child_plane_invariant=" << childPlaneInvariant
               << " topology_count_invariant=" << topologyCountInvariant
               << " surface_area_mm2=" << meshArea/(mm*mm)
               << " volume_mm3=" << std::abs(signedMeshVolume)/(mm*mm*mm)
               << " reference_surface_area_mm2="
               << referenceMeshArea/(mm*mm)
               << " reference_volume_mm3="
               << std::abs(signedReferenceMeshVolume)/(mm*mm*mm)
               << " area_delta_mm2=" << areaDelta/(mm*mm)
               << " signed_volume_delta_mm3="
               << signedVolumeDelta/(mm*mm*mm)
               << " area_tolerance_mm2=" << areaTolerance/(mm*mm)
               << " volume_tolerance_mm3="
               << volumeTolerance/(mm*mm*mm)
               << " metric_invariant=" << metricInvariant
               << G4endl;
        if(!topologyCountInvariant) {
          G4ExceptionDescription message;
          message << "Conforming tessellated refinement count mismatch: "
                  << "name=" << name
                  << " refinement=" << activeRefinement
                  << " facets=" << meshFacetCount
                  << " expected_facets=" << expectedFacetCount
                  << " reference_facets=" << referenceFacetCount
                  << " expected_reference_facets=" << baseFacetCount
                  << " canonical_edges="
                  << canonicalEdgePointCache.size()
                  << " expected_canonical_edges="
                  << expectedCanonicalEdges
                  << " uniform_parents=" << uniformlyRefinedParentCount
                  << " uniform_children=" << uniformChildFacetCount
                  << " boundary_fan_parents=" << boundaryFanParentCount
                  << " boundary_fan_children=" << boundaryFanChildFacetCount;
          G4Exception("DetectorConstruction::BuildFanoutPreview",
                      "HODO_V21_MESH_TOPOLOGY_COUNT",
                      FatalException,message);
        }
        if(!childPlaneInvariant) {
          G4ExceptionDescription message;
          message << "Conforming tessellated refinement left an original "
                  << "triangle plane beyond the numerical bound: "
                  << "name=" << name
                  << " refinement=" << activeRefinement
                  << " max_residual_nm=" << maxChildPlaneResidual/nm
                  << " tolerance_nm=" << childPlaneTolerance/nm;
          G4Exception("DetectorConstruction::BuildFanoutPreview",
                      "HODO_V21_MESH_CHILD_PLANE",
                      FatalException,message);
        }
        if(defects!=0) {
          G4ExceptionDescription message;
          message << "Tessellated solid " << name
                  << " failed CheckStructure with " << defects
                  << " defects after conforming refinement.";
          G4Exception("DetectorConstruction::BuildFanoutPreview",
                      "HODO_V21_MESH_STRUCTURE",FatalException,message);
        }
        if(!metricInvariant) {
          G4ExceptionDescription message;
          message << "Tessellated surface-preserving refinement changed "
                  << "area or signed volume beyond the floating-point bound: "
                  << "name=" << name
                  << " refinement=" << activeRefinement
                  << " area_delta_mm2=" << areaDelta/(mm*mm)
                  << " area_tolerance_mm2=" << areaTolerance/(mm*mm)
                  << " volume_delta_mm3="
                  << signedVolumeDelta/(mm*mm*mm)
                  << " volume_tolerance_mm3="
                  << volumeTolerance/(mm*mm*mm);
          G4Exception("DetectorConstruction::BuildFanoutPreview",
                      "HODO_V21_MESH_METRIC_INVARIANT",
                      FatalException,message);
        }
        return swept;
      };

      // Native-solid OAT candidates for the H1 round fibre.  Both candidates
      // retain the exact sampled centreline, layer radii, endpoint planes,
      // placement, materials, and optical coefficients.  They differ only in
      // how the piecewise centreline is represented between samples.
      const auto segmentRotation = [&](const G4ThreeVector& tangent) {
        G4ThreeVector localX=verticalAxis.cross(tangent);
        if(localX.mag2()<1.e-18)
          localX=axis==0 ? G4ThreeVector(0.,1.,0.)
                         : G4ThreeVector(-1.,0.,0.);
        localX=localX.unit();
        const G4ThreeVector localY=tangent.cross(localX).unit();
        return G4RotationMatrix(localX,localY,tangent);
      };
      const auto auditNativeEndpoints = [&](
          const G4String& name, G4MultiUnion* solid,
          const std::vector<G4ThreeVector>& tangents,
          G4double sampleRadius) {
        const G4double surfaceTolerance =
            G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
        const G4double probe=std::max(
            100.*surfaceTolerance,10.*fH1RoundMiterJointOverlap);
        const auto startRadial=
            sampleRadius*segmentRotation(tangents.front()).colX();
        const auto endRadial=
            sampleRadius*segmentRotation(tangents.back()).colX();
        const auto startOutside=
            centres.front()+startRadial-probe*tangents.front();
        const auto startInside=
            centres.front()+startRadial+probe*tangents.front();
        const auto endInside=
            centres.back()+endRadial-probe*tangents.back();
        const auto endOutside=
            centres.back()+endRadial+probe*tangents.back();
        if(solid->Inside(startOutside)!=kOutside ||
           solid->Inside(startInside)==kOutside ||
           solid->Inside(endInside)==kOutside ||
           solid->Inside(endOutside)!=kOutside) {
          G4ExceptionDescription message;
          message << "Native round sweep failed signed flat-end probes: "
                  << name << " probe_nm=" << probe/nm
                  << " states=(" << solid->Inside(startOutside) << ","
                  << solid->Inside(startInside) << ","
                  << solid->Inside(endInside) << ","
                  << solid->Inside(endOutside) << ")";
          G4Exception("DetectorConstruction::BuildFanoutPreview",
                      "HODO_V21_NATIVE_ENDCAP",FatalException,message);
        }
      };
      const auto makeNativeMiter = [&](const G4String& name,
                                       G4double innerRadius,
                                       G4double outerRadius,
                                       G4bool displacedWorldNodes) {
        std::vector<G4ThreeVector> tangents;
        tangents.reserve(centres.size()-1);
        for(std::size_t segment=0;segment+1<centres.size();++segment) {
          const auto delta=centres[segment+1]-centres[segment];
          if(!std::isfinite(delta.mag2())||delta.mag2()<=0.) {
            G4ExceptionDescription message;
            message << "Invalid native-miter centreline segment: name="
                    << name << " segment=" << segment
                    << " delta=" << delta;
            G4Exception("DetectorConstruction::BuildFanoutPreview",
                        "HODO_V21_NATIVE_SEGMENT",FatalException,message);
          }
          tangents.push_back(delta.unit());
        }
        std::vector<G4ThreeVector> jointNormals(centres.size());
        jointNormals.front()=tangents.front();
        jointNormals.back()=tangents.back();
        for(std::size_t joint=1;joint+1<centres.size();++joint) {
          const auto sum=tangents[joint-1]+tangents[joint];
          if(sum.mag2()<1.e-12) {
            G4ExceptionDescription message;
            message << "Native-miter consecutive segments are antiparallel: "
                    << "name=" << name << " joint=" << joint;
            G4Exception("DetectorConstruction::BuildFanoutPreview",
                        "HODO_V21_NATIVE_MITER",FatalException,message);
          }
          jointNormals[joint]=sum.unit();
        }
        // G4MultiUnion 11.4.1 combines touching-node normals in each node's
        // local coordinates.  Recursively choose the otherwise-free twist of
        // every axisymmetric segment so a shared miter normal has the same
        // local tuple in both neighbours.  This is a numerical invariant, not
        // a physical rotation of the circular fibre.
        std::vector<G4RotationMatrix> rotations;
        rotations.reserve(tangents.size());
        rotations.push_back(segmentRotation(tangents.front()));
        for(std::size_t segment=1;segment<tangents.size();++segment) {
          const auto& zAxis=tangents[segment];
          const auto previousLocal=
              rotations.back().inverse()*jointNormals[segment];
          const G4double previousRho=std::hypot(
              previousLocal.x(),previousLocal.y());
          G4ThreeVector xAxis;
          G4ThreeVector yAxis;
          if(previousRho>1.e-8) {
            const G4double currentZ=
                jointNormals[segment].dot(zAxis);
            const auto p=(jointNormals[segment]-
                          currentZ*zAxis).unit();
            const auto q=zAxis.cross(p).unit();
            const G4double alpha=previousLocal.x()/previousRho;
            const G4double beta=previousLocal.y()/previousRho;
            xAxis=(alpha*p-beta*q).unit();
            yAxis=(beta*p+alpha*q).unit();
          } else {
            xAxis=rotations.back().colX()-
                rotations.back().colX().dot(zAxis)*zAxis;
            if(xAxis.mag2()<1.e-18) {
              const auto fallback=segmentRotation(zAxis);
              xAxis=fallback.colX();
            } else {
              xAxis=xAxis.unit();
            }
            yAxis=zAxis.cross(xAxis).unit();
          }
          rotations.emplace_back(xAxis,yAxis,zAxis);
          const auto currentLocal=
              rotations.back().inverse()*jointNormals[segment];
          const G4double tupleDefect=(currentLocal-previousLocal).mag();
          if(tupleDefect>1.e-8) {
            G4ExceptionDescription message;
            message << "Native-miter shared-normal local tuple mismatch: "
                    << "name=" << name << " joint=" << segment
                    << " defect=" << tupleDefect
                    << " previous=" << previousLocal
                    << " current=" << currentLocal;
            G4Exception("DetectorConstruction::BuildFanoutPreview",
                        "HODO_V21_NATIVE_MITER_TWIST",FatalException,message);
          }
        }
        auto* multi=new G4MultiUnion(name);
        for(std::size_t segment=0;segment<tangents.size();++segment) {
          const auto& tangent=tangents[segment];
          const G4double lowExtension=
              segment==0 ? 0. : fH1RoundMiterJointOverlap;
          const G4double highExtension=
              segment+1==tangents.size()
                  ? 0. : fH1RoundMiterJointOverlap;
          const auto lowAnchor=centres[segment]-lowExtension*tangent;
          const auto highAnchor=
              centres[segment+1]+highExtension*tangent;
          const auto midpoint=.5*(lowAnchor+highAnchor);
          const G4double halfLength=.5*(highAnchor-lowAnchor).mag();
          const auto& rotation=rotations[segment];
          const auto inverse=rotation.inverse();
          const auto lowNormalLocal=
              inverse*(segment==0 ? -tangent : -jointNormals[segment]);
          const auto highNormalLocal=inverse*(
              segment+1==tangents.size()
                  ? tangent : jointNormals[segment+1]);
          if(lowNormalLocal.z()>=0. || highNormalLocal.z()<=0.) {
            G4ExceptionDescription message;
            message << "Invalid local G4CutTubs normals: name=" << name
                    << " segment=" << segment
                    << " low=" << lowNormalLocal
                    << " high=" << highNormalLocal;
            G4Exception("DetectorConstruction::BuildFanoutPreview",
                        "HODO_V21_NATIVE_MITER_NORMAL",FatalException,message);
          }
          G4VSolid* cut=nullptr;
          const G4double lowTransverse=std::hypot(
              lowNormalLocal.x(),lowNormalLocal.y());
          const G4double highTransverse=std::hypot(
              highNormalLocal.x(),highNormalLocal.y());
          if(lowTransverse<1.e-14&&highTransverse<1.e-14) {
            // Avoid G4CutTubs' own GeomSolids1001 "use G4Tubs" warning for an
            // exactly straight run.  The two solids are identical here.
            cut=new G4Tubs(
                name+"_Segment_"+std::to_string(segment),
                innerRadius,outerRadius,halfLength,0.,2.*CLHEP::pi);
          } else {
            cut=new G4CutTubs(
                name+"_Segment_"+std::to_string(segment),
                innerRadius,outerRadius,halfLength,0.,2.*CLHEP::pi,
                lowNormalLocal,highNormalLocal);
          }
          const G4Transform3D transform(rotation,midpoint);
          if(displacedWorldNodes) {
            // Geant4 11.4.1 G4MultiUnion compares touching-node normals
            // without rotating them out of each node frame.  Wrapping the
            // analytic primitive makes SurfaceNormal return world-frame
            // normals; identity AddNode then avoids that implementation
            // defect.  This is an explicit OAT because G4DisplacedSolid is
            // documented as an internal Boolean helper.
            auto* displaced=new G4DisplacedSolid(
                name+"_WorldSegment_"+std::to_string(segment),
                cut,transform);
            multi->AddNode(
                *displaced,G4Transform3D(G4RotationMatrix(),{}));
          } else {
            multi->AddNode(*cut,transform);
          }
        }
        multi->Voxelize();
        multi->SetAccurateSafety(true);
        if(multi->GetNumberOfSolids()!=
           static_cast<G4int>(tangents.size())) {
          G4Exception("DetectorConstruction::BuildFanoutPreview",
                      "HODO_V21_NATIVE_NODE_COUNT",FatalException,
                      "Native-miter node-count audit failed.");
        }
        const G4double endpointSampleRadius=
            innerRadius==0. ? .5*outerRadius
                            : .5*(innerRadius+outerRadius);
        std::vector<G4double> seamRadii;
        if(innerRadius==0.) {
          seamRadii={0.,.5*outerRadius,.9*outerRadius};
        } else {
          const auto band=outerRadius-innerRadius;
          seamRadii={
              innerRadius+.1*band,
              .5*(innerRadius+outerRadius),
              outerRadius-.1*band};
        }
        G4int seamProbeCount=0;
        for(std::size_t joint=1;joint+1<centres.size();++joint) {
          const auto e1=jointNormals[joint].orthogonal().unit();
          const auto e2=jointNormals[joint].cross(e1).unit();
          for(const auto sampleRadius:seamRadii) {
           const G4int phaseClasses=sampleRadius==0. ? 1 : 2;
           const G4int azimuths=sampleRadius==0. ? 1 : 64;
           for(G4int phaseClass=0;
               phaseClass<phaseClasses;++phaseClass) {
            for(G4int azimuth=0;azimuth<azimuths;++azimuth) {
              const G4double phi=2.*CLHEP::pi*
                  (azimuth+.5*phaseClass)/64.;
              const auto inPlane=
                  std::cos(phi)*e1+std::sin(phi)*e2;
              const G4double projectedRadius=std::sqrt(std::max(
                  0.,1.-std::pow(
                      inPlane.dot(tangents[joint-1]),2)));
              if(projectedRadius<1.e-12) {
                G4Exception(
                    "DetectorConstruction::BuildFanoutPreview",
                    "HODO_V21_NATIVE_CONNECTIVITY",FatalException,
                    "Degenerate native-miter joint probe direction.");
              }
              const auto jointProbe=centres[joint]+
                  sampleRadius/projectedRadius*inPlane;
              G4bool seamInside=multi->Inside(jointProbe)==kInside;
              for(const auto delta:{10.*nm,100.*nm}) {
                seamInside=seamInside &&
                    multi->Inside(
                        jointProbe-delta*jointNormals[joint])==kInside &&
                    multi->Inside(
                        jointProbe+delta*jointNormals[joint])==kInside;
              }
              ++seamProbeCount;
              if(!seamInside) {
                G4ExceptionDescription message;
                message << "Native-miter union has a disconnected or "
                        << "artificial internal seam: name=" << name
                        << " joint=" << joint
                        << " radius_mm=" << sampleRadius/mm
                        << " phase_class=" << phaseClass
                        << " azimuth=" << azimuth
                        << " overlap_nm="
                        << fH1RoundMiterJointOverlap/nm
                        << " states=(joint:"
                        << multi->Inside(jointProbe)
                        << ",minus10:"
                        << multi->Inside(
                            jointProbe-10.*nm*jointNormals[joint])
                        << ",plus10:"
                        << multi->Inside(
                            jointProbe+10.*nm*jointNormals[joint])
                        << ",minus100:"
                        << multi->Inside(
                            jointProbe-100.*nm*jointNormals[joint])
                        << ",plus100:"
                        << multi->Inside(
                            jointProbe+100.*nm*jointNormals[joint])
                        << ")";
                G4Exception(
                    "DetectorConstruction::BuildFanoutPreview",
                    "HODO_V21_NATIVE_CONNECTIVITY",
                    FatalException,message);
              }
            }
           }
          }
        }
        auditNativeEndpoints(
            name,multi,tangents,endpointSampleRadius);
        G4double centrelineLength=0.;
        for(std::size_t segment=0;segment+1<centres.size();++segment)
          centrelineLength+=(centres[segment+1]-centres[segment]).mag();
        const G4double crossArea=CLHEP::pi*
            (outerRadius*outerRadius-innerRadius*innerRadius);
        const G4double zeroOverlapVolume=crossArea*centrelineLength;
        const G4double zeroOverlapArea=
            2.*CLHEP::pi*(outerRadius+innerRadius)*centrelineLength+
            2.*crossArea;
        const G4double overlapNodeSumUpper=crossArea*(
            2.*fH1RoundMiterJointOverlap*(tangents.size()-1));
        G4cout << "[v21 native sweep audit] name=" << name
               << " representation=native_miter_union"
               << " cut_tubs_nodes=" << tangents.size()
               << " joint_overlap_nm=" << fH1RoundMiterJointOverlap/nm
               << " inner_radius_mm=" << innerRadius/mm
               << " outer_radius_mm=" << outerRadius/mm
               << " node_frame="
               << (displacedWorldNodes ? "displaced_world" : "public_addnode")
               << " centreline_length_mm=" << centrelineLength/mm
               << " zero_overlap_volume_mm3="
               << zeroOverlapVolume/(mm*mm*mm)
               << " zero_overlap_surface_area_mm2="
               << zeroOverlapArea/(mm*mm)
               << " overlap_node_sum_volume_upper_mm3="
               << overlapNodeSumUpper/(mm*mm*mm)
               << " seam_probe_count=" << seamProbeCount
               << " seam_expected_state=kInside"
               << " seam_offsets_nm=10,100"
               << " endpoint_caps=flat" << G4endl;
        return multi;
      };
      struct BalancedMiterResult {
        G4VSolid* solid=nullptr;
        G4Transform3D localToWorld;
        G4int leafCount=0;
        G4int depth=0;
      };
      const auto makeNativeBalancedMiter = [&](
          const G4String& name,G4double innerRadius,
          G4double outerRadius) {
        std::vector<G4ThreeVector> tangents;
        tangents.reserve(centres.size()-1);
        for(std::size_t segment=0;segment+1<centres.size();++segment) {
          const auto delta=centres[segment+1]-centres[segment];
          if(!std::isfinite(delta.mag2())||delta.mag2()<=0.) {
            G4ExceptionDescription message;
            message << "Invalid balanced-miter centreline segment: name="
                    << name << " segment=" << segment
                    << " delta=" << delta;
            G4Exception("DetectorConstruction::BuildFanoutPreview",
                        "HODO_V21_BALANCED_SEGMENT",
                        FatalException,message);
          }
          tangents.push_back(delta.unit());
        }
        std::vector<G4ThreeVector> jointNormals(centres.size());
        jointNormals.front()=tangents.front();
        jointNormals.back()=tangents.back();
        for(std::size_t joint=1;joint+1<centres.size();++joint) {
          const auto sum=tangents[joint-1]+tangents[joint];
          if(sum.mag2()<1.e-12) {
            G4ExceptionDescription message;
            message << "Balanced-miter consecutive segments are "
                    << "antiparallel: name=" << name
                    << " joint=" << joint;
            G4Exception("DetectorConstruction::BuildFanoutPreview",
                        "HODO_V21_BALANCED_MITER",
                        FatalException,message);
          }
          jointNormals[joint]=sum.unit();
        }
        struct TreeNode {
          G4VSolid* solid=nullptr;
          G4Transform3D localToWorld;
          G4int minimumLeafDepth=0;
          G4int maximumLeafDepth=0;
        };
        std::vector<TreeNode> leaves;
        leaves.reserve(tangents.size());
        for(std::size_t segment=0;segment<tangents.size();++segment) {
          const auto rotation=segmentRotation(tangents[segment]);
          const auto inverse=rotation.inverse();
          const auto lowNormalLocal=inverse*(
              segment==0 ? -tangents[segment]
                         : -jointNormals[segment]);
          const auto highNormalLocal=inverse*(
              segment+1==tangents.size()
                  ? tangents[segment] : jointNormals[segment+1]);
          if(lowNormalLocal.z()>=0.||highNormalLocal.z()<=0.) {
            G4ExceptionDescription message;
            message << "Invalid balanced G4CutTubs normals: name="
                    << name << " segment=" << segment
                    << " low=" << lowNormalLocal
                    << " high=" << highNormalLocal;
            G4Exception("DetectorConstruction::BuildFanoutPreview",
                        "HODO_V21_BALANCED_NORMAL",
                        FatalException,message);
          }
          const auto midpoint=.5*(
              centres[segment]+centres[segment+1]);
          const G4double halfLength=.5*(
              centres[segment+1]-centres[segment]).mag();
          G4VSolid* leaf=nullptr;
          const G4double lowTransverse=std::hypot(
              lowNormalLocal.x(),lowNormalLocal.y());
          const G4double highTransverse=std::hypot(
              highNormalLocal.x(),highNormalLocal.y());
          const auto leafName=
              name+"_BalancedLeaf_"+std::to_string(segment);
          if(lowTransverse<1.e-14&&highTransverse<1.e-14) {
            leaf=new G4Tubs(
                leafName,innerRadius,outerRadius,halfLength,
                0.,2.*CLHEP::pi);
          } else {
            leaf=new G4CutTubs(
                leafName,innerRadius,outerRadius,halfLength,
                0.,2.*CLHEP::pi,lowNormalLocal,highNormalLocal);
          }
          leaves.push_back(
              {leaf,G4Transform3D(rotation,midpoint),0,0});
        }
        G4int nodeOrdinal=0;
        const auto buildBalancedRange=[&](
            auto&& self,std::size_t first,
            std::size_t last) -> TreeNode {
          if(last-first==1)return leaves[first];
          const auto middle=first+(last-first)/2;
          const auto left=self(self,first,middle);
          const auto right=self(self,middle,last);
          const auto rightToLeft=
              left.localToWorld.inverse()*right.localToWorld;
          const G4bool finalNode=
              first==0&&last==leaves.size();
          const G4String nodeName=finalNode
              ? name
              : G4String(name+"_BalancedNode_"+
                         std::to_string(nodeOrdinal++));
          return {
              new G4UnionSolid(
                  nodeName,left.solid,right.solid,rightToLeft),
              left.localToWorld,
              std::min(
                  left.minimumLeafDepth,
                  right.minimumLeafDepth)+1,
              std::max(
                  left.maximumLeafDepth,
                  right.maximumLeafDepth)+1};
        };
        const auto root=buildBalancedRange(
            buildBalancedRange,0,leaves.size());
        if(!root.solid||
           root.minimumLeafDepth!=5||
           root.maximumLeafDepth!=6||
           tangents.size()!=49) {
          G4ExceptionDescription message;
          message << "Balanced Boolean tree invariant failed: name="
                  << name << " leaves=" << tangents.size()
                  << " minimum_leaf_depth="
                  << root.minimumLeafDepth
                  << " maximum_leaf_depth="
                  << root.maximumLeafDepth;
          G4Exception("DetectorConstruction::BuildFanoutPreview",
                      "HODO_V21_BALANCED_TREE",
                      FatalException,message);
        }
        const auto toLocal=[&](const G4ThreeVector& worldPoint) {
          return root.localToWorld.getRotation().inverse()*
              (worldPoint-root.localToWorld.getTranslation());
        };
        const auto inside=[&](const G4ThreeVector& worldPoint) {
          return root.solid->Inside(toLocal(worldPoint));
        };

        std::vector<G4double> seamRadii;
        if(innerRadius==0.) {
          seamRadii={0.,.5*outerRadius,.9*outerRadius};
        } else {
          const auto band=outerRadius-innerRadius;
          seamRadii={
              innerRadius+.1*band,
              .5*(innerRadius+outerRadius),
              outerRadius-.1*band};
        }
        G4int seamProbeCount=0;
        for(std::size_t joint=1;joint+1<centres.size();++joint) {
          const auto e1=jointNormals[joint].orthogonal().unit();
          const auto e2=jointNormals[joint].cross(e1).unit();
          for(const auto sampleRadius:seamRadii) {
            const G4int phaseClasses=sampleRadius==0. ? 1 : 2;
            const G4int azimuths=sampleRadius==0. ? 1 : 64;
            for(G4int phaseClass=0;
                phaseClass<phaseClasses;++phaseClass) {
              for(G4int azimuth=0;azimuth<azimuths;++azimuth) {
                const G4double phi=2.*CLHEP::pi*
                    (azimuth+.5*phaseClass)/64.;
                const auto inPlane=
                    std::cos(phi)*e1+std::sin(phi)*e2;
                const G4double projected=std::sqrt(std::max(
                    0.,1.-std::pow(
                        inPlane.dot(tangents[joint-1]),2)));
                if(projected<1.e-12) {
                  G4Exception(
                      "DetectorConstruction::BuildFanoutPreview",
                      "HODO_V21_BALANCED_SEAM",FatalException,
                      "Degenerate balanced seam probe.");
                }
                const auto jointProbe=centres[joint]+
                    sampleRadius/projected*inPlane;
                G4bool seamInside=inside(jointProbe)==kInside;
                for(const auto delta:{10.*nm,100.*nm}) {
                  seamInside=seamInside&&
                      inside(jointProbe-
                             delta*jointNormals[joint])==kInside&&
                      inside(jointProbe+
                             delta*jointNormals[joint])==kInside;
                }
                ++seamProbeCount;
                if(!seamInside) {
                  G4ExceptionDescription message;
                  message << "Balanced union artificial/internal seam: "
                          << "name=" << name
                          << " joint=" << joint
                          << " radius_mm=" << sampleRadius/mm
                          << " phase_class=" << phaseClass
                          << " azimuth=" << azimuth
                          << " states=(joint:"
                          << inside(jointProbe)
                          << ",minus10:"
                          << inside(jointProbe-
                              10.*nm*jointNormals[joint])
                          << ",plus10:"
                          << inside(jointProbe+
                              10.*nm*jointNormals[joint])
                          << ",minus100:"
                          << inside(jointProbe-
                              100.*nm*jointNormals[joint])
                          << ",plus100:"
                          << inside(jointProbe+
                              100.*nm*jointNormals[joint])
                          << ")";
                  G4Exception(
                      "DetectorConstruction::BuildFanoutPreview",
                      "HODO_V21_BALANCED_SEAM",
                      FatalException,message);
                }
              }
            }
          }
        }

        const G4double endpointRadius=
            innerRadius==0. ? .5*outerRadius
                            : .5*(innerRadius+outerRadius);
        const G4double endpointProbe=std::max(
            100.*G4GeometryTolerance::GetInstance()->
                GetSurfaceTolerance(),100.*nm);
        const auto startRadial=
            endpointRadius*segmentRotation(
                tangents.front()).colX();
        const auto endRadial=
            endpointRadius*segmentRotation(
                tangents.back()).colX();
        if(inside(centres.front()+startRadial-
                  endpointProbe*tangents.front())!=kOutside||
           inside(centres.front()+startRadial+
                  endpointProbe*tangents.front())!=kInside||
           inside(centres.back()+endRadial-
                  endpointProbe*tangents.back())!=kInside||
           inside(centres.back()+endRadial+
                  endpointProbe*tangents.back())!=kOutside) {
          G4ExceptionDescription message;
          message << "Balanced union signed endpoint audit failed: name="
                  << name << " probe_nm=" << endpointProbe/nm;
          G4Exception("DetectorConstruction::BuildFanoutPreview",
                      "HODO_V21_BALANCED_ENDCAP",
                      FatalException,message);
        }

        G4double centrelineLength=0.;
        for(std::size_t segment=0;segment+1<centres.size();++segment)
          centrelineLength+=(
              centres[segment+1]-centres[segment]).mag();
        const G4double crossArea=CLHEP::pi*
            (outerRadius*outerRadius-innerRadius*innerRadius);
        G4cout << "[v21 native sweep audit] name=" << name
               << " representation=native_balanced_union"
               << " public_api=G4UnionSolid"
               << " leaves=" << tangents.size()
               << " minimum_leaf_depth="
               << root.minimumLeafDepth
               << " maximum_leaf_depth="
               << root.maximumLeafDepth
               << " inner_radius_mm=" << innerRadius/mm
               << " outer_radius_mm=" << outerRadius/mm
               << " centreline_length_mm=" << centrelineLength/mm
               << " zero_overlap_volume_mm3="
               << crossArea*centrelineLength/(mm*mm*mm)
               << " seam_probe_count=" << seamProbeCount
               << " seam_expected_state=kInside"
               << " seam_offsets_nm=10,100"
               << " endpoint_caps=flat" << G4endl;
        return BalancedMiterResult{
            root.solid,root.localToWorld,
            static_cast<G4int>(tangents.size()),
            root.maximumLeafDepth};
      };
      const auto makeNativeCapsule = [&](const G4String& name,
                                         G4double radius) {
        std::vector<G4ThreeVector> tangents;
        tangents.reserve(centres.size()-1);
        auto* multi=new G4MultiUnion(name);
        for(std::size_t segment=0;segment+1<centres.size();++segment) {
          const auto delta=centres[segment+1]-centres[segment];
          if(!std::isfinite(delta.mag2())||delta.mag2()<=0.) {
            G4ExceptionDescription message;
            message << "Invalid native-capsule centreline segment: name="
                    << name << " segment=" << segment
                    << " delta=" << delta;
            G4Exception("DetectorConstruction::BuildFanoutPreview",
                        "HODO_V21_NATIVE_SEGMENT",FatalException,message);
          }
          const auto tangent=delta.unit();
          tangents.push_back(tangent);
          auto* tube=new G4Tubs(
              name+"_Segment_"+std::to_string(segment),
              0.,radius,.5*delta.mag(),0.,2.*CLHEP::pi);
          multi->AddNode(
              *tube,G4Transform3D(
                  segmentRotation(tangent),
                  .5*(centres[segment]+centres[segment+1])));
        }
        for(std::size_t joint=1;joint+1<centres.size();++joint) {
          auto* orb=new G4Orb(
              name+"_Joint_"+std::to_string(joint),radius);
          multi->AddNode(
              *orb,G4Transform3D(G4RotationMatrix(),centres[joint]));
        }
        multi->Voxelize();
        multi->SetAccurateSafety(true);
        const auto expectedNodes=
            static_cast<G4int>(2*tangents.size()-1);
        if(multi->GetNumberOfSolids()!=expectedNodes) {
          G4Exception("DetectorConstruction::BuildFanoutPreview",
                      "HODO_V21_NATIVE_NODE_COUNT",FatalException,
                      "Native-capsule node-count audit failed.");
        }
        for(std::size_t joint=1;joint+1<centres.size();++joint)
          if(multi->Inside(centres[joint])==kOutside) {
            G4ExceptionDescription message;
            message << "Native-capsule union is disconnected at joint: name="
                    << name << " joint=" << joint;
            G4Exception("DetectorConstruction::BuildFanoutPreview",
                        "HODO_V21_NATIVE_CONNECTIVITY",FatalException,message);
        }
        auditNativeEndpoints(name,multi,tangents,0.);
        G4double centrelineLength=0.;
        for(std::size_t segment=0;segment+1<centres.size();++segment)
          centrelineLength+=(centres[segment+1]-centres[segment]).mag();
        const G4double tubeReferenceVolume=
            CLHEP::pi*radius*radius*centrelineLength;
        const G4double tubeReferenceArea=
            2.*CLHEP::pi*radius*centrelineLength+
            2.*CLHEP::pi*radius*radius;
        G4cout << "[v21 native sweep audit] name=" << name
               << " representation=native_capsule_union"
               << " cylinder_nodes=" << tangents.size()
               << " sphere_nodes=" << tangents.size()-1
               << " radius_mm=" << radius/mm
               << " centreline_length_mm=" << centrelineLength/mm
               << " tube_reference_volume_mm3="
               << tubeReferenceVolume/(mm*mm*mm)
               << " tube_reference_surface_area_mm2="
               << tubeReferenceArea/(mm*mm)
               << " endpoint_caps=flat" << G4endl;
        return multi;
      };
      const auto makeLayer = [&](const G4String& name, G4double halfWidth,
                                 G4double layerShift,
                                 G4bool h2OuterBoundaryFix=false)
          -> G4VSolid* {
        if(!round ||
           fH1RoundSweepRepresentation=="tessellated")
          return makeSwept(
              name,halfWidth,layerShift,h2OuterBoundaryFix);
        if(fH1RoundSweepRepresentation=="native_miter_subtraction_union")
          return makeNativeMiter(name,0.,halfWidth,false);
        if(fH1RoundSweepRepresentation=="native_capsule_union")
          return makeNativeCapsule(name,halfWidth);
        G4Exception("DetectorConstruction::BuildFanoutPreview",
                    "HODO_V21_ROUND_REPRESENTATION",FatalException,
                    "Unreachable round sweep representation.");
        return nullptr;
      };

      constexpr G4double nominalOuterSize=.4995*mm;
      const G4bool singleClad=
          roundFibre ? GetH1SingleClad() : GetH2SingleClad();
      const G4double effectiveFibreOverlapClearance =
          roundFibre ? std::min(fFibreOverlapClearance, 0.01498*mm)
                     : fFibreOverlapClearance;
      const G4double outerSize =
          nominalOuterSize-effectiveFibreOverlapClearance;
      const G4double innerSize=roundFibre
          ? (singleClad ? outerSize : .4845*mm) : .4845*mm;
      const G4double coreSize=roundFibre
          ? (singleClad ? .4795*mm : .4695*mm) : .4695*mm;
      if(outerSize<=coreSize ||
         (!singleClad&&outerSize<=innerSize)) {
        G4ExceptionDescription message;
        message << "fibreOverlapClearance removes the cladding shell: "
                << "prefix=" << prefix
                << " channel=" << channel
                << " outer_mm=" << outerSize/mm
                << " inner_mm=" << innerSize/mm
                << " core_mm=" << coreSize/mm
                << " requested_clearance_mm=" << fFibreOverlapClearance/mm
                << " effective_clearance_mm="
                << effectiveFibreOverlapClearance/mm
                << " single_clad=" << singleClad;
        G4Exception("DetectorConstruction::BuildFanoutPreview",
                    "HODO_FIBRE_OVERLAP_CLEARANCE",
                    FatalException,message);
      }
      const auto segmentDistance = [](
          const G4ThreeVector& p1,const G4ThreeVector& q1,
          const G4ThreeVector& p2,const G4ThreeVector& q2) {
        const auto u=q1-p1;
        const auto v=q2-p2;
        const auto w=p1-p2;
        const G4double a=u.dot(u);
        const G4double b=u.dot(v);
        const G4double c=v.dot(v);
        const G4double d=u.dot(w);
        const G4double e=v.dot(w);
        const G4double small=1.e-24*mm*mm;
        G4double sD=a*c-b*b;
        G4double tD=sD;
        G4double sN=0.;
        G4double tN=0.;
        if (a<=small&&c<=small) return (p1-p2).mag();
        if (a<=small) {
          sN=0.; sD=1.; tN=e; tD=c;
        } else if (c<=small) {
          tN=0.; tD=1.; sN=-d; sD=a;
        } else {
          if (sD<small) {
            sN=0.; sD=1.; tN=e; tD=c;
          } else {
            sN=b*e-c*d;
            tN=a*e-b*d;
            if (sN<0.) { sN=0.; tN=e; tD=c; }
            else if (sN>sD) { sN=sD; tN=e+b; tD=c; }
          }
        }
        if (tN<0.) {
          tN=0.;
          if (-d<0.) sN=0.;
          else if (-d>a) sN=sD;
          else { sN=-d; sD=a; }
        } else if (tN>tD) {
          tN=tD;
          if (-d+b<0.) sN=0.;
          else if (-d+b>a) sN=sD;
          else { sN=-d+b; sD=a; }
        }
        const G4double sc=std::abs(sN)<small ? 0. : sN/sD;
        const G4double tc=std::abs(tN)<small ? 0. : tN/tD;
        return (w+sc*u-tc*v).mag();
      };
      const G4bool deterministicSquareFanoutCheck =
          !roundFibre && fFibreOverlapClearance>0.;
      if (deterministicSquareFanoutCheck) {
        G4double minimumGap=DBL_MAX;
        for (const auto& previous : priorSquareCentrelines) {
          for (std::size_t i=0;i+1<centres.size();++i) {
            for (std::size_t j=0;j+1<previous.size();++j) {
              const G4double centreDistance=segmentDistance(
                  centres[i],centres[i+1],previous[j],previous[j+1]);
              minimumGap=std::min(
                  minimumGap,centreDistance-2.*outerSize);
            }
          }
        }
        const G4double requiredGap=std::max(
            G4GeometryTolerance::GetInstance()->GetSurfaceTolerance(),
            1.*nm);
        if (minimumGap<=requiredGap) {
          G4ExceptionDescription message;
          message << "Deterministic square fan-out centreline clearance "
                  << "failed: prefix=" << prefix
                  << " channel=" << channel
                  << " minimum_gap_mm=" << minimumGap/mm
                  << " required_gap_mm=" << requiredGap/mm
                  << " outer_size_mm=" << outerSize/mm
                  << " requested_clearance_mm="
                  << fFibreOverlapClearance/mm
                  << " effective_clearance_mm="
                  << effectiveFibreOverlapClearance/mm;
          G4Exception("DetectorConstruction::BuildFanoutPreview",
                      "HODO_V21_SQUARE_FANOUT_CLEARANCE",
                      FatalException,message);
        }
        if (channel==GetPlaneFibreCount(plane)-1) {
          G4cout << "[v21 square fanout deterministic clearance] prefix="
                 << prefix << " plane=" << plane
                 << " minimum_gap_mm=" << minimumGap/mm
                 << " pSurfChk=replaced_for_square_fanout"
                 << G4endl;
        }
      }
      const G4bool balancedNativeMiter=
          roundFibre&&fH1RoundSweepRepresentation=="native_balanced_union";
      const G4bool directNativeMiter=
          roundFibre&&(fH1RoundSweepRepresentation=="native_miter_union"||
                  fH1RoundSweepRepresentation==
                      "native_miter_public_union"||
                  balancedNativeMiter);
      const G4bool displacedMiterNodes=
          fH1RoundSweepRepresentation=="native_miter_union";
      G4Transform3D roundSweepLocalToWorld;
      G4VSolid* outerSolid=nullptr;
      G4VSolid* coreSolid=nullptr;
      G4VSolid* innerSolid=nullptr;
      if(balancedNativeMiter) {
        const auto suffix=std::to_string(channel);
        const auto outerName=
            prefix+"OuterShellSolid_"+suffix;
        const auto coreName=prefix+"CoreSolid_"+suffix;
        const auto outerResult=makeNativeBalancedMiter(
            outerName,singleClad?coreSize:innerSize,outerSize);
        const auto coreResult=makeNativeBalancedMiter(
            coreName,0.,coreSize);
        const auto transformResidual=[&](
            const G4Transform3D& first,
            const G4Transform3D& second) {
          G4double residual=(
              first.getTranslation()-second.getTranslation()).mag();
          for(const auto axisVector:{
                  G4ThreeVector(1.,0.,0.),
                  G4ThreeVector(0.,1.,0.),
                  G4ThreeVector(0.,0.,1.)})
            residual=std::max(
                residual,
                (first.getRotation()*axisVector-
                 second.getRotation()*axisVector).mag()*mm);
          return residual;
        };
        if(transformResidual(
               outerResult.localToWorld,
               coreResult.localToWorld)>
           4096.*std::numeric_limits<G4double>::epsilon()*mm) {
          G4Exception(
              "DetectorConstruction::BuildFanoutPreview",
              "HODO_V21_BALANCED_FRAME",FatalException,
              "Balanced core and cladding roots do not share one frame.");
        }
        roundSweepLocalToWorld=coreResult.localToWorld;
        outerSolid=outerResult.solid;
        coreSolid=coreResult.solid;
        const auto registerProbeTransform=[&](
            const G4String& solidName) {
          const auto inserted=fNavigatorSolidLocalToWorld.emplace(
              solidName,roundSweepLocalToWorld);
          if(!inserted.second) {
            G4ExceptionDescription message;
            message << "Duplicate balanced solid probe transform: "
                    << solidName;
            G4Exception(
                "DetectorConstruction::BuildFanoutPreview",
                "HODO_V21_BALANCED_NAME",
                FatalException,message);
          }
        };
        registerProbeTransform(outerName);
        registerProbeTransform(coreName);
        if(!singleClad) {
          const auto innerName=
              prefix+"InnerShellSolid_"+suffix;
          const auto innerResult=makeNativeBalancedMiter(
              innerName,coreSize,innerSize);
          if(transformResidual(
                 innerResult.localToWorld,
                 roundSweepLocalToWorld)>
             4096.*std::numeric_limits<G4double>::epsilon()*mm) {
            G4Exception(
                "DetectorConstruction::BuildFanoutPreview",
                "HODO_V21_BALANCED_FRAME",FatalException,
                "Balanced inner-cladding root frame mismatch.");
          }
          innerSolid=innerResult.solid;
          registerProbeTransform(innerName);
        }
      } else if(directNativeMiter) {
        const auto suffix=std::to_string(channel);
        outerSolid=makeNativeMiter(
            prefix+"OuterShellSolid_"+suffix,
            singleClad?coreSize:innerSize,outerSize,
            displacedMiterNodes);
        coreSolid=makeNativeMiter(
            prefix+"CoreSolid_"+suffix,0.,coreSize,
            displacedMiterNodes);
        if(!singleClad)
          innerSolid=makeNativeMiter(
              prefix+"InnerShellSolid_"+suffix,coreSize,innerSize,
              displacedMiterNodes);
      } else {
        outerSolid=makeLayer(
            prefix+"OuterSolid_"+std::to_string(channel),outerSize,0.,
            !roundFibre);
        coreSolid=makeLayer(
            prefix+"CoreSolid_"+std::to_string(channel),
            coreSize,outerSize-coreSize,false);
        if(!singleClad)
          innerSolid=makeLayer(
              prefix+"InnerSolid_"+std::to_string(channel),innerSize,
              .5*(outerSize-innerSize),false);
      }
      if(!roundFibre&&axis==0&&channel==9) {
        auto* corrected=dynamic_cast<
            H2OuterBoundaryTessellatedSolid*>(outerSolid);
        auto* tessellated=dynamic_cast<G4TessellatedSolid*>(outerSolid);
        if(!corrected||!tessellated||
           tessellated->GetNumberOfFacets()!=400||
           tessellated->GetVoxels().GetCountOfVoxels()<=1) {
          G4ExceptionDescription message;
          message << "Gate-3 H2 outer frozen-ray solid contract failed: "
                  << "derived=" << (corrected!=nullptr)
                  << " facets="
                  << (tessellated
                          ? tessellated->GetNumberOfFacets() : -1)
                  << " actual_voxels="
                  << (tessellated
                          ? tessellated->GetVoxels().GetCountOfVoxels()
                          : -1);
          G4Exception("DetectorConstruction::BuildFanoutPreview",
                      "HODO_GATE3_H2_OUTER_BOUNDARY_CORRECTION_FROZEN_RAY",
                      FatalException,message);
        }
        struct FrozenRay {
          G4double longitudinal;
          G4double sideInset;
          G4double normalInset;
          G4ThreeVector localDirection;
        };
        const std::array<FrozenRay,2> frozenRays{{
            {0.2482968209821247*mm,18.423145e-9*mm,
             .771408e-9*mm,
             {1.92759535827335e-05,-.9999999972683715,
              -7.135610921920089e-05}},
            {0.2502501692693591*mm,13.093095e-9*mm,
             .771501e-9*mm,
             {2.21644573541356e-05,-.9999999997460738,
              4.072948288929316e-06}}}};
        const auto& activeFrame=stationFrames.front();
        // Keep frozen probes anchored to the active centre, not the new tail midpoint.
        const auto activeCentre=fFibreFarEndExtension==0.
            ? .5*(centres[0]+centres[1])
            : centres[1]-activeHalfLength*activeFrame.tangent;
        const G4double lookAhead=std::max(
            128.*G4GeometryTolerance::GetInstance()->GetSurfaceTolerance(),
            1.*nm);
        for(std::size_t ray=0;ray<frozenRays.size();++ray) {
          const auto& value=frozenRays[ray];
          const auto worldPoint=activeCentre+
              value.longitudinal*activeFrame.tangent+
              (-outerSize+value.sideInset)*activeFrame.side+
              (-outerSize+value.normalInset)*activeFrame.normal;
          const auto leavingDirection=(
              value.localDirection.x()*activeFrame.side+
              value.localDirection.y()*activeFrame.normal+
              value.localDirection.z()*activeFrame.tangent).unit();
          const auto enteringDirection=-leavingDirection;
          const auto surfaceState=tessellated->Inside(worldPoint);
          const auto leavingLookAheadState=
              tessellated->Inside(worldPoint+lookAhead*leavingDirection);
          const auto enteringLookAheadState=
              tessellated->Inside(worldPoint+lookAhead*enteringDirection);
          const G4double baseLeavingDistance=
              tessellated->G4TessellatedSolid::DistanceToIn(
                  worldPoint,leavingDirection);
          const G4double correctedLeavingDistance=
              corrected->DistanceToIn(worldPoint,leavingDirection);
          const G4double correctedEnteringDistance=
              corrected->DistanceToIn(worldPoint,enteringDirection);
          const G4double leavingDistanceToOut=
              tessellated->DistanceToOut(
                  worldPoint,leavingDirection,false,nullptr,nullptr);
          const G4bool correctedLeaving=
              correctedLeavingDistance>=kInfinity;
          if(surfaceState!=kSurface||
             leavingLookAheadState!=kOutside||
             enteringLookAheadState==kOutside||
             baseLeavingDistance!=0.||!correctedLeaving||
             correctedEnteringDistance!=0.||leavingDistanceToOut!=0.) {
            G4ExceptionDescription message;
            message << "Gate-3 H2 outer frozen-ray directional contract failed: "
                    << "ray=" << ray
                    << " states=(surface:" << surfaceState
                    << ",leaving_lookahead:" << leavingLookAheadState
                    << ",entering_lookahead:" << enteringLookAheadState
                    << ") base_DIn_leave_mm="
                    << baseLeavingDistance/mm
                    << " corrected_DIn_leave_mm="
                    << correctedLeavingDistance/mm
                    << " corrected_DIn_enter_mm="
                    << correctedEnteringDistance/mm
                    << " DOut_leave_mm=" << leavingDistanceToOut/mm
                    << " lookahead_nm=" << lookAhead/nm
                    << " voxels="
                    << tessellated->GetVoxels().GetCountOfVoxels();
            G4Exception("DetectorConstruction::BuildFanoutPreview",
                        "HODO_GATE3_H2_OUTER_BOUNDARY_CORRECTION_FROZEN_RAY",
                        FatalException,message);
          }
          ++h2OuterFrozenRayAuditCount;
          G4cout << "[gate3 h2 outer boundary correction frozen ray] PASS ray="
                 << ray
                 << " base_DIn_leave_mm=" << baseLeavingDistance/mm
                 << " corrected_DIn_leave_mm="
                 << correctedLeavingDistance/mm
                 << " corrected_DIn_enter_mm="
                 << correctedEnteringDistance/mm
                 << " DOut_leave_mm=" << leavingDistanceToOut/mm
                 << " lookahead_nm=" << lookAhead/nm
                 << " voxels="
                 << tessellated->GetVoxels().GetCountOfVoxels()
                 << G4endl;
        }
      }
      auto* core=new G4LogicalVolume(
          coreSolid, roundFibre ? fSCSF78 : fBCF20,
          prefix+"CoreLV_"+std::to_string(channel));
      ApplyFibreUserLimits(core);
      const G4bool siblingShells=
          roundFibre ? fH1FanoutHierarchy=="sibling_shells"
                     : fH2FanoutHierarchy=="sibling_shells";
      G4VSolid* outerPlacedSolid=outerSolid;
      G4VSolid* innerPlacedSolid=innerSolid;
      if(siblingShells&&!directNativeMiter) {
        outerPlacedSolid=new G4SubtractionSolid(
            prefix+"OuterShellSolid_"+std::to_string(channel),
            outerSolid,singleClad?coreSolid:innerSolid);
        if(!singleClad)
          innerPlacedSolid=new G4SubtractionSolid(
              prefix+"InnerShellSolid_"+std::to_string(channel),
              innerSolid,coreSolid);
      }
      const auto roundSolidLocalPoint=[&](
          const G4ThreeVector& worldPoint) {
        if(!balancedNativeMiter)return worldPoint;
        return roundSweepLocalToWorld.getRotation().inverse()*
            (worldPoint-roundSweepLocalToWorld.getTranslation());
      };
      const auto roundSolidInside=[&](
          const G4VSolid* solid,const G4ThreeVector& worldPoint) {
        return solid->Inside(roundSolidLocalPoint(worldPoint));
      };
      const auto roundSolidNormal=[&](
          const G4VSolid* solid,const G4ThreeVector& worldPoint) {
        const auto localNormal=
            solid->SurfaceNormal(roundSolidLocalPoint(worldPoint));
        return balancedNativeMiter
            ? roundSweepLocalToWorld.getRotation()*localNormal
            : localNormal;
      };
      if(roundFibre&&fH1RoundSweepRepresentation!="tessellated") {
        // Deterministic material-partition probes at every segment midpoint.
        // These catch a disconnected union or a core/shell Boolean inversion
        // before the navigator sees an event.
        for(std::size_t segment=0;segment+1<centres.size();++segment) {
          const auto tangent=
              (centres[segment+1]-centres[segment]).unit();
          auto radial=verticalAxis.cross(tangent);
          if(radial.mag2()<1.e-18)
            radial=axis==0 ? G4ThreeVector(0.,1.,0.)
                           : G4ThreeVector(-1.,0.,0.);
          radial=radial.unit();
          const auto midpoint=.5*(centres[segment]+centres[segment+1]);
          const G4double outerInnerSize=
              singleClad ? coreSize : innerSize;
          const auto outerProbe=
              midpoint+.5*(outerInnerSize+outerSize)*radial;
          const auto coreProbe=midpoint;
          G4bool partitionOk=
              roundSolidInside(coreSolid,coreProbe)!=kOutside &&
              roundSolidInside(outerPlacedSolid,coreProbe)==kOutside &&
              roundSolidInside(coreSolid,outerProbe)==kOutside &&
              roundSolidInside(outerPlacedSolid,outerProbe)!=kOutside;
          if(!singleClad) {
            const auto innerProbe=
                midpoint+.5*(coreSize+innerSize)*radial;
            partitionOk=partitionOk &&
                roundSolidInside(innerPlacedSolid,coreProbe)==kOutside &&
                roundSolidInside(innerPlacedSolid,innerProbe)!=kOutside &&
                roundSolidInside(coreSolid,innerProbe)==kOutside &&
                roundSolidInside(outerPlacedSolid,innerProbe)==kOutside &&
                roundSolidInside(innerPlacedSolid,outerProbe)==kOutside;
          }
          if(!partitionOk) {
            G4ExceptionDescription message;
            message << "Native round material-partition audit failed: prefix="
                    << prefix << " channel=" << channel
                    << " segment=" << segment
                    << " representation=" << fH1RoundSweepRepresentation
                    << " single_clad=" << singleClad;
            G4Exception("DetectorConstruction::BuildFanoutPreview",
                        "HODO_V21_NATIVE_PARTITION",FatalException,message);
          }
        }
        if(directNativeMiter) {
          const G4double coreBandProbe=
              .25*std::min(coreSize,innerSize-coreSize);
          const G4double outerBandProbe=
              .25*(outerSize-innerSize);
          for(std::size_t joint=1;joint+1<centres.size();++joint) {
            const auto incoming=
                (centres[joint]-centres[joint-1]).unit();
            const auto outgoing=
                (centres[joint+1]-centres[joint]).unit();
            const auto miterNormal=(incoming+outgoing).unit();
            const auto e1=miterNormal.orthogonal().unit();
            const auto e2=miterNormal.cross(e1).unit();
            for(G4int phaseClass=0;phaseClass<2;++phaseClass) {
             for(G4int azimuth=0;azimuth<64;++azimuth) {
              const G4double phi=2.*CLHEP::pi*
                  (azimuth+.5*phaseClass)/64.;
              const auto inPlane=
                  std::cos(phi)*e1+std::sin(phi)*e2;
              const G4double projected=std::sqrt(std::max(
                  0.,1.-std::pow(inPlane.dot(incoming),2)));
              if(projected<1.e-12) {
                G4Exception("DetectorConstruction::BuildFanoutPreview",
                            "HODO_V21_NATIVE_INTERFACE",FatalException,
                            "Degenerate direct-annular interface probe.");
              }
              const auto atRadius=[&](G4double radius) {
                return centres[joint]+radius/projected*inPlane;
              };
              auto* coreNeighbour=
                  singleClad ? outerPlacedSolid : innerPlacedSolid;
              const auto coreBoundary=atRadius(coreSize);
              const auto coreInterior=atRadius(coreSize-coreBandProbe);
              const auto neighbourInterior=
                  atRadius(coreSize+coreBandProbe);
              const auto coreBoundaryState=
                  roundSolidInside(coreSolid,coreBoundary);
              const auto neighbourBoundaryState=
                  roundSolidInside(coreNeighbour,coreBoundary);
              G4double coreNormalResidual=kInfinity;
              G4bool coreNormalOk=false;
              if(coreBoundaryState==kSurface&&
                 neighbourBoundaryState==kSurface) {
                const auto coreNormal=
                    roundSolidNormal(coreSolid,coreBoundary);
                const auto neighbourNormal=
                    roundSolidNormal(coreNeighbour,coreBoundary);
                const G4bool finiteNormals=
                    std::isfinite(coreNormal.mag2())&&
                    std::isfinite(neighbourNormal.mag2());
                coreNormalResidual=
                    finiteNormals
                        ? (coreNormal+neighbourNormal).mag()
                        : kInfinity;
                coreNormalOk=finiteNormals&&
                    std::abs(coreNormal.mag2()-1.)<1.e-10&&
                    std::abs(neighbourNormal.mag2()-1.)<1.e-10&&
                    coreNormalResidual<1.e-10;
              }
              G4bool interfaceOk=
                  coreBoundaryState==kSurface &&
                  neighbourBoundaryState==kSurface &&
                  coreNormalOk &&
                  roundSolidInside(coreSolid,coreInterior)==kInside &&
                  roundSolidInside(
                      coreNeighbour,coreInterior)==kOutside &&
                  roundSolidInside(
                      coreSolid,neighbourInterior)==kOutside &&
                  roundSolidInside(
                      coreNeighbour,neighbourInterior)==kInside;
              G4double outerNormalResidual=0.;
              if(!singleClad) {
                const auto outerBoundary=atRadius(innerSize);
                const auto innerInterior=
                    atRadius(innerSize-outerBandProbe);
                const auto outerInterior=
                    atRadius(innerSize+outerBandProbe);
                const auto innerBoundaryState=
                    roundSolidInside(
                        innerPlacedSolid,outerBoundary);
                const auto outerBoundaryState=
                    roundSolidInside(
                        outerPlacedSolid,outerBoundary);
                G4bool outerNormalOk=false;
                if(innerBoundaryState==kSurface&&
                   outerBoundaryState==kSurface) {
                  const auto innerNormal=roundSolidNormal(
                      innerPlacedSolid,outerBoundary);
                  const auto outerNormal=roundSolidNormal(
                      outerPlacedSolid,outerBoundary);
                  const G4bool finiteNormals=
                      std::isfinite(innerNormal.mag2())&&
                      std::isfinite(outerNormal.mag2());
                  outerNormalResidual=
                      finiteNormals
                          ? (innerNormal+outerNormal).mag()
                          : kInfinity;
                  outerNormalOk=finiteNormals&&
                      std::abs(innerNormal.mag2()-1.)<1.e-10&&
                      std::abs(outerNormal.mag2()-1.)<1.e-10&&
                      outerNormalResidual<1.e-10;
                }
                interfaceOk=interfaceOk &&
                    innerBoundaryState==kSurface &&
                    outerBoundaryState==kSurface &&
                    outerNormalOk &&
                    roundSolidInside(
                        innerPlacedSolid,innerInterior)==kInside &&
                    roundSolidInside(
                        outerPlacedSolid,innerInterior)==kOutside &&
                    roundSolidInside(
                        innerPlacedSolid,outerInterior)==kOutside &&
                    roundSolidInside(
                        outerPlacedSolid,outerInterior)==kInside;
              }
              if(!interfaceOk) {
                G4ExceptionDescription message;
                message << "Direct-annular radial interface audit failed: "
                        << "prefix=" << prefix << " channel=" << channel
                        << " joint=" << joint
                        << " phase_class=" << phaseClass
                        << " azimuth=" << azimuth
                        << " single_clad=" << singleClad
                        << " core_boundary="
                        << coreBoundaryState
                        << " neighbour_boundary="
                        << neighbourBoundaryState
                        << " core_normal_sum=" << coreNormalResidual
                        << " outer_normal_sum=" << outerNormalResidual
                        << " core_at_core_interior="
                        << roundSolidInside(coreSolid,coreInterior)
                        << " neighbour_at_core_interior="
                        << roundSolidInside(
                            coreNeighbour,coreInterior)
                        << " core_at_neighbour_interior="
                        << roundSolidInside(
                            coreSolid,neighbourInterior)
                        << " neighbour_at_neighbour_interior="
                        << roundSolidInside(
                            coreNeighbour,neighbourInterior);
                G4Exception("DetectorConstruction::BuildFanoutPreview",
                            "HODO_V21_NATIVE_INTERFACE",
                            FatalException,message);
              }
             }
            }
          }
        }
      }
      auto* outer=new G4LogicalVolume(
          outerPlacedSolid, singleClad ? fInnerCladding : fOuterCladding,
          prefix+"OuterLV_"+std::to_string(channel));
      ApplyFibreUserLimits(outer);
      G4VPhysicalVolume* corePV=nullptr;
      G4LogicalVolume* inner=nullptr;
      G4VPhysicalVolume* innerPV=nullptr;
      if (singleClad && !siblingShells) {
        corePV=new G4PVPlacement(nullptr, {}, core, prefix+"Core", outer,
                                 false, 0, false);
      } else if (!singleClad && !siblingShells) {
        inner=new G4LogicalVolume(
            innerPlacedSolid, fInnerCladding,
            prefix+"InnerLV_"+std::to_string(channel));
        ApplyFibreUserLimits(inner);
        innerPV=new G4PVPlacement(nullptr, {}, inner, prefix+"InnerClad", outer,
                                 false, 0, false);
        corePV=new G4PVPlacement(nullptr, {}, core, prefix+"Core", inner,
                                 false, 0, false);
        inner->SetVisAttributes(G4VisAttributes::GetInvisible());
      }
      new G4LogicalSkinSurface(
          prefix+"CoreSurface_"+std::to_string(channel), core, fFiberSurface);
      auto* coreVis=new G4VisAttributes(
          G4Colour(colour.GetRed(),colour.GetGreen(),colour.GetBlue(),.90));
      coreVis->SetForceSolid(true);
      core->SetVisAttributes(coreVis);
      outer->SetVisAttributes(
          new G4VisAttributes(G4Colour(.7,.9,1.,.15)));
      const G4int localCopy=axis*32+channel;
      const auto physicalName=
          round ? (axis==0 ? "H1_X" : "H1_Y")
                : (axis==0 ? "H2_X" : "H2_Y");
      G4VPhysicalVolume* outerPV=nullptr;
      const G4bool fibrePlacementSurfaceCheck = false;
      if(balancedNativeMiter) {
        outerPV=new G4PVPlacement(
            roundSweepLocalToWorld,outer,physicalName,
            world,false,localCopy,fibrePlacementSurfaceCheck);
      } else {
        outerPV=new G4PVPlacement(
            nullptr,{},outer,physicalName,
            world,false,localCopy,fibrePlacementSurfaceCheck);
      }
      if(siblingShells) {
        if(!singleClad) {
          inner=new G4LogicalVolume(
              innerPlacedSolid,fInnerCladding,
              prefix+"InnerLV_"+std::to_string(channel));
          ApplyFibreUserLimits(inner);
          inner->SetVisAttributes(G4VisAttributes::GetInvisible());
          if(balancedNativeMiter) {
            innerPV=new G4PVPlacement(
                roundSweepLocalToWorld,inner,prefix+"InnerClad",
                world,false,localCopy,fibrePlacementSurfaceCheck);
          } else {
            innerPV=new G4PVPlacement(
                nullptr,{},inner,prefix+"InnerClad",world,
                false,localCopy,fibrePlacementSurfaceCheck);
          }
        }
        if(balancedNativeMiter) {
          corePV=new G4PVPlacement(
              roundSweepLocalToWorld,core,prefix+"Core",
              world,false,localCopy,fibrePlacementSurfaceCheck);
        } else {
          corePV=new G4PVPlacement(
              nullptr,{},core,prefix+"Core",world,
              false,localCopy,fibrePlacementSurfaceCheck);
        }
      }

      const G4double greaseHalfLength=.5*fGreaseGap;
      constexpr G4double windowHalfLength=.05*mm;
      constexpr G4double sensorHalfLength=.005*mm;
      const G4double activeHalf=.5*fSiPMActiveSize;
      auto* greaseSolid=new G4Box(
          prefix+"GreaseSolid_"+std::to_string(channel),
          activeHalf,activeHalf,greaseHalfLength);
      auto* windowSolid=new G4Box(
          prefix+"WindowSolid_"+std::to_string(channel),
          activeHalf,activeHalf,windowHalfLength);
      auto* sensorSolid=new G4Box(
          prefix+"SensorSolid_"+std::to_string(channel),
          activeHalf,activeHalf,sensorHalfLength);
      auto* grease=new G4LogicalVolume(
          greaseSolid,fGrease,prefix+"GreaseLV_"+std::to_string(channel));
      auto* window=new G4LogicalVolume(
          windowSolid,fSiPMWindow,prefix+"WindowLV_"+std::to_string(channel));
      auto* sensor=new G4LogicalVolume(
          sensorSolid,fSilicon,prefix+"SensorLV_"+std::to_string(channel));
      grease->SetVisAttributes(
          new G4VisAttributes(G4Colour(1.,1.,.3,.35)));
      window->SetVisAttributes(
          new G4VisAttributes(G4Colour(.8,.8,.8,.25)));
      sensor->SetVisAttributes(
          new G4VisAttributes(G4Colour(.15,.15,.15)));
      const auto& terminalFrame=stationFrames.back();
      const auto& unit=terminalFrame.tangent;
      const G4ThreeVector offset=fSiPMOffsetU*terminalFrame.side+
                                 fSiPMOffsetV*terminalFrame.normal;
      const G4ThreeVector endpoint=centres.back();
      fFanoutEndpoint[transmissionChannel]=endpoint;
      fFanoutTerminalDirection[transmissionChannel]=unit;
      fFanoutEndpointReady[transmissionChannel]=true;
      const G4RotationMatrix terminalRotation(
          terminalFrame.side,terminalFrame.normal,terminalFrame.tangent);
      const G4Transform3D greaseTransform(
          terminalRotation,endpoint+offset+unit*greaseHalfLength);
      const G4Transform3D windowTransform(
          terminalRotation,
          endpoint+offset+unit*(2.*greaseHalfLength+windowHalfLength));
      const G4Transform3D sensorTransform(
          terminalRotation,
          endpoint+offset+unit*(2.*greaseHalfLength+2.*windowHalfLength+
                                sensorHalfLength));

      // Fail closed if any independently expressed interface is not the same
      // signed terminal plane to floating-point precision.  These invariants
      // are deliberately much tighter than the configured Geant4 surface
      // tolerance; the latter must not hide an internally inconsistent frame.
      const G4double lengthScale=std::max(1.*mm,endpoint.mag());
      const G4double contactTolerance=
          4096.*std::numeric_limits<G4double>::epsilon()*lengthScale;
      const auto failContact=[&](const G4String& detail,G4double residual) {
        G4ExceptionDescription message;
        message << "Fan-out terminal contact invariant failed: prefix="
                << prefix << " channel=" << channel
                << " detail=" << detail
                << " residual_mm=" << residual/mm
                << " allowed_mm=" << contactTolerance/mm;
        G4Exception("DetectorConstruction::BuildFanoutPreview",
                    "HODO_V21_TERMINAL_CONTACT",FatalException,message);
      };
      const G4ThreeVector coreFace=endpoint;
      const G4ThreeVector greaseEntrance=
          greaseTransform.getTranslation()-unit*greaseHalfLength;
      const G4ThreeVector greaseExit=
          greaseTransform.getTranslation()+unit*greaseHalfLength;
      const G4ThreeVector windowEntrance=
          windowTransform.getTranslation()-unit*windowHalfLength;
      const G4ThreeVector windowExit=
          windowTransform.getTranslation()+unit*windowHalfLength;
      const G4ThreeVector sensorEntrance=
          sensorTransform.getTranslation()-unit*sensorHalfLength;
      // A nonzero configured SiPM U/V offset is intentionally allowed.  The
      // invariant here therefore compares the signed axial planes, while the
      // deterministic probes below prove that the fibre face remains covered
      // by grease at the core centre.
      const G4double coreGreasePlane=
          std::abs((greaseEntrance-coreFace).dot(unit));
      const G4double greaseWindowContact=(greaseExit-windowEntrance).mag();
      const G4double windowSensorContact=(windowExit-sensorEntrance).mag();
      if (coreGreasePlane > contactTolerance)
        failContact("core_to_grease_signed_plane",coreGreasePlane);
      if (greaseWindowContact > contactTolerance)
        failContact("grease_to_window",greaseWindowContact);
      if (windowSensorContact > contactTolerance)
        failContact("window_to_sipm",windowSensorContact);

      // Deterministic signed probes certify the expected material ordering on
      // both sides of all three contacts.  No warning is suppressed: an
      // ambiguous surface or reversed/gapped interface is a fatal geometry
      // construction error.
      const G4double probeDistance=std::max(
          64.*G4GeometryTolerance::GetInstance()->GetSurfaceTolerance(),
          1.e-6*mm);
      const G4Transform3D identityTransform;
      const auto& corePlacementTransform=
          balancedNativeMiter
              ? roundSweepLocalToWorld : identityTransform;
      const auto classify=[](const G4VSolid* solid,
                             const G4Transform3D& transform,
                             const G4ThreeVector& worldPoint) {
        const auto localPoint=transform.getRotation().inverse()*
            (worldPoint-transform.getTranslation());
        return solid->Inside(localPoint);
      };
      const auto requireClass=[&](const G4String& detail,
                                  const G4VSolid* solid,
                                  const G4Transform3D& transform,
                                  const G4ThreeVector& point,
                                  EInside expected) {
        const auto observed=classify(solid,transform,point);
        if (observed != expected) {
          G4ExceptionDescription message;
          message << "Fan-out signed contact probe failed: prefix=" << prefix
                  << " channel=" << channel << " detail=" << detail
                  << " expected=" << static_cast<G4int>(expected)
                  << " observed=" << static_cast<G4int>(observed)
                  << " probe_distance_mm=" << probeDistance/mm;
          G4Exception("DetectorConstruction::BuildFanoutPreview",
                      "HODO_V21_TERMINAL_PROBE",FatalException,message);
        }
      };
      requireClass("core_before_cap",coreSolid,corePlacementTransform,
                   coreFace-unit*probeDistance,kInside);
      requireClass("core_on_cap",coreSolid,corePlacementTransform,
                   coreFace,kSurface);
      requireClass("core_after_cap",coreSolid,corePlacementTransform,
                   coreFace+unit*probeDistance,kOutside);
      requireClass("grease_before_entry",greaseSolid,greaseTransform,
                   coreFace-unit*probeDistance,kOutside);
      requireClass("grease_on_entry",greaseSolid,greaseTransform,
                   coreFace,kSurface);
      requireClass("grease_after_entry",greaseSolid,greaseTransform,
                   coreFace+unit*probeDistance,kInside);
      const auto requireTwoSidedContact=
          [&](const G4String& label,const G4ThreeVector& contact,
              const G4VSolid* upstream,const G4Transform3D& upstreamTransform,
              const G4VSolid* downstream,
              const G4Transform3D& downstreamTransform) {
        requireClass(label+"_upstream_inside",upstream,upstreamTransform,
                     contact-unit*probeDistance,kInside);
        requireClass(label+"_upstream_surface",upstream,upstreamTransform,
                     contact,kSurface);
        requireClass(label+"_upstream_outside",upstream,upstreamTransform,
                     contact+unit*probeDistance,kOutside);
        requireClass(label+"_downstream_outside",downstream,
                     downstreamTransform,contact-unit*probeDistance,kOutside);
        requireClass(label+"_downstream_surface",downstream,
                     downstreamTransform,contact,kSurface);
        requireClass(label+"_downstream_inside",downstream,
                     downstreamTransform,contact+unit*probeDistance,kInside);
      };
      requireTwoSidedContact("grease_window",greaseExit,greaseSolid,
                             greaseTransform,windowSolid,windowTransform);
      requireTwoSidedContact("window_sipm",windowExit,windowSolid,
                             windowTransform,sensorSolid,sensorTransform);
      ++terminalFrameAuditCount;

      const G4double stackLength=
          2.*greaseHalfLength+2.*windowHalfLength+2.*sensorHalfLength;
      const CouplingEnvelope envelope{
          plane,channel,prefix,
          endpoint+offset+unit*(.5*stackLength),
          {terminalFrame.side,terminalFrame.normal,terminalFrame.tangent},
          {activeHalf,activeHalf,.5*stackLength}};
      const auto separatingGap=[](const CouplingEnvelope& first,
                                  const CouplingEnvelope& second) {
        G4double maximumGap=-DBL_MAX;
        const auto evaluateAxis=[&](G4ThreeVector axis) {
          if (axis.mag2() < 1.e-24) return;
          axis=axis.unit();
          G4double firstRadius=0.;
          G4double secondRadius=0.;
          for (G4int component=0;component<3;++component) {
            firstRadius+=first.halfLengths[component]*
                std::abs(axis.dot(first.axes[component]));
            secondRadius+=second.halfLengths[component]*
                std::abs(axis.dot(second.axes[component]));
          }
          const G4double gap=
              std::abs(axis.dot(second.centre-first.centre))-
              firstRadius-secondRadius;
          maximumGap=std::max(maximumGap,gap);
        };
        for (const auto& axis:first.axes) evaluateAxis(axis);
        for (const auto& axis:second.axes) evaluateAxis(axis);
        for (const auto& firstAxis:first.axes)
          for (const auto& secondAxis:second.axes)
            evaluateAxis(firstAxis.cross(secondAxis));
        return maximumGap;
      };
      for (const auto& placed:couplingEnvelopes) {
        const auto clearance=separatingGap(placed,envelope);
        if (clearance <=
            G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()) {
          G4ExceptionDescription message;
          message << "Deterministic coupling-stack clearance failed: "
                  << placed.prefix << " channel=" << placed.channel
                  << " versus " << envelope.prefix
                  << " channel=" << envelope.channel
                  << " separating_gap_mm=" << clearance/mm;
          G4Exception("DetectorConstruction::BuildFanoutPreview",
                      "HODO_V21_COUPLING_CLEARANCE",FatalException,message);
        }
        auto& minimum=placed.plane==envelope.plane
            ? minimumSamePlaneClearance
            : minimumDifferentPlaneClearance;
        minimum=std::min(minimum,clearance);
      }
      couplingEnvelopes.push_back(envelope);

      // The generic stochastic placement checker cannot distinguish the
      // intentional coplanar fibre-cap/grease contact from penetration at
      // tessellated cap fan edges.  It produced a non-reproducible 11.8 nm
      // warning after the frames were made exactly coincident.  Only these
      // three terminal placements use pSurfChk=false; the all-channel signed
      // interface probes and deterministic 15-axis stack-clearance audit
      // above replace that exact scope, while every other placement and the
      // global runtime failure regex remain unchanged.
      auto* greasePV=new G4PVPlacement(greaseTransform,grease,
          plane < 2 ? "H1_Grease" : "H2_Grease",
          world,false,localCopy,false);
      new G4LogicalBorderSurface(
          prefix+"CoreEndSurface_"+std::to_string(channel),
          corePV,greasePV,fEndFaceSurface);
      new G4LogicalBorderSurface(
          prefix+"CladEndSurface_"+std::to_string(channel),
          outerPV,greasePV,fEndFaceSurface);
      auto* windowPV=new G4PVPlacement(windowTransform,window,
          "SiPMWindow",world,false,localCopy,false);
      const G4int globalCopy=round ? localCopy : 64+localCopy;
      auto* sensorPV=new G4PVPlacement(sensorTransform,sensor,
          "SiPM",world,false,globalCopy,false);
      new G4LogicalBorderSurface(
          prefix+"SiPMSurface_"+std::to_string(channel),
          windowPV,sensorPV,fSiPMSurface);
      if (deterministicSquareFanoutCheck)
        priorSquareCentrelines.push_back(centres);
    }
  };

  // Keep the two complete hodoscopes on distinct Z bands in the preview.
  // Their 4x8 endpoint grids span 12 mm, so +/-9 mm leaves a visible gap
  // while remaining between the trigger paddles at +/-25 mm.
  buildAssembly(true, 0, fHodoscopeFanoutZ[0],
                "H1XFanout", G4Colour(.10,1.,.30));
  buildAssembly(true, 1, fHodoscopeFanoutZ[0],
                "H1YFanout", G4Colour(.10,1.,.30));
  buildAssembly(false, 0, fHodoscopeFanoutZ[1],
                "H2XFanout", G4Colour(.15,.45,1.));
  buildAssembly(false, 1, fHodoscopeFanoutZ[1],
                "H2YFanout", G4Colour(.15,.45,1.));
  if(!fH2PhysicalRound &&
     (h2OuterBoundaryAuditCount!=2*GetPlaneFibreCount(2)||
      h2OuterBoundaryCloneAuditCount!=1||
      h2OuterFrozenRayAuditCount!=2)) {
    G4ExceptionDescription message;
    message << "Gate-3 H2 outer boundary audit count failed: channels="
            << h2OuterBoundaryAuditCount
            << " clones=" << h2OuterBoundaryCloneAuditCount
            << " frozen_rays=" << h2OuterFrozenRayAuditCount;
    G4Exception("DetectorConstruction::BuildFanoutPreview",
                "HODO_GATE3_H2_OUTER_BOUNDARY_CORRECTION_COUNT",
                FatalException,message);
  }
  G4cout << "[gate3 h2 outer boundary correction audit] PASS channels="
         << h2OuterBoundaryAuditCount
         << " clones=" << h2OuterBoundaryCloneAuditCount
         << " frozen_rays=" << h2OuterFrozenRayAuditCount
         << " h2_physical_round=" << fH2PhysicalRound
         << " facets_per_channel=400"
         << " voxelization=baseline"
         << " fallback=base_zero_after_8_surface_lookaheads"
         << G4endl;
  const auto physicalFibreCount=64+2*GetPlaneFibreCount(2);
  if (terminalFrameAuditCount != physicalFibreCount) {
    G4ExceptionDescription message;
    message << "Expected " << physicalFibreCount << " audited fan-out terminal frames, got "
            << terminalFrameAuditCount;
    G4Exception("DetectorConstruction::BuildFanoutPreview",
                "HODO_V21_TERMINAL_AUDIT_COUNT",FatalException,message);
  }
  G4cout << "[v21 terminal frame audit] PASS channels="
         << terminalFrameAuditCount
         << " shared_frame=core_cap,grease,window,sipm"
         << " signed_contact_probes=PASS" << G4endl;
  if (couplingEnvelopes.size() != static_cast<std::size_t>(physicalFibreCount) ||
      !std::isfinite(minimumSamePlaneClearance) ||
      !std::isfinite(minimumDifferentPlaneClearance)) {
    G4ExceptionDescription message;
    message << "Incomplete deterministic coupling clearance audit: envelopes="
            << couplingEnvelopes.size()
            << " same_plane_min_mm=" << minimumSamePlaneClearance/mm
            << " different_plane_min_mm="
            << minimumDifferentPlaneClearance/mm;
    G4Exception("DetectorConstruction::BuildFanoutPreview",
                "HODO_V21_COUPLING_CLEARANCE_COUNT",FatalException,message);
  }
  G4cout << "[v21 coupling clearance audit] PASS envelopes="
         << couplingEnvelopes.size()
         << " same_plane_min_mm=" << minimumSamePlaneClearance/mm
         << " different_plane_min_mm=" << minimumDifferentPlaneClearance/mm
         << " method=deterministic_15_axis_OBB_SAT"
         << " generic_terminal_pSurfChk=replaced" << G4endl;
  const auto printTransmission=[this](G4int first,const G4String& label) {
    G4double minT=1.,maxT=0.,sumT=0.,minR=DBL_MAX,maxR=0.;
    G4int count=0;
    for(G4int channel=first;channel<first+64;++channel) {
      if(!HasPhysicalFibre(channel/32,channel%32))continue;
      ++count;
      minT=std::min(minT,fFanoutChannelTransmission[channel]);
      maxT=std::max(maxT,fFanoutChannelTransmission[channel]);
      sumT+=fFanoutChannelTransmission[channel];
      if(fFanoutMinRadius[channel]>0.) {
        minR=std::min(minR,fFanoutMinRadius[channel]);
        maxR=std::max(maxR,fFanoutMinRadius[channel]);
      }
    }
    G4cout<<"[fanout channel model] "<<label
          <<" transmission min/mean/max="<<minT<<"/"<<sumT/count<<"/"<<maxT
          <<" radius_mm min/max="<<minR/mm<<"/"<<maxR/mm<<G4endl;
  };
  printTransmission(0,"H1");
  printTransmission(64,"H2");
  G4cout << "[fanout optical] REPLACEMENT hodoscope: continuous core/cladding "
            "32 mm active region + 70 mm 48-station smootherstep fan-out + "
            "0.10 mm grease gap and 1.3x1.3 mm S14160-1310PS optical faces "
            "at assumed 4 mm pitch; "
            "holder omitted and legacy 165 mm geometry omitted."
         << G4endl;
}

void DetectorConstruction::BuildFanoutPreviewLite(G4LogicalVolume* world) {
  constexpr G4double straightExitLength = 20.*mm;
  constexpr G4double fanoutLength = 100.*mm;
  constexpr G4double pitch = 4.*mm;
  constexpr G4double activeHalfLength = 16.*mm;
  constexpr G4double roundRadius = .24*mm;
  constexpr G4double squareHalf = .24*mm;
  {
    // Photo-like viewer, not physics geometry:
    //   - central X/Y hodoscope ribbons cross in a compact active region
    //   - each plane leaves the active region as a broad ribbon
    //   - only after that does it gently spread toward a black readout bar
    //   - the opposite side is shown as a short taped/far-end strip
    //
    // The previous lightweight viewer drew individual 3-D swept fibres and
    // looked like a comb of floating sticks.  That was faithful to neither
    // the photographs nor visual inspection needs, so this top-level viewer
    // deliberately communicates the assembled ribbon layout instead.
    constexpr G4double bandHalf = 14.*mm;
    constexpr G4double ribbonThickness = .18*mm;
    const G4double farEndLength = fFibreFarEndExtension>0.
        ? fFibreFarEndExtension : 34.*mm;
    constexpr G4double readoutBarHalf = 17.*mm;
    constexpr G4double readoutBarThick = 1.2*mm;
    constexpr G4double tapeHalfLength = 4.*mm;

    auto* h1xVis = new G4VisAttributes(G4Colour(.10, 1.00, .35, .34));
    auto* h1yVis = new G4VisAttributes(G4Colour(.45, 1.00, .45, .34));
    auto* h2xVis = new G4VisAttributes(G4Colour(.15, .55, 1.00, .34));
    auto* h2yVis = new G4VisAttributes(G4Colour(.35, .80, 1.00, .34));
    auto* guideVis = new G4VisAttributes(G4Colour(.96, .96, .78, .86));
    auto* boardVis = new G4VisAttributes(G4Colour(.02, .02, .02, .86));
    auto* tapeVis = new G4VisAttributes(G4Colour(1.00, .72, .05, .70));
    for(auto* vis : {h1xVis,h1yVis,h2xVis,h2yVis,guideVis,boardVis,tapeVis})
      vis->SetForceSolid(true);

    const auto rotationAlong = [](const G4ThreeVector& direction) {
      const auto tangent = direction.unit();
      const G4ThreeVector reference =
          std::abs(tangent.z()) < .90 ? G4ThreeVector(0.,0.,1.)
                                      : G4ThreeVector(1.,0.,0.);
      const auto side = reference.cross(tangent).unit();
      const auto normal = tangent.cross(side).unit();
      return new G4RotationMatrix(side, normal, tangent);
    };
    const auto smoother = [](G4double value) {
      const G4double u = std::clamp(value,0.,1.);
      return u*u*u*(10.-15.*u+6.*u*u);
    };
    const auto planePoint =
        [this](G4int plane, G4int axis, G4double longitudinal,
               G4double transverse, G4double z) {
          const auto nominal = axis == 0
              ? G4ThreeVector(longitudinal, transverse, z)
              : G4ThreeVector(transverse, longitudinal, z);
          return PlaneNominalToWorldPosition(plane, nominal);
        };
    const auto placeBoxAlong =
        [&](const G4String& name, G4Material* material, G4VisAttributes* vis,
            const G4ThreeVector& a, const G4ThreeVector& b,
            G4double halfWidth, G4double halfThickness, G4int copyNo) {
          const auto delta = b-a;
          const auto length = delta.mag();
          if(length <= 1.e-6*mm) return;
          auto* solid = new G4Box(name+"Solid", halfWidth, halfThickness,
                                  .5*length);
          auto* logical = new G4LogicalVolume(solid, material, name+"LV");
          logical->SetVisAttributes(vis);
          new G4PVPlacement(rotationAlong(delta), .5*(a+b), logical, name,
                            world, false, copyNo, false);
        };
    const auto placeBoard =
        [&](const G4String& name, G4int plane, G4int axis,
            G4double longitudinal, G4double z, G4bool readout) {
          auto* solid = axis == 0
              ? new G4Box(name+"Solid", readoutBarThick, readoutBarHalf,
                          3.5*mm)
              : new G4Box(name+"Solid", readoutBarHalf, readoutBarThick,
                          3.5*mm);
          auto* logical = new G4LogicalVolume(
              solid, readout ? fSilicon : fTriggerScint, name+"LV");
          logical->SetVisAttributes(readout ? boardVis : tapeVis);
          const auto centre = planePoint(plane,axis,longitudinal,0.,z);
          new G4PVPlacement(nullptr, centre, logical, name, world, false,
                            plane+(readout?100:200), false);
        };
    const auto placeGuide =
        [&](const G4String& name, G4int plane, G4int axis,
            G4double offset, G4double sourceZ, G4double targetZ,
            G4double targetOffset, G4int copyBase) {
          std::vector<G4ThreeVector> point;
          point.push_back(planePoint(plane,axis,-activeHalfLength-farEndLength,
                                     offset,sourceZ));
          point.push_back(planePoint(plane,axis,-activeHalfLength,
                                     offset,sourceZ));
          point.push_back(planePoint(plane,axis,activeHalfLength,
                                     offset,sourceZ));
          point.push_back(planePoint(plane,axis,
                                     activeHalfLength+straightExitLength,
                                     offset,sourceZ));
          for(G4int step=1; step<=fFanoutPreviewLiteSegments; ++step) {
            const G4double t =
                static_cast<G4double>(step)/fFanoutPreviewLiteSegments;
            const G4double smooth = step == fFanoutPreviewLiteSegments
                ? 1. : smoother(t);
            const G4double longitudinal =
                activeHalfLength+straightExitLength+fanoutLength*t;
            const G4double transverse = offset+(targetOffset-offset)*smooth;
            const G4double z = sourceZ+(targetZ-sourceZ)*smooth;
            point.push_back(planePoint(plane,axis,longitudinal,transverse,z));
          }
          for(std::size_t i=1;i<point.size();++i)
            placeBoxAlong(name+"Guide"+std::to_string(i), fSCSF78, guideVis,
                          point[i-1], point[i], .06*mm, .06*mm,
                          copyBase+static_cast<G4int>(i));
        };
    const auto buildRibbonPlane =
        [&](G4bool round, G4int axis, const G4String& prefix,
            G4VisAttributes* ribbonVis) {
          const G4int plane = (round ? 0 : 2)+axis;
          const G4bool roundFibre = round || (!round && fH2PhysicalRound);
          const G4double detectorZ =
              round ? fHodoscopeFanoutZ[0] : fHodoscopeFanoutZ[1];
          const G4double sourceZ = GetValidationSurfaceZ(plane+2);
          const G4double targetZ =
              round ? detectorZ+axis*2.2*mm
                    : (roundFibre ? detectorZ-axis*2.2*mm
                                  : detectorZ-axis*1.*mm);
          const G4Material* constMaterial = roundFibre ? fSCSF78 : fBCF20;
          auto* material = const_cast<G4Material*>(constMaterial);

          // Far/taped side, compact active crossing, straight exit, and
          // tapered fan-out ribbon.  Boxes are deliberately thin in local-Z
          // so the photo-like top-view structure is readable.
          placeBoxAlong(prefix+"FarRibbon", material, ribbonVis,
                        planePoint(plane,axis,-activeHalfLength-farEndLength,
                                   0.,sourceZ),
                        planePoint(plane,axis,-activeHalfLength,0.,sourceZ),
                        bandHalf, ribbonThickness, plane*1000+1);
          placeBoxAlong(prefix+"ActiveRibbon", material, ribbonVis,
                        planePoint(plane,axis,-activeHalfLength,0.,sourceZ),
                        planePoint(plane,axis,activeHalfLength,0.,sourceZ),
                        bandHalf, ribbonThickness, plane*1000+2);
          placeBoxAlong(prefix+"StraightExitRibbon", material, ribbonVis,
                        planePoint(plane,axis,activeHalfLength,0.,sourceZ),
                        planePoint(plane,axis,
                                   activeHalfLength+straightExitLength,
                                   0.,sourceZ),
                        bandHalf, ribbonThickness, plane*1000+3);
          G4ThreeVector previous =
              planePoint(plane,axis,activeHalfLength+straightExitLength,
                         0.,sourceZ);
          for(G4int step=1; step<=fFanoutPreviewLiteSegments; ++step) {
            const G4double t =
                static_cast<G4double>(step)/fFanoutPreviewLiteSegments;
            const G4double smooth = step == fFanoutPreviewLiteSegments
                ? 1. : smoother(t);
            const G4double longitudinal =
                activeHalfLength+straightExitLength+fanoutLength*t;
            const G4double z = sourceZ+(targetZ-sourceZ)*smooth;
            const auto current = planePoint(plane,axis,longitudinal,0.,z);
            const G4double width = bandHalf+(16.*mm-bandHalf)*smooth;
            placeBoxAlong(prefix+"FanoutRibbon"+std::to_string(step),
                          material, ribbonVis, previous, current,
                          width, ribbonThickness, plane*1000+10+step);
            previous=current;
          }
          placeBoard(prefix+"ReadoutBar", plane, axis,
                     activeHalfLength+straightExitLength+fanoutLength+
                     3.*mm, targetZ, true);
          placeBoard(prefix+"TapeEnd", plane, axis,
                     -activeHalfLength-farEndLength-tapeHalfLength,
                     sourceZ, false);

          // Sparse guide lines make the ribbon look like a fibre bundle
          // without turning the GUI into the earlier comb/forest.
          const std::array<G4double,5> offsets{{
              -11.*mm,-5.5*mm,0.*mm,5.5*mm,11.*mm}};
          for(std::size_t index=0; index<offsets.size(); ++index) {
            const G4double targetOffset =
                offsets[index]*1.12;
            placeGuide(prefix, plane, axis, offsets[index], sourceZ, targetZ,
                       targetOffset, plane*10000+
                       static_cast<G4int>(index)*100);
          }
        };

    buildRibbonPlane(true, 0, "H1XPhoto", h1xVis);
    buildRibbonPlane(true, 1, "H1YPhoto", h1yVis);
    buildRibbonPlane(false, 0, "H2XPhoto", h2xVis);
    buildRibbonPlane(false, 1, "H2YPhoto", h2yVis);

    G4cout << "[fanoutPreviewLite] photo-like ribbon viewer: central cross + "
              "straight exits + tapered readout ribbons; visualization only."
           << G4endl;
    return;
  }

  auto* h1Vis = new G4VisAttributes(G4Colour(.1, 1., .35, .92));
  auto* h1YVis = new G4VisAttributes(G4Colour(.45, 1., .45, .92));
  auto* h2Vis = new G4VisAttributes(G4Colour(.15, .55, 1., .92));
  auto* h2YVis = new G4VisAttributes(G4Colour(.35, .80, 1., .92));
  auto* boardVis = new G4VisAttributes(G4Colour(.12, .12, .12, .65));
  auto* sipmVis = new G4VisAttributes(G4Colour(1., .92, .15, .85));
  auto* activeVis = new G4VisAttributes(G4Colour(1., 1., 1., .24));
  for(auto* vis : {h1Vis, h1YVis, h2Vis, h2YVis, boardVis, sipmVis}) {
    vis->SetForceSolid(true);
  }
  activeVis->SetForceWireframe(true);

  const auto shouldDrawChannel = [this](G4int channel) {
    return channel == 0 || channel == 31 ||
           channel%fFanoutPreviewLiteChannelStride == 0;
  };
  const auto smoother = [](G4double value) {
    const G4double u = std::clamp(value, 0., 1.);
    return u*u*u*(10.-15.*u+6.*u*u);
  };
  const auto buildRotation = [](const G4ThreeVector& direction) {
    const auto tangent = direction.unit();
    const G4ThreeVector reference =
        std::abs(tangent.z()) < .90 ? G4ThreeVector(0.,0.,1.)
                                    : G4ThreeVector(1.,0.,0.);
    const auto side = reference.cross(tangent).unit();
    const auto normal = tangent.cross(side).unit();
    return new G4RotationMatrix(side, normal, tangent);
  };

  const auto placeSegment =
      [&](const G4String& prefix, G4int plane, G4int channel, G4int segment,
          G4bool round, const G4ThreeVector& a, const G4ThreeVector& b,
          G4VisAttributes* vis) {
        const auto delta = b-a;
        const auto length = delta.mag();
        if(length <= 1.e-6*mm) return;
        const auto centre = .5*(a+b);
        G4LogicalVolume* logical = nullptr;
        if(round) {
          auto* solid = new G4Tubs(
              prefix+"LiteSegSolid_"+std::to_string(channel)+"_"+
                  std::to_string(segment),
              0., roundRadius, .5*length, 0., 360.*deg);
          logical = new G4LogicalVolume(
              solid, fSCSF78,
              prefix+"LiteSegLV_"+std::to_string(channel)+"_"+
                  std::to_string(segment));
        } else {
          auto* solid = new G4Box(
              prefix+"LiteSegSolid_"+std::to_string(channel)+"_"+
                  std::to_string(segment),
              squareHalf, squareHalf, .5*length);
          logical = new G4LogicalVolume(
              solid, fBCF20,
              prefix+"LiteSegLV_"+std::to_string(channel)+"_"+
                  std::to_string(segment));
        }
        logical->SetVisAttributes(vis);
        new G4PVPlacement(
            buildRotation(delta), centre, logical,
            prefix+"LiteSegment", world, false, plane*1000+channel*64+segment,
            false);
      };

  const auto placeMarker =
      [&](const G4String& prefix, G4int plane, G4int channel,
          const G4ThreeVector& centre, G4bool round) {
        G4LogicalVolume* logical = nullptr;
        if(round) {
          auto* solid = new G4Tubs(prefix+"LiteSiPMSolid_"+std::to_string(channel),
                                   0., .45*mm, .12*mm, 0., 360.*deg);
          logical = new G4LogicalVolume(
              solid, fSilicon,
              prefix+"LiteSiPMLV_"+std::to_string(channel));
        } else {
          auto* solid = new G4Box(prefix+"LiteSiPMSolid_"+std::to_string(channel),
                                  .45*mm, .45*mm, .12*mm);
          logical = new G4LogicalVolume(
              solid, fSilicon,
              prefix+"LiteSiPMLV_"+std::to_string(channel));
        }
        logical->SetVisAttributes(sipmVis);
        new G4PVPlacement(nullptr, centre, logical, prefix+"LiteSiPM",
                          world, false, plane*1000+channel, false);
      };

  const auto placePlaneBox =
      [&](const G4String& name, G4int plane, G4int axis, G4bool board) {
        const G4double zCentre =
            board ? GetPlaneSourceZ(plane, 15) : GetValidationSurfaceZ(plane+2);
        const G4double longitudinal =
            board ? activeHalfLength+straightExitLength+fanoutLength
                  : 0.;
        const G4double transverseHalf = board ? 16.*mm : 17.*mm;
        const G4double zHalf = board ? 8.5*mm : .18*mm;
        const G4double thin = board ? .35*mm : 17.*mm;
        auto* solid = board
            ? new G4Box(name+"Solid", thin, transverseHalf, zHalf)
            : new G4Box(name+"Solid", 17.*mm, 17.*mm, zHalf);
        auto* logical = new G4LogicalVolume(
            solid, board ? fSilicon : fTriggerScint, name+"LV");
        logical->SetVisAttributes(board ? boardVis : activeVis);
        G4ThreeVector centre(longitudinal, 0., zCentre);
        centre = PlaneNominalToWorldPosition(plane, centre);
        new G4PVPlacement(nullptr, centre, logical, name, world, false, plane,
                          false);
      };

  const auto buildPlane =
      [&](G4bool round, G4int axis, const G4String& prefix,
          G4VisAttributes* vis) {
        const G4int plane = (round ? 0 : 2)+axis;
        const G4bool roundFibre = round || (!round && fH2PhysicalRound);
        const G4double detectorZ =
            round ? fHodoscopeFanoutZ[0] : fHodoscopeFanoutZ[1];
        const G4double planeZ = round ? detectorZ+axis*2.2*mm
                                      : (roundFibre ? detectorZ-axis*2.2*mm
                                                    : detectorZ-axis*1.*mm);
        placePlaneBox(prefix+"ActiveEnvelope", plane, axis, false);
        placePlaneBox(prefix+"SiPMBoard", plane, axis, true);

        for(G4int channel=0; channel<32; ++channel) {
          if(!shouldDrawChannel(channel)) continue;
          const G4int row = GetPlaneFanoutRow(plane, channel);
          const G4int col = GetPlaneFanoutColumn(plane, channel);
          const G4double sourceTransverse =
              GetPlaneSourceTransverse(plane, channel);
          const G4double sourceZ = GetPlaneSourceZ(plane, channel);
          const G4double targetTransverse = (col-3.5)*pitch;
          const G4double targetZ = planeZ+(row-1.5)*pitch;

          std::vector<G4ThreeVector> centres;
          centres.reserve(fFanoutPreviewLiteSegments+2);
          // In the lightweight viewer, the compact detector region is shown
          // by the active envelope.  Draw individual representative fibres
          // only after they leave that region; otherwise the GUI becomes a
          // forest of floating sticks and no longer communicates the geometry.
          centres.push_back(
              G4ThreeVector(activeHalfLength, sourceTransverse, sourceZ));
          centres.push_back(
              G4ThreeVector(activeHalfLength+straightExitLength,
                            sourceTransverse, sourceZ));
          for(G4int step=1; step<=fFanoutPreviewLiteSegments; ++step) {
            const G4double t =
                static_cast<G4double>(step)/fFanoutPreviewLiteSegments;
            const G4double smooth = step == fFanoutPreviewLiteSegments
                ? 1. : smoother(t);
            const G4double longitudinal =
                activeHalfLength+straightExitLength+fanoutLength*t;
            const G4double transverse =
                sourceTransverse+(targetTransverse-sourceTransverse)*smooth;
            const G4double curveZ = sourceZ+(targetZ-sourceZ)*smooth;
            centres.push_back(
                G4ThreeVector(longitudinal, transverse, curveZ));
          }

          for(auto& centre : centres)
            centre = PlaneNominalToWorldPosition(plane, centre);

          const G4int globalChannel=(round ? 0 : 64)+axis*32+channel;
          fFanoutEndpoint[globalChannel]=centres.back();
          fFanoutTerminalDirection[globalChannel]=
              (centres.back()-centres[centres.size()-2]).unit();
          fFanoutEndpointReady[globalChannel]=true;

          for(std::size_t segment=1; segment<centres.size(); ++segment) {
            placeSegment(prefix, plane, channel,
                         static_cast<G4int>(segment-1), roundFibre,
                         centres[segment-1], centres[segment], vis);
          }
          placeMarker(prefix, plane, channel, centres.back(), roundFibre);
        }
      };

  buildPlane(true, 0, "H1X", h1Vis);
  buildPlane(true, 1, "H1Y", h1YVis);
  buildPlane(false, 0, "H2X", h2Vis);
  buildPlane(false, 1, "H2Y", h2YVis);

  G4cout << "[fanoutPreviewLite] visualization-only geometry: stride="
         << fFanoutPreviewLiteChannelStride
         << " segments=" << fFanoutPreviewLiteSegments
         << " active straight half-length=" << activeHalfLength/mm
         << " mm, straight exit=" << straightExitLength/mm
         << " mm, fan-out length=" << fanoutLength/mm << " mm."
         << G4endl;
}

void DetectorConstruction::BuildTrigger(G4LogicalVolume* world, G4int trigger) {
  if(trigger<0||trigger>=4||!GetTriggerEnabled(trigger))return;
  auto* solid = new G4Box("TriggerSolid", GetTriggerHalfLength(trigger),
                          GetTriggerHalfWidth(trigger),
                          GetTriggerHalfThickness(trigger));
  const auto suffix = std::to_string(trigger+1);
  auto* logical = new G4LogicalVolume(solid, fTriggerScint,
                                      "Trigger"+suffix+"LV");
  new G4LogicalSkinSurface("Trigger"+suffix+"BlackWrap", logical,
                           fBlackWrapSurface);
  auto* rotation = new G4RotationMatrix;
  rotation->rotateZ(GetTriggerRotationAngle(trigger));
  const auto position=GetTriggerCentre(trigger);
  new G4PVPlacement(rotation, position, logical, "Trigger"+suffix,
                    world, false, 0, true);
  const std::array<G4Colour,4> colours{{
      G4Colour(1., .2, .2, .45),
      G4Colour(.9, .5, .1, .45),
      G4Colour(.2, .8, 1., .40),
      G4Colour(.3, 1., .3, .40)}};
  auto* vis = new G4VisAttributes(colours[trigger]);
  vis->SetForceSolid(true);
  logical->SetVisAttributes(vis);
}

void DetectorConstruction::PlaceUpstreamMaterials(G4LogicalVolume* world) {
  auto* nist = G4NistManager::Instance();
  auto* concrete = nist->FindOrBuildMaterial("G4_CONCRETE");
  auto* slabSolid = new G4Box("ConcreteSlabSolid", 1.2*m, 1.2*m,
                              fConcreteThickness/2.);
  auto* slab = new G4LogicalVolume(slabSolid, concrete, "ConcreteSlabLV");
  auto* slabVis = new G4VisAttributes(G4Colour(.45, .45, .45, .20));
  slabVis->SetForceSolid(true);
  slab->SetVisAttributes(slabVis);
  for (G4int i=0; i<3; ++i) {
    const auto z = (50. + i*25.)*cm;
    new G4PVPlacement(nullptr, {0.,0.,z}, slab, "ConcreteFloor", world,
                      false, i, true);
  }

  // Composition/thickness are placeholders until the real fabric is measured.
  auto* polymer = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  auto* clothSolid = new G4Box("BlackoutSheetSolid", 6.*cm, 6.*cm,
                               fBlackoutThickness/2.);
  auto* cloth = new G4LogicalVolume(clothSolid, polymer, "BlackoutSheetLV");
  cloth->SetVisAttributes(new G4VisAttributes(G4Colour(.05,.05,.05,.65)));
  new G4PVPlacement(nullptr, {0.,0.,35.*mm}, cloth, "BlackoutSheet", world,
                    false, 0, true);
}

void DetectorConstruction::BuildRoundHodoscope(G4LogicalVolume* world, G4double z) {
  constexpr G4double halfLength = 82.5*mm;
  // The measured nominal diameter is 1 mm. A 1 um construction clearance
  // removes floating-point overlaps between nominally touching neighbours.
  constexpr G4double outerRadius = 0.4995*mm;
  const G4bool singleClad = GetH1SingleClad();
  const G4double innerRadius = singleClad ? outerRadius : 0.4845*mm;
  const G4double coreRadius = singleClad ? 0.4795*mm : 0.4695*mm;
  constexpr G4double greaseHalfLength = 0.05*mm;
  constexpr G4double windowHalfLength = 0.05*mm;
  constexpr G4double sensorHalfLength = 0.005*mm;
  constexpr G4double mirrorHalfLength = 0.01*mm;

  auto* outerSolid = new G4Tubs("H1OuterCladSolid", 0., outerRadius, halfLength, 0., 360.*deg);
  auto* innerSolid = new G4Tubs("H1InnerCladSolid", 0., innerRadius, halfLength, 0., 360.*deg);
  auto* coreSolid = new G4Tubs("H1CoreSolid", 0., coreRadius, halfLength, 0., 360.*deg);
  auto* outer = new G4LogicalVolume(outerSolid,
      singleClad ? fInnerCladding : fOuterCladding, "H1OuterCladLV");
  auto* inner = new G4LogicalVolume(innerSolid, fInnerCladding, "H1InnerCladLV");
  auto* core = new G4LogicalVolume(coreSolid, fSCSF78, "H1CoreLV");
  ApplyFibreUserLimits(outer);
  ApplyFibreUserLimits(inner);
  ApplyFibreUserLimits(core);
  if (singleClad) {
    new G4PVPlacement(nullptr, {}, core, "H1Core", outer, false, 0, false);
  } else {
    new G4PVPlacement(nullptr, {}, inner, "H1InnerClad", outer, false, 0, false);
    new G4PVPlacement(nullptr, {}, core, "H1Core", inner, false, 0, false);
  }
  new G4LogicalSkinSurface("H1CoreSurface", core, fFiberSurface);

  auto* greaseSolid = new G4Tubs("H1GreaseSolid", 0., outerRadius, greaseHalfLength, 0., 360.*deg);
  // The real S14160-1310PS active area is 1.3 x 1.3 mm, but the unmeasured fan-out is not explicitly bent
  // in this geometry. Model only the optical column behind the 1 mm fibre to
  // avoid overlapping neighbouring packages at the detector-plane pitch.
  auto* windowSolid = new G4Tubs("H1SiPMWindowSolid", 0., outerRadius, windowHalfLength, 0., 360.*deg);
  auto* sensorSolid = new G4Tubs("H1SiPMSolid", 0., outerRadius, sensorHalfLength, 0., 360.*deg);
  auto* mirrorSolid = new G4Tubs("H1MirrorSolid", 0., outerRadius, mirrorHalfLength, 0., 360.*deg);
  auto* grease = new G4LogicalVolume(greaseSolid, fGrease, "H1GreaseLV");
  auto* window = new G4LogicalVolume(windowSolid, fSiPMWindow, "H1SiPMWindowLV");
  auto* sensor = new G4LogicalVolume(sensorSolid, fSilicon, "H1SiPMLV");
  auto* mirror = new G4LogicalVolume(mirrorSolid, fSilicon, "H1MirrorLV");
  new G4LogicalSkinSurface("H1MirrorSkin", mirror, fMirrorSurface);

  auto* coreVis = new G4VisAttributes(G4Colour(.1, 1., .35));
  coreVis->SetForceSolid(true); core->SetVisAttributes(coreVis);
  outer->SetVisAttributes(new G4VisAttributes(G4Colour(.7, .9, 1., .18)));
  grease->SetVisAttributes(new G4VisAttributes(G4Colour(1., 1., .3, .35)));
  window->SetVisAttributes(new G4VisAttributes(G4Colour(.8, .8, .8, .25)));
  sensor->SetVisAttributes(new G4VisAttributes(G4Colour(.15, .15, .15)));

  constexpr G4double rowDz = 0.8660254*mm;
  for (G4int axis = 0; axis < 2; ++axis) {
    const G4int plane = axis;
    auto* rotation = new G4RotationMatrix;
    if (axis == 0) rotation->rotateY(90.*deg);
    else rotation->rotateX(90.*deg);
    if (GetPlaneYaw(plane) != 0.) {
      G4RotationMatrix inverseYaw;
      inverseYaw.rotateZ(-GetPlaneYaw(plane));
      *rotation = (*rotation) * inverseYaw;
    }
    const G4ThreeVector nominalUnit =
        axis == 0 ? G4ThreeVector(1.,0.,0.) : G4ThreeVector(0.,1.,0.);
    const auto unit = PlaneNominalToWorldDirection(plane, nominalUnit);
    for (G4int channel = 0; channel < 32; ++channel) {
      const G4int row = channel/16;
      const G4int index = channel%16;
      const G4double transverse = (index-7.5)*mm + row*.5*mm;
      const G4double localZ = z + row*rowDz + axis*2.2*mm;
      const G4ThreeVector nominalCentre =
          axis == 0 ? G4ThreeVector(0., transverse, localZ)
                    : G4ThreeVector(transverse, 0., localZ);
      const auto centre = PlaneNominalToWorldPosition(plane, nominalCentre);
      const G4int copy = axis*32 + channel;
      new G4PVPlacement(rotation, centre, outer, axis == 0 ? "H1_X" : "H1_Y",
                        world, false, copy, true);
      if (!fFanoutPreview) {
      new G4PVPlacement(rotation, centre + unit*(halfLength+greaseHalfLength), grease,
                        "H1_Grease", world, false, copy, true);
      auto* windowPV = new G4PVPlacement(rotation,
                        centre + unit*(halfLength+2.*greaseHalfLength+windowHalfLength),
                        window, "SiPMWindow", world, false, copy, true);
      auto* sensorPV = new G4PVPlacement(rotation,
                        centre + unit*(halfLength+2.*greaseHalfLength+2.*windowHalfLength+sensorHalfLength),
                        sensor, "SiPM", world, false, copy, true);
      new G4LogicalBorderSurface("H1GreaseToSiPM_"+std::to_string(copy),
                                 windowPV, sensorPV, fSiPMSurface);
      }
      if (fFarEndReflectivity > 0.)
        new G4PVPlacement(rotation, centre - unit*(halfLength+mirrorHalfLength), mirror,
                          "H1_FarMirror", world, false, copy, true);
    }
  }
}

void DetectorConstruction::BuildRoundHodoscopeH2(G4LogicalVolume* world, G4double z) {
  constexpr G4double halfLength = 82.5*mm;
  constexpr G4double outerRadius = 0.4995*mm;
  const G4bool singleClad = GetH2SingleClad();
  const G4double innerRadius = singleClad ? outerRadius : 0.4845*mm;
  const G4double coreRadius = singleClad ? 0.4795*mm : 0.4695*mm;
  constexpr G4double greaseHalfLength = 0.05*mm;
  constexpr G4double windowHalfLength = 0.05*mm;
  constexpr G4double sensorHalfLength = 0.005*mm;
  constexpr G4double mirrorHalfLength = 0.01*mm;

  auto* outerSolid = new G4Tubs("H2RoundOuterCladSolid", 0., outerRadius, halfLength, 0., 360.*deg);
  auto* innerSolid = new G4Tubs("H2RoundInnerCladSolid", 0., innerRadius, halfLength, 0., 360.*deg);
  auto* coreSolid = new G4Tubs("H2RoundCoreSolid", 0., coreRadius, halfLength, 0., 360.*deg);
  auto* outer = new G4LogicalVolume(outerSolid,
      singleClad ? fInnerCladding : fOuterCladding, "H2RoundOuterCladLV");
  auto* inner = new G4LogicalVolume(innerSolid, fInnerCladding, "H2RoundInnerCladLV");
  auto* core = new G4LogicalVolume(coreSolid, fSCSF78, "H2RoundCoreLV");
  ApplyFibreUserLimits(outer);
  ApplyFibreUserLimits(inner);
  ApplyFibreUserLimits(core);
  if (singleClad) {
    new G4PVPlacement(nullptr, {}, core, "H2Core", outer, false, 0, false);
  } else {
    new G4PVPlacement(nullptr, {}, inner, "H2InnerClad", outer, false, 0, false);
    new G4PVPlacement(nullptr, {}, core, "H2Core", inner, false, 0, false);
  }
  new G4LogicalSkinSurface("H2RoundCoreSurface", core, fFiberSurface);

  auto* greaseSolid = new G4Tubs("H2RoundGreaseSolid", 0., outerRadius, greaseHalfLength, 0., 360.*deg);
  auto* windowSolid = new G4Tubs("H2RoundSiPMWindowSolid", 0., outerRadius, windowHalfLength, 0., 360.*deg);
  auto* sensorSolid = new G4Tubs("H2RoundSiPMSolid", 0., outerRadius, sensorHalfLength, 0., 360.*deg);
  auto* mirrorSolid = new G4Tubs("H2RoundMirrorSolid", 0., outerRadius, mirrorHalfLength, 0., 360.*deg);
  auto* grease = new G4LogicalVolume(greaseSolid, fGrease, "H2RoundGreaseLV");
  auto* window = new G4LogicalVolume(windowSolid, fSiPMWindow, "H2RoundSiPMWindowLV");
  auto* sensor = new G4LogicalVolume(sensorSolid, fSilicon, "H2RoundSiPMLV");
  auto* mirror = new G4LogicalVolume(mirrorSolid, fSilicon, "H2RoundMirrorLV");
  new G4LogicalSkinSurface("H2RoundMirrorSkin", mirror, fMirrorSurface);

  auto* coreVis = new G4VisAttributes(G4Colour(.1, 1., .35));
  coreVis->SetForceSolid(true); core->SetVisAttributes(coreVis);
  outer->SetVisAttributes(new G4VisAttributes(G4Colour(.7, .9, 1., .18)));
  grease->SetVisAttributes(new G4VisAttributes(G4Colour(1., 1., .3, .35)));
  window->SetVisAttributes(new G4VisAttributes(G4Colour(.8, .8, .8, .25)));
  sensor->SetVisAttributes(new G4VisAttributes(G4Colour(.15, .15, .15)));

  constexpr G4double rowDz = 0.8660254*mm;
  for (G4int axis = 0; axis < 2; ++axis) {
    const G4int plane = 2 + axis;
    auto* rotation = new G4RotationMatrix;
    if (axis == 0) rotation->rotateY(90.*deg);
    else rotation->rotateX(90.*deg);
    if (GetPlaneYaw(plane) != 0.) {
      G4RotationMatrix inverseYaw;
      inverseYaw.rotateZ(-GetPlaneYaw(plane));
      *rotation = (*rotation) * inverseYaw;
    }
    const G4ThreeVector nominalUnit =
        axis == 0 ? G4ThreeVector(1.,0.,0.) : G4ThreeVector(0.,1.,0.);
    const auto unit = PlaneNominalToWorldDirection(plane, nominalUnit);
    for (G4int channel = 0; channel < 32; ++channel) {
      const G4int row = channel/16;
      const G4int index = channel%16;
      const G4double transverse = (index-7.5)*mm + row*.5*mm;
      const G4double localZ = z - row*rowDz - axis*2.2*mm;
      const G4ThreeVector nominalCentre =
          axis == 0 ? G4ThreeVector(0., transverse, localZ)
                    : G4ThreeVector(transverse, 0., localZ);
      const auto centre = PlaneNominalToWorldPosition(plane, nominalCentre);
      const G4int localCopy = axis*32 + channel;
      const G4int globalCopy = 64 + localCopy;
      new G4PVPlacement(rotation, centre, outer, axis == 0 ? "H2_X" : "H2_Y",
                        world, false, localCopy, true);
      if (!fFanoutPreview) {
        new G4PVPlacement(rotation, centre + unit*(halfLength+greaseHalfLength), grease,
                          "H2_Grease", world, false, localCopy, true);
        auto* windowPV = new G4PVPlacement(rotation,
                          centre + unit*(halfLength+2.*greaseHalfLength+windowHalfLength),
                          window, "SiPMWindow", world, false, localCopy, true);
        auto* sensorPV = new G4PVPlacement(rotation,
                          centre + unit*(halfLength+2.*greaseHalfLength+2.*windowHalfLength+sensorHalfLength),
                          sensor, "SiPM", world, false, globalCopy, true);
        new G4LogicalBorderSurface("H2RoundGreaseToSiPM_"+std::to_string(localCopy),
                                   windowPV, sensorPV, fSiPMSurface);
      }
      if (fFarEndReflectivity > 0.)
        new G4PVPlacement(rotation, centre - unit*(halfLength+mirrorHalfLength), mirror,
                          "H2_FarMirror", world, false, localCopy, true);
    }
  }
}

void DetectorConstruction::BuildSquareHodoscope(G4LogicalVolume* world, G4double z) {
  constexpr G4double halfLength = 82.5*mm;
  constexpr G4double outerHalf = 0.4995*mm;
  constexpr G4double innerHalf = 0.4845*mm;
  constexpr G4double coreHalf = 0.4695*mm;
  constexpr G4double greaseHalfLength = 0.05*mm;
  constexpr G4double windowHalfLength = 0.05*mm;
  constexpr G4double sensorHalfLength = 0.005*mm;
  constexpr G4double mirrorHalfLength = 0.01*mm;
  const G4bool singleClad = GetH2SingleClad();

  auto* outerSolid = new G4Box("H2OuterCladSolid", outerHalf, outerHalf, halfLength);
  auto* innerSolid = new G4Box("H2InnerCladSolid", innerHalf, innerHalf, halfLength);
  auto* coreSolid = new G4Box("H2CoreSolid", coreHalf, coreHalf, halfLength);
  auto* outer = new G4LogicalVolume(
      outerSolid, singleClad ? fInnerCladding : fOuterCladding, "H2OuterCladLV");
  auto* inner = new G4LogicalVolume(innerSolid, fInnerCladding, "H2InnerCladLV");
  auto* core = new G4LogicalVolume(coreSolid, fBCF20, "H2CoreLV");
  ApplyFibreUserLimits(outer);
  ApplyFibreUserLimits(inner);
  ApplyFibreUserLimits(core);
  if (singleClad) {
    new G4PVPlacement(nullptr, {}, core, "H2Core", outer, false, 0, false);
  } else {
    new G4PVPlacement(nullptr, {}, inner, "H2InnerClad", outer, false, 0, false);
    new G4PVPlacement(nullptr, {}, core, "H2Core", inner, false, 0, false);
  }
  new G4LogicalSkinSurface("H2CoreSurface", core, fFiberSurface);

  auto* greaseSolid = new G4Box("H2GreaseSolid", outerHalf, outerHalf, greaseHalfLength);
  auto* windowSolid = new G4Box("H2SiPMWindowSolid", outerHalf, outerHalf, windowHalfLength);
  auto* sensorSolid = new G4Box("H2SiPMSolid", outerHalf, outerHalf, sensorHalfLength);
  auto* mirrorSolid = new G4Box("H2MirrorSolid", outerHalf, outerHalf, mirrorHalfLength);
  auto* grease = new G4LogicalVolume(greaseSolid, fGrease, "H2GreaseLV");
  auto* window = new G4LogicalVolume(windowSolid, fSiPMWindow, "H2SiPMWindowLV");
  auto* sensor = new G4LogicalVolume(sensorSolid, fSilicon, "H2SiPMLV");
  auto* mirror = new G4LogicalVolume(mirrorSolid, fSilicon, "H2MirrorLV");
  new G4LogicalSkinSurface("H2MirrorSkin", mirror, fMirrorSurface);

  core->SetVisAttributes(new G4VisAttributes(G4Colour(.15, .55, 1.)));
  outer->SetVisAttributes(new G4VisAttributes(G4Colour(.7, .9, 1., .18)));
  sensor->SetVisAttributes(new G4VisAttributes(G4Colour(.15, .15, .15)));

  for (G4int axis = 0; axis < 2; ++axis) {
    const G4int plane = 2 + axis;
    auto* rotation = new G4RotationMatrix;
    if (axis == 0) rotation->rotateY(90.*deg);
    else rotation->rotateX(90.*deg);
    if (GetPlaneYaw(plane) != 0.) {
      G4RotationMatrix inverseYaw;
      inverseYaw.rotateZ(-GetPlaneYaw(plane));
      *rotation = (*rotation) * inverseYaw;
    }
    const G4ThreeVector nominalUnit =
        axis == 0 ? G4ThreeVector(1.,0.,0.) : G4ThreeVector(0.,1.,0.);
    const auto unit = PlaneNominalToWorldDirection(plane, nominalUnit);
    for (G4int channel = 0; channel < 32; ++channel) {
      const G4double transverse = (channel-15.5)*mm;
      const G4double planeZ = z - axis*1.*mm;
      const G4ThreeVector nominalCentre =
          axis == 0 ? G4ThreeVector(0., transverse, planeZ)
                    : G4ThreeVector(transverse, 0., planeZ);
      const auto centre = PlaneNominalToWorldPosition(plane, nominalCentre);
      const G4int localCopy = axis*32 + channel;
      const G4int globalCopy = 64 + localCopy;
      new G4PVPlacement(rotation, centre, outer, axis == 0 ? "H2_X" : "H2_Y",
                        world, false, localCopy, true);
      if (!fFanoutPreview) {
      new G4PVPlacement(rotation, centre + unit*(halfLength+greaseHalfLength), grease,
                        "H2_Grease", world, false, localCopy, true);
      auto* windowPV = new G4PVPlacement(rotation,
                        centre + unit*(halfLength+2.*greaseHalfLength+windowHalfLength),
                        window, "SiPMWindow", world, false, localCopy, true);
      auto* sensorPV = new G4PVPlacement(rotation,
                        centre + unit*(halfLength+2.*greaseHalfLength+2.*windowHalfLength+sensorHalfLength),
                        sensor, "SiPM", world, false, globalCopy, true);
      new G4LogicalBorderSurface("H2GreaseToSiPM_"+std::to_string(localCopy),
                                 windowPV, sensorPV, fSiPMSurface);
      }
      if (fFarEndReflectivity > 0.)
        new G4PVPlacement(rotation, centre - unit*(halfLength+mirrorHalfLength), mirror,
                          "H2_FarMirror", world, false, localCopy, true);
    }
  }
}
