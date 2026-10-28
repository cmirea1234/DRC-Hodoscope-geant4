#pragma once

#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "globals.hh"
#include <array>
#include <map>

class G4GenericMessenger;
class G4LogicalVolume;
class G4Material;
class G4OpticalSurface;
class G4UserLimits;
class G4VPhysicalVolume;

class DetectorConstruction final : public G4VUserDetectorConstruction {
 public:
  // The Stage-2 source table is a frozen validation fixture.  Keep the row
  // contract available to the generator and ROOT provenance writer without
  // duplicating a magic number in those two consumers.
  static constexpr G4int kGeometryProbeExpectedRows = 533;
  // Stage-3 uses the same event-indexed, fail-closed source-table pattern as
  // Stage-2, but with charged 4 GeV muon probes and a separately frozen row
  // contract.
  static constexpr G4int kGenerationProbeExpectedRows = 196960;

  DetectorConstruction();
  ~DetectorConstruction() override;
  G4VPhysicalVolume* Construct() override;

  G4double GetHitThreshold() const { return fHitThreshold; }
  G4bool GetLightStraightFibres() const { return fLightStraightFibres; }
  G4bool GetLightStraightReadout() const { return fLightStraightReadout; }
  G4double GetGreaseGap() const { return fGreaseGap; }
  G4bool GetCosmicBeam() const { return fCosmicBeam; }
  G4long GetCosmicSeed() const { return fCosmicSeed; }
  G4int GetReplayEventIdOffset() const { return fReplayEventIdOffset; }
  G4int GetTransportSeedEventOffset() const {
    return fTransportSeedEventOffset;
  }
  G4bool GetFixedPhotonTest() const { return fFixedPhotonTest; }
  const G4String& GetFixedPhotonPlane() const { return fFixedPhotonPlane; }
  const G4String& GetFixedPhotonAngularMode() const {
    return fFixedPhotonAngularMode;
  }
  const G4String& GetFixedPhotonSpectrum() const { return fFixedPhotonSpectrum; }
  const G4String& GetFixedPhotonSource() const { return fFixedPhotonSource; }
  G4int GetFixedPhotonChannel() const { return fFixedPhotonChannel; }
  G4int GetFixedPhotonChannelForEvent(G4int eventId) const {
    const G4int plane=fFixedPhotonPlane=="h1x"?0:fFixedPhotonPlane=="h1y"?1:
                      fFixedPhotonPlane=="h2x"?2:3;
    return fFixedPhotonChannel < 0 ? eventId%GetPlaneFibreCount(plane) : fFixedPhotonChannel;
  }
  G4double GetFixedPhotonWavelength() const { return fFixedPhotonWavelength; }
  G4double GetFixedPhotonPolarAngle() const { return fFixedPhotonPolarAngle; }
  G4double GetFixedPhotonAzimuth() const { return fFixedPhotonAzimuth; }
  G4double GetFixedPhotonLongitudinalOffset() const {
    return fFixedPhotonLongitudinalOffset;
  }
  G4bool GetGate4OpticalAudit() const { return fGate4OpticalAudit; }
  G4bool GetGate4BoundaryTrace() const { return fGate4BoundaryTrace; }
  G4bool GetGate4SpectrumOnlyKill() const { return fGate4SpectrumOnlyKill; }
  G4bool GetRecordSpatialClosureNtuples() const {
    return fRecordSpatialClosureNtuples;
  }
  G4bool GetPathOnlyMuonTruth() const {
    return fPathOnlyMuonTruth;
  }
  G4bool GetRecordGeometryControlNtuples() const {
    return fRecordGeometryControlNtuples;
  }
  G4bool GetGeometryProbeMode() const { return fGeometryProbeMode; }
  const G4String& GetGeometryProbeFile() const { return fGeometryProbeFile; }
  G4int GetGeometryProbeRowOffset() const { return fGeometryProbeRowOffset; }
  G4bool GetGenerationProbeMode() const { return fGenerationProbeMode; }
  const G4String& GetGenerationProbeFile() const {
    return fGenerationProbeFile;
  }
  const G4String& GetGenerationProbeSha256() const {
    return fGenerationProbeSha256;
  }
  G4int GetGenerationProbeRowOffset() const {
    return fGenerationProbeRowOffset;
  }
  G4bool GetRecordGenerationNtuples() const {
    return fRecordGenerationNtuples;
  }
  G4double GetFibreMaxStep() const { return fFibreMaxStep; }
  G4bool GetGenerationOnlyKill() const { return fGenerationOnlyKill; }
  G4double GetOpticalPhotonMaxTime() const { return fOpticalPhotonMaxTime; }
  G4bool GetSkipOpticalUntilT1T2() const { return fSkipOpticalUntilT1T2; }
  G4bool GetGate6ElectronicsProbeMode() const {
    return fGate6ElectronicsProbeMode;
  }
  G4int GetGate6InjectedChannelForEvent(G4int eventId) const {
    return fGate6InjectedChannel < 0 ? eventId%128 : fGate6InjectedChannel;
  }
  G4int GetGate6InjectedPe() const { return fGate6InjectedPe; }
  G4double GetGate6InjectedStartTime() const {
    return fGate6InjectedStartTime;
  }
  G4double GetGate6InjectedSpacing() const {
    return fGate6InjectedSpacing;
  }
  G4ThreeVector GetFanoutEndpoint(G4int channel) const {
    return channel >= 0 && channel < 128 ? fFanoutEndpoint[channel]
                                         : G4ThreeVector();
  }
  G4ThreeVector GetFanoutTerminalDirection(G4int channel) const {
    return channel >= 0 && channel < 128 ? fFanoutTerminalDirection[channel]
                                         : G4ThreeVector();
  }
  G4bool HasFanoutEndpoint(G4int channel) const {
    return channel >= 0 && channel < 128 && fFanoutEndpointReady[channel];
  }
  G4double GetMaxZenith() const { return fMaxZenith; }
  G4double GetCosmicTargetHalfX() const { return fCosmicTargetHalfX; }
  G4double GetCosmicTargetHalfY() const { return fCosmicTargetHalfY; }
  const G4String& GetCosmicPositionModel() const { return fCosmicPositionModel; }
  G4double GetCosmicSourceHalfX() const { return fCosmicSourceHalfX; }
  G4double GetCosmicSourceHalfY() const { return fCosmicSourceHalfY; }
  G4int GetUnselectedEventStride() const { return fUnselectedEventStride; }
  G4bool GetSkipUntriggeredDigitization() const { return fSkipUntriggeredDigitization; }
  const G4String& GetCosmicEnergySpectrum() const { return fCosmicEnergySpectrum; }
  const G4String& GetCosmicAngularDistribution() const {
    return fCosmicAngularDistribution;
  }
  G4double GetCosmicZenithPower() const { return fCosmicZenithPower; }
  G4double GetCosmicMinEnergy() const { return fCosmicMinEnergy; }
  G4double GetCosmicMaxEnergy() const { return fCosmicMaxEnergy; }
  G4int GetTargetCoincidences() const { return fTargetCoincidences; }
  const G4String& GetTriggerCoincidenceMode() const {
    return fTriggerCoincidenceMode;
  }
  G4double GetH1AdcThreshold() const { return fH1AdcThreshold; }
  G4double GetH2AdcThreshold() const { return fH2AdcThreshold; }
  G4double GetCrossTalkProbability() const { return fCrossTalkProbability; }
  G4double GetAfterPulseProbability() const { return fAfterPulseProbability; }
  G4double GetDarkCountRate() const { return fDarkCountRate; }
  G4double GetReadoutGate() const { return fReadoutGate; }
  G4double GetElectronicsNoise() const { return fElectronicsNoise; }
  G4bool GetUseMeasuredNoise() const { return fUseMeasuredNoise; }
  G4double GetAdcNoiseScale() const { return fAdcNoiseScale; }
  G4int GetPedestalWindowSamples() const { return fPedestalWindowSamples; }
  G4int GetPeakWindowSamples() const { return fPeakWindowSamples; }
  G4bool GetUseWaveformDigitization() const { return fUseWaveformDigitization; }
  G4double GetAdcSamplePeriod() const { return fAdcSamplePeriod; }
  G4double GetSinglePeFwhmSamples() const { return fSinglePeFwhmSamples; }
  G4double GetAdcPeakWindowDelay() const { return fAdcPeakWindowDelay; }
  G4double GetXYCoincidenceWindow() const { return fXYCoincidenceWindow; }
  G4bool GetFanoutPreview() const { return fFanoutPreview; }
  G4double GetFibreFarEndExtension() const { return fFibreFarEndExtension; }
  G4bool GetActivePlaneOverlap() const { return fActivePlaneOverlap; }
  std::array<G4double,2> GetActiveLongitudinalBounds(G4int plane) const;
  G4double GetFanoutOuterSize(G4int plane) const {
    return .4995*mm-((plane<2||fH2PhysicalRound)
        ? std::min(fFibreOverlapClearance,.01498*mm):fFibreOverlapClearance);
  }
  const G4String& GetSinglePeTemplateFile(G4int group) const {
    if(group<0||group>3)return fSinglePeTemplateFile[0];
    return fSinglePeTemplateFile[group];
  }
  G4bool GetRecordPhotonNtuple() const { return fRecordPhotonNtuple; }
  G4bool GetRecordChannelNtuple() const { return fRecordChannelNtuple; }
  G4bool GetRecordSelectedDiagnosticsOnly() const { return fRecordSelectedDiagnosticsOnly; }
  G4bool GetRecordPhotonAuditNtuple() const {
    return fRecordPhotonAuditNtuple;
  }
  G4bool GetRecordPhotonAuditRows() const {
    return fRecordPhotonAuditRows;
  }
  G4bool GetRecordActiveCoreChordBranches() const {
    return fRecordActiveCoreChordBranches;
  }
  G4int GetMicrocellCount() const {
    if (fMicrocellCount > 0) return fMicrocellCount;
    return 16663;
  }
  G4double GetBiasVoltage() const { return fBiasVoltage; }
  G4bool GetUseBuilding() const { return fUseBuilding; }
  G4double GetH1AttenuationLength() const { return fH1AttenuationLength; }
  G4double GetH2AttenuationLength() const { return fH2AttenuationLength; }
  G4double GetH1ScintillationYield() const { return fH1ScintillationYield; }
  G4double GetH2ScintillationYield() const { return fH2ScintillationYield; }
  const G4Material* GetH1CoreMaterial() const { return fSCSF78; }
  const G4Material* GetH2CoreMaterial() const { return fBCF20; }
  G4bool GetH1SingleClad() const { return fH1CladdingModel == "single"; }
  G4bool GetH2SingleClad() const { return fH2CladdingModel == "single"; }
  G4bool GetH1Uses22PointEmission() const { return fH1EmissionModel == "scsf78_22"; }
  G4int GetH1RoundSides() const { return fH1RoundSides; }
  G4int GetH1ActiveFacetRefinement() const {
    return fH1ActiveFacetRefinement;
  }
  G4double GetH1RoundPolygonPhase() const;
  const G4String& GetH1RoundSweepRepresentation() const {
    return fH1RoundSweepRepresentation;
  }
  G4int GetH1RoundSweepRepresentationCode() const {
    if (fH1RoundSweepRepresentation == "tessellated") return 0;
    if (fH1RoundSweepRepresentation == "native_miter_union") return 1;
    if (fH1RoundSweepRepresentation == "native_capsule_union") return 2;
    if (fH1RoundSweepRepresentation ==
        "native_miter_subtraction_union") return 3;
    if (fH1RoundSweepRepresentation ==
        "native_miter_public_union") return 4;
    if (fH1RoundSweepRepresentation ==
        "native_balanced_union") return 5;
    return -1;
  }
  G4double GetH1RoundMiterJointOverlap() const {
    return fH1RoundMiterJointOverlap;
  }
  const G4String& GetH1FanoutHierarchy() const {
    return fH1FanoutHierarchy;
  }
  G4int GetH1FanoutHierarchyCode() const {
    return fH1FanoutHierarchy == "sibling_shells" ? 1 : 0;
  }
  const G4String& GetH2FanoutHierarchy() const {
    return fH2FanoutHierarchy;
  }
  G4int GetH2FanoutHierarchyCode() const {
    return fH2FanoutHierarchy == "sibling_shells" ? 1 : 0;
  }
  const G4String& GetPlaneLayout(G4int plane) const {
    return plane >= 0 && plane < 4 ? fPlaneLayout[plane] : fPlaneLayout[0];
  }
  G4int GetPlaneLayoutCode(G4int plane) const {
    const auto& layout = GetPlaneLayout(plane);
    if (layout == "single1x32") return 0;
    if (layout == "2x16_aligned") return 1;
    if (layout == "2x16_staggered") return 2;
    return -1;
  }
  const G4String& GetPlaneChannelOrder(G4int plane) const {
    return plane >= 0 && plane < 4
        ? fPlaneChannelOrder[plane] : fPlaneChannelOrder[0];
  }
  G4int GetPlaneChannelOrderCode(G4int plane) const {
    const auto& order = GetPlaneChannelOrder(plane);
    if (order == "ascending") return 0;
    if (order == "serpentine") return 1;
    if (order == "descending") return 2;
    return -1;
  }
  G4double GetPlaneRowZSpacing(G4int plane) const {
    return plane >= 0 && plane < 4 ? fPlaneRowZSpacing[plane] : 0.;
  }
  G4double GetPlaneSourceTransverse(G4int plane, G4int channel) const;
  G4double GetPlaneSourceZ(G4int plane, G4int channel) const;
  G4double GetValidationSurfaceZ(G4int surface) const;
  G4double GetTriggerZ(G4int trigger) const {
    return trigger >= 0 && trigger < 4 ? fTriggerZ[trigger] : 0.;
  }
  G4bool GetTriggerEnabled(G4int trigger) const {
    return trigger >= 0 && trigger < 4 ? fTriggerEnabled[trigger] : false;
  }
  G4double GetHodoscopeFanoutZ(G4int module) const {
    return module >= 0 && module < 2 ? fHodoscopeFanoutZ[module] : 0.;
  }
  G4bool GetVirtualRoundModelEnabled() const {
    return fVirtualRoundModelEnabled;
  }
  G4bool GetVirtualRoundH3Enabled() const {
    return fVirtualRoundH3Enabled;
  }
  G4bool GetH2PhysicalRound() const { return fH2PhysicalRound; }
  G4int GetH2SquareFibreCount() const { return fH2SquareFibreCount; }
  G4int GetPlaneFibreCount(G4int plane) const {
    return plane>=2 && plane<4 && !fH2PhysicalRound ? fH2SquareFibreCount : 32;
  }
  G4bool HasPhysicalFibre(G4int plane,G4int channel) const {
    return plane>=0 && plane<4 && channel>=0 && channel<GetPlaneFibreCount(plane);
  }
  G4double GetVirtualRoundH3Z() const { return fVirtualRoundH3Z; }
  G4double GetVirtualRoundH3OffsetX() const {
    return fVirtualRoundH3OffsetX;
  }
  G4double GetVirtualRoundH3OffsetY() const {
    return fVirtualRoundH3OffsetY;
  }
  G4double GetVirtualRoundH3Yaw() const {
    return fVirtualRoundH3Yaw;
  }
  G4int GetPlaneFanoutRow(G4int plane, G4int channel) const;
  G4int GetPlaneFanoutColumn(G4int plane, G4int channel) const;
  G4bool GetMixedLayoutPattern() const;
  void SetLayoutPreset(G4String preset);
  const G4String& GetDiagnosticRootFile() const { return fDiagnosticRootFile; }
  G4double GetPDEScale() const { return fPDEScale; }
  G4double GetCouplingEfficiency() const { return fCouplingEfficiency; }
  G4double GetBendTransmission() const { return fBendTransmission; }
  G4double GetEndFaceEfficiency() const { return fEndFaceEfficiency; }
  G4double GetFanoutChannelTransmission(G4int channel) const {
    return channel>=0&&channel<128 ? fFanoutChannelTransmission[channel] : 1.;
  }
  G4double GetFanoutMinRadius(G4int channel) const {
    return channel>=0&&channel<128 ? fFanoutMinRadius[channel] : 0.;
  }
  G4double GetScintillationYieldScale() const { return fScintillationYieldScale; }
  G4double GetEffectiveLightScaleH1() const { return fEffectiveLightScaleH1; }
  G4double GetEffectiveLightScaleH2() const { return fEffectiveLightScaleH2; }
  G4double GetH1BirksConstant() const { return fH1BirksConstant; }
  G4double GetH2BirksConstant() const { return fH2BirksConstant; }
  G4double GetTriggerBirksConstant() const { return fTriggerBirksConstant; }
  G4double GetTriggerScintillationYield() const { return fTriggerScintillationYield; }
  G4double GetStressLossPerM() const { return fStressLossPerM/m; }
  G4double GetCrackLossPerM() const { return fCrackLossPerM/m; }
  G4double GetPackingLossPerM() const { return fPackingLossPerM/m; }
  G4double GetAdhesiveLossPerM() const { return fAdhesiveLossPerM/m; }
  G4double GetH1CladdingDamagePerM() const { return fH1CladdingDamagePerM/m; }
  G4double GetH2CladdingDamagePerM() const { return fH2CladdingDamagePerM/m; }
  G4double GetAdhesiveIntrusionLength() const { return fAdhesiveIntrusionLength; }
  G4double GetTriggerYaw(G4int trigger) const {
    return trigger >= 0 && trigger < 4 ? fTriggerYaw[trigger] : 0.;
  }
  G4double GetTriggerHalfLength() const { return 0.5*fTriggerLength; }
  G4double GetTriggerHalfWidth() const { return 0.5*fTriggerWidth; }
  G4double GetTriggerHalfThickness() const { return 0.5*fTriggerThickness; }
  G4double GetTriggerHalfLength(G4int trigger) const {
    if(trigger==2)return 0.5*fTrigger3Length;
    if(trigger==3)return 0.5*fTrigger4Length;
    return GetTriggerHalfLength();
  }
  G4double GetTriggerHalfWidth(G4int trigger) const {
    if(trigger==2)return 0.5*fTrigger3Width;
    if(trigger==3)return 0.5*fTrigger4Width;
    return GetTriggerHalfWidth();
  }
  G4double GetTriggerHalfThickness(G4int trigger) const {
    if(trigger==2)return 0.5*fTrigger3Thickness;
    if(trigger==3)return 0.5*fTrigger4Thickness;
    return GetTriggerHalfThickness();
  }
  G4double GetTriggerOffsetX(G4int trigger) const {
    return trigger >= 0 && trigger < 4 ? fTriggerOffsetX[trigger] : 0.;
  }
  G4double GetTriggerOffsetY(G4int trigger) const {
    return trigger >= 0 && trigger < 4 ? fTriggerOffsetY[trigger] : 0.;
  }
  const G4String& GetAlignmentReference() const {
    return fAlignmentReference;
  }
  G4int GetAlignmentReferenceCode() const {
    return fAlignmentReference == "ideal_axis" ? 1 : 0;
  }
  G4bool GetRecordAlignmentStudyNtuple() const {
    return fRecordAlignmentStudyNtuple;
  }
  G4ThreeVector GetTriggerCentre(G4int trigger) const;
  G4double GetTriggerRotationAngle(G4int trigger) const;
  G4ThreeVector TriggerWorldToLocalPosition(
      G4int trigger,const G4ThreeVector& position) const;
  G4ThreeVector TriggerWorldToLocalDirection(
      G4int trigger,const G4ThreeVector& direction) const;
  G4double GetPlaneYaw(G4int plane) const {
    return plane >= 0 && plane < 4 ? fPlaneYaw[plane] : 0.;
  }
  G4bool GetH2ModuleRotate180() const { return fH2ModuleRotate180; }
  G4double GetHodoscopeOffsetX(G4int module) const {
    return module >= 0 && module < 2 ? fHodoscopeOffsetX[module] : 0.;
  }
  G4double GetHodoscopeOffsetY(G4int module) const {
    return module >= 0 && module < 2 ? fHodoscopeOffsetY[module] : 0.;
  }
  // V15 alignment helpers map between the frozen V14 nominal plane frame and
  // the yawed world frame.  A zero yaw returns the input without arithmetic so
  // the nominal random stream and geometry coordinates remain bit-identical.
  G4ThreeVector PlaneNominalToWorldPosition(
      G4int plane, const G4ThreeVector& position) const;
  G4ThreeVector PlaneWorldToNominalPosition(
      G4int plane, const G4ThreeVector& position) const;
  G4ThreeVector PlaneNominalToWorldDirection(
      G4int plane, const G4ThreeVector& direction) const;
  G4ThreeVector PlaneWorldToNominalDirection(
      G4int plane, const G4ThreeVector& direction) const;
  G4ThreeVector PlaneNominalPivot(G4int plane) const;
  G4int GetGainBoard(G4int plane) const {
    if (plane < 0 || plane >= 4) return plane + 1;
    return fGainBoard[plane];
  }

 private:
  void DefineMaterials();
  void DefineSurfaces();
  void BuildTrigger(G4LogicalVolume*, G4int);
  void BuildRoundHodoscope(G4LogicalVolume*, G4double);
  void BuildRoundHodoscopeH2(G4LogicalVolume*, G4double);
  void BuildSquareHodoscope(G4LogicalVolume*, G4double);
  void BuildFanoutPreview(G4LogicalVolume*);
  void BuildLightStraightFibres(G4LogicalVolume*);
  void BuildFanoutPreviewLite(G4LogicalVolume*);
  void PlaceUpstreamMaterials(G4LogicalVolume*);
  void RunNavigatorSolidProbe(G4String);
  void ApplyFibreUserLimits(G4LogicalVolume*) const;

  G4double fHitThreshold = 0.05*MeV;
  G4bool fUseBuilding = false;
  G4double fConcreteThickness = 15.*cm;
  G4double fBlackoutThickness = 0.2*mm;
  G4bool fBlackoutCoverAllFibres = false;
  G4bool fCosmicBeam = false;
  G4long fCosmicSeed = 310713;
  G4int fReplayEventIdOffset = 0;
  // Diagnostic replay control. The production value is exactly zero.
  // A positive offset lets a one-event macro reproduce the transport RNG of
  // a later event without simulating all preceding events.
  G4int fTransportSeedEventOffset = 0;
  G4bool fFixedPhotonTest = false;
  G4String fFixedPhotonPlane = "h1x";
  G4String fFixedPhotonAngularMode = "isotropic";
  G4String fFixedPhotonSpectrum = "mono";
  // Gate-1 diagnostic source location. "active" preserves the predecessor
  // fixture; "terminal" injects 1 um before the modelled fibre end and tests
  // the configured fibre-end -> grease -> window -> SiPM channel ownership.
  G4String fFixedPhotonSource = "active";
  G4int fFixedPhotonChannel = 8;
  G4double fFixedPhotonWavelength = 450.*nm;
  // Gate-4 fixed-direction and straight-fibre source controls. Defaults keep
  // every predecessor fixed-photon run unchanged.
  G4double fFixedPhotonPolarAngle = 0.*deg;
  G4double fFixedPhotonAzimuth = 0.*deg;
  G4double fFixedPhotonLongitudinalOffset = 0.*mm;
  G4bool fGate4OpticalAudit = false;
  G4bool fGate4BoundaryTrace = false;
  G4bool fGate4SpectrumOnlyKill = false;
  G4String fCosmicEnergySpectrum = "modified_gaisser";
  // Sea-level muon zenith model.  The default preserves the predecessor
  // hard-coded p(cos(theta)) proportional to cos^2(theta).  Use
  // flat_plane_cos2 for horizontal-generation-plane sampling of an
  // I(theta) proportional to cos^2(theta) flux, i.e. p(cos(theta)) ∝ cos^3.
  // modified_gaisser_joint samples E and cos(theta) together from
  // dI/dE(E,theta)*cos(theta) on the horizontal generation plane.
  G4String fCosmicAngularDistribution = "cos2";
  G4double fCosmicZenithPower = 2.;
  G4double fCosmicMinEnergy = 0.2*GeV;
  G4double fCosmicMaxEnergy = 1000.*GeV;
  // Isolated light-model candidate; never enabled in the preserved original.
  G4bool fLightDirectShell = false;
  G4int fLightMaxVoxels = 4000;
  G4int fLightRoundSegments = 47;
  G4bool fLightStraightFanout = false;
  G4bool fLightStraightFibres = false;
  G4bool fLightStraightReadout = false;
  G4bool fFanoutPreview = false;
  // Explicit opt-in preserves historical macros; revised setup uses 40 mm.
  G4double fFibreFarEndExtension = 0.*mm;
  G4bool fActivePlaneOverlap = false;
  // Keep the original fan-out geometry but use linear interpolation from
  // detector-side packed fibres to the SiPM-side grid.  Default false
  // preserves the validated smootherstep path.
  G4bool fFanoutLinearFanout = false;
  // Visualization-only fan-out preview. It keeps the same detector/trigger
  // coordinates but draws only representative, simple solids so remote GUI
  // checks remain responsive. It is not a physics-production geometry.
  G4bool fFanoutPreviewLite = false;
  G4int fFanoutPreviewLiteChannelStride = 8;
  G4int fFanoutPreviewLiteSegments = 8;
  // Optional geometry-construction clearance for dense fan-out overlap
  // diagnostics. Default zero preserves validated predecessor geometry.
  // Final-validation geometry uses a tiny envelope clearance to prevent
  // Geant4 overlap/navigation warnings at packed active/fan-out seams.  This
  // is a geometric representation tolerance, not an empirical optical loss.
  // H1 round fibres cap this value internally so the single-clad shell remains
  // navigable; H2 square fibres use the full value for deterministic clearance.
  G4double fFibreOverlapClearance = 0.015*mm;
  G4double fMaxZenith = 45.*deg;
  G4double fCosmicTargetHalfX = 16.*mm;
  G4double fCosmicTargetHalfY = 16.*mm;
  G4String fCosmicPositionModel = "target_backproject";
  G4double fCosmicSourceHalfX = 2200.*mm;
  G4double fCosmicSourceHalfY = 2200.*mm;
  G4int fUnselectedEventStride = 1;
  G4bool fSkipUntriggeredDigitization = false;
  G4int fTargetCoincidences = 0;
  G4String fTriggerCoincidenceMode = "t1t2";
  // Trigger scintillator local full dimensions. G4Box placement converts
  // these to half extents. The legacy trigger was 30 x 10 x 5 mm full size;
  // geometry-acceptance studies can set 60 x 10 x 5 mm explicitly.
  G4double fTriggerLength = 30.*mm;
  G4double fTriggerWidth = 10.*mm;
  G4double fTriggerThickness = 5.*mm;
  G4double fTrigger3Length = 30.*mm;
  G4double fTrigger3Width = 10.*mm;
  G4double fTrigger3Thickness = 5.*mm;
  G4double fTrigger4Length = 30.*mm;
  G4double fTrigger4Width = 10.*mm;
  G4double fTrigger4Thickness = 5.*mm;
  G4double fH1AdcThreshold = 120.;
  G4double fH2AdcThreshold = 120.;
  G4double fFiberSurfaceRoughness = 0.;
  G4double fEndFaceRoughness = 0.;
  G4double fFarEndReflectivity = 0.;
  // V15 retains the V14 fibre -> grease -> SiPM-window interfaces explicitly.
  // Keep this residual multiplier at unity unless an independently measured
  // non-geometrical contact loss is supplied.
  G4double fCouplingEfficiency = 1.0;
  G4double fGreaseGap = 0.10*mm;
  G4double fSiPMActiveSize = 1.30*mm;
  G4double fSiPMOffsetU = 0.*mm;
  G4double fSiPMOffsetV = 0.*mm;
  G4double fBendTransmission = 1.0;
  G4double fFanoutReferenceTransmissionH1 = 1.0;
  G4double fFanoutReferenceTransmissionH2 = 1.0;
  G4double fFanoutReferenceRadius = 30.*mm;
  G4double fFanoutCurvaturePower = 1.0;
  G4double fFanoutChannelRms = 0.;
  G4double fPDEScale = 1.0;
  G4double fScintillationYieldScale = 1.0;
  // Post-transport/PDE effective light scales.  These are not assigned to a
  // specific aging/fan-out/coupling mechanism; they represent unresolved
  // assembled-detector light budget only when explicitly scanned.
  G4double fEffectiveLightScaleH1 = 1.0;
  G4double fEffectiveLightScaleH2 = 1.0;
  // Generic polystyrene reference, not a fibre-specific measurement.
  // Command values are numeric coefficients in mm/MeV because the Geant4 UI
  // unit table does not define this compound unit.
  G4double fH1BirksConstant = 0.126;
  G4double fH2BirksConstant = 0.126;
  G4double fTriggerBirksConstant = 0.126;
  G4double fTriggerScintillationYield = 10000.;
  G4double fMicrobendAmplitude = 0.*mm;
  G4double fMicrobendPeriod = 10.*mm;
  G4double fBendOvality = 0.;
  G4double fCladdingAsymmetry = 0.;
  G4double fStressLossPerM = 0.;
  G4double fCrackLossPerM = 0.;
  G4double fPackingLossPerM = 0.;
  G4double fAdhesiveIntrusionLength = 0.*mm;
  G4double fAdhesiveLossPerM = 0.;
  G4double fH1CladdingDamagePerM = 0.;
  G4double fH2CladdingDamagePerM = 0.;
  G4double fEndFaceEfficiency = 1.;
  G4double fCrossTalkProbability = 0.01;
  G4double fAfterPulseProbability = 0.;
  G4double fDarkCountRate = 120.e3/second;
  G4double fReadoutGate = 100.*ns;
  G4double fElectronicsNoise = 0.25;
  G4bool fUseMeasuredNoise = false;
  G4double fAdcNoiseScale = 1.;
  G4int fPedestalWindowSamples = 100;
  G4int fPeakWindowSamples = 400;
  G4bool fUseWaveformDigitization = true;
  G4double fAdcSamplePeriod = .2*ns;
  G4double fSinglePeFwhmSamples = 62.;
  G4double fAdcPeakWindowDelay = -20.*ns;
  G4double fXYCoincidenceWindow = 10.*ns;
  G4String fSinglePeTemplateFile[4] = {
    "single_pe_template_group1.csv", "single_pe_template_group2.csv",
    "single_pe_template_group3.csv", "single_pe_template_group4.csv"};
  G4bool fRecordPhotonNtuple = true;
  // Optional long-format electronics diagnostics (128 rows/event).  Keep off
  // for normal optical production to avoid unnecessarily large ROOT files.
  G4bool fRecordChannelNtuple = false;
  // Output-only filter: keep the complete event tree and physics/RNG unchanged.
  G4bool fRecordSelectedDiagnosticsOnly = false;
  // V15 diagnostic-only lineage record. It must not alter geometry, optical
  // properties, random-number draws, or the canonical 227-column event tree.
  G4bool fRecordPhotonAuditNtuple = false;
  // Keep the full lineage/closure calculation above active while optionally
  // suppressing only the large per-photon ROOT rows.  Default true preserves
  // every existing full-audit macro that enables fRecordPhotonAuditNtuple.
  G4bool fRecordPhotonAuditRows = true;
  // Exact V16 port: when enabled, append four primary-only chord columns
  // clipped to the nominal [-16,+16] mm active straight. Keep default false
  // so the V18-D compatibility regression retains the canonical 227-column
  // schema.
  G4bool fRecordActiveCoreChordBranches = false;
  // Passive Gate-1 relational diagnostics. This option only records state
  // already produced by transport and must never draw a random number or
  // mutate a track.
  G4bool fRecordSpatialClosureNtuples = false;
  // Compact structural-acceptance mode: record primary muon trigger/fibre
  // truth while intentionally suppressing photon/channel/audit rows.
  G4bool fPathOnlyMuonTruth = false;
  // Stage-2 geometry-control extension. It is legal only together with the
  // Gate-1 spatial rows and appends trigger-volume truth without changing
  // the established Gate-1 schema when disabled.
  G4bool fRecordGeometryControlNtuples = false;
  // Deterministic geantino source table used only by geometry controls.  The
  // normal primary generator remains the default and does not read this file.
  G4bool fGeometryProbeMode = false;
  G4String fGeometryProbeFile = "../data/validation/gate2_geometry_probe_cases.tsv";
  G4int fGeometryProbeRowOffset = 0;
  // Stage-3 controls are default-OFF and PreInit-only.  The expected digest
  // deliberately has no built-in value: an enabled validation macro must
  // provide the digest of its sealed source table explicitly.
  G4bool fGenerationProbeMode = false;
  G4String fGenerationProbeFile =
      "../data/validation/gate3_generation_probe_cases.tsv";
  G4String fGenerationProbeSha256 = "";
  G4int fGenerationProbeRowOffset = 0;
  G4bool fRecordGenerationNtuples = false;
  G4double fFibreMaxStep = 0.*mm;
  G4bool fGenerationOnlyKill = false;
  // Production safety valve for long-lived optical photons.  The default
  // value preserves all established validation macros; production macros may
  // set a positive value when late photons are outside the digitizer window.
  G4double fOpticalPhotonMaxTime = 0.*ns;
  // Production optimization for T1&&T2-denominator studies: optical photons
  // are stacked until charged-particle tracking is done, so non-triggered
  // events can skip optical transport without changing conditional results.
  G4bool fSkipOpticalUntilT1T2 = false;
  // Gate-6 default-OFF digitizer fixture. It injects exact primary p.e.
  // arrival times after event-state reset and never participates in normal
  // optical transport.
  G4bool fGate6ElectronicsProbeMode = false;
  G4int fGate6InjectedChannel = -1;
  G4int fGate6InjectedPe = 1;
  G4double fGate6InjectedStartTime = 0.*ns;
  G4double fGate6InjectedSpacing = 0.*ns;
  G4int fMicrocellCount = 0;
  G4double fBiasVoltage = 45.;
  G4String fSiPMModel = "S14160-1310PS";
  G4double fH1AttenuationLength = 4.67*m;
  G4double fH2AttenuationLength = 3.5*m;
  G4double fH1ScintillationYield = 14000.;
  G4double fH2ScintillationYield = 8000.;
  G4String fH1CladdingModel = "single";
  G4String fH2CladdingModel = "single";
  G4String fH1EmissionModel = "scsf78_22";
  // V21 production geometry. Compatibility aliases retain the historical H1
  // names, but the approved native balanced Boolean construction is now the
  // canonical process default.
  G4int fH1RoundSides = 64;
  G4int fH1ActiveFacetRefinement = 1;
  G4String fH1FanoutHierarchy = "sibling_shells";
  // Preserve the predecessor square-fibre construction unless a macro
  // explicitly selects the navigation-safe sibling-shell representation.
  // Production V21 macros record and require that explicit selection.
  G4String fH2FanoutHierarchy = "nested";
  // Independent predecessor geometry and short full-optical gates approved
  // this representation with refinement=1 and exact touching common planes.
  G4String fH1RoundSweepRepresentation = "native_balanced_union";
  // Non-zero overlap remains a sensitivity only, never an implicit physical
  // loss or a production default.
  G4double fH1RoundMiterJointOverlap = 0.*nm;
  // Fail-closed direct-solid replay controls.  A macro must provide an exact
  // condition label and row count; the frozen mode characterises the original
  // tessellated failure while candidate_gate requires a valid outside-to-
  // inside transition independent of the old captured surface point.
  G4String fNavigatorSolidProbeMode = "candidate_gate";
  G4String fNavigatorSolidProbeExpectedCondition = "";
  G4int fNavigatorSolidProbeExpectedRows = 0;
  // Named direct-solid probes are expressed in the historical common
  // world-like frame. Public balanced Boolean roots instead live in their
  // leftmost primitive frame, so preserve the exact rigid transform used by
  // their sole physical placement.
  std::map<G4String,G4Transform3D> fNavigatorSolidLocalToWorld;
  // Process default is the V18-D-compatible plane layout. The unconfirmed
  // mixed photo hypothesis is available only through layoutPreset or explicit
  // plane commands; it must never become an implicit construction fact.
  // Plane order: H1X, H1Y, H2X, H2Y.
  std::array<G4String,4> fPlaneLayout{{
      "2x16_staggered", "2x16_staggered",
      "single1x32", "single1x32"}};
  // ASSUMED channel numbering. For two-row layouts, ascending means
  // channels 0..15 in row 0 and 16..31 in row 1, both left-to-right.
  std::array<G4String,4> fPlaneChannelOrder{{
      "ascending", "ascending", "ascending", "ascending"}};
  // ASSUMED row-Z separation. Values for single-row layouts are provenance
  // only and do not move any fibre.
  std::array<G4double,4> fPlaneRowZSpacing{{
      0.8660254*mm, 0.8660254*mm, 1.*mm, 1.*mm}};
  G4String fDiagnosticRootFile = "../outputs/optical_v21_diagnostics";
  // Legacy trigger-yaw command convention: the placed angle is -command for
  // Trigger 1 and 90 degrees-command for Trigger 2.  Preserve that physical
  // placement while documenting the sign explicitly at the UI boundary.
  std::array<G4double,4> fTriggerYaw{{0., 0., 0., 0.}};
  // Centre Z positions for structural-acceptance geometry. Defaults preserve
  // the predecessor compact stack. Structural studies may widen the stack
  // explicitly by macro without changing validated default runs.
  std::array<G4bool,4> fTriggerEnabled{{true, true, false, false}};
  std::array<G4double,4> fTriggerZ{{25.*mm, -25.*mm, 50.*mm, -50.*mm}};
  std::array<G4double,2> fHodoscopeFanoutZ{{9.*mm, -9.*mm}};
  // Alignment-study extensions.  The default legacy reference and zero
  // offsets preserve the approved predecessor geometry exactly.  The
  // ideal_axis reference places both trigger physical centres on global Z.
  G4String fAlignmentReference = "legacy";
  std::array<G4double,4> fTriggerOffsetX{{0., 0., 0., 0.}};
  std::array<G4double,4> fTriggerOffsetY{{0., 0., 0., 0.}};
  // Analysis-only path oracle for comparing all-round hodoscope layouts
  // without changing the validated physical H1/H2 geometry.  It evaluates
  // whether the primary muon line crosses round X/Y active-core planes using
  // the current H1 round pitch/layout.  It does not create Geant4 volumes,
  // optical photons, channels, or digitization.
  G4bool fVirtualRoundModelEnabled = false;
  G4bool fVirtualRoundH3Enabled = false;
  G4double fVirtualRoundH3Z = 0.*mm;
  G4double fVirtualRoundH3OffsetX = 0.*mm;
  G4double fVirtualRoundH3OffsetY = 0.*mm;
  G4double fVirtualRoundH3Yaw = 0.*deg;
  // Physical all-round geometry option. Default false preserves the validated
  // mixed H1-round/H2-square detector. When true, H2 is constructed from
  // round SCSF-78 fibres with H2_X/H2_Y plane ownership and H2 copy-numbering.
  G4bool fH2PhysicalRound = false;
  G4int fH2SquareFibreCount = 32; // legacy default; current hardware macro explicitly uses 29
  G4bool fRecordAlignmentStudyNtuple = false;
  // Plane order: H1X, H1Y, H2X, H2Y.
  std::array<G4double,4> fPlaneYaw{{0., 0., 0., 0.}};
  // V16-compatible rigid module translations. Index 0=H1, 1=H2. Zero has an
  // arithmetic-free transform path so the V18-D default remains bit-exact.
  std::array<G4double,2> fHodoscopeOffsetX{{0., 0.}};
  std::array<G4double,2> fHodoscopeOffsetY{{0., 0.}};
  // Working mapping: [H1X,H1Y,H2X,H2Y] -> hodoscope/groups [3,1,2,4].
  // Plane order: H1X(round X), H1Y(round Y), H2X(square X), H2Y(square Y).
  // Hodoscope 4 = square Y is photo-confirmed; the other three are provisional
  // assignments inferred from the same physical layout.
  G4int fGainBoard[4] = {3, 1, 2, 4};
  // Structural layout option: keep the H2 fibre/fan-out geometry unchanged
  // internally, but place the whole H2 module with a rigid 180 degree
  // rotation about global Z. This is intentionally separate from the V15
  // +/-15 degree yaw-delta controls.
  G4bool fH2ModuleRotate180 = false;
  std::array<G4double,128> fFanoutChannelTransmission{};
  std::array<G4double,128> fFanoutMinRadius{};
  std::array<G4ThreeVector,128> fFanoutEndpoint{};
  std::array<G4ThreeVector,128> fFanoutTerminalDirection{};
  std::array<G4bool,128> fFanoutEndpointReady{};

  G4GenericMessenger* fMessenger = nullptr;
  G4Material* fSCSF78 = nullptr;
  G4Material* fBCF20 = nullptr;
  G4Material* fInnerCladding = nullptr;
  G4Material* fOuterCladding = nullptr;
  G4Material* fGrease = nullptr;
  G4Material* fSiPMWindow = nullptr;
  G4Material* fSilicon = nullptr;
  G4Material* fTriggerScint = nullptr;
  G4OpticalSurface* fSiPMSurface = nullptr;
  G4OpticalSurface* fMirrorSurface = nullptr;
  G4OpticalSurface* fFiberSurface = nullptr;
  G4OpticalSurface* fEndFaceSurface = nullptr;
  G4OpticalSurface* fBlackWrapSurface = nullptr;
  // One immutable limit object is shared by every fibre core/cladding
  // logical volume.  It is never allocated or attached for the default
  // zero-step-limit configuration.
  G4UserLimits* fFibreUserLimits = nullptr;
};
