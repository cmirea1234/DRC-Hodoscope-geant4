#include "DetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Exception.hh"
#include <algorithm>

void DetectorConstruction::BuildLightStraightFibres(G4LogicalVolume* world) {
  if(fUseBuilding || fLightStraightFanout ||
     (!fLightStraightReadout && (!fPathOnlyMuonTruth || fUseWaveformDigitization || fFixedPhotonTest)) ||
     (fLightStraightReadout && (fSiPMOffsetU!=0. || fSiPMOffsetV!=0.)))
    G4Exception("BuildLightStraightFibres","LIGHT_STRAIGHT_FIBRES_CONFIG",FatalException,
                "Straight model requires no building/bridge; bare mode requires path-only; equivalent readout requires centred coupling.");
  const std::array<G4String,4> prefix{{"H1XFanout","H1YFanout","H2XFanout","H2YFanout"}};
  const std::array<G4String,4> names{{"H1_X","H1_Y","H2_X","H2_Y"}};
  for(int plane=0;plane<4;++plane) {
    const bool round=plane<2;
    const bool single=round?GetH1SingleClad():GetH2SingleClad();
    const auto bounds=GetActiveLongitudinalBounds(plane);
    const auto first=fFibreFarEndExtension>0.?bounds[0]-fFibreFarEndExtension:-16.*mm;
    // Keep the previous longitudinal extent, not the curved arc length.
    constexpr G4double last=136.*mm;
    const auto half=.5*(last-first), middle=.5*(last+first);
    const auto outerSize=GetFanoutOuterSize(plane);
    const auto innerSize=single?outerSize:.4845*mm;
    const auto coreSize=round?(single?.4795*mm:.4695*mm):.4695*mm;
    const auto makeSolid=[&](const G4String& suffix,G4double radius)->G4VSolid* {
      if(round) return new G4Tubs(prefix[plane]+suffix,0.,radius,half,0.,360.*deg);
      return new G4Box(prefix[plane]+suffix,radius,radius,half);
    };
    auto* outerSolid=makeSolid("StraightOuterSolid",outerSize);
    auto* coreSolid=makeSolid("StraightCoreSolid",coreSize);
    auto* innerSolid=single?nullptr:makeSolid("StraightInnerSolid",innerSize);
    const auto nominalUnit=plane%2==0?G4ThreeVector(1,0,0):G4ThreeVector(0,1,0);
    const auto unit=PlaneNominalToWorldDirection(plane,nominalUnit).unit();
    const auto side=G4ThreeVector(0,0,1).cross(unit).unit();
    const auto normal=unit.cross(side).unit();
    const G4RotationMatrix rotation(side,normal,unit);
    for(int channel=0;channel<GetPlaneFibreCount(plane);++channel) {
      const auto suffix="_"+std::to_string(channel);
      auto* outer=new G4LogicalVolume(outerSolid,single?fInnerCladding:fOuterCladding,prefix[plane]+"OuterLV"+suffix);
      auto* core=new G4LogicalVolume(coreSolid,round?fSCSF78:fBCF20,prefix[plane]+"CoreLV"+suffix);
      auto* mother=outer;
      G4VPhysicalVolume* innerPV=nullptr;
      if(!single) {
        auto* inner=new G4LogicalVolume(innerSolid,fInnerCladding,prefix[plane]+"InnerLV"+suffix);
        ApplyFibreUserLimits(inner);
        innerPV=new G4PVPlacement(nullptr,{},inner,prefix[plane]+"InnerClad",outer,false,0,false);
        inner->SetVisAttributes(G4VisAttributes::GetInvisible());
        mother=inner;
      }
      ApplyFibreUserLimits(outer); ApplyFibreUserLimits(core);
      auto* corePV=new G4PVPlacement(nullptr,{},core,prefix[plane]+"Core",mother,false,0,false);
      new G4LogicalSkinSurface(prefix[plane]+"CoreSurface"+suffix,core,fFiberSurface);
      auto* vis=new G4VisAttributes(round?G4Colour(.1,1.,.3):G4Colour(.15,.45,1.));
      vis->SetForceSolid(true); core->SetVisAttributes(vis);
      outer->SetVisAttributes(new G4VisAttributes(G4Colour(.7,.9,1.,.15)));
      const auto transverse=GetPlaneSourceTransverse(plane,channel);
      const auto z=GetPlaneSourceZ(plane,channel);
      const auto nominal=plane%2==0?G4ThreeVector(middle,transverse,z):G4ThreeVector(transverse,middle,z);
      const auto centre=PlaneNominalToWorldPosition(plane,nominal);
      const int local=(plane%2)*32+channel;
      auto* outerPV=new G4PVPlacement(G4Transform3D(rotation,centre),outer,names[plane],world,false,local,false);
      const int global=plane*32+channel;
      // Compatibility getter always denotes the fibre cap, not the sensor centre.
      fFanoutEndpoint[global]=centre+half*unit;
      fFanoutTerminalDirection[global]=unit;
      fFanoutEndpointReady[global]=true;
      if(fLightStraightReadout) {
        // Deliberately footprint-matched equivalent receiver, NOT the real
        // 1.3 mm square SiPM. Preserve materials, axial layer thicknesses,
        // end-face surface and wavelength-dependent PDE of the curved model.
        // Smaller aperture/shorter straight length are explicit control changes.
        const auto layer=[&](const G4String& tag,G4Material* material,G4double halfZ) {
          G4VSolid* solid=round
            ? static_cast<G4VSolid*>(new G4Tubs(prefix[plane]+tag+"Solid"+suffix,0.,outerSize,halfZ,0.,360.*deg))
            : static_cast<G4VSolid*>(new G4Box(prefix[plane]+tag+"Solid"+suffix,outerSize,outerSize,halfZ));
          return new G4LogicalVolume(solid,material,prefix[plane]+tag+"LV"+suffix);
        };
        auto* grease=layer("Grease",fGrease,.5*fGreaseGap);
        auto* window=layer("SiPMWindow",fSiPMWindow,.05*mm);
        auto* sensor=layer("SiPM",fSilicon,.005*mm);
        grease->SetVisAttributes(new G4VisAttributes(G4Colour(1.,1.,.3,.35)));
        window->SetVisAttributes(new G4VisAttributes(G4Colour(.8,.8,.8,.25)));
        sensor->SetVisAttributes(new G4VisAttributes(G4Colour(1.,.85,.05)));
        const auto end=fFanoutEndpoint[global];
        auto* greasePV=new G4PVPlacement(G4Transform3D(rotation,end+unit*(.5*fGreaseGap)),grease,
            round?"H1_Grease":"H2_Grease",world,false,local,false);
        auto* windowPV=new G4PVPlacement(G4Transform3D(rotation,end+unit*(fGreaseGap+.05*mm)),window,
            "SiPMWindow",world,false,local,false);
        auto* sensorPV=new G4PVPlacement(G4Transform3D(rotation,end+unit*(fGreaseGap+.105*mm)),sensor,
            "SiPM",world,false,global,false);
        new G4LogicalBorderSurface(prefix[plane]+"StraightCoreEnd"+suffix,corePV,greasePV,fEndFaceSurface);
        new G4LogicalBorderSurface(prefix[plane]+"StraightOuterEnd"+suffix,outerPV,greasePV,fEndFaceSurface);
        if(innerPV)
          new G4LogicalBorderSurface(prefix[plane]+"StraightInnerEnd"+suffix,innerPV,greasePV,fEndFaceSurface);
        new G4LogicalBorderSurface(prefix[plane]+"StraightSiPM"+suffix,windowPV,sensorPV,fSiPMSurface);
        G4cout<<"LIGHT_EQUIVALENT_READOUT channel="<<global<<" plane="<<plane
              <<" footprint="<<(round?"circle":"square")<<" outer_half_mm="<<outerSize/mm
              <<" grease_mm="<<fGreaseGap/mm<<" window_mm=0.1 sensor_mm=0.01"
              <<" nominal_chip_mm="<<fSiPMActiveSize/mm<<" physical_chip_geometry=false"<<G4endl;
      }
    }
    G4cout<<"LIGHT_PARALLEL_FIBRES plane="<<plane<<" count="<<GetPlaneFibreCount(plane)
          <<" first_mm="<<first/mm<<" last_mm="<<last/mm<<" length_mm="<<2*half/mm
          <<" core_half_mm="<<coreSize/mm<<" outer_half_mm="<<outerSize/mm<<G4endl;
  }
  G4cout<<"LIGHT_PARALLEL_FIBRES_COMPLETE fanout=0 sipm="<<fLightStraightReadout
        <<" grease="<<fLightStraightReadout<<" window="<<fLightStraightReadout<<" building=0"
        <<" scope="<<(fLightStraightReadout?"equivalent_readout_not_physical_chip":"geometry_path_edep_not_detection_efficiency")<<G4endl;
}
