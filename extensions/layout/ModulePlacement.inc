void DetectorConstruction::RegisterModuleVolume(int module,G4VPhysicalVolume* volume) {
  if(volume && volume->GetMotherLogical() && volume->GetMotherLogical()->GetName()=="WorldLV")
    fModuleVolumes.emplace_back(module,volume);
}
void DetectorConstruction::ApplyModuleLayout(G4VPhysicalVolume* worldPV) {
  if(!fModuleLayout.configured)return;
  if(!fFanoutPreview||fFanoutPreviewLite||fUseBuilding||fTriggerEnabled[2]||fTriggerEnabled[3]||fVirtualRoundModelEnabled)
    throw std::runtime_error("Module poses require baseline full fanout, no building/extra triggers/virtual model");
  if(fHodoscopeFanoutZ!=std::array<G4double,2>{{20*mm,-20*mm}} ||
     (LegacyGetTriggerCentre(0)-G4ThreeVector(-25*mm,0,80*mm)).mag()>1e-12*mm ||
     (LegacyGetTriggerCentre(1)-G4ThreeVector(0,-25*mm,-75*mm)).mag()>1e-12*mm)
    throw std::runtime_error("Use the original baseline macro; change module poses instead of legacy Z/trigger offsets");
  if(fModuleLayout.changed()&&fCosmicBeam)
    throw std::runtime_error("Changed module poses currently require an explicit beam or fixed-photon source");
  std::array<int,4> counts{};
  std::ofstream audit("layout_placements.jsonl");audit<<std::setprecision(17);
  for(const auto& item:fModuleVolumes){
    const int id=item.first;auto* pv=item.second;const auto& pose=fModuleLayout.poses[id];
    ++counts[id];
    const auto oldT=pv->GetObjectTranslation();const auto oldR=pv->GetObjectRotationValue();
    if(!pose.identity){
      pv->SetTranslation(pose.point(oldT));
      pv->SetRotation(new G4RotationMatrix((pose.rotation*oldR).inverse()));
    }
    G4ThreeVector low,high;pv->GetLogicalVolume()->GetSolid()->BoundingLimits(low,high);
    for(int mask=0;mask<8;++mask){
      const G4ThreeVector corner((mask&1)?high.x():low.x(),(mask&2)?high.y():low.y(),(mask&4)?high.z():low.z());
      const auto point=pv->GetObjectRotationValue()*corner+pv->GetObjectTranslation();
      for(int axis=0;axis<3;++axis)if(std::abs(point[axis])>=fModuleLayout.worldHalf[axis])
        throw std::runtime_error("Module leaves world: enlarge world_half_size_mm or change pose");
    }
    const auto t=pv->GetObjectTranslation();const auto r=pv->GetObjectRotationValue();
    audit<<"{\"module\":"<<id<<",\"name\":\""<<pv->GetName()<<"\",\"copy\":"<<pv->GetCopyNo()<<",\"position_mm\":["<<t.x()/mm<<","<<t.y()/mm<<","<<t.z()/mm<<"],\"rotation\":[";
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)audit<<(i||j?",":"")<<r(i,j);
    audit<<"]}\n";
  }
  if(counts!=std::array<int,4>{{1,1,320,232}})
    throw std::runtime_error("Module membership count mismatch; expected trigger/core/cladding/readout assemblies");
  for(int channel=0;channel<128;++channel)if(fFanoutEndpointReady[channel]){
    const auto& pose=fModuleLayout.poses[2+channel/64];
    fFanoutEndpoint[channel]=pose.point(fFanoutEndpoint[channel]);
    fFanoutTerminalDirection[channel]=pose.direction(fFanoutTerminalDirection[channel]);
  }
  fModuleLayout.applied=true;
  G4Navigator navigator;navigator.SetWorldVolume(worldPV);
  int checks=0;
  for(int plane=0;plane<4;++plane)for(int channel=0;channel<GetPlaneFibreCount(plane);++channel){
    const auto t=GetPlaneSourceTransverse(plane,channel),z=GetPlaneSourceZ(plane,channel);
    const auto nominal=plane%2==0?G4ThreeVector(0,t,z):G4ThreeVector(t,0,z);
    const auto world=PlaneNominalToWorldPosition(plane,nominal);
    if((PlaneWorldToNominalPosition(plane,world)-nominal).mag()>1e-8*mm)throw std::runtime_error("Pose inverse mismatch");
    auto* pv=navigator.LocateGlobalPointAndSetup(world,nullptr,false);
    if(!pv||pv->GetLogicalVolume()->GetMaterial()!=(plane<2?fSCSF78:fBCF20))throw std::runtime_error("Moved core navigator check failed");
    ++checks;
  }
  for(const auto& item:fModuleVolumes){auto* expected=item.second;
    if(expected->GetName()!="SiPM" && item.first>=2)continue;
    auto* found=navigator.LocateGlobalPointAndSetup(expected->GetObjectTranslation(),nullptr,false);
    if(found!=expected)throw std::runtime_error("Moved sensor/trigger navigator identity mismatch");++checks;
  }
  G4cout<<"MODULE_POSE_NAVIGATION_PASS checks="<<checks<<G4endl;
}
G4double DetectorConstruction::ValidationSurfaceDistance(int surface,const G4ThreeVector& p) const {
  if(surface<2)return TriggerWorldToLocalPosition(surface,p).z();
  const int plane=surface-2;return PlaneWorldToNominalPosition(plane,p).z()-GetValidationSurfaceZ(surface);
}
