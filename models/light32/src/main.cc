#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "V21GeometryTolerance.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4RunManagerFactory.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4VModularPhysicsList.hh"

namespace {
class OpticalPhysicsNoPhotonScintillation final : public G4OpticalPhysics {
 public:
  void ConstructProcess() override {
    G4OpticalPhysics::ConstructProcess();
    auto* manager=G4OpticalPhoton::Definition()->GetProcessManager();
    if(!manager)return;
    auto* processes=manager->GetProcessList();
    for(G4int index=processes->size()-1;index>=0;--index) {
      auto* process=(*processes)[index];
      if(process&&process->GetProcessName()=="Scintillation")
        manager->RemoveProcess(process);
    }
  }
};
}

int main(int argc, char** argv) {
  // Geant4 requires this call before the first solid is created. Keep it
  // before UI, run-manager, and detector construction so the configured
  // world-scale tolerance is unambiguous and fail-closed.
  hodoscope::v21::ConfigureGeometryToleranceBeforeGeometry();

  const G4bool interactive = argc == 1 || G4String(argv[1]) == "--gui";
  auto* ui = interactive ? new G4UIExecutive(argc, argv) : nullptr;
  // V21 deliberately uses one serial Geant4 kernel per process. The official
  // geometry-tolerance singleton is thread-local, while the public setter
  // cannot be called in an MT worker after the shared solid store is built.
  // Campaign-level parallelism is therefore multiple independent processes.
  auto* runManager =
      G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);
  G4cout << "V21_RUN_MANAGER mode=serial process_parallelism=external"
         << G4endl;
  auto* detector = new DetectorConstruction;
  runManager->SetUserInitialization(detector);

  // This detector study needs precision electromagnetic muon transport,
  // decay, scintillation, and optical photon transport.  FTFP_BERT also
  // enabled rare muon-nuclear/hypernuclear secondaries; long runs crashed in
  // G4HyperH4::Definition/G4Track destruction on Geant4 11.4.1.  Keep those
  // rare hadronic channels out of the production model until that upstream
  // failure is independently resolved.
  auto* physics = new G4VModularPhysicsList;
  physics->SetVerboseLevel(0);
  physics->RegisterPhysics(new G4EmStandardPhysics_option4);
  physics->RegisterPhysics(new G4DecayPhysics);
  physics->RegisterPhysics(new G4StepLimiterPhysics);
  // Preserve the approved optical process registration exactly.  Even an
  // apparently redundant process-table cleanup belongs to the bit-exact
  // software fixture until an isolated regression proves it removable.
  physics->RegisterPhysics(new OpticalPhysicsNoPhotonScintillation);
  runManager->SetUserInitialization(physics);
  runManager->SetUserInitialization(new ActionInitialization(detector));

  auto* vis = new G4VisExecutive("quiet");
  vis->Initialize();
  auto* uiManager = G4UImanager::GetUIpointer();
  G4int commandStatus = 0;
  if (ui) {
    const G4String macro = argc >= 3 ? argv[2] : "gui.mac";
    commandStatus = uiManager->ApplyCommand("/control/execute " + macro);
    if (commandStatus == 0) {
      ui->SessionStart();
    } else {
      G4cerr << "V21 startup macro failed: " << macro
             << " (G4 command status " << commandStatus << ")" << G4endl;
    }
    delete ui;
  } else {
    commandStatus =
      uiManager->ApplyCommand(G4String("/control/execute ") + argv[1]);
    if (commandStatus != 0) {
      G4cerr << "V21 batch macro failed: " << argv[1]
             << " (G4 command status " << commandStatus << ")" << G4endl;
    }
  }
  delete vis;
  delete runManager;
  return commandStatus == 0 ? 0 : 20;
}
