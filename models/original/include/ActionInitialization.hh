#pragma once

#include "G4VUserActionInitialization.hh"

class DetectorConstruction;

class ActionInitialization final : public G4VUserActionInitialization {
 public:
  explicit ActionInitialization(const DetectorConstruction* detector) : fDetector(detector) {}
  void Build() const override;
  void BuildForMaster() const override;
 private:
  const DetectorConstruction* fDetector;
};

