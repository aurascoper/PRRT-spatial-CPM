// Study 1 reference/kernel arm: Lu-177 (ENSDF via G4RadioactiveDecay)
// transport in homogeneous water; voxel tally.
// Adapted from Geant4 rdecay02 (v11.4.2) canonical patterns:
//   - ion primary via G4IonTable::GetIonTable()->GetIon(Z, A, 0.)
//   - BiasedRDPhysics (G4RadioactiveDecay wrapper) from rdecay02's PhysicsList
//   - G4NuclideTable thresholds + Auger cascade as in rdecay02
//   - Run-class accumulation/merge per G4MTRun pattern
// Modes:
//   field   : world = ROI (n*pitch)^3; sources sampled ∝ activity map
//   kernel  : world = kernel box (kbox*pitch)^3; single source at center voxel
// Both modes write <out>.bin (float64 MeV/voxel, C-order) + <out>.json.

#ifndef REF_APP_HH
#define REF_APP_HH

#include "G4Run.hh"
#include "G4ParticleGun.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserRunAction.hh"
#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include <vector>
#include <cstdint>

struct Geom {
    long n = 0;               // voxels per axis
    double pitch_um = 0.0;
    long kbox = 0;            // kernel-box voxels per axis (kernel mode)
    std::vector<double> activity;   // n^3, C-order (field mode)
};

extern Geom g_geom;
extern G4String g_mode;       // "field" | "kernel"
extern G4String g_outPrefix;

class StudyRun : public G4Run {
public:
    StudyRun() : G4Run() { fEdep.assign(GetNumberOfEventToBeProcessed() > 0 ? 0 : 0, 0.0); }
    void Setup(size_t nvox) { fEdep.assign(nvox, 0.0); }
    void AddEdep(size_t idx, double e) { fEdep[idx] += e; }
    void RecordDecay() { fDecays++; }
    const std::vector<double>& Edep() const { return fEdep; }
    G4long Decays() const { return fDecays; }
    void Merge(const G4Run* aRun) override {
        const StudyRun* local = static_cast<const StudyRun*>(aRun);
        if (fEdep.size() != local->fEdep.size()) fEdep = local->fEdep;  // first merge sizes
        for (size_t i = 0; i < fEdep.size(); ++i) fEdep[i] += local->fEdep[i];
        fDecays += local->fDecays;
        G4Run::Merge(aRun);
    }
private:
    std::vector<double> fEdep;
    G4long fDecays = 0;
};

class StudyDetector : public G4VUserDetectorConstruction {
public:
    G4VPhysicalVolume* Construct() override;
};

class StudyPhysics : public G4VModularPhysicsList {
public:
    StudyPhysics();
    void ConstructParticle() override;
    void ConstructProcess() override;
    void SetCuts() override;
};

class StudyPrimary : public G4VUserPrimaryGeneratorAction {
public:
    StudyPrimary();
    void GeneratePrimaries(G4Event* anEvent) override;
private:
    G4ParticleGun* fGun = nullptr;
    std::vector<double> fCum;   // activity CDF (field mode)
    double fTotal = 0.0;
    G4int fIonZ = 71, fIonA = 177;
};

class StudyRunAction : public G4UserRunAction {
public:
    G4Run* GenerateRun() override { return new StudyRun; }  // non-const in G4 11.4
    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;
};

class StudyStepping : public G4UserSteppingAction {
public:
    void UserSteppingAction(const G4Step* step) override;
};

class StudyActions : public G4VUserActionInitialization {
public:
    void Build() const override;
    void BuildForMaster() const override;
};

#endif