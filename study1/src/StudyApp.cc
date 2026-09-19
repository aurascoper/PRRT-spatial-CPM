#include "StudyApp.hh"
#include "G4Box.hh"
#include "G4NistManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4Geantino.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmParameters.hh"
#include "G4NuclideTable.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4RadioactiveDecay.hh"
#include "G4PhysicsListHelper.hh"
#include "G4GenericIon.hh"
#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclearLevelData.hh"
#include "G4Radioactivation.hh"
#include <fstream>
#include <algorithm>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

Geom g_geom;
G4String g_mode = "field";
G4String g_outPrefix = "out";

// ---------------- detector ----------------
G4VPhysicalVolume* StudyDetector::Construct() {
    auto* nist = G4NistManager::Instance();
    auto* water = nist->FindOrBuildMaterial("G4_WATER");
    long n = (g_mode == "kernel") ? g_geom.kbox : g_geom.n;
    double half = 0.5 * n * g_geom.pitch_um * CLHEP::um;
    auto* box = new G4Box("world", half, half, half);
    auto* log = new G4LogicalVolume(box, water, "world");
    return new G4PVPlacement(nullptr, {}, log, "world", nullptr, false, 0);
}

// ---------------- physics ----------------
StudyPhysics::StudyPhysics() : G4VModularPhysicsList() {
    G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(0.1 * CLHEP::picosecond);
    G4NuclideTable::GetInstance()->SetLevelTolerance(1.0 * CLHEP::eV);

    RegisterPhysics(new G4EmStandardPhysics());
    auto* param = G4EmParameters::Instance();
    param->SetAugerCascade(true);
    // fine stepping for 10-40 um voxels: dRoverRange 0.1, final 1 um
    param->SetStepFunction(0.1, 1 * CLHEP::um);

    RegisterPhysics(new G4DecayPhysics());
}

void StudyPhysics::ConstructProcess() {
    // canonical rdecay02 pattern: register the RDM process on G4GenericIon.
    // G4Radioactivation = analogue (unbiased) RDM; use it for the reference
    // arm (no biasing anywhere in this study).
    G4VModularPhysicsList::ConstructProcess();
    G4LossTableManager* man = G4LossTableManager::Instance();
    G4VAtomDeexcitation* ad = man->AtomDeexcitation();
    if (!ad) {
        G4EmParameters::Instance()->SetAugerCascade(true);
        ad = new G4UAtomicDeexcitation();
        man->SetAtomDeexcitation(ad);
        ad->InitialiseAtomicDeexcitation();
    }
    G4DeexPrecoParameters* deex = G4NuclearLevelData::GetInstance()->GetParameters();
    deex->SetStoreICLevelData(true);
    deex->SetMaxLifeTime(G4NuclideTable::GetInstance()->GetThresholdOfHalfLife() / std::log(2.));
    G4PhysicsListHelper::GetPhysicsListHelper()->RegisterProcess(
        new G4Radioactivation(), G4GenericIon::GenericIon());
}

void StudyPhysics::ConstructParticle() {
    G4BosonConstructor b; b.ConstructParticle();
    G4LeptonConstructor l; l.ConstructParticle();
    G4MesonConstructor m; m.ConstructParticle();
    G4BaryonConstructor ba; ba.ConstructParticle();
    G4IonConstructor io; io.ConstructParticle();
    G4ShortLivedConstructor s; s.ConstructParticle();
}

void StudyPhysics::SetCuts() {
    // 1 um production cuts: fine enough for 10 um voxels
    SetCutValue(1 * CLHEP::um, "e-");
    SetCutValue(1 * CLHEP::um, "e+");
    SetCutValue(1 * CLHEP::um, "gamma");
    SetCutValue(1 * CLHEP::um, "proton");
}

// ---------------- primary ----------------
StudyPrimary::StudyPrimary() {
    fGun = new G4ParticleGun(1);
    fGun->SetParticleEnergy(0 * CLHEP::eV);
    fGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    if (g_mode == "field") {
        double total = 0.0;
        fCum.reserve(g_geom.activity.size());
        for (double a : g_geom.activity) { total += a; fCum.push_back(total); }
        fTotal = total;
    }
}

void StudyPrimary::GeneratePrimaries(G4Event* anEvent) {
    G4ParticleDefinition* ion =
        G4IonTable::GetIonTable()->GetIon(fIonZ, fIonA, 0.0);  // Lu-177 g.s.
    fGun->SetParticleDefinition(ion);
    fGun->SetParticleCharge(0. * CLHEP::eplus);

    long n = (g_mode == "kernel") ? g_geom.kbox : g_geom.n;
    double pitch = g_geom.pitch_um * CLHEP::um;
    double half = 0.5 * n * pitch;
    G4double x, y, z;
    if (g_mode == "kernel") {
        // single source at the center voxel's center. For odd kbox the
        // center voxel center is the world origin (0,0,0); all kernel boxes
        // in this study are odd (361/181/91/51/47).
        x = y = z = 0.0;
    } else {
        double u = G4UniformRand() * fTotal;
        size_t idx = std::lower_bound(fCum.begin(), fCum.end(), u) - fCum.begin();
        idx = std::min(idx, fCum.size() - 1);
        long k = idx % n; long tmp = idx / n;
        long j = tmp % n; long i = tmp / n;
        // voxel (i,j,k) center in world coordinates (world spans
        // [-half, +half]; voxel 0 center = 0.5*pitch - half):
        x = (i + 0.5) * pitch - half;
        y = (j + 0.5) * pitch - half;
        z = (k + 0.5) * pitch - half;
    }
    // ORIGIN BUG FIX (protocol v1.4): v1.0-v1.3 wrote (x-half+half, ...) ==
    // (x, ...), placing sources in the positive octant [0, 2*half]^3 so
    // ~7/8 of field primaries were born OUTSIDE the world and transported
    // nothing, and the kernel source sat at a box corner capturing one
    // octant of the dose spread. Caught by energy conservation: per-decay
    // deposit was 0.01643 MeV vs 0.1352 MeV smoke anchor = 12.2% ~= 1/8.
    fGun->SetParticlePosition(G4ThreeVector(x, y, z));
    fGun->GeneratePrimaryVertex(anEvent);
    reinterpret_cast<StudyRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun())->RecordDecay();
}

// ---------------- run actions ----------------
void StudyRunAction::BeginOfRunAction(const G4Run* aRun) {
    StudyRun* run = const_cast<StudyRun*>(static_cast<const StudyRun*>(aRun));
    long n = (g_mode == "kernel") ? g_geom.kbox : g_geom.n;
    run->Setup((size_t)n * n * n);
}

void StudyRunAction::EndOfRunAction(const G4Run* aRun) {
    const StudyRun* run = static_cast<const StudyRun*>(aRun);
    // master-end only: write file once
    if (!G4Threading::IsMasterThread()) return;
    long n = (g_mode == "kernel") ? g_geom.kbox : g_geom.n;
    std::string bin = g_outPrefix + ".bin";
    std::ofstream f(bin, std::ios::binary);
    if (!f.is_open()) {
        G4cerr << "[study] FATAL: cannot open output " << bin
               << " (directory missing?) — REFUSING to report success" << G4endl;
        G4Exception("StudyRunAction", "study01", FatalException,
                    "output file could not be opened");
    }
    f.write(reinterpret_cast<const char*>(run->Edep().data()),
            std::streamsize(run->Edep().size() * sizeof(double)));
    f.flush();
    if (!f.good()) {
        G4Exception("StudyRunAction", "study02", FatalException,
                    "output file write failed");
    }
    json meta;
    meta["mode"] = g_mode;
    meta["n"] = n;
    meta["pitch_um"] = g_geom.pitch_um;
    meta["decays_simulated"] = run->Decays();
    meta["units"] = "MeV per voxel, C-order";
    meta["geant4_version"] = "11.4.2";
    meta["physics"] = "G4EmStandardPhysics + G4Decay + G4RadioactiveDecay (ENSDF, RDM 6.1.2)";
    meta["cuts_um"] = {{"e", 1}, {"gamma", 1}};
    meta["step_function"] = "dRoverRange 0.1, final 1 um";
    std::ofstream mj(g_outPrefix + ".json");
    if (!mj.is_open()) {
        G4Exception("StudyRunAction", "study03", FatalException,
                    "meta output file could not be opened");
    }
    mj << meta.dump(1) << "\n";
    mj.flush();
    if (!mj.good()) {
        G4Exception("StudyRunAction", "study04", FatalException,
                    "meta write failed");
    }
    G4cout << "[study] wrote " << bin << "  decays=" << run->Decays() << G4endl;
}

// ---------------- stepping ----------------
void StudyStepping::UserSteppingAction(const G4Step* step) {
    G4double e = step->GetTotalEnergyDeposit();
    if (e <= 0.) return;
    StudyRun* run = static_cast<StudyRun*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    long n = (g_mode == "kernel") ? g_geom.kbox : g_geom.n;
    double pitch = g_geom.pitch_um * CLHEP::um;
    double half = 0.5 * n * pitch;
    const G4ThreeVector& pos = step->GetPreStepPoint()->GetPosition();
    long i = (long)std::floor((pos.x() + half) / pitch);
    long j = (long)std::floor((pos.y() + half) / pitch);
    long k = (long)std::floor((pos.z() + half) / pitch);
    if (i < 0 || j < 0 || k < 0 || i >= n || j >= n || k >= n) return;
    run->AddEdep((size_t(i) * n + j) * n + k, e);
}

// ---------------- action init ----------------
void StudyActions::Build() const {
    SetUserAction(new StudyPrimary);
    SetUserAction(new StudyRunAction);
    SetUserAction(new StudyStepping);
}
void StudyActions::BuildForMaster() const {
    SetUserAction(new StudyRunAction);
}

// ---------------- main ----------------
#include "G4RunManagerFactory.hh"
#include "G4Threading.hh"
#include <cstring>

int main(int argc, char** argv) {
    if (argc < 6) {
        G4cout << "usage: study1_app mode(field|kernel) pitch_um geom_bin_or_'-' out_prefix nEvents [threads] [kbox]\n";
        return 1;
    }
    g_mode = argv[1];
    g_geom.pitch_um = std::stod(argv[2]);
    long nEvents = std::stol(argv[5]);
    int threads = (argc > 6) ? std::stoi(argv[6]) : 1;
    long kbox = (argc > 7) ? std::stol(argv[7]) : 0;

    G4Random::setTheEngine(new CLHEP::RanecuEngine);
    // v1.2: per-run seed from argv (independent runs MUST NOT share a seed,
    // else h1 is a strict subset of h2 and the statistical ladder is
    // vacuous — proven 2026-09-16: h2/h1 = 3.9988 with r1<=r2 at every
    // nonzero voxel). Reproducibility = same seed reproduces same run.
    long seed = (argc > 8) ? std::stol(argv[8]) : 20260916;
    G4Random::setTheSeed(seed);

    if (g_mode == "field") {
        // geom bin: int32 cell_id[n^3] then float64 activity[n^3]; n from sidecar json
        std::string gp = argv[3];
        std::ifstream gj(gp + "_meta.json");   // geometry.py writes <prefix>_meta.json
        json gm = json::parse(gj);
        g_geom.n = gm["n"];
        size_t nv = (size_t)g_geom.n * g_geom.n * g_geom.n;
        std::vector<int32_t> cid(nv);
        g_geom.activity.resize(nv);
        std::ifstream gb(gp + ".bin", std::ios::binary);
        gb.read(reinterpret_cast<char*>(cid.data()), std::streamsize(nv * sizeof(int32_t)));
        gb.read(reinterpret_cast<char*>(g_geom.activity.data()), std::streamsize(nv * sizeof(double)));
    } else {
        g_geom.kbox = kbox;
    }
    g_outPrefix = argv[4];

    auto* run = G4RunManagerFactory::CreateRunManager(G4RunManagerType::MT);
    run->SetNumberOfThreads(threads);
    run->SetUserInitialization(new StudyDetector);
    run->SetUserInitialization(new StudyPhysics);
    run->SetUserInitialization(new StudyActions);
    run->Initialize();
    run->BeamOn(nEvents);
    delete run;
    return 0;
}