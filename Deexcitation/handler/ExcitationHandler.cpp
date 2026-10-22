//
// Created by Artem Novikov on 17.05.2023.
//

#include "ExcitationHandler.h"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Alpha.hh>
#include <G4BaryonConstructor.hh>
#include <G4BosonConstructor.hh>
#include <G4Deuteron.hh>
#include <G4Electron.hh>
#include <G4Evaporation.hh>
#include <G4FermiBreakUpAN.hh>
#include <G4Gamma.hh>
#include <G4GenericIon.hh>
#include <G4He3.hh>
#include <G4IonConstructor.hh>
#include <G4Ions.hh>
#include <G4LeptonConstructor.hh>
#include <G4LorentzVector.hh>
#include <G4MesonConstructor.hh>
#include <G4Neutron.hh>
#include <G4NistManager.hh>
#include <G4ParticleTable.hh>
#include <G4PhotonEvaporation.hh>
#include <G4PionMinus.hh>
#include <G4PionPlus.hh>
#include <G4PionZero.hh>
#include <G4ProcessManager.hh>
#include <G4Proton.hh>
#include <G4StatMF.hh>
#include <G4StateManager.hh>
#include <G4Triton.hh>

#include <algorithm>
#include <string>

namespace {
  constexpr size_t kEvaporationIterationThreshold = 1e3;

  const std::string error_no_model = "no model was applied, check conditions";

  class FermiBreakUpWrapper : public G4FermiBreakUpAN {
   public:
    using G4FermiBreakUpAN::G4FermiBreakUpAN;

    void BreakFragment(G4FragmentVector* results, G4Fragment* the_nucleus) override {
      auto* deletable_fragment = new G4Fragment(*the_nucleus);
      auto old_size = results->size();
      G4FermiBreakUpAN::BreakFragment(results, deletable_fragment);
      if (old_size == results->size()) {
        // NOLINTNEXTLINE(clang-analyzer-cplusplus.NewDelete)
        delete deletable_fragment;
      }
    }
  };

  class DataCleaner {
   public:
    explicit DataCleaner(G4FragmentVector& results) : results_(results) {}

    ~DataCleaner() {
      for (auto* ptr : results_) {
        delete ptr;
      }
    }

   private:
    G4FragmentVector& results_;
  };

  void ClearSingularResults(G4FragmentVector& results, G4Fragment* initial) {
    for (auto* fragment_ptr : results) {
      if (fragment_ptr != initial) {
        delete fragment_ptr;
      }
    }
  }

  constexpr G4int HashParticle(G4int atomic_mass, G4int charge) { return atomic_mass * 1000 + charge; }

  G4ParticleDefinition* SpecialParticleDefinition(const G4Fragment& fragment) {
    switch (HashParticle(fragment.GetA_asInt(), fragment.GetZ_asInt())) {
      case HashParticle(0, 0): {
        return G4Gamma::GammaDefinition();
      }

      case HashParticle(0, -1): {
        return G4Electron::ElectronDefinition();
      }

      case HashParticle(-1, 1): {
        return G4PionPlus::PionPlus();
      }

      case HashParticle(-1, -1): {
        return G4PionMinus::PionMinus();
      }

      case HashParticle(-1, 0): {
        return G4PionZero::PionZero();
      }

      case HashParticle(1, 0): {
        return G4Neutron::NeutronDefinition();
      }

      case HashParticle(1, 1): {
        return G4Proton::ProtonDefinition();
      }

      case HashParticle(2, 1): {
        return G4Deuteron::DeuteronDefinition();
      }

      case HashParticle(3, 1): {
        return G4Triton::TritonDefinition();
      }

      case HashParticle(3, 2): {
        return G4He3::He3Definition();
      }

      case HashParticle(4, 2): {
        return G4Alpha::AlphaDefinition();
      }

      default: {
        return nullptr;
      }
    }
  }

  void EvaporationError(const G4Fragment& fragment, const G4Fragment& current_fragment, size_t iteration_count) {
    G4ExceptionDescription ed;
    ed << "Infinite loop in the de-excitation module: " << iteration_count << " iterations \n"
       << "      Initial fragment: \n"
       << fragment << "\n      Current fragment: \n"
       << current_fragment;
    G4Exception("ExcitationHandler::BreakItUp", "", FatalException, ed, "Stop execution");
  }
}  // namespace

ExcitationHandler::ExcitationHandler()
    : multi_fragmentation_model_(DefaultMultiFragmentation())
    , fermi_break_up_model_(DefaultFermiBreakUp())
    , evaporation_model_(DefaultEvaporation())
    , photon_evaporation_model_(DefaultPhotonEvaporation())
    , neutron_decay_model_(DefaultNeutronDecay())
    , multi_fragmentation_condition_(DefaultMultiFragmentationCondition())
    , fermi_condition_(DefaultFermiBreakUpCondition())
    , photon_evaporation_condition_(DefaultPhotonEvaporationCondition())
    , evaporation_condition_(DefaultEvaporationCondition())
    , neutron_decay_condition_(DefaultNeutronDecayCondition()) {
  evaporation_model_->SetFermiBreakUp(fermi_break_up_model_.get());
  evaporation_model_->SetPhotonEvaporation(photon_evaporation_model_.get());

  G4BosonConstructor::ConstructParticle();
  G4LeptonConstructor::ConstructParticle();
  G4MesonConstructor::ConstructParticle();
  G4BaryonConstructor::ConstructParticle();
  G4IonConstructor::ConstructParticle();

  auto* gion = G4GenericIon::GenericIon();
  auto* manager = new G4ProcessManager(gion);
  manager->SetVerboseLevel(0);
  gion->SetProcessManager(manager);

  G4StateManager::GetStateManager()->SetNewState(G4State_Init);  // To let create ions
  auto* particle_table = G4ParticleTable::GetParticleTable();
  auto* ion_table = particle_table->GetIonTable();
  particle_table->SetReadiness();
  ion_table->CreateAllIon();
  ion_table->CreateAllIsomer();
}

ExcitationHandler::~ExcitationHandler() {
  // NOLINTNEXTLINE(bugprone-unused-return-value)
  photon_evaporation_model_.release();  // otherwise, SegFault in evaporation destructor
}

std::vector<G4ReactionProduct> ExcitationHandler::BreakItUp(const G4Fragment& fragment) {
  auto* nist = G4NistManager::Instance();
  G4FragmentVector results;
  const auto cleaner = DataCleaner(results);
  FragmentQueue evaporation_queue;
  FragmentQueue photon_evaporation_queue;

  // In case A <= 1 the fragment will not perform any nucleon emission
  auto initial_fragment_ptr = std::make_unique<G4Fragment>(fragment);
  if (neutron_decay_condition_(fragment)) {
    ApplyPureNeutronDecay(std::move(initial_fragment_ptr), results);
  } else if (IsStable(fragment, nist)) {
    results.push_back(initial_fragment_ptr.release());
  } else {
    if (multi_fragmentation_condition_(fragment)) {
      ApplyMultiFragmentation(std::move(initial_fragment_ptr), results, evaporation_queue);
    } else {
      evaporation_queue.emplace(std::move(initial_fragment_ptr));
    }

    for (size_t iteration_count = 0; !evaporation_queue.empty(); ++iteration_count) {
      auto fragment_ptr = std::move(evaporation_queue.front());
      evaporation_queue.pop();

      // infinite loop check
      if (iteration_count == kEvaporationIterationThreshold) {
        EvaporationError(fragment, *fragment_ptr, iteration_count);
        return {};
        // process is dead
      }

      // NeutronDecay part
      if (neutron_decay_condition_(*fragment_ptr)) {
        ApplyPureNeutronDecay(std::move(fragment_ptr), results);
        continue;
      }

      // FermiBreakUp part
      if (fermi_condition_(*fragment_ptr)) {
        ApplyFermiBreakUp(std::move(fragment_ptr), results, photon_evaporation_queue);
        continue;
      }

      // Evaporation part
      if (evaporation_condition_(*fragment_ptr)) {
        ApplyEvaporation(std::move(fragment_ptr), results, evaporation_queue);
        continue;
      }

      throw std::runtime_error(error_no_model);
    }

    // Photon Evaporation part
    while (!photon_evaporation_queue.empty()) {
      auto fragment_ptr = std::move(photon_evaporation_queue.front());
      photon_evaporation_queue.pop();

      if (photon_evaporation_condition_(*fragment_ptr)) {
        ApplyPhotonEvaporation(std::move(fragment_ptr), results);
        continue;
      }

      throw std::runtime_error(error_no_model);
    }
  }

  // Geant4 randomly throws neutron fragments, here's a bandaid for it

  G4FragmentVector buffer;
  for (auto it = results.begin(); it != results.end();) {
    if (neutron_decay_condition_(**it)) {
      buffer.push_back(*it);
      it = results.erase(it);
    } else {
      ++it;
    }
  }

  for (const auto& fragment_to_decay : buffer) {
    ApplyPureNeutronDecay(std::unique_ptr<G4Fragment>(fragment_to_decay), results);
  }

  auto reaction_products = ConvertResults(results);

  return reaction_products;
}

void ExcitationHandler::NeutronDecay::BreakFragment(G4FragmentVector& results, const G4Fragment& fragment) {
  if (fragment.GetZ_asInt() != 0) {
    throw std::runtime_error("only Z = 0 particles can be decayed by NeutronDecay, but got: A = " +
                             std::to_string(fragment.GetA_asInt()) + ", Z = " + std::to_string(fragment.GetZ_asInt()));
  }

  if (fragment.GetA_asInt() == 1) {
    return;
  }

  auto masses = std::vector<G4double>(fragment.GetA_asInt(), CLHEP::neutron_mass_c2);
  auto momentum = fragment.GetMomentum();
  const auto mass_excess = momentum.m() - CLHEP::neutron_mass_c2 * fragment.GetA_asInt();
  if (mass_excess < 10. * CLHEP::eV) {
    momentum.setE(momentum.e() + 10. * CLHEP::eV - mass_excess);
  }

  const auto particles_momentum = phase_space_decay_.CalculateDecay(momentum, masses);
  if (particles_momentum.empty()) {
    std::stringstream ss;
    ss << "NeutronDecay is unable to break particle with "
       << "A = " << fragment.GetA_asInt() << ", Z = " << fragment.GetZ_asInt() << ", P = " << momentum;
    throw std::runtime_error(ss.str());
  }

  const auto boost = momentum.boostVector();
  for (const auto& particle_momentum : particles_momentum) {
    results.emplace_back(new G4Fragment(1, 0, CLHEP::boostOf(particle_momentum, boost)));
  }
}

std::unique_ptr<G4VMultiFragmentation> ExcitationHandler::DefaultMultiFragmentation() {
  return std::make_unique<G4StatMF>();
}

std::unique_ptr<G4VFermiBreakUp> ExcitationHandler::DefaultFermiBreakUp() {
  auto model = std::make_unique<FermiBreakUpWrapper>();
  model->Initialise();
  return model;
}

std::unique_ptr<G4VEvaporation> ExcitationHandler::DefaultEvaporation() {
  auto evaporation = std::make_unique<G4Evaporation>();
  return evaporation;
}

std::unique_ptr<G4VEvaporationChannel> ExcitationHandler::DefaultPhotonEvaporation() {
  return std::make_unique<G4PhotonEvaporation>();
}

std::unique_ptr<ExcitationHandler::NeutronDecay> ExcitationHandler::DefaultNeutronDecay() {
  return std::make_unique<NeutronDecay>();
}

ExcitationHandler::Condition ExcitationHandler::DefaultMultiFragmentationCondition() {
  return [](const G4Fragment& fragment) -> bool {
    constexpr G4int kMaxAtomicMass = 19;
    constexpr G4int kMaxCharge = 9;
    constexpr G4double kLowerBoundTransitionMf = 3 * CLHEP::MeV;
    constexpr G4double kUpperBoundTransitionMf = 5 * CLHEP::MeV;

    const auto atomic_mass = fragment.GetA_asInt();
    const auto charge = fragment.GetZ_asInt();

    if (atomic_mass < kMaxAtomicMass && charge < kMaxCharge) {
      return false;
    }

    const auto exitation_energy = fragment.GetExcitationEnergy();

    const auto scale = 1. / (2. * (kUpperBoundTransitionMf - kLowerBoundTransitionMf));
    const auto energy_offset = (kUpperBoundTransitionMf + kLowerBoundTransitionMf) / 2.;
    const auto transition_prob = 0.5 * std::tanh((exitation_energy / atomic_mass - energy_offset) / scale) + 0.5;

    const auto random = G4RandFlat::shoot();

    if (exitation_energy < kLowerBoundTransitionMf * atomic_mass) {
      return false;
    }

    if (random < transition_prob && exitation_energy < kUpperBoundTransitionMf * atomic_mass) {
      return true;
    }

    if (random > transition_prob && exitation_energy < kUpperBoundTransitionMf * atomic_mass) {
      return false;
    }

    if (exitation_energy > kUpperBoundTransitionMf * atomic_mass) {
      return true;
    }

    return false;
  };
}

ExcitationHandler::Condition ExcitationHandler::DefaultFermiBreakUpCondition() {
  return
      [](const G4Fragment& fragment) -> bool { return fragment.GetZ_asInt() < MAX_Z && fragment.GetA_asInt() < MAX_A; };
}

ExcitationHandler::Condition ExcitationHandler::DefaultEvaporationCondition() {
  return [](const G4Fragment&) { return true; };
}

ExcitationHandler::Condition ExcitationHandler::DefaultPhotonEvaporationCondition() {
  return [](const G4Fragment&) { return true; };
}

ExcitationHandler::Condition ExcitationHandler::DefaultNeutronDecayCondition() {
  return [](const G4Fragment& fragment) { return fragment.GetA_asInt() > 1 && fragment.GetZ_asInt() == 0; };
}

bool ExcitationHandler::IsGroundState(const G4Fragment& fragment) const {
  return fragment.GetExcitationEnergy() < stable_threshold_;
}

bool ExcitationHandler::IsStable(const G4Fragment& fragment, const G4NistManager* nist) const {
  return fragment.GetA_asInt() <= 1 ||
         (IsGroundState(fragment) && nist->GetIsotopeAbundance(fragment.GetZ_asInt(), fragment.GetA_asInt()) > 0);
}

void ExcitationHandler::ApplyMultiFragmentation(std::unique_ptr<G4Fragment>&& fragment, G4FragmentVector& results,
                                                FragmentQueue& next_stage) {
  auto fragments = std::unique_ptr<G4FragmentVector>(multi_fragmentation_model_->BreakItUp(*fragment));
  if (fragments == nullptr || fragments->size() <= 1) {
    ClearSingularResults(*fragments, fragment.get());
    next_stage.emplace(fragment.release());
    return;
  }

  GroupFragments(std::move(*fragments), results, next_stage);
}

void ExcitationHandler::ApplyFermiBreakUp(std::unique_ptr<G4Fragment>&& fragment, G4FragmentVector& results,
                                          FragmentQueue& next_stage) {
  G4FragmentVector fragments;
  fermi_break_up_model_->BreakFragment(&fragments, fragment.get());

  if (fragments.size() <= 1) {
    ClearSingularResults(fragments, fragment.get());
    next_stage.emplace(fragment.release());
    return;
  }

  GroupFragments(std::move(fragments), results, next_stage);
}

void ExcitationHandler::ApplyEvaporation(std::unique_ptr<G4Fragment>&& fragment, G4FragmentVector& results,
                                         FragmentQueue& next_stage) {
  G4FragmentVector fragments;
  evaporation_model_->BreakFragment(&fragments, fragment.get());

  if (!fragments.empty() && fragments.back() == fragment.get()) {
    fragments.pop_back();
  }

  if (fragments.empty()) {
    results.emplace_back(fragment.release());
    return;
  }

  fragments.push_back(fragment.release());
  GroupFragments(std::move(fragments), results, next_stage);
}

void ExcitationHandler::ApplyPhotonEvaporation(std::unique_ptr<G4Fragment>&& fragment, G4FragmentVector& results) {
  // photon de-excitation only for hot fragments
  if (!IsGroundState(*fragment)) {
    G4FragmentVector fragments;

    photon_evaporation_model_->BreakUpChain(&fragments, fragment.get());

    for (auto* fragment_ptr : fragments) {
      results.emplace_back(fragment_ptr);
    }
  }

  // primary fragment is kept
  results.emplace_back(fragment.release());
}

void ExcitationHandler::ApplyPureNeutronDecay(std::unique_ptr<G4Fragment>&& fragment, G4FragmentVector& results) {
  size_t old_size = results.size();
  neutron_decay_model_->BreakFragment(results, *fragment);

  if (old_size == results.size()) {
    results.emplace_back(fragment.release());
  }
}

void ExcitationHandler::GroupFragments(G4FragmentVector&& fragments, G4FragmentVector& results,
                                       FragmentQueue& next_stage) {
  auto* nist = G4NistManager::Instance();

  for (auto* fragment_ptr : fragments) {
    // gamma, p, n or stable nuclei
    if (neutron_decay_condition_(*fragment_ptr)) {
      ApplyPureNeutronDecay(std::unique_ptr<G4Fragment>(fragment_ptr), results);
    } else if (IsStable(*fragment_ptr, nist)) {
      results.emplace_back(fragment_ptr);
    } else {
      next_stage.emplace(fragment_ptr);
    }
  }
}

std::vector<G4ReactionProduct> ExcitationHandler::ConvertResults(const G4FragmentVector& results) {
  std::vector<G4ReactionProduct> reaction_products;
  reaction_products.reserve(results.size());
  auto* ion_table = G4ParticleTable::GetParticleTable()->GetIonTable();

  for (const auto& fragment_ptr : results) {
    auto* fragment_definition = SpecialParticleDefinition(*fragment_ptr);
    if (fragment_definition == nullptr) {
      auto excitation_energy = fragment_ptr->GetExcitationEnergy();
      auto level = fragment_ptr->GetFloatingLevelNumber();
      if (IsGroundState(*fragment_ptr)) {
        excitation_energy = 0;
        level = 0;
      }
      fragment_definition = ion_table->GetIon(fragment_ptr->GetZ_asInt(), fragment_ptr->GetA_asInt(), excitation_energy,
                                              G4Ions::FloatLevelBase(level));
    }
    // fragment wasn't found, ground state is created
    if (fragment_definition == nullptr) {
      fragment_definition = ion_table->GetIon(fragment_ptr->GetZ_asInt(), fragment_ptr->GetA_asInt(), 0, noFloat, 0);
      if (fragment_definition == nullptr) {
        throw std::runtime_error("ion table isn't created");
      }
      G4double ion_mass = fragment_definition->GetPDGMass();
      if (fragment_ptr->GetMomentum().e() <= ion_mass) {
        fragment_ptr->SetMomentum(G4LorentzVector(ion_mass));
      } else {
        auto momentum = fragment_ptr->GetMomentum();
        G4double momentum_modulus = std::sqrt((momentum.e() - ion_mass) * (momentum.e() + ion_mass));
        momentum.setVect(momentum.vect().unit() * momentum_modulus);
        fragment_ptr->SetMomentum(momentum);
      }
    }

    reaction_products.emplace_back(fragment_definition);
    reaction_products.back().SetMomentum(fragment_ptr->GetMomentum().vect());
    reaction_products.back().SetTotalEnergy((fragment_ptr->GetMomentum()).e());
    reaction_products.back().SetFormationTime(fragment_ptr->GetCreationTime());
  }

  return reaction_products;
}
