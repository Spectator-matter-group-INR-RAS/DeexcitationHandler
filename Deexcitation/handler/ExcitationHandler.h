#ifndef DEEXCITATION_HANDLER_EXCITATIONHANDLER_H_
#define DEEXCITATION_HANDLER_EXCITATIONHANDLER_H_

#include <G4FermiPhaseDecay.hh>
#include <G4Fragment.hh>
#include <G4NistManager.hh>
#include <G4ReactionProduct.hh>
#include <G4VEvaporation.hh>
#include <G4VFermiBreakUp.hh>
#include <G4VMultiFragmentation.hh>

#include <functional>
#include <memory>
#include <queue>
#include <utility>
#include <vector>

class ExcitationHandler {
 private:
  using FragmentQueue = std::queue<std::unique_ptr<G4Fragment>>;

 public:
  class NeutronDecay {
   public:
    NeutronDecay() = default;

    void BreakFragment(G4FragmentVector& results, const G4Fragment& fragment);

   private:
    G4FermiPhaseDecay phase_space_decay_;
  };

  using Condition = std::function<bool(const G4Fragment&)>;

  ExcitationHandler();

  ExcitationHandler(const ExcitationHandler&) = delete;

  ExcitationHandler(ExcitationHandler&&) = default;

  ~ExcitationHandler();

  ExcitationHandler& operator=(const ExcitationHandler&) = delete;

  ExcitationHandler& operator=(ExcitationHandler&&) = default;

  std::vector<G4ReactionProduct> BreakItUp(const G4Fragment& fragment);

  // parameters setters
  ExcitationHandler& SetMultiFragmentation(
      std::unique_ptr<G4VMultiFragmentation>&& model = DefaultMultiFragmentation()) {
    multi_fragmentation_model_ = std::move(model);
    return *this;
  }

  ExcitationHandler& SetFermiBreakUp(std::unique_ptr<G4VFermiBreakUp>&& model = DefaultFermiBreakUp(),
                                     bool set_models = true) {
    fermi_break_up_model_ = std::move(model);
    if (set_models) {
      evaporation_model_->SetFermiBreakUp(fermi_break_up_model_.get());
    }
    return *this;
  }

  ExcitationHandler& SetEvaporation(std::unique_ptr<G4VEvaporation>&& model = DefaultEvaporation(),
                                    bool set_models = true) {
    evaporation_model_ = std::move(model);
    if (set_models) {
      evaporation_model_->SetPhotonEvaporation(photon_evaporation_model_.get());
      evaporation_model_->SetFermiBreakUp(fermi_break_up_model_.get());
    }
    return *this;
  }

  ExcitationHandler& SetPhotonEvaporation(std::unique_ptr<G4VEvaporationChannel>&& model = DefaultPhotonEvaporation(),
                                          bool set_models = true) {
    photon_evaporation_model_ = std::move(model);
    if (set_models) {
      evaporation_model_->SetPhotonEvaporation(photon_evaporation_model_.get());
    }
    return *this;
  }

  ExcitationHandler& SetNeutronDecay(std::unique_ptr<NeutronDecay>&& model = DefaultNeutronDecay()) {
    neutron_decay_model_ = std::move(model);
    return *this;
  }

  template <class F>
  ExcitationHandler& SetMultiFragmentationCondition(F&& func) {
    multi_fragmentation_condition_ = std::forward<F>(func);
    return *this;
  }

  // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
  ExcitationHandler& SetMultiFragmentationCondition() {
    return SetMultiFragmentationCondition(DefaultMultiFragmentationCondition());
  }

  template <class F>
  ExcitationHandler& SetFermiBreakUpCondition(F&& func) {
    fermi_condition_ = std::forward<F>(func);
    return *this;
  }

  // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
  ExcitationHandler& SetFermiBreakUpCondition() { return SetFermiBreakUpCondition(DefaultFermiBreakUpCondition()); }

  template <class F>
  ExcitationHandler& SetEvaporationCondition(F&& func) {
    evaporation_condition_ = std::forward<F>(func);
    return *this;
  }

  // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
  ExcitationHandler& SetEvaporationCondition() { return SetEvaporationCondition(DefaultEvaporationCondition()); }

  template <class F>
  ExcitationHandler& SetPhotonEvaporationCondition(F&& func) {
    photon_evaporation_condition_ = std::forward<F>(func);
    return *this;
  }

  // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
  ExcitationHandler& SetPhotonEvaporationCondition() {
    return SetPhotonEvaporationCondition(DefaultPhotonEvaporationCondition());
  }

  template <class F>
  ExcitationHandler& SetNeutronDecayCondition(F&& func) {
    neutron_decay_condition_ = std::forward<F>(func);
    return *this;
  }

  // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
  ExcitationHandler& SetNeutronDecayCondition() { return SetNeutronDecayCondition(DefaultNeutronDecayCondition()); }

  ExcitationHandler& SetStableThreshold(double threshold) {
    stable_threshold_ = threshold;
    return *this;
  }

  // parameters getters
  std::unique_ptr<NeutronDecay>& GetNeutronDecay() { return neutron_decay_model_; }

  const std::unique_ptr<NeutronDecay>& GetNeutronDecay() const { return neutron_decay_model_; }

  std::unique_ptr<G4VMultiFragmentation>& GetMultiFragmentation() { return multi_fragmentation_model_; }

  const std::unique_ptr<G4VMultiFragmentation>& GetMultiFragmentation() const { return multi_fragmentation_model_; }

  std::unique_ptr<G4VFermiBreakUp>& GetFermiBreakUp() { return fermi_break_up_model_; }

  const std::unique_ptr<G4VFermiBreakUp>& GetFermiBreakUp() const { return fermi_break_up_model_; }

  std::unique_ptr<G4VEvaporation>& GetEvaporation() { return evaporation_model_; }

  const std::unique_ptr<G4VEvaporation>& GetEvaporation() const { return evaporation_model_; }

  Condition& GetMultiFragmentationCondition() { return multi_fragmentation_condition_; }

  const Condition& GetMultiFragmentationCondition() const { return multi_fragmentation_condition_; }

  Condition& GetFermiBreakUpCondition() { return fermi_condition_; }

  const Condition& GetFermiBreakUpCondition() const { return fermi_condition_; }

  Condition& GetEvaporationCondition() { return evaporation_condition_; }

  const Condition& GetEvaporationCondition() const { return evaporation_condition_; }

  Condition& GetPhotonEvaporationCondition() { return photon_evaporation_condition_; }

  const Condition& GetPhotonEvaporationCondition() const { return photon_evaporation_condition_; }

  Condition& GetNeutronDecayCondition() { return neutron_decay_condition_; }

  const Condition& GetNeutronDecayCondition() const { return neutron_decay_condition_; }

  double GetStableThreshold() const { return stable_threshold_; }

 protected:
  // default models and conditions
  static std::unique_ptr<G4VMultiFragmentation> DefaultMultiFragmentation();

  static std::unique_ptr<G4VFermiBreakUp> DefaultFermiBreakUp();

  static std::unique_ptr<G4VEvaporation> DefaultEvaporation();

  static std::unique_ptr<G4VEvaporationChannel> DefaultPhotonEvaporation();

  static std::unique_ptr<NeutronDecay> DefaultNeutronDecay();

  static Condition DefaultMultiFragmentationCondition();

  static Condition DefaultFermiBreakUpCondition();

  static Condition DefaultEvaporationCondition();

  static Condition DefaultPhotonEvaporationCondition();

  static Condition DefaultNeutronDecayCondition();

  bool IsGroundState(const G4Fragment& fragment) const;

  bool IsStable(const G4Fragment& fragment, const G4NistManager* nist) const;

  void ApplyMultiFragmentation(std::unique_ptr<G4Fragment>&& fragment, G4FragmentVector& results,
                               FragmentQueue& next_stage);

  void ApplyFermiBreakUp(std::unique_ptr<G4Fragment>&& fragment, G4FragmentVector& results, FragmentQueue& next_stage);

  void ApplyEvaporation(std::unique_ptr<G4Fragment>&& fragment, G4FragmentVector& results, FragmentQueue& next_stage);

  void ApplyPhotonEvaporation(std::unique_ptr<G4Fragment>&& fragment, G4FragmentVector& results);

  void ApplyPureNeutronDecay(std::unique_ptr<G4Fragment>&& fragment, G4FragmentVector& results);

  void GroupFragments(G4FragmentVector&& fragments, G4FragmentVector& results, FragmentQueue& next_stage);

  std::vector<G4ReactionProduct> ConvertResults(const G4FragmentVector& results);

  std::unique_ptr<G4VMultiFragmentation> multi_fragmentation_model_;
  std::unique_ptr<G4VFermiBreakUp> fermi_break_up_model_;
  std::unique_ptr<G4VEvaporation> evaporation_model_;
  std::unique_ptr<G4VEvaporationChannel> photon_evaporation_model_;
  std::unique_ptr<NeutronDecay> neutron_decay_model_;

  Condition multi_fragmentation_condition_;
  Condition fermi_condition_;
  Condition photon_evaporation_condition_;
  Condition evaporation_condition_;
  Condition neutron_decay_condition_;

  double stable_threshold_ = 0.;
};

#endif  // DEEXCITATION_HANDLER_EXCITATIONHANDLER_H_
