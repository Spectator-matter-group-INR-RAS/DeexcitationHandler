#include "Deexcitation/G4HandlerFactory.h"

#include "Deexcitation/G4HandlerConverter.h"
#include "Deexcitation/handler/ExcitationHandler.h"

#include <CLHEP/Units/SystemOfUnits.h>
#include <COLA.hh>
#include <G4FermiDataTypes.hh>
#include <G4Fragment.hh>
#include <G4Types.hh>
#include <Randomize.hh>

#include <cmath>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>

namespace {

  double StodWithFactor(const std::string& value) {
    const double num = std::stod(value);
    if (value.ends_with("eV")) {
      return num * CLHEP::eV;
    }
    if (value.ends_with("MeV")) {
      return num * CLHEP::MeV;
    }
    if (value.ends_with("keV")) {
      return num * CLHEP::keV;
    }
    if (value.ends_with("GeV")) {
      return num * CLHEP::GeV;
    }
    return num;
  }

  struct Config {
    explicit Config(const std::unordered_map<std::string, std::string>& params) {
      if (const auto it = params.find("A"); it != params.end()) {
        atomic_mass = std::stoi(it->second);
      }
      if (const auto it = params.find("Z"); it != params.end()) {
        charge = std::stoi(it->second);
      }
      if (const auto it = params.find("lowerMfThreshold"); it != params.end()) {
        lower_mf_threshold = StodWithFactor(it->second);
      }
      if (const auto it = params.find("upperMfThreshold"); it != params.end()) {
        upper_mf_threshold = StodWithFactor(it->second);
      }
      if (const auto it = params.find("stableThreshold"); it != params.end()) {
        stable_threshold = StodWithFactor(it->second);
      }
    }

    std::optional<int> atomic_mass;
    std::optional<int> charge;
    std::optional<double> stable_threshold;
    std::optional<double> lower_mf_threshold;
    std::optional<double> upper_mf_threshold;
  };

}  // namespace

namespace cola {

  std::unique_ptr<VFilter> G4HandlerFactory::Create(const std::unordered_map<std::string, std::string>& meta_data) {
    const Config config(meta_data);

    auto model = std::make_unique<ExcitationHandler>();

    if (config.stable_threshold.has_value()) {
      model->SetStableThreshold(*config.stable_threshold);
    }

    model->SetFermiBreakUpCondition(
        [max_a = config.atomic_mass.value_or(static_cast<int>(MAX_A)),
         max_z = config.charge.value_or(static_cast<int>(MAX_Z))](const G4Fragment& fragment) -> bool {
          return fragment.GetZ_asInt() < max_z && fragment.GetA_asInt() < max_a;
        });

    model->SetMultiFragmentationCondition(
        [max_a = config.atomic_mass.value_or(static_cast<int>(MAX_A)),
         max_z = config.charge.value_or(static_cast<int>(MAX_Z)),
         lower_bound_transition_mf = config.lower_mf_threshold.value_or(3 * CLHEP::MeV),
         upper_bound_transition_mf =
             config.upper_mf_threshold.value_or(5 * CLHEP::MeV)](const G4Fragment& fragment) -> bool {
          const auto a = fragment.GetA_asInt();
          const auto z = fragment.GetZ_asInt();
          const auto ex = fragment.GetExcitationEnergy();
          if (a < max_a && z < max_z) {
            return false;
          }
          const G4double a_e = 1 / (2. * (upper_bound_transition_mf - lower_bound_transition_mf));
          const G4double e0 = (upper_bound_transition_mf + lower_bound_transition_mf) / 2.;
          const G4double weight = G4RandFlat::shoot();
          const G4double trans_f = 0.5 * std::tanh((ex / a - e0) / a_e) + 0.5;

          if (ex < lower_bound_transition_mf * a) {
            return false;
          }
          if (weight < trans_f && ex < upper_bound_transition_mf * a) {
            return true;
          }
          if (weight > trans_f && ex < upper_bound_transition_mf * a) {
            return false;
          }
          if (ex > upper_bound_transition_mf * a) {
            return true;
          }
          return false;
        });

    return std::make_unique<G4HandlerConverter>(std::move(model));
  }

}  // namespace cola
