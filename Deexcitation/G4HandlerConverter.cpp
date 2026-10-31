// NOLINTNEXTLINE(misc-include-cleaner)
#include "Deexcitation/G4HandlerConverter.h"

#include "Deexcitation/handler/ExcitationHandler.h"

#include <EventData.hh>
#include <G4Fragment.hh>
#include <G4LorentzVector.hh>
#include <G4ReactionProduct.hh>
#include <G4Types.hh>

#include <cstdint>
#include <memory>
#include <utility>

namespace {

  G4Fragment ColaToG4(const cola::Particle& particle) {
    const auto [atomic_mass, charge] = particle.GetAZ();

    const G4LorentzVector momentum(particle.momentum.x, particle.momentum.y, particle.momentum.z, particle.momentum.e);
    return G4Fragment(static_cast<G4int>(atomic_mass), static_cast<G4int>(charge),
                      G4LorentzVector(0., 0., 0., momentum.mag()));
  }

  cola::Particle G4ToCola(const G4ReactionProduct& fragment) {
    return cola::Particle{
        .position = cola::LorentzVector{.e = 0., .x = 0., .y = 0., .z = 0.},
        .momentum =
            cola::LorentzVector{
                .e = fragment.GetTotalEnergy(),
                .x = fragment.GetMomentum().x(),
                .y = fragment.GetMomentum().y(),
                .z = fragment.GetMomentum().z(),
            },
        .pdg_code = cola::AZToPdg({static_cast<uint16_t>(fragment.GetDefinition()->GetAtomicMass()),
                                   static_cast<uint16_t>(fragment.GetDefinition()->GetAtomicNumber())}),
        .p_class = cola::ParticleClass::kProduced,
    };
  }

}  // namespace

namespace cola {

  G4HandlerConverter::G4HandlerConverter(std::unique_ptr<ExcitationHandler>&& model) : model_(std::move(model)) {}

  std::unique_ptr<EventData> G4HandlerConverter::operator()(std::unique_ptr<EventData>&& data) {
    EventParticles results;
    for (const auto& particle : data->particles) {
      const auto particle_class = particle.p_class;
      if (particle_class == ParticleClass::kSpectatorA || particle_class == ParticleClass::kSpectatorB) {
        if (particle.GetAZ().first == 1) {
          results.emplace_back(particle);
          results.back().position = {};
          continue;
        }
        auto model_result = model_->BreakItUp(ColaToG4(particle));
        const G4LorentzVector initial_momentum(particle.momentum.x, particle.momentum.y, particle.momentum.z,
                                               particle.momentum.e);
        const auto boost = initial_momentum.boostVector();
        for (auto& product : model_result) {
          G4LorentzVector momentum(product.GetMomentum(), product.GetTotalEnergy());
          momentum.boost(boost);
          product.SetMomentum(momentum.vect());
          product.SetTotalEnergy(momentum.e());
        }

        for (const auto& fragment : model_result) {
          results.emplace_back(G4ToCola(fragment));
          results.back().p_class = particle_class;
        }
      } else {
        results.push_back(particle);
      }
    }

    data->particles = std::move(results);
    return std::move(data);
  }

}  // namespace cola
