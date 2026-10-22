//
// Created by Artem Novikov on 19.07.2023.
//

#include "Deexcitation/handler/ExcitationHandler.h"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4NucleiProperties.hh>
#include <Randomize.hh>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

namespace {
  struct Fragments {
    G4FragmentVector values;
    ~Fragments() {
      for (auto* fragment : values) {
        delete fragment;
      }
    }
  };

  void ExpectMomentum(const G4LorentzVector& actual, const G4LorentzVector& expected) {
    const auto tolerance = 1.e-8 * std::max(1., std::abs(expected.e()));
    EXPECT_NEAR(actual.e(), expected.e(), tolerance);
    EXPECT_NEAR(actual.x(), expected.x(), tolerance);
    EXPECT_NEAR(actual.y(), expected.y(), tolerance);
    EXPECT_NEAR(actual.z(), expected.z(), tolerance);
  }

  void ExpectNeutronDecay(G4int mass, const G4LorentzVector& momentum) {
    CLHEP::HepRandom::setTheSeed(42);
    const G4Fragment source(mass, 0, momentum);
    ExcitationHandler::NeutronDecay model;
    Fragments products;
    model.BreakFragment(products.values, source);
    ASSERT_EQ(products.values.size(), static_cast<size_t>(mass));
    G4LorentzVector total;
    for (const auto* product : products.values) {
      EXPECT_EQ(product->GetA_asInt(), 1);
      EXPECT_EQ(product->GetZ_asInt(), 0);
      total += product->GetMomentum();
    }
    ExpectMomentum(total, momentum);
  }

  struct Configuration {
    G4int mass;
    G4int charge;
    G4double excitation;
  };

  class ConfigurationsFixture : public ::testing::TestWithParam<Configuration> {
   protected:
    void SetUp() override {
      // Disable SMM because of the known conservation bugs in Geant4 StatMF.
      model.SetMultiFragmentationCondition([](const G4Fragment&) { return false; });
    }

    G4Fragment MakeFragment(const G4ThreeVector& beta) const {
      const auto& configuration = GetParam();
      G4LorentzVector momentum(
          0., 0., 0.,
          G4NucleiProperties::GetNuclearMass(configuration.mass, configuration.charge) + configuration.excitation);
      momentum.boost(beta);
      return G4Fragment(configuration.mass, configuration.charge, momentum);
    }

    ExcitationHandler model;
  };

  // Excitations measured at the handler input in event 1, seed 27182.
  INSTANTIATE_TEST_SUITE_P(ConservationTests, ConfigurationsFixture,
                           ::testing::Values(Configuration{29, 14, 125.61367721807255},
                                             Configuration{50, 20, 216.57522910003172},
                                             Configuration{77, 37, 351.32331353027257}));
}  // namespace

TEST(NeutronDecayConservation, RestFrame) {
  ExpectNeutronDecay(2, G4LorentzVector(0., 0., 0., 2. * CLHEP::neutron_mass_c2 + 10.));
}

TEST(NeutronDecayConservation, MovingParent) {
  G4LorentzVector momentum(0., 0., 0., 3. * CLHEP::neutron_mass_c2 + 20.);
  momentum.boost(G4ThreeVector(0.3, -0.2, 0.7));
  ExpectNeutronDecay(3, momentum);
}

TEST(NeutronDecayConservation, RecordedPbPbFragment) {
  // GMST output in event 482, seed 27182, sqrt(sNN) = 5020 GeV.
  ExpectNeutronDecay(3,
                     G4LorentzVector(-2559.5258298818876, 3452.9688877832536, -17767573.153734092, 17767573.899202898));
}

TEST_P(ConfigurationsFixture, MassAndChargeConservation) {
  for (const auto beta : {G4ThreeVector(), G4ThreeVector(0.3, -0.2, 0.7)}) {
    SCOPED_TRACE(::testing::Message() << "beta=" << beta);
    CLHEP::HepRandom::setTheSeed(42);
    const auto products = model.BreakItUp(MakeFragment(beta));
    ASSERT_FALSE(products.empty());
    G4int total_mass = 0;
    G4int total_charge = 0;
    for (const auto& product : products) {
      total_mass += product.GetDefinition()->GetBaryonNumber();
      total_charge += static_cast<G4int>(std::lround(product.GetDefinition()->GetPDGCharge()));
    }
    EXPECT_EQ(total_mass, GetParam().mass);
    EXPECT_EQ(total_charge, GetParam().charge);
  }
}

TEST_P(ConfigurationsFixture, Vector4Conservation) {
  for (const auto beta : {G4ThreeVector(), G4ThreeVector(0.3, -0.2, 0.7)}) {
    SCOPED_TRACE(::testing::Message() << "beta=" << beta);
    CLHEP::HepRandom::setTheSeed(42);
    const auto fragment = MakeFragment(beta);
    const auto products = model.BreakItUp(fragment);
    ASSERT_FALSE(products.empty());
    G4LorentzVector total_momentum;
    for (const auto& product : products) {
      total_momentum += G4LorentzVector(product.GetMomentum(), product.GetTotalEnergy());
    }
    ExpectMomentum(total_momentum, fragment.GetMomentum());
  }
}
