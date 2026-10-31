#include "Deexcitation/G4HandlerConverter.h"
#include "Deexcitation/G4HandlerFactory.h"
#include "Deexcitation/handler/ExcitationHandler.h"

#include <CLHEP/Units/SystemOfUnits.h>
#include <COLA.hh>
#include <G4LorentzVector.hh>
#include <G4NucleiProperties.hh>
#include <Randomize.hh>
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <sstream>

namespace {

  constexpr auto kPipelineConfigXml = R"(<?xml version="1.0" encoding="UTF-8" ?>
<program>
    <generator name="generator"/>
    <converter name="G4DeexcitationHandler"/>
    <writer name="writer"/>
</program>
)";

  class TestGenerator final : public cola::VGenerator {
   public:
    explicit TestGenerator(cola::EventParticles particles) : particles_(std::move(particles)) {}

    std::unique_ptr<cola::EventData> operator()() override {
      auto data = std::make_unique<cola::EventData>();
      data->particles = particles_;
      return data;
    }

   private:
    cola::EventParticles particles_;
  };

  class TestGeneratorFactory final : public cola::VGeneratorFactory {
   public:
    explicit TestGeneratorFactory(cola::EventParticles particles) : particles_(std::move(particles)) {}

    std::unique_ptr<cola::VFilter> Create(const std::unordered_map<std::string, std::string>& /*meta_data*/) override {
      return std::make_unique<TestGenerator>(particles_);
    }

    const std::string& GetFilterName() const override {
      static const std::string name{"generator"};
      return name;
    }

   private:
    cola::EventParticles particles_;
  };

  class TestWriter final : public cola::VWriter {
   public:
    explicit TestWriter(std::shared_ptr<std::vector<std::unique_ptr<cola::EventData>>> sink) : sink_(std::move(sink)) {}

    void operator()(std::unique_ptr<cola::EventData>&& data) override { sink_->emplace_back(std::move(data)); }

   private:
    std::shared_ptr<std::vector<std::unique_ptr<cola::EventData>>> sink_;
  };

  class TestWriterFactory final : public cola::VWriterFactory {
   public:
    explicit TestWriterFactory(std::shared_ptr<std::vector<std::unique_ptr<cola::EventData>>> sink)
        : sink_(std::move(sink)) {}

    std::unique_ptr<cola::VFilter> Create(const std::unordered_map<std::string, std::string>& /*meta_data*/) override {
      return std::make_unique<TestWriter>(sink_);
    }

    const std::string& GetFilterName() const override {
      static const std::string name{"writer"};
      return name;
    }

   private:
    std::shared_ptr<std::vector<std::unique_ptr<cola::EventData>>> sink_;
  };

  void RegisterTestPipeline(cola::MetaProcessor& meta_processor,
                            std::shared_ptr<std::vector<std::unique_ptr<cola::EventData>>> sink,
                            cola::EventParticles particles) {
    meta_processor.Register(std::make_unique<TestGeneratorFactory>(std::move(particles)));
    meta_processor.Register(std::make_unique<cola::G4HandlerFactory>());
    meta_processor.Register(std::make_unique<TestWriterFactory>(std::move(sink)));
  }

  cola::Particle MakeSpectator(std::uint16_t mass, std::uint16_t charge, const G4LorentzVector& momentum) {
    return cola::Particle{
        .position = {},
        .momentum = {.e = momentum.e(), .x = momentum.x(), .y = momentum.y(), .z = momentum.z()},
        .pdg_code = cola::AZToPdg({mass, charge}),
        .p_class = cola::ParticleClass::kSpectatorA,
    };
  }

  cola::EventParticles ConvertSpectator(const cola::Particle& particle) {
    auto model = std::make_unique<ExcitationHandler>();
    model->SetMultiFragmentationCondition([](const G4Fragment&) { return false; });
    cola::G4HandlerConverter converter(std::move(model));
    auto event = std::make_unique<cola::EventData>();
    event->particles.push_back(particle);
    CLHEP::HepRandom::setTheSeed(42);
    return converter(std::move(event))->particles;
  }

}  // namespace

TEST(TestModule, TestG4Handler) {
  auto sink = std::make_shared<std::vector<std::unique_ptr<cola::EventData>>>();

  {
    const cola::Particle particle{
        .position = cola::LorentzVector{},
        .momentum =
            cola::LorentzVector{
                .e = std::sqrt(3 * std::pow(100 * CLHEP::MeV, 2) + std::pow(4 * 938 * CLHEP::MeV, 2)),
                .x = 100 * CLHEP::MeV,
                .y = 100 * CLHEP::MeV,
                .z = 100 * CLHEP::MeV,
            },
        .pdg_code = cola::AZToPdg({4, 2}),
        .p_class = cola::ParticleClass::kSpectatorA,
    };
    sink->clear();
    cola::MetaProcessor mp;
    RegisterTestPipeline(mp, sink, {particle});
    std::istringstream xml(kPipelineConfigXml);
    cola::ColaRunManager manager(mp.Parse(xml));
    manager.Run(1);

    ASSERT_EQ(sink->size(), 1u);
    EXPECT_GE((*sink)[0]->particles.size(), 1u);
  }

  {
    const cola::Particle particle{
        .position = cola::LorentzVector{},
        .momentum =
            cola::LorentzVector{
                .e = std::sqrt(3 * std::pow(100 * CLHEP::MeV, 2) + std::pow(5 * 938 * CLHEP::MeV, 2)),
                .x = 100 * CLHEP::MeV,
                .y = 100 * CLHEP::MeV,
                .z = 100 * CLHEP::MeV,
            },
        .pdg_code = cola::AZToPdg({5, 3}),
        .p_class = cola::ParticleClass::kSpectatorA,
    };
    sink->clear();
    cola::MetaProcessor mp;
    RegisterTestPipeline(mp, sink, {particle});
    std::istringstream xml(kPipelineConfigXml);
    cola::ColaRunManager manager(mp.Parse(xml));
    manager.Run(1);

    ASSERT_EQ(sink->size(), 1u);
    EXPECT_GE((*sink)[0]->particles.size(), 1u);
  }
}

TEST(ModuleExport, LoadCOLAModuleExposesSingleG4DeexcitationFactory) {
  auto module = std::unique_ptr<cola::VModule>(LoadCOLAModule());
  ASSERT_NE(module, nullptr);

  const auto filters = module->GetModuleFilters();
  ASSERT_EQ(filters.size(), 1u);
  ASSERT_TRUE(filters.contains("G4DeexcitationHandler"));

  const auto* factory = filters.at("G4DeexcitationHandler").get();
  ASSERT_NE(factory, nullptr);
  EXPECT_EQ(factory->GetFilterName(), "G4DeexcitationHandler");
  EXPECT_EQ(factory->GetFilterType(), cola::FilterType::kConverter);
}

TEST(TestModule, BoostedDeexcitationMatchesRestFrame) {
  struct Configuration {
    std::uint16_t mass;
    std::uint16_t charge;
    double excitation;
  };
  const std::array configurations{Configuration{29, 14, 125.61367721807255}, Configuration{50, 20, 216.57522910003172},
                                  Configuration{77, 37, 351.32331353027257}};
  constexpr double gamma = 2660.;
  const G4ThreeVector boost(0., 0., std::sqrt(1. - 1. / (gamma * gamma)));
  for (const auto& configuration : configurations) {
    SCOPED_TRACE(::testing::Message() << "A=" << configuration.mass << ", Z=" << configuration.charge);
    G4LorentzVector initial(
        0., 0., 0.,
        G4NucleiProperties::GetNuclearMass(configuration.mass, configuration.charge) + configuration.excitation);
    initial.boost(boost);
    const auto rest_particle =
        MakeSpectator(configuration.mass, configuration.charge, G4LorentzVector(0., 0., 0., initial.mag()));
    const auto lab_particle = MakeSpectator(configuration.mass, configuration.charge, initial);

    ::testing::internal::CaptureStdout();
    const auto rest_products = ConvertSpectator(rest_particle);
    const auto lab_products = ConvertSpectator(lab_particle);
    const auto messages = ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(messages.find("Fermi001"), std::string::npos);
    EXPECT_EQ(messages.find("G4Fragment::CalculateExcitationEnergy"), std::string::npos);
    ASSERT_FALSE(lab_products.empty());
    ASSERT_EQ(lab_products.size(), rest_products.size());

    G4LorentzVector total;
    int total_mass = 0;
    int total_charge = 0;
    for (std::size_t i = 0; i < lab_products.size(); ++i) {
      const auto& actual = lab_products[i];
      const auto& rest = rest_products[i];
      G4LorentzVector expected(rest.momentum.x, rest.momentum.y, rest.momentum.z, rest.momentum.e);
      expected.boost(initial.boostVector());
      const auto tolerance = 1.e-8 * std::max(1., std::abs(expected.e()));
      EXPECT_EQ(actual.pdg_code, rest.pdg_code);
      EXPECT_EQ(actual.p_class, lab_particle.p_class);
      EXPECT_NEAR(actual.momentum.e, expected.e(), tolerance);
      EXPECT_NEAR(actual.momentum.x, expected.x(), tolerance);
      EXPECT_NEAR(actual.momentum.y, expected.y(), tolerance);
      EXPECT_NEAR(actual.momentum.z, expected.z(), tolerance);
      const auto [mass, charge] = actual.GetAZ();
      total_mass += mass;
      total_charge += charge;
      total += G4LorentzVector(actual.momentum.x, actual.momentum.y, actual.momentum.z, actual.momentum.e);
    }
    EXPECT_EQ(total_mass, configuration.mass);
    EXPECT_EQ(total_charge, configuration.charge);
    const auto tolerance = 1.e-8 * initial.e();
    EXPECT_NEAR(total.e(), initial.e(), tolerance);
    EXPECT_NEAR(total.x(), initial.x(), tolerance);
    EXPECT_NEAR(total.y(), initial.y(), tolerance);
    EXPECT_NEAR(total.z(), initial.z(), tolerance);
  }
}

TEST(TestModule, SingleNucleonsPreserveMomentumWithoutWarnings) {
  // Recorded neutron: Pb-Pb, 5020 GeV, seed 42, event 3.
  const auto neutron = MakeSpectator(
      1, 0, G4LorentzVector(1354.4213866527077, -688.7791168955431, 8940463.5233726669, 8940463.7018677164));
  G4LorentzVector proton_momentum(0., 0., 0., G4NucleiProperties::GetNuclearMass(1, 1));
  proton_momentum.boost(G4ThreeVector(0., 0., std::sqrt(1. - 1. / (2660. * 2660.))));
  auto proton = MakeSpectator(1, 1, proton_momentum);
  proton.p_class = cola::ParticleClass::kSpectatorB;
  for (auto particle : {neutron, proton}) {
    SCOPED_TRACE(::testing::Message() << "PDG=" << particle.pdg_code);
    particle.position = {.e = 4., .x = 1., .y = 2., .z = 3.};
    ::testing::internal::CaptureStdout();
    const auto products = ConvertSpectator(particle);
    const auto messages = ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(messages.find("G4Fragment::CalculateExcitationEnergy"), std::string::npos);
    ASSERT_EQ(products.size(), 1u);
    const auto& product = products.front();
    EXPECT_EQ(product.pdg_code, particle.pdg_code);
    EXPECT_EQ(product.p_class, particle.p_class);
    EXPECT_EQ(product.momentum.e, particle.momentum.e);
    EXPECT_EQ(product.momentum.x, particle.momentum.x);
    EXPECT_EQ(product.momentum.y, particle.momentum.y);
    EXPECT_EQ(product.momentum.z, particle.momentum.z);
    EXPECT_EQ(product.position.e, 0.);
    EXPECT_EQ(product.position.x, 0.);
    EXPECT_EQ(product.position.y, 0.);
    EXPECT_EQ(product.position.z, 0.);
  }
}
