#ifndef DEEXCITATION_G4HANDLERFACTORY_H_
#define DEEXCITATION_G4HANDLERFACTORY_H_

#include <COLA.hh>

#include <memory>
#include <string>
#include <unordered_map>

namespace cola {

  class G4HandlerFactory final : public VConverterFactory {
   public:
    std::unique_ptr<VFilter> Create(const std::unordered_map<std::string, std::string>& meta_data) override;

    const std::string& GetFilterName() const override {
      static const std::string name{"G4DeexcitationHandler"};
      return name;
    }
  };

}  // namespace cola

#endif  // DEEXCITATION_G4HANDLERFACTORY_H_
