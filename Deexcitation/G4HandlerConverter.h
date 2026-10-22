#ifndef DEEXCITATION_G4HANDLERCONVERTER_H_
#define DEEXCITATION_G4HANDLERCONVERTER_H_

#include <COLA.hh>
#include <EventData.hh>

#include <memory>

class ExcitationHandler;

namespace cola {
  class G4HandlerConverter final : public cola::VConverter {
   public:
    explicit G4HandlerConverter(std::unique_ptr<ExcitationHandler>&& model);

    std::unique_ptr<cola::EventData> operator()(std::unique_ptr<cola::EventData>&& data) final;

   private:
    std::unique_ptr<ExcitationHandler> model_;
  };
}  // namespace cola

#endif  // DEEXCITATION_G4HANDLERCONVERTER_H_
