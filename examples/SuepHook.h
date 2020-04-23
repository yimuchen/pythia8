#ifndef SuepHook_h
#define SuepHook_h

#include "Pythia8/Pythia.h"

#include "suep_shower.h"

#include <memory>

class SuepHook : public Pythia8::UserHooks {
public:
  SuepHook(int idMediator, int idDark, float temperature);
  ~SuepHook() {}

  bool initAfterBeams() override;

  bool canVetoProcessLevel() override { return true; }
  bool doVetoProcessLevel(Pythia8::Event& event) override;

protected:
  int idMediator_, idDark_;
  float temperature_, mMediator_, mDark_;
  std::unique_ptr<Suep_shower> suep_shower_;
};

#endif
