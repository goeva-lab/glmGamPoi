#ifndef TATAMI_HELP_H
#define TATAMI_HELP_H

#include "tatami/tatami.hpp"

// Oracle sub-class which just always fetches a compile-time constant index
template <int IndexValue> class ConstIndexOracle final : public tatami::Oracle<int> {
public:
  ConstIndexOracle(const int length) : my_length(sanisizer::cast<tatami::PredictionIndex>(length)) {}
  tatami::PredictionIndex total() const { return my_length; }
  int get(const tatami::PredictionIndex i) const { return IndexValue; }

private:
  tatami::PredictionIndex my_length;
};

#endif