#ifndef TATAMI_HELP_H
#define TATAMI_HELP_H

#include "tatami/tatami.hpp"

// Oracle sub-class which just always fetches a compile-time constant index
template <int IndexValue> class ConstIndexOracle final : public tatami::Oracle<int> {
private:
  tatami::PredictionIndex length_;

public:
  ConstIndexOracle(const int length) : length_(sanisizer::cast<tatami::PredictionIndex>(length)) {}
  tatami::PredictionIndex total() const { return length_; }
  int get(const tatami::PredictionIndex i) const { return IndexValue; }
};

#endif