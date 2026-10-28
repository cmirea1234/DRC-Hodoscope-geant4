#ifndef HODOSCOPE_HIT_CLASSIFIER_HH
#define HODOSCOPE_HIT_CLASSIFIER_HH

#include <array>
#include <cstddef>
#include <limits>

namespace hodoscope::hit {

inline constexpr std::size_t kChannelsPerPlane = 32;
inline constexpr std::size_t kPlaneCount = 4;
inline constexpr std::size_t kChannelCount =
    kChannelsPerPlane * kPlaneCount;

// Run15105 raw-array conventions. All ranges are half-open.
inline constexpr std::size_t kPedestalFirst = 1;
inline constexpr std::size_t kPedestalLast = 101;
inline constexpr std::size_t kPeakFirst = 200;
inline constexpr std::size_t kPeakLast = 600;
inline constexpr std::size_t kPeakSamples = kPeakLast - kPeakFirst;

enum Plane : std::size_t {
  kRoundX = 0,   // MID 1
  kRoundY = 1,   // MID 9
  kSquareX = 2,  // MID 15
  kSquareY = 3   // MID 13
};

struct PlaneDecision {
  std::array<bool, kPlaneCount> plane{};
  bool round = false;
  bool square = false;
};

struct Decision : PlaneDecision {
  bool trigger1 = false;
  bool trigger2 = false;
  bool eligible = false;
};

template <class Values, class Thresholds>
PlaneDecision ClassifyPlanes(const Values& values,
                             const Thresholds& thresholds) {
  PlaneDecision result;
  for (std::size_t plane = 0; plane < kPlaneCount; ++plane) {
    const auto first = plane * kChannelsPerPlane;
    const auto last = first + kChannelsPerPlane;
    for (auto channel = first; channel < last; ++channel) {
      // The transferred Park/Yewon classifier uses strict greater-than.
      if (values[channel] > thresholds[channel]) {
        result.plane[plane] = true;
        break;
      }
    }
  }
  result.round = result.plane[kRoundX] && result.plane[kRoundY];
  result.square = result.plane[kSquareX] && result.plane[kSquareY];
  return result;
}

template <class Values, class Thresholds, class Scalar>
Decision Classify(Scalar trigger1Value, Scalar trigger2Value,
                  Scalar triggerThreshold, const Values& values,
                  const Thresholds& thresholds) {
  const auto planes = ClassifyPlanes(values, thresholds);
  Decision result;
  result.plane = planes.plane;
  result.round = planes.round;
  result.square = planes.square;
  // Trigger scintillator hits use inclusive energy-deposit thresholds.
  result.trigger1 = trigger1Value >= triggerThreshold;
  result.trigger2 = trigger2Value >= triggerThreshold;
  result.eligible = result.trigger1 && result.trigger2;
  return result;
}

template <class Samples>
double MeanHalfOpen(const Samples& samples, std::size_t first,
                    std::size_t last) {
  double sum = 0.;
  for (auto index = first; index < last; ++index) sum += samples[index];
  return sum / static_cast<double>(last - first);
}

template <class Samples>
double NegativePulsePeak(const Samples& samples, std::size_t pedestalFirst,
                         std::size_t pedestalLast, std::size_t peakFirst,
                         std::size_t peakLast) {
  const double pedestal =
      MeanHalfOpen(samples, pedestalFirst, pedestalLast);
  double minimum = std::numeric_limits<double>::infinity();
  for (auto index = peakFirst; index < peakLast; ++index) {
    const double value = samples[index];
    if (value < minimum) minimum = value;
  }
  return pedestal - minimum;
}

template <class Samples>
double PositivePulsePeak(const Samples& samples, std::size_t first,
                         std::size_t last) {
  double maximum = -std::numeric_limits<double>::infinity();
  for (auto index = first; index < last; ++index) {
    const double value = samples[index];
    if (value > maximum) maximum = value;
  }
  return maximum;
}

}  // namespace hodoscope::hit

#endif
