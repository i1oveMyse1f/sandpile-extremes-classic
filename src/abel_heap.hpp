//#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx")
#pragma GCC optimize 03
//#pragma GCC optimize("unroll-loops")

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

std::vector<int> read_sequence(std::ifstream &in) {
  std::vector<int> sequence;
  int s;
  while (in >> s) {
    sequence.push_back(s);
  }
  return sequence;
}

std::pair<std::vector<double>, std::vector<double>>
read_conditional_probability(std::ifstream &in) {
  double y, prob;

  string bin_middles_str, conditional_probability_str;
  std::getline(in, bin_middles_str);
  std::getline(in, conditional_probability_str);

  stringstream bin_middles_in(bin_middles_str);
  stringstream conditional_probability_in(bin_middles_str);
  vector<double> bin_middle, conditional_probability;

  while (bin_middles_in >> y) {
    bin_middle.push_back(y);
  }

  while (conditional_probability_in >> prob) {
    conditional_probability.push_back(prob);
  }

  return {bin_middle, conditional_probability};
}

std::map<int, int> calculate_distribution(const std::vector<int> &sequence) {
  std::map<int, int> distribution;
  for (int x : sequence) {
    ++distribution[x];
  }
  return distribution;
}

std::vector<bool> calculate_x(const std::vector<int> &s, int ETA) {
  int n = s.size();
  std::vector<bool> x(n);
  for (int i = 0; i < n; ++i)
    x[i] = s[i] > ETA;
  return x;
}

std::vector<double> calculate_y(const std::vector<int> &s, double A) {
  int n = s.size();
  std::vector<double> y(n);
  for (int i = 1; i < n; ++i)
    y[i] = A * y[i - 1] + A * s[i - 1];
  return y;
}

std::vector<double> calculate_bins_bounds(const std::vector<double> &y,
                                          int K_BINS) {
  auto min_y = *std::min_element(y.begin(), y.end());
  auto max_y = *std::max_element(y.begin(), y.end());

  std::vector<double> bins_bounds(K_BINS + 1);
  bins_bounds[0] = min_y;
  bins_bounds[K_BINS] = max_y;

  double bin_size = (max_y - min_y) / K_BINS;
  for (int i = 1; i < K_BINS; ++i) {
    bins_bounds[i] = min_y + i * bin_size;
  }

  return bins_bounds;
}

std::vector<double> calculate_conditional_probability_by_bins(
    const std::vector<bool> &x, const std::vector<double> &y,
    const std::vector<double> &bins_bounds) {
  int n = x.size();
  int K_BINS = int(bins_bounds.size()) - 1;
  std::vector<std::vector<int>> count_in_bin(K_BINS, std::vector<int>(2));

  for (int i = 0; i < n; ++i) {
    int n_bin = lower_bound(bins_bounds.begin() + 1, bins_bounds.end(), y[i]) -
                bins_bounds.begin() - 1;
    ++count_in_bin[n_bin][x[i]];
  }

  std::vector<double> conditional_probability_X1y(K_BINS);
  for (int i = 0; i < K_BINS; ++i) {
    conditional_probability_X1y[i] =
        count_in_bin[i][1] * 1. / (count_in_bin[i][0] + count_in_bin[i][1]);
  }
  return conditional_probability_X1y;
}

template <class T> int find_nearest(const vector<T> &a, T x) {
  if (a.size() == 1) {
    return 0;
  }

  int pos = std::lower_bound(a.begin(), a.end(), x) - a.begin();
  if (pos == a.size()) {
    --pos;
  }

  if (pos + 1 != a.size() && std::abs(a[pos + 1] - x) < std::abs(x - a[pos])) {
    pos = pos + 1;
  }

  return pos;
}

vector<double>
calculate_predictions(const vector<double> &y,
                      const vector<double> &bins_middle,
                      const vector<double> &conditional_probability) {
  vector<double> predictions;
  predictions.reserve(y.size());

  for (double yi : y) {
    predictions.push_back(
        conditional_probability[find_nearest(bins_middle, yi)]);
  }

  return predictions;
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
roc_curve(const std::vector<bool> &y_true, const std::vector<double> &y_pred) {
  int n = y_true.size();
  std::vector<pair<double, bool>> y_pred_true(n);
  for (int i = 0; i < n; ++i) {
    y_pred_true[i] = {y_pred[i], y_true[i]};
  }
  std::sort(y_pred_true.begin(), y_pred_true.end());

  std::vector<double> tpr;
  std::vector<double> fpr;
  std::vector<double> thresholds;

  int tp = 0, fp = 0, tn = 0, fn = 0;
  for (int i = 0; i < n; ++i) {
    if (y_true[i] == 0)
      tn++;
    else
      fn++;
  }

  for (int i = 0; i < n;) {
    int j = i;
    while (j < n && y_pred_true[i].first == y_pred_true[j].first) {
      if (y_pred_true[j].second == 1) {
        tp++;
        fn--;
      } else {
        fp++;
        tn--;
      }

      ++j;
    }

    fpr.push_back(fp + tn == 0 ? 0 : (double)fp / (fp + tn));
    tpr.push_back(tp + fn == 0 ? 0 : (double)tp / (tp + fn));
    thresholds.push_back(y_pred_true[i].first);

    i = j;
  }

  return {fpr, tpr, thresholds};
}

void print_conditional_probability(
    ostream &out, const vector<double> &bins_bounds,
    const vector<double> &conditional_probability) {
  for (int i = 0; i < int(bins_bounds.size()) - 1; ++i) {
    double bin_middle = (bins_bounds[i + 1] + bins_bounds[i]) / 2;
    out << bin_middle << ' ';
  }
  out << '\n';
  for (int i = 0; i < conditional_probability.size(); ++i) {
    out << conditional_probability[i] << ' ';
  }
  out << '\n';
}

void print_roc_curve(ostream &out, const vector<double> &fpr,
                     const vector<double> &tpr,
                     const vector<double> &thresholds) {
  for (auto x : fpr) {
    out << std::fixed << std::setprecision(5) << x << ' ';
  }
  out << '\n';

  for (auto x : tpr) {
    out << std::fixed << std::setprecision(5) << x << ' ';
  }
  out << '\n';

  for (auto x : thresholds) {
    out << std::fixed << std::setprecision(5) << x << ' ';
  }
  out << '\n';
}
