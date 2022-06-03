#include "abel_heap.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 7) {
        cerr << "Expected 6 arguments: INPUT SEQUENCE, INPUT CONDITIONAL PROB, OUTPUT ROC CURVE, L, ETA, A\n";
        exit(1);
    }

    std::ifstream in_sequence(argv[1]);
    std::ifstream in_conditional_probability(argv[2]);
    std::ofstream out(argv[3]);
    auto L = strtol(argv[4], nullptr, 10);
    auto ETA = strtol(argv[5], nullptr, 10);
    auto A = strtod(argv[6], nullptr);
    
    auto TRUNCATE_BEGIN = L * L * 10;

    auto s = read_sequence(in_sequence);
    auto x = calculate_x(s, ETA);
    auto y = calculate_y(s, A);

    s.clear();
    x.erase(x.begin(), x.begin() + TRUNCATE_BEGIN);
    y.erase(y.begin(), y.begin() + TRUNCATE_BEGIN);
    
    auto [bin_middles, conditional_probability] = read_conditional_probability(in_conditional_probability);
    auto predictions = calculate_predictions(y, bin_middles, conditional_probability);
    auto [fpr, tpr, thresholds] = roc_curve(x, predictions);

    print_roc_curve(out, fpr, tpr, thresholds);

    return 0;
}
