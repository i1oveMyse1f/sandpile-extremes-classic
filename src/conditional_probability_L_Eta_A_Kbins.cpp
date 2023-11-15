#include "abel_heap.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 7) {
        cerr << "Expected 6 arguments: INPUT SEQUENCE, OUTPUT CONDITIONAL PROB, L, ETA, A, K_BINS\n";
        exit(1);
    }

    std::ifstream in_sequence(argv[1]);
    std::ofstream out(argv[2]);
    auto L = strtol(argv[3], nullptr, 10);
    auto ETA = strtol(argv[4], nullptr, 10);
    auto A = strtod(argv[5], nullptr);
    auto K_BINS = strtol(argv[6], nullptr, 10);

    auto TRUNCATE_BEGIN = L * L * 10;

    auto s = read_sequence(in_sequence);
    auto x = calculate_x(s, ETA);
    auto y = calculate_y(s, A);

    s.clear();
    x.erase(x.begin(), x.begin() + TRUNCATE_BEGIN);
    y.erase(y.begin(), y.begin() + TRUNCATE_BEGIN);

    auto bins_bounds = calculate_bins_bounds(y, K_BINS);
    auto condition_probability = calculate_conditional_probability_by_bins(x, y, bins_bounds);

    print_conditional_probability(out, bins_bounds, condition_probability);

    return 0;
}
