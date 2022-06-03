#include "abel_heap.hpp"

using namespace std;

void print_distribution(std::ofstream& out, const map<int, int>& distribution) {
    for (auto [size, count] : distribution) {
        out << size << ' ' << count << '\n';
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Expected 2 arguments: INPUT, OUTPUT\n";
        exit(1);
    }

    std::ifstream in(argv[1]);
    std::ofstream out(argv[2]);

    auto sequence = read_sequence(in);
    auto distribution = calculate_distribution(sequence);
    print_distribution(out, distribution);

    return 0;
}
