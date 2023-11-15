#include "abel_heap.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Expected 1 arguments: INPUT SEQUENCE\n";
        exit(1);
    }

    std::ifstream in(argv[1]);

    auto sequence = read_sequence(in);
    double eta;
    while (std::cin >> eta) {
        int count = std::count_if(sequence.begin(), sequence.end(), [&eta](int x){
            return x > eta;
        });
        std::cout << count << ' ';
    }

    std::cout << sequence.size() << ' ';

    return 0;
}
