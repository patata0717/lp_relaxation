#include <iostream>
#include <string>
#include "hasCandy.hpp"
using namespace candy;

int main(int argc, char** argv) {
    // Optional: allow setting candy via command line: ./a.out x1 x2 x3
    if (argc == 4) {
        setCandy(std::stoi(argv[1]), std::stoi(argv[2]), std::stoi(argv[3]));
    } else {
        setCandy(3,5,2); // default example; change as you like
    }

    std::cout << "Candy is hidden in 10x10x10 (0..9 each axis).\n";
    std::cout << "Enter six integers: a1 a2 a3 b1 b2 b3 (inclusive ranges), or 'q' to quit.\n";

    std::cout << std::boolalpha;
    for (;;) {
        std::cout << "> ";
        std::string t1;
        if (!(std::cin >> t1)) break;           // EOF
        if (t1 == "q" || t1 == "Q") break;

        int a1;
        try { a1 = std::stoi(t1); }
        catch (...) { std::cout << "invalid a1\n"; continue; }

        int a2, a3, b1, b2, b3;
        if (!(std::cin >> a2 >> a3 >> b1 >> b2 >> b3)) break;

        auto inRange = [](int v){ return 0 <= v && v < N; };
        if (!inRange(a1) || !inRange(a2) || !inRange(a3) ||
            !inRange(b1) || !inRange(b2) || !inRange(b3)) {
            std::cout << "out of range (must be 0..9)\n";
            continue;
        }

        bool inside = hasCandy(a1,a2,a3, b1,b2,b3); // α/β order doesn't matter
        std::cout << (inside ? "yes" : "no") << "\n";
    }
    return 0;
}
