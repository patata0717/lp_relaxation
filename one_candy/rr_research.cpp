#include <iostream>
#include <string>
#include "hasCandy.hpp"
using namespace candy;

enum Axis { X1 = 0, X2 = 1, X3 = 2 };

// Build the query in canonical order (x1, x2, x3), split only the chosen axis (lower half)
// If answer is "yes", keep lower half; if "no", keep upper half.
static void query_and_update(Axis ax,
                             int& L1, int& U1,
                             int& L2, int& U2,
                             int& L3, int& U3)
{
    int a1 = L1, b1 = U1;
    int a2 = L2, b2 = U2;
    int a3 = L3, b3 = U3;

    if (ax == X1) {
        int m = (L1 + U1) / 2;
        b1 = m;                                   // test lower half on x1
        bool yes = hasCandy(a1,a2,a3, b1,b2,b3);
        std::cout << "hasCandy("<<a1<<","<<a2<<","<<a3<<", "
                  << b1 <<","<<b2<<","<<b3<<") -> " << (yes?"yes":"no") << "\n";
        if (yes) U1 = m; else L1 = m + 1;
    } else if (ax == X2) {
        int m = (L2 + U2) / 2;
        b2 = m;                                   // test lower half on x2
        bool yes = hasCandy(a1,a2,a3, b1,b2,b3);
        std::cout << "hasCandy("<<a1<<","<<a2<<","<<a3<<", "
                  << b1 <<","<<b2<<","<<b3<<") -> " << (yes?"yes":"no") << "\n";
        if (yes) U2 = m; else L2 = m + 1;
    } else {
        int m = (L3 + U3) / 2;
        b3 = m;                                   // test lower half on x3
        bool yes = hasCandy(a1,a2,a3, b1,b2,b3);
        std::cout << "hasCandy("<<a1<<","<<a2<<","<<a3<<", "
                  << b1 <<","<<b2<<","<<b3<<") -> " << (yes?"yes":"no") << "\n";
        if (yes) U3 = m; else L3 = m + 1;
    }
}

// Pick the next axis in round-robin order, skipping any axis that's already fixed
static Axis next_axis(Axis cur, int L1,int U1,int L2,int U2,int L3,int U3) {
    for (int k = 1; k <= 3; ++k) {
        Axis cand = static_cast<Axis>((static_cast<int>(cur) + k) % 3);
        if ((cand == X1 && L1 < U1) ||
            (cand == X2 && L2 < U2) ||
            (cand == X3 && L3 < U3)) return cand;
    }
    return cur; // all fixed
}

int main(int argc, char** argv) {
    // Set your hidden candy position here (or pass via argv like "3 5 2")
    if (argc >= 4) setCandy(std::stoi(argv[1]), std::stoi(argv[2]), std::stoi(argv[3]));
    else           setCandy(3,5,2);

    // Choose starting axis via argument current_branch_variable: x1/x2/x3 or 0/1/2
    Axis ax = X1;
    if (argc >= 5) {
        std::string s = argv[4];
        if (s=="x1"||s=="X1"||s=="0") ax = X1;
        else if (s=="x2"||s=="X2"||s=="1") ax = X2;
        else if (s=="x3"||s=="X3"||s=="2") ax = X3;
    }

    int L1=0,U1=9, L2=0,U2=9, L3=0,U3=9;
    int steps = 0;

    while (L1<U1 || L2<U2 || L3<U3) {
        // if chosen axis is already fixed, skip to the next available
        if ((ax==X1 && L1==U1) || (ax==X2 && L2==U2) || (ax==X3 && L3==U3)) {
            ax = next_axis(ax, L1,U1, L2,U2, L3,U3);
            continue;
        }
        ++steps;
        query_and_update(ax, L1,U1, L2,U2, L3,U3);
        ax = next_axis(ax, L1,U1, L2,U2, L3,U3);
    }

    std::cout << "Found ("<<L1<<","<<L2<<","<<L3<<") in " << steps << " queries\n";
    return 0;
}
