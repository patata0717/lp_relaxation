#include <algorithm>
#include <iostream>
#include <stdexcept>

namespace candy {
constexpr int N = 10;

struct Candy { int x1=0, x2=0, x3=0; };
inline Candy g;

// 0-based coordinates (0..9)
inline void setCandy(int x1, int x2, int x3) {
    auto ok = [](int v){ return 0 <= v && v < N; };
    if (!ok(x1) || !ok(x2) || !ok(x3)) throw std::out_of_range("setCandy: 0..9 only");
    g.x1 = x1; g.x2 = x2; g.x3 = x3;
}

// hasCandy(α1,α2,α3, β1,β2,β3) with inclusive ranges; α/β order doesn’t matter
inline bool hasCandy(int a1, int a2, int a3, int b1, int b2, int b3) {
    int lo1 = std::min(a1,b1), hi1 = std::max(a1,b1);
    int lo2 = std::min(a2,b2), hi2 = std::max(a2,b2);
    int lo3 = std::min(a3,b3), hi3 = std::max(a3,b3);
    return (lo1 <= g.x1 && g.x1 <= hi1) &&
           (lo2 <= g.x2 && g.x2 <= hi2) &&
           (lo3 <= g.x3 && g.x3 <= hi3);
}
} // namespace candy