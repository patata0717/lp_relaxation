#include <bits/stdc++.h>
using namespace std;

constexpr int N = 10;                 // 0..9 on each axis
struct Candy { int x1, x2, x3, val; }; // val in [1..100]

// ---------------- Oracle (simulation for this experiment) ----------------
static vector<Candy> gCandies;     // exactly 3 candies, distinct points
static int g_queries = 0;

// hasCandy: returns true if any candy exists in inclusive box
bool hasCandy(int a1, int a2, int a3, int b1, int b2, int b3) {
    ++g_queries;
    int lo1=min(a1,b1), hi1=max(a1,b1);
    int lo2=min(a2,b2), hi2=max(a2,b2);
    int lo3=min(a3,b3), hi3=max(a3,b3);
    bool ans=false;
    for (const auto& c : gCandies) {
        if (lo1<=c.x1 && c.x1<=hi1 &&
            lo2<=c.x2 && c.x2<=hi2 &&
            lo3<=c.x3 && c.x3<=hi3) { ans=true; break; }
    }
    cout << "hasCandy("<<lo1<<","<<lo2<<","<<lo3<<", "
                     <<hi1<<","<<hi2<<","<<hi3<<") -> " << (ans?"yes":"no") << "\n";
    return ans;
}

// Fathom: only valid at a single grid point. If a candy sits here, return {true, value}
pair<bool,int> fathom_point(int x1,int x2,int x3){
    for (const auto& c : gCandies)
        if (c.x1==x1 && c.x2==x2 && c.x3==x3) return {true, c.val};
    return {false, 0};
}

// ---------------- Branch & Search ----------------
enum Axis { X1=0, X2=1, X3=2 };

struct Box {
    int L1=0,U1=9, L2=0,U2=9, L3=0,U3=9;
    Axis nextAxis = X1;  // which axis to try to split next
};

struct Best {
    int value = 0;              // LB: largest fathomed value so far
    int x1=-1, x2=-1, x3=-1;    // where LB was found
    int foundCount = 0;         // how many distinct candy points we have fathomed
};

// choose next axis in round-robin, skipping already-fixed axes
static Axis pick_axis(Axis cur, const Box& b) {
    for (int k=0;k<3;++k) {
        Axis a = static_cast<Axis>((static_cast<int>(cur)+k)%3);
        if ((a==X1 && b.L1<b.U1) || (a==X2 && b.L2<b.U2) || (a==X3 && b.L3<b.U3))
            return a;
    }
    return cur; // all fixed
}

// Recurse: evaluate current box (caller guarantees parent hadCandy==true)
// Strategy: if singleton -> fathom; else split by chosen axis, evaluate both children,
// recurse only into the children that contain at least one candy.
void search_box(const Box& b, Best& best) {
    // Early exit if we already know the maximum possible outcome
    if (best.foundCount >= 3 || best.value == 100) return;

    // If box is a single point, fathom it
    if (b.L1==b.U1 && b.L2==b.U2 && b.L3==b.U3) {
        auto [has,val] = fathom_point(b.L1,b.L2,b.L3);
        if (has) {
            ++best.foundCount;
            if (val > best.value) {
                best.value = val;
                best.x1=b.L1; best.x2=b.L2; best.x3=b.L3;
                cout << "FATHOM ("<<b.L1<<","<<b.L2<<","<<b.L3<<") value="<<val
                     << "  -> LB="<<best.value<<" (found="<<best.foundCount<<")\n";
            } else {
                cout << "FATHOM ("<<b.L1<<","<<b.L2<<","<<b.L3<<") value="<<val
                     << "  (LB stays "<<best.value<<", found="<<best.foundCount<<")\n";
            }
        }
        return;
    }

    // pick axis to split (round-robin, skipping fixed axes)
    Axis ax = pick_axis(b.nextAxis, b);
    // Build children boxes
    Box low=b, upp=b;
    if (ax==X1) {
        int m=(b.L1+b.U1)/2;
        low.U1 = m;
        upp.L1 = m+1;
    } else if (ax==X2) {
        int m=(b.L2+b.U2)/2;
        low.U2 = m;
        upp.L2 = m+1;
    } else { // X3
        int m=(b.L3+b.U3)/2;
        low.U3 = m;
        upp.L3 = m+1;
    }
    // Next axis for both children (round-robin)
    low.nextAxis = static_cast<Axis>((static_cast<int>(ax)+1)%3);
    upp.nextAxis = low.nextAxis;

    // Evaluate children with hasCandy; recurse only where true
    bool lowHas = hasCandy(low.L1,low.L2,low.L3, low.U1,low.U2,low.U3);
    if (lowHas) search_box(low, best);

    bool uppHas = hasCandy(upp.L1,upp.L2,upp.L3, upp.U1,upp.U2,upp.U3);
    if (uppHas) search_box(upp, best);
}

// ---------------- Demo harness ----------------
int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Default candies (coords and values). You can override via argv:
    // ./a.out  (uses defaults)
    // ./a.out x1 x2 x3 v1  y1 y2 y3 v2  z1 z2 z3 v3  [start_axis]
    vector<Candy> defaults = { {3,5,2,88}, {8,1,9,42}, {0,9,4,95} };

    if (argc >= 13) {
        gCandies = {
            {stoi(argv[1]), stoi(argv[2]), stoi(argv[3]), stoi(argv[4])},
            {stoi(argv[5]), stoi(argv[6]), stoi(argv[7]), stoi(argv[8])},
            {stoi(argv[9]), stoi(argv[10]), stoi(argv[11]), stoi(argv[12])}
        };
    } else {
        gCandies = defaults;
    }

    // Starting axis (optional 14th arg): x1/x2/x3 or 0/1/2
    Box root;
    if (argc >= 14) {
        string s = argv[13];
        if (s=="x1"||s=="X1"||s=="0") root.nextAxis=X1;
        else if (s=="x2"||s=="X2"||s=="1") root.nextAxis=X2;
        else if (s=="x3"||s=="X3"||s=="2") root.nextAxis=X3;
    }

    cout << "Candies:\n";
    for (auto c : gCandies)
        cout << "  ("<<c.x1<<","<<c.x2<<","<<c.x3<<") val="<<c.val<<"\n";
    cout << "Search begins...\n";

    // Root box is the whole cube; we know it contains candy, so skip testing root
    Best best;
    search_box(root, best);

    cout << "\n=== Summary ===\n";
    cout << "Queries: " << g_queries << "\n";
    cout << "Fathomed candies: " << best.foundCount << " / 3\n";
    if (best.value>0)
        cout << "Best candy: value="<<best.value<<" at ("<<best.x1<<","<<best.x2<<","<<best.x3<<")\n";
    else
        cout << "No candy fathomed (unexpected)\n";
    if (best.value==100)
        cout << "Early stop reason: LB reached 100 (global max).\n";
    return 0;
}
