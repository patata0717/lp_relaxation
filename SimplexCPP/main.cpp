#include <vector>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <cctype>
#include "Simplex.h"            // your patched header (with feasibility checks)
#include "include/fraction.h"

using namespace std;

// ===== Globals used by LP_relaxation (set them before the first call) =====
static int        g_next_axis_rr = 0;             // round-robin pointer
static fraction   g_LB = 0;   // incumbent objective (best integer)
static vector<int> g_best_x;                      // incumbent integer solution


// Returns an S_RESULTS that contains:
//   - x1..xn : optimal original variables (after adding alpha)
//   - y1..yn : shifted decisions (helpful for debugging)
//   - Pmax   : true objective (already shifted back)
//   - slacks (including ub1..ubn) from the shifted problem
S_RESULTS Call_simplex(const MATRIX& problem, int num_vars, const vector<int>& a, const vector<int>& b){
    // ---- Basic validation ----
    if ((int)a.size() != num_vars || (int)b.size() != num_vars)
        throw runtime_error("a,b must have length num_vars");
    for (int j = 0; j < num_vars; ++j) {
        if (b[j] < a[j]) {
            ostringstream msg; msg << "Invalid box: b["<<j<<"] < a["<<j<<"]";
            throw runtime_error(msg.str());
        }
    }
    if (problem.empty()) throw runtime_error("Empty matrix");

    // Find a non-P row to deduce sizes
    int first_con = -1;
    for (int i = 0; i < (int)problem.size(); ++i)
        if (problem[i].Key() != "P") { first_con = i; break; }
    if (first_con < 0) throw runtime_error("No constraint rows found");

    const ROW& sample   = problem[first_con].Row();
    const int orig_cols = (int)sample.size();
    const int orig_rhs  = orig_cols - 1;
    const int orig_slacks = orig_cols - 1 - num_vars;
    if (orig_slacks < 0) throw runtime_error("num_vars exceeds column count");

    // Grab objective coefficients c from P row (P stores -c)
    const LINE* prow_ptr = nullptr;
    for (int i = 0; i < (int)problem.size(); ++i)
        if (problem[i].Key() == "P") { prow_ptr = &problem[i]; break; }
    if (!prow_ptr) throw runtime_error("P row not found");
    const ROW& prow = prow_ptr->Row();

    vector<fraction> alpha(num_vars), beta(num_vars), u(num_vars), c(num_vars);
    fraction c_dot_a(0);
    for (int j = 0; j < num_vars; ++j) {
        alpha[j] = fraction(a[j]);
        beta[j]  = fraction(b[j]);
        u[j]     = beta[j] - alpha[j];
        c[j]     = fraction(0) - prow[j];   // because prow has -c
        c_dot_a += c[j] * alpha[j];
    }

    // New matrix sizing
    const int new_slacks = orig_slacks + num_vars;     // add one UB slack per variable
    const int new_cols   = num_vars + new_slacks + 1;  // y + slacks + RHS
    MATRIX shifted;

    // ---- Shift each constraint row: A y + s = (b - A a) ----
    for (int i = 0; i < (int)problem.size(); ++i) {
        if (problem[i].Key() == "P") continue;

        const ROW& r = problem[i].Row();
        ROW row(new_cols, fraction(0));

        // copy A (first num_vars)
        for (int j = 0; j < num_vars; ++j) row[j] = r[j];

        // copy original slacks into same block
        for (int k = 0; k < orig_slacks; ++k) row[num_vars + k] = r[num_vars + k];

        // RHS shift: b' = b - A*alpha
        fraction rhs = r[orig_rhs];
        for (int j = 0; j < num_vars; ++j) rhs -= r[j] * alpha[j];
        row[new_cols - 1] = rhs;

        shifted.push_back( LINE(problem[i].Key(), row) );
    }

    // ---- Add upper-bound rows: y_j + t_j = u_j ----
    for (int j = 0; j < num_vars; ++j) {
        ROW row(new_cols, fraction(0));
        row[j] = fraction(1);                               // y_j
        row[num_vars + orig_slacks + j] = fraction(1);      // t_j
        row[new_cols - 1] = u[j];                           // RHS = β-α

        ostringstream key; key << "ub" << (j+1);
        shifted.push_back( LINE(key.str(), row) );
    }

    // ---- Objective row: keep -c on y, zeros on slacks & RHS ----
    {
        ROW row(new_cols, fraction(0));
        for (int j = 0; j < num_vars; ++j) row[j] = fraction(0) - c[j]; // store -c again
        shifted.push_back( LINE("P", row) );
    }

    // ---- Solve on (y,s,t)
    S_RESULTS res_shifted = Simplex::SolveEq(shifted, num_vars);

    // ---- Build final results in x-coordinates ----
    S_RESULTS out;

    // a) x* = y* + a  (and also return y*)
    for (int j = 0; j < num_vars; ++j) {
        fraction yj = res_shifted.GetValue(string("x") + to_string(j+1)); // these are y_j
        fraction xj = yj + alpha[j];
        out.Add( S_RESULT( string("x") + to_string(j+1), xj ) );
        out.Add( S_RESULT( string("y") + to_string(j+1), yj ) );
    }

    // b) True objective: P(x) = P(y) + c·a
    fraction Py = res_shifted.GetValue("Pmax");
    out.Add( S_RESULT("Pmax", Py + c_dot_a) );

    // c) Pass through other slacks/labels (skip x1..xn; we've added x/y already)
    for (const auto& it : res_shifted.GetAll()) {
        const string& k = it.Key();
        if (k == "Pmax") continue;
        bool is_x = (k.size()>=2 && k[0]=='x' && isdigit((unsigned char)k[1]));
        if (!is_x) out.Add(it);
    }

    return out;
}




// ----- exact helpers: no floating point, no <cmath> -----

// return max m in [lo,hi] s.t. m <= x  (assumes LP respects the box, so x∈[lo,hi])
static inline int floor_in_box(const fraction& x, int lo, int hi) {
    int L = lo, R = hi, ans = lo - 1;
    while (L <= R) {
        int m = (L + R) >> 1;
        fraction d = x - m;                   // needs operator-(int)
        if (d < fraction(0)) {                // needs operator<(fraction)
            R = m - 1;
        } else {
            ans = m;
            L = m + 1;
        }
    }
    if (ans < lo) ans = lo - 1;
    if (ans > hi) ans = hi;
    return ans;
}

static inline bool is_integer_in_box(const fraction& x, int lo, int hi) {
    int f = floor_in_box(x, lo, hi);
    return (x - f) == fraction(0);            // needs operator==(fraction)
}

// pick next non-integer variable in round-robin order, skipping fixed a[j]==b[j]
static int pick_next_var_rr(const S_RESULTS& sol,
                            const vector<int>& a, const vector<int>& b) {
    const int N = (int)a.size();
    for (int t = 0; t < N; ++t) {
        int j = (g_next_axis_rr + t) % N;
        if (a[j] >= b[j]) continue; // already fixed
        const fraction xj = sol.GetValue(string("x") + to_string(j + 1));
        if (!is_integer_in_box(xj, a[j], b[j])) return j;
    }
    return -1; // all integer or fixed
}

void LP_relaxation(const MATRIX& problem,
                   S_RESULTS& ans,
                   vector<int>& a,
                   vector<int>& b)
{
    // ---- Lightweight logging helpers (function-static) ----
    static long long s_node_id = 0;
    static int       s_depth   = 0;

    const long long node_id = ++s_node_id;
    const int N = (int)a.size();

    auto indent = [&](){ return std::string(s_depth * 2, ' '); };
    auto print_box = [&](const vector<int>& A, const vector<int>& B){
        std::cout << indent() << "[node " << node_id << "] box=[";
        for (int i = 0; i < (int)A.size(); ++i) {
            if (i) std::cout << " × ";
            std::cout << A[i] << ".." << B[i];
        }
        std::cout << "]\n";
    };

    // ---- Log: entering node ----
    print_box(a, b);

    // ---- Solve LP on current box ----
    bool feasible = true;
    try {
        ans = Call_simplex(problem, N, a, b);
    } catch (const std::exception&) {
        feasible = false;
    }

    if (!feasible) {
        std::cout << indent() << "[node " << node_id << "] prune: infeasible\n";
        return; // prune THIS node; parent will continue with sibling
    }

    const fraction UB = ans.GetValue("Pmax");
    std::cout << indent() << "[node " << node_id << "] ";
    for (int j = 0; j < (int)a.size(); ++j) {
        if (j) std::cout << ", ";
        std::cout << "x" << (j+1) << "=" << ans.GetValue(string("x") + to_string(j+1));
    }
    std::cout << ", UB=" << UB << "  LB=" << g_LB << "\n";
    if (UB <= g_LB) {
        std::cout << indent() << "[node " << node_id << "] prune: UB <= LB\n";
        return;
    }

    // ---- Check integrality ----
    bool all_int = true;
    vector<int> x_int(N, 0);
    for (int j = 0; j < N; ++j) {
        const fraction xj = ans.GetValue(string("x") + to_string(j + 1));
        if (!is_integer_in_box(xj, a[j], b[j])) { all_int = false; break; }
        x_int[j] = floor_in_box(xj, a[j], b[j]); // exact when integer
    }

    if (all_int) {
        std::cout << indent() << "[node " << node_id << "] FATHOM: x*=(";
        for (int j = 0; j < N; ++j) {
            if (j) std::cout << ",";
            std::cout << x_int[j];
        }
        std::cout << ")  P=" << UB;
        if (UB > g_LB) std::cout << "  -> new LB";
        std::cout << "\n";

        if (UB > g_LB) { g_LB = UB; g_best_x = x_int; }
        return;
    }

    // ---- Branch (round-robin among non-integers) ----
    const int k = pick_next_var_rr(ans, a, b);
    if (k < 0) {
        std::cout << indent() << "[node " << node_id << "] prune: no non-integer var (safety)\n";
        return;
    }

    const fraction xk = ans.GetValue(string("x") + to_string(k + 1));
    const int f = floor_in_box(xk, a[k], b[k]);
    const int c = f + 1; // ceil, since xk is fractional in (f, f+1)

    std::cout << indent()
              << "[node " << node_id << "] BRANCH on x" << (k + 1)
              << " = " << xk << "  -> left: <= " << f << " | right: >= " << c << "\n";

    // ---- Recurse: left child ----
    {
        vector<int> a_left = a, b_left = b;
        b_left[k] = std::min(b_left[k], f);
        if (b_left[k] >= a_left[k]) {
            ++s_depth;
            LP_relaxation(problem, ans, a_left, b_left);
            --s_depth;
        } else {
            std::cout << indent() << "  (skip left: empty box)\n";
        }
    }

    // ---- Recurse: right child ----
    {
        vector<int> a_right = a, b_right = b;
        a_right[k] = std::max(a_right[k], c);
        if (b_right[k] >= a_right[k]) {
            ++s_depth;
            LP_relaxation(problem, ans, a_right, b_right);
            --s_depth;
        } else {
            std::cout << indent() << "  (skip right: empty box)\n";
        }
    }
}



int main() {

    // Max z = 3x1 + 2x2 + 4x3
    // s.t.
    //   2x1 + 1x2 + 1x3 <= 9
    //   1x1 + 3x2 + 2x3 <= 13
    //   x1, x2, x3 >= 0 (and boxed by a..b)


    // 20 items, capacity = 300
    // weights = [12,13,12,11,22,17,20,25,11,29,17,18,13,30,25,22,28,25,13,22]
    // values  = [29,28,20,19, 7, 6,15,33,32,29,17,39,29,16,39,12,17,11,17,37]
    
    MATRIX problem = {
        // capacity constraint: sum(w_i * x_i) + s1 = 300
        LINE("w1", {
            12,13,12,11,22,17,20,25,11,29,
            17,18,13,30,25,22,28,25,13,22,
            1,   // s1
            300  // RHS
        }),
        // objective (maximize): store -values on x's, 0 on slack, 0 RHS
        LINE("P", {
            -29,-28,-20,-19,-7,-6,-15,-33,-32,-29,
            -17,-39,-29,-16,-39,-12,-17,-11,-17,-37,
             0,  // s1
             0   // RHS
        })
    };


    // Box α < x < β
    vector<int> a = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, };
    vector<int> b = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, };

    S_RESULTS ans = Call_simplex(problem, 20, a, b);

    for (int j=0; j < a.size(); ++j)
        cout << "x" << (j+1) << ": " << ans.GetValue("x"+to_string(j+1)) << "\n";
    cout << "Pmax: " << ans.GetValue("Pmax")<< "\n";

    // Box α < x < β
    S_RESULTS ans_root;
    LP_relaxation(problem, ans_root, a, b);
    
    // Optional: report incumbent
    cout << "LB (best integer Pmax): " << g_LB << "\n";
    for (int j=0; j < a.size(); ++j)
        cout << "x" << (j+1) << " = " << g_best_x[j] << "\n";

    return 0;
}
