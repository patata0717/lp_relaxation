#include <vector>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <cctype>
#include "Simplex.h"            // your patched header (with feasibility checks)
#include "include/fraction.h"

using namespace std;

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

void LP_relaxation(const MATRIX& problem, S_RESULTS& ans, vector<int>& a_right, vector<int>& b_right) {

    pick_axis(curr axis);

    vector<int> a_left, b_left;
    ans = Call_simplex(problem, N_var, a_left, b_left);
    if (fathomed) {
        update LB;
        return;
    }
    LP_relaxation(problem, ans, a_right, b_right);

    vector<int> a_right, b_right;
    ans = Call_simplex(problem, N_var, a_right, b_right);
    if (fathomed) {
        update LB;
        return;
    }
    LP_relaxation(problem, ans, a_right, b_right);
}

int main() {

    // Max z = 5x1 + 8x2
    // s.t. x1 + x2 <= 6
    //      5x1 + 9x2 <= 45
    //      x1, x2 >= 0
    MATRIX problem = {
        LINE("w1", { 1, 1, 1, 0, 6 }),    // x1 + x2 + s1 = 6
        LINE("w2", { 5, 9, 0, 1, 45 }),   // 5x1 + 9x2 + s2 = 45
        LINE("P",  { -5, -8, 0, 0, 0 })  // maximize 5x1 + 8x2  -> store -coeffs
    };

    // Box α < x < β
    vector<int> a = { 0, 5 };
    vector<int> b = { 1, 10 };

    S_RESULTS ans = Call_simplex(problem, 2, a, b);

    cout << "x1: "   << ans.GetValue("x1")  << "\n";
    cout << "x2: "   << ans.GetValue("x2")  << "\n";
    cout << "Pmax: " << ans.GetValue("Pmax")<< "\n";
    return 0;
}
