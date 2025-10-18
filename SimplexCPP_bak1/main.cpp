#include <iostream>
#include "Simplex.h"

using namespace std;

double LB = 0;

/*
bool fathomed(S_RESULTS &result) {
    // if all variables are integers
}

void Search_region(axis, region) {
    Call_simplex();
    if (fathomed(result)) {
        // update LB is result > curr LB
        //
        return;
    }

}

Call_simplex(vector<int> a, vector<int> b) {
}

*/
int main()
{
	MATRIX problem;
	S_RESULTS result;

	// --- Fractional optimum (non-grid point) ---
    // Max z = x1 + x2
    // s.t. 2x1 + x2 <= 4
    //      x1 + 2x2 <= 4
    //      x1, x2 >= 0
    problem = {
        LINE("w1", { 2, 1, 1, 0, 4 }),   // 2x1 + x2 + s1 = 4
        LINE("w2", { 1, 2, 0, 1, 4 }),   // x1 + 2x2 + s2 = 4
        LINE("P",  { -1, -1, 0, 0, 0 })  // maximize x1 + x2  -> store -coeffs
    };
    result = Simplex::SolveEq(problem, 2);
    
    // Expect x1 = 4/3, x2 = 4/3, slacks = 0, Pmax = 8/3
    cout << "Fractional optimum test:" << endl;
    cout << "x1: " << result.GetValue("x1") << endl;
    cout << "x2: " << result.GetValue("x2") << endl;
    cout << "w1: " << result.GetValue("w1") << endl;
    cout << "w2: " << result.GetValue("w2") << endl;
    cout << "Pmax: " << result.GetValue("Pmax") << endl;

    // --- Fractional optimum (non-grid point) ---
    // Max z = 5x1 + 8x2
    // s.t. x1 + x2 <= 6
    //      5x1 + 9x2 <= 45
    //      x1, x2 >= 0
    problem = {
        LINE("w1", { 1, 1, 1, 0, 6 }),    // x1 + x2 + s1 = 6
        LINE("w2", { 5, 9, 0, 1, 45 }),   // 5x1 + 9x2 + s2 = 45
        LINE("P",  { -5, -8, 0, 0, 0 })  // maximize 5x1 + 8x2  -> store -coeffs
    };
    result = Simplex::SolveEq(problem, 2);

    // Expect x1 = 9/4, x2 = 15/4, slacks = 0, Pmax = 165/4
    cout << "\nFractional optimum test:" << endl;
    cout << "x1: " << result.GetValue("x1") << endl;
    cout << "x2: " << result.GetValue("x2") << endl;
    cout << "w1: " << result.GetValue("w1") << endl;
    cout << "w2: " << result.GetValue("w2") << endl;
    cout << "Pmax: " << result.GetValue("Pmax") << endl;


    cout << "\nLP with box constraints via variable shift n" << endl;
    // ---------- your original LP (shifted) ----------
    // Max: 5 x1 + 8 x2  
    // s.t. x1 + x2 ≤ 6
    //      5x1 + 9x2 ≤ 45
    //      x1, x2 >= 0
    // and search only inside α<x<β (box constraints)

    // x = y + α:
    // Row1:  y1 +  y2 + s1 = 6  - (a1 + a2)
    // Row2: 5y1 + 9y2 + s2 = 45 - (5a1 + 9a2)
    // B1:    y1 + t1 = b1 - a1
    // B2:    y2 + t2 = b2 - a2
    //        y1, y2 >= 0

    // a = lower bounds; b = upper bounds in x
    const int a1 = 0, a2 = 5;
    const int b1 = 1, b2 = 10;
    const int c1 = 5,  c2 = 8; // objective
    
    // Columns: [y1, y2, s1, s2, t1, t2, RHS]
    MATRIX prob = {
        LINE("w1",  { 1, 1, 1, 0, 0, 0,  6 - (a1 + a2) }),   // y1 + y2 ≤ 3
        LINE("w2",  { 5, 9, 0, 1, 0, 0, 45 - (5*a1 + 9*a2) }), // 5y1 + 9y2 ≤ 30
        LINE("w3", { 1, 0, 0, 0, 1, 0,  b1 - a1 }),        // y1 ≤ 7
        LINE("w4", { 0, 1, 0, 0, 0, 1,  b2 - a2 }),        // y2 ≤ 10
        LINE("P",   { -c1, -c2, 0, 0, 0, 0, 0 })
    };
    
    S_RESULTS res = Simplex::SolveEq(prob, 2);
    
    // Read decisions (first two columns)
    fraction y1 = res.GetValue("x1");
    fraction y2 = res.GetValue("x2");
    
    // Read slacks by row name (optional)
    fraction s1 = res.GetValue("w1");
    fraction s2 = res.GetValue("w2");
    fraction t1 = res.GetValue("w3");
    fraction t2 = res.GetValue("w4");
    
    // Reconstruct x and objective
    fraction x1 = y1 + a1;
    fraction x2 = y2 + a2;
    fraction Pmax = res.GetValue("Pmax") + c1 * a1 + c2 * a2;
    
    cout << "x*: " << x1 << ", " << x2 << "\n";
    cout << "Pmax: " << Pmax << "\n";
    cout << "Slacks: s1=" << s1 << ", s2=" << s2 << ", t1=" << t1 << ", t2=" << t2 << "\n";

    cout << "\nLP with box constraints via variable shift n" << endl;
    // ---------- your original LP (shifted) ----------
    // Max: x1 + x2  
    // s.t. 2x1 - 2x2 ≤ -1
    //      -8x1 + 10x2 ≤ 13
    //      x1, x2 >= 0
    // and search only inside α<x<β (box constraints)

    // x = y + α:
    // Row1:  2y1 - 2y2 + s1 = -1 - 2a1 + 2a2
    // Row2:  -8y1 + 10y2 + s2 = 13 + 8a1 - 10a2
    // B1:    y1 + t1 = b1 - a1
    // B2:    y2 + t2 = b2 - a2
    // Obj:   max y1 + y2
    //        y1, y2 >= 0

    // a = lower bounds; b = upper bounds in x
    const int a1 = 2,  a2 = 0;
    const int b1 = 2,  b2 = 2;
    const int c1 = 1,  c2 = 1; // objective
    
    // Columns: [y1, y2, s1, s2, t1, t2, RHS]
    MATRIX prob = {
        LINE("w1", { 2, -2, 1, 0, 0, 0,  -1 - 2*a1 + 2*a2 }),
        LINE("w2", { -8, 10, 0, 1, 0, 0, 13 + 8*a1 - 10*a2 }),
        LINE("w3", { 1, 0, 0, 0, 1, 0,  b1 - a1 }),        
        LINE("w4", { 0, 1, 0, 0, 0, 1,  b2 - a2 }),        
        LINE("P",  { -c1, -c2, 0, 0, 0, 0, 0 })
    };
    
    S_RESULTS res = Simplex::SolveEq(prob, 2);
    
    // Read decisions (first two columns)
    fraction y1 = res.GetValue("x1");
    fraction y2 = res.GetValue("x2");
    
    // Read slacks by row name (optional)
    fraction s1 = res.GetValue("w1");
    fraction s2 = res.GetValue("w2");
    fraction t1 = res.GetValue("w3");
    fraction t2 = res.GetValue("w4");
    
    // Reconstruct x and objective
    fraction x1 = y1 + a1;
    fraction x2 = y2 + a2;
    fraction Pmax = res.GetValue("Pmax") + c1 * a1 + c2 * a2;
    
    cout << "x: " << x1 << ", " << x2 << "\n";
    cout << "Pmax: " << Pmax << "\n";
    cout << "Slacks: s1=" << s1 << ", s2=" << s2 << ", t1=" << t1 << ", t2=" << t2 << "\n";

    return 0;
}

