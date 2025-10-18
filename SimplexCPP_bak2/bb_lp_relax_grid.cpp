#include <bits/stdc++.h>
#include "Simplex.h"
using namespace std;

struct LP { int n; vector<vector<int>> A; vector<int> b; vector<int> c; };
struct Box { vector<int> L, U; int nextAxis = 0; };
struct Best { fraction LB = fraction(-1000000,1); vector<int> xbest; };

static bool all_fixed(const Box& B){ for(int j=0;j<(int)B.L.size();++j) if(B.L[j]!=B.U[j]) return false; return true; }
static int pick_axis_rr(const Box& B){
    for(int k=0;k<(int)B.L.size();++k){ int a=(B.nextAxis+k)%B.L.size(); if(B.L[a]<B.U[a]) return a; }
    return B.nextAxis;
}

static LINE make_row_shifted(const string& name, const vector<int>& coeffs,
                             const vector<int>& L, int rhs, int n, int slackIdx, int totalSlacks){
    const int W = n + totalSlacks + 1;
    ROW row(W, fraction(0,1));
    long long shift = 0;
    for(int j=0;j<n;++j){ row[j]=fraction(coeffs[j],1); shift += (long long)coeffs[j]*L[j]; }
    row[n+slackIdx]=fraction(1,1);
    row[W-1]=fraction(rhs - (int)shift, 1);
    return LINE(name,row);
}
static LINE make_bound_row_up(const string& name,int var,int up,int n,int slackIdx,int totalSlacks){
    const int W = n + totalSlacks + 1;
    ROW row(W, fraction(0,1));
    row[var]=fraction(1,1);
    row[n+slackIdx]=fraction(1,1);
    row[W-1]=fraction(up,1);
    return LINE(name,row);
}

// Build tableau in shifted vars x' = x - L (x' >= 0), only upper bounds x' <= U-L
static MATRIX build_tableau_shift(const LP& P, const Box& B, fraction& obj_const){
    const int n=P.n, mCore=(int)P.A.size(), mBounds=n, m=mCore+mBounds, W=n+m+1;
    MATRIX M; M.reserve(m+1); int slack=0;
    for(int i=0;i<mCore;++i) M.push_back(make_row_shifted("w"+to_string(++slack), P.A[i], B.L, P.b[i], n, slack-1, m));
    for(int j=0;j<n;++j){ int up=B.U[j]-B.L[j]; M.push_back(make_bound_row_up("w"+to_string(++slack), j, up, n, slack-1, m)); }
    ROW prow(W, fraction(0,1)); for(int j=0;j<n;++j) prow[j]=fraction(-P.c[j],1); M.push_back(LINE("P",prow));
    long long oc=0; for(int j=0;j<n;++j) oc += (long long)P.c[j]*B.L[j]; obj_const=fraction((int)oc,1);
    return M;
}

static tuple<bool,fraction,vector<fraction>> solve_box(const LP& P, const Box& B){
    try{
        fraction objc(0,1);
        MATRIX M = build_tableau_shift(P,B,objc);
        // trace like your run:
        cerr << "P, B = ("; for(int j=0;j<P.n;++j){ if(j) cerr<<","; cerr<<P.c[j]; } cerr << "), (";
        for(int j=0;j<P.n;++j){ if(j) cerr<<","; cerr<<B.L[j]; }
        cerr << ")..(";
        for(int j=0;j<P.n;++j){ if(j) cerr<<","; cerr<<B.U[j]; }
        cerr << ")\n";

        S_RESULTS res = Simplex::SolveEq(M,P.n);
        fraction ub = res.GetValue("Pmax") + objc;
        vector<fraction> xprime(P.n);
        for(int j=0;j<P.n;++j) xprime[j] = res.GetValue("x"+to_string(j+1));
        return {true, ub, xprime};
    }catch(const exception&){ return {false, fraction(0,1), {}}; }
}

// Try to interpret an LP solution x' as an integer grid point in original vars
static bool extract_integer_point(const vector<fraction>& xprime, const Box& B, vector<int>& xint){
    const int n = (int)xprime.size();
    xint.assign(n,0);
    for(int j=0;j<n;++j){
        bool hit=false;
        for(int k=B.L[j]; k<=B.U[j]; ++k){
            // x = L + x' ; x' integral iff x equals some integer in [L,U]
            // We only see x' values; check x' == k-L
            if (xprime[j] == (k - B.L[j])) { xint[j]=k; hit=true; break; }
        }
        if(!hit) return false;
    }
    return true;
}

// Direct feasibility and objective at a single grid point (no LP call)
static bool point_feasible_and_val(const LP& P, const Box& B, fraction& val){
    for(int i=0;i<(int)P.A.size();++i){
        long long lhs=0; for(int j=0;j<P.n;++j) lhs += (long long)P.A[i][j]*B.L[j];
        if(lhs > P.b[i]) return false;
    }
    long long obj=0; for(int j=0;j<P.n;++j) obj += (long long)P.c[j]*B.L[j];
    val = fraction((int)obj,1);
    return true;
}

static void search(const LP& P, const Box& B, Best& best){
    // Leaf: evaluate directly (fixes the [0,0] -> 14 issue)
    if(all_fixed(B)){
        fraction val(0,1);
        if(point_feasible_and_val(P,B,val) && val > best.LB){ best.LB=val; best.xbest=B.L; }
        return;
    }

    auto [feas, ub, xprime] = solve_box(P,B);
    if(!feas) return;                 // infeasible → fathom
    if(!(ub > best.LB)) return;       // prune

    // Fathom if LP optimum is integer and within the box
    vector<int> xint;
    if(!xprime.empty() && extract_integer_point(xprime, B, xint)){
        if(ub > best.LB){ best.LB = ub; best.xbest = xint; }
        return;
    }

    // Split round-robin at midpoint
    int ax = pick_axis_rr(B);
    int mid = (B.L[ax] + B.U[ax]) / 2;
    Box low=B, upp=B;
    low.U[ax]=mid;  upp.L[ax]=mid+1;
    low.nextAxis=(ax+1)%P.n; upp.nextAxis=low.nextAxis;

    search(P, low, best);
    search(P, upp, best);
}

// Demo LPs
static LP demo2(){ LP P; P.n=2; P.c={1,2}; P.A={{2,1},{1,2},{3,2}}; P.b={14,14,20}; return P; }
static LP demo3(){ LP P; P.n=3; P.c={6,5,4}; P.A={{1,1,1},{2,1,0},{0,2,1}}; P.b={12,14,14}; return P; }

static void run_case(const LP& P, const string& name){
    Box root; root.L.assign(P.n,0); root.U.assign(P.n,9); root.nextAxis=0;
    Best best; search(P, root, best);
    cout << "\n=== " << name << " ===\n";
    cout << "LB = " << best.LB << "\n";
    if(!best.xbest.empty()){ cout << "x* (grid) = ["; for(int j=0;j<P.n;++j){ if(j) cout<<","; cout<<best.xbest[j]; } cout<<"]\n"; }
    else cout << "No grid point fathomed.\n";
}

int main(){
    run_case(demo2(), "2-variable case");
    run_case(demo3(), "3-variable case");
    return 0;
}
