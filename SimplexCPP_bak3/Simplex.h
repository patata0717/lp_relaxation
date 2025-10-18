#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include "include/fraction.h"

using namespace std;

// ---------- Basic types ----------
typedef vector<fraction> ROW;

class LINE
{
private:
    string _key;
    ROW _row;

public:
    LINE(const string& key = "", const ROW& row = {})
        : _key(key), _row(row)
    {}

    string Key() const { return _key; }
    string& Key() { return _key; }

    ROW Row() const { return _row; }
    ROW& Row() { return _row; }
};

typedef vector<LINE> MATRIX;

class S_RESULT
{
private:
    string _key;
    fraction _value;

public:
    S_RESULT(const string &key = "", const fraction &value = fraction(0))
        : _key(key), _value(value)
    {}

    string Key() const { return _key; }
    fraction Value() const { return _value; }

    bool IsNull() const { return _key.length() == 0; }
};

class S_RESULTS
{
private:
    vector<string> keys;
    vector<S_RESULT> results;

public:
    S_RESULTS(const vector<S_RESULT>& res = {})
    {
        for (int i = 0; i < (int)res.size(); i++)
        {
            const S_RESULT& val = res[i];
            if (val.IsNull()) continue;
            keys.push_back(val.Key());
            results.push_back(val);
        }
    }

    void Add(const S_RESULT& value)
    {
        if (value.IsNull()) return;
        keys.push_back(value.Key());
        results.push_back(value);
    }

    S_RESULT GetResult(const string& key) const
    {
        for (int i = 0; i < (int)keys.size(); i++)
            if (key == keys[i]) return results[i];
        throw runtime_error("Key does not exist");
    }

    fraction GetValue(const string& key) const
    {
        for (int i = 0; i < (int)keys.size(); i++)
            if (key == keys[i]) return results[i].Value();
        throw runtime_error("Key does not exist");
    }

    vector<S_RESULT> GetAll() const { return results; }
    int Length() const { return (int)keys.size(); }
};

// --- Add this helper near your other helpers in Simplex class ---
static inline void CheckFeasibilityOrThrow(const MATRIX& arr, const S_RESULTS& ans, int no_variables)
{
    if (arr.empty()) return;
    const int ncols   = (int)arr[0].Row().size();
    const int rhs_col = ncols - 1;

    // Gather variable values (decisions x1.. and slacks w1..)
    vector<fraction> val(ncols - 1, fraction(0));

    // Decision variables must be >= 0
    for (int j = 0; j < no_variables; ++j) {
        ostringstream s; s << "x" << (j + 1);
        const fraction xj = ans.GetValue(s.str());
        if (xj < fraction(0)) {
            ostringstream msg; msg << "Infeasible: " << s.str() << " < 0";
            throw runtime_error(msg.str());
        }
        val[j] = xj;
    }

    // Slacks w1..wK must be >= 0 (K = (#cols - 1 RHS) - no_variables)
    const int slack_count = max(0, ncols - 1 - no_variables);
    for (int k = 0; k < slack_count; ++k) {
        ostringstream s; s << "w" << (k + 1);
        const fraction wk = ans.GetValue(s.str());
        if (wk < fraction(0)) {
            ostringstream msg; msg << "Infeasible: " << s.str() << " < 0";
            throw runtime_error(msg.str());
        }
        val[no_variables + k] = wk;
    }

    // Each constraint row must satisfy exact equality: sum_j a_ij * v_j == b_i
    for (int i = 0; i < (int)arr.size(); ++i) {
        if (arr[i].Key() == "P") continue;
        const ROW& row = arr[i].Row();
        fraction lhs(0);
        for (int j = 0; j < rhs_col; ++j) lhs += row[j] * val[j];
        if (lhs != row[rhs_col]) {
            ostringstream msg;
            msg << "Infeasible: row '" << arr[i].Key() << "' violated (LHS=" << lhs
                << ", RHS=" << row[rhs_col] << ")";
            throw runtime_error(msg.str());
        }
    }
}

// --- Add this helper ---
static inline void GuardImmediateInfeasibility(const MATRIX& arr, int no_variables) {
    if (arr.empty()) return;
    const int rhs_col = (int)arr[0].Row().size() - 1;

    for (int i = 0; i < (int)arr.size(); ++i) {
        if (arr[i].Key() == "P") continue;
        const ROW& r = arr[i].Row();
        // If RHS is negative and all coefficients (excluding RHS) are >= 0,
        // then with y >= 0 the minimum LHS is 0 -> infeasible immediately.
        bool all_nonneg = true;
        for (int j = 0; j < rhs_col; ++j) {
            if (r[j] < fraction(0)) { all_nonneg = false; break; }
        }
        if (all_nonneg && r[rhs_col] < fraction(0)) {
            throw runtime_error("Infeasible: a row has nonnegative coefficients but negative RHS.");
        }
    }
}


// ---------- Simplex core ----------
class Simplex
{
public:
    static S_RESULTS SolveEq(const MATRIX& array, const int no_variables)
    {
		GuardImmediateInfeasibility(array, no_variables);  // <--- add this line
        bool optimal;
        int col = 0;

        bool artificial_var = HasArtificialVariable(array, no_variables);

        MATRIX arr = SimplexEngine(array);                 // one pivot
        while (ShouldIterate(arr))                         // drive negatives in P row to nonnegative
            arr = SimplexEngine(arr);

        if (artificial_var)
        {
            optimal = IsOptimal(arr, no_variables, col);   // ensure first no_variables cost coeffs are non-positive
            while (!optimal)
            {
                arr = SimplexEngine(arr, optimal, col);
                optimal = IsOptimal(arr, no_variables, col);
            }
        }

        // Build results by inspecting the final tableau
        S_RESULTS ans;
        for (int i = 0; i < (int)arr.size(); i++)
        {
            S_RESULT res = ResolveResult(arr, i, no_variables);
            ans.Add(res);
        }

        // Ensure x1..x_no_variables exist (fill missing with 0)
        for (int j = 0; j < no_variables; j++)
        {
            ostringstream s; s << "x" << (j + 1);
            try { (void)ans.GetValue(s.str()); }
            catch (const exception&) { ans.Add(S_RESULT(s.str(), fraction(0))); }
        }

        // Ensure all slack names w1..wK exist (fill missing with 0)
        // K = (#cols - 1 RHS) - no_variables
        int slack_count = 0;
        if (!arr.empty())
        {
            int ncols = (int)arr[0].Row().size();
            slack_count = max(0, ncols - 1 - no_variables);
        }
        for (int k = 0; k < slack_count; k++)
        {
            ostringstream s; s << "w" << (k + 1);
            try { (void)ans.GetValue(s.str()); }
            catch (const exception&) { ans.Add(S_RESULT(s.str(), fraction(0))); }
        }

		CheckFeasibilityOrThrow(arr, ans, no_variables);

        return ans;
    }

private:
    // One simplex iteration. If 'optimal' is false, force the entering column = optimal_col.
    static MATRIX SimplexEngine(const MATRIX& _arr, const bool optimal = true, const int optimal_col = 0)
    {
        MATRIX arr = _arr;

        // Locate the objective (P) row
        int key_row = -1;
        for (int i = 0; i < (int)arr.size(); i++)
        {
            if (arr[i].Key() == "P") { key_row = i; break; }
        }
        if (key_row < 0) throw runtime_error("P row not found");

        // Choose entering column
        int key_col = 0;
        const ROW& prow = arr[key_row].Row();
        const int rhs_col = (int)prow.size() - 1;

        if (optimal)
        {
            // pick the most negative coefficient in P row (excluding RHS)
            fraction min_value = fraction(0);
            for (int j = 0; j < rhs_col; j++)
            {
                const fraction& c = prow[j];
                if (c < min_value)
                {
                    min_value = c;
                    key_col = j;
                }
            }
            // If all c_j >= 0, tableau is already optimal for this phase; still do a no-op pivot guard below
        }
        else
        {
            key_col = optimal_col; // provided by IsOptimal
        }

        // Minimum nonnegative ratio test (allow 0); only rows with positive divisor
        int pivot_row = -1;
        fraction best_ratio;            // uninitialized sentinel
        bool found = false;
        
        for (int i = 0; i < (int)arr.size(); i++)
        {
            if (i == key_row) continue;                 // skip P row
            const ROW& r = arr[i].Row();
            fraction a = r[key_col];
            if (a <= fraction(0)) continue;             // primal simplex: only a > 0
        
            const int rhs_col = (int)r.size() - 1;
            fraction ratio = r[rhs_col] / a;            // can be 0 (degenerate)
            if (!found || ratio < best_ratio)
            {
                best_ratio = ratio;
                pivot_row = i;
                found = true;
            }
        }
        
        if (!found)
            throw runtime_error("Unbounded or infeasible: no valid pivot row found."); 

        // Normalize pivot row
        fraction piv = arr[pivot_row].Row()[key_col];
        for (int j = 0; j < (int)arr[pivot_row].Row().size(); j++)
            arr[pivot_row].Row()[j] /= piv;

        // Eliminate pivot column entries from other rows
        const ROW main_row = arr[pivot_row].Row();
        for (int i = 0; i < (int)arr.size(); i++)
        {
            if (i == pivot_row) continue;
            ROW& r = arr[i].Row();
            fraction coeff = r[key_col];
            if (coeff == fraction(0)) continue;
            for (int j = 0; j < (int)r.size(); j++)
                r[j] -= main_row[j] * coeff;
        }

        return arr;
    }

    // After Phase 0, enforce non-positivity of reduced costs on the first no_variables columns
    static inline bool IsOptimal(const MATRIX &arr, const int no_variables, int &col)
    {
        col = 0;
        for (int r = 0; r < (int)arr.size(); r++)
        {
            if (arr[r].Key() == "P")
            {
                const ROW& row = arr[r].Row();
                for (int j = 0; j < no_variables; j++)
                {
                    if (row[j] > fraction(0))
                    {
                        col = j;
                        return false; // not optimal
                    }
                }
                break;
            }
        }
        return true; // optimal for this phase
    }

    // Continue iterating while any P-row coefficient (excluding RHS) is negative
    static inline bool ShouldIterate(const MATRIX &arr)
    {
        for (int r = 0; r < (int)arr.size(); r++)
        {
            if (arr[r].Key() == "P")
            {
                const ROW& row = arr[r].Row();
                const int rhs_col = (int)row.size() - 1;
                for (int j = 0; j < rhs_col; j++)
                {
                    if (row[j] < fraction(0)) return true;
                }
                break;
            }
        }
        return false;
    }

    // Heuristic: if any non-P row has negative RHS or any slack column has a negative coeff, we’ll need Phase I
    static inline bool HasArtificialVariable(const MATRIX& arr, const int no_variables)
    {
        if (arr.empty()) return false;
        const int rhs_col = (int)arr[0].Row().size() - 1;

        for (int i = 0; i < (int)arr.size(); i++)
        {
            if (arr[i].Key() == "P") continue;
            const ROW& row = arr[i].Row();

            // negative RHS → infeasible BFS → needs Phase I
            if (row[rhs_col] < fraction(0)) return true;

            // any negative coefficient in slack region → could imply surplus/artificial (not fully implemented)
            for (int j = no_variables; j < rhs_col; j++)
                if (row[j] < fraction(0)) return true;
        }
        return false;
    }

    // Robust basic-variable detection:
    // A row defines a basic var only if there exists a column j whose column is an identity vector
    static inline S_RESULT ResolveResult(const MATRIX &arr, int row_index, int no_variables)
    {
        const LINE& line = arr[row_index];
        const ROW& row  = line.Row();
        const int ncols = (int)row.size();
        const int rhs_col = ncols - 1;

        // Objective row: return Pmax
        if (line.Key() == "P")
            return S_RESULT("Pmax", row[rhs_col]);

        // Search for an identity column in this row
        int pos = -1;
        for (int j = 0; j < rhs_col; ++j)
        {
            if (row[j] != fraction(1)) continue;

            bool identity = true;
            for (int k = 0; k < (int)arr.size(); ++k)
            {
                if (k == row_index) continue;
                if (arr[k].Row()[j] != fraction(0)) { identity = false; break; }
            }
            if (identity) { pos = j; break; }
        }

        if (pos < 0) return S_RESULT(); // no basic variable detected in this row

        ostringstream name;
        if (pos < no_variables)
            name << "x" << (pos + 1);
        else
            name << "w" << (pos - no_variables + 1);

        return S_RESULT(name.str(), row[rhs_col]);
    }
};
