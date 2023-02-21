#pragma once
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <complex>

#define MAX_N 1001
#define MAX_M 1001

typedef long long lld;
typedef unsigned long long llu;
using namespace std;

class Simplex{
    typedef vector<vector<double>> matrix;
    private:
        int n, m;
        double A[MAX_M][MAX_N], b[MAX_M], c[MAX_N], v;
        int N[MAX_N], B[MAX_M]; // nonbasic & basic

        void pivot(int x, int y);

        // Run a single iteration of the simplex algorithm.
        // Returns: 0 if OK, 1 if STOP, -1 if UNBOUNDED
        int iterate_simplex();

        // (Possibly) converts the LP into a slack form with a feasible basic solution.
        // Returns 0 if OK, -1 if INFEASIBLE
        int initialise_simplex();

        // Runs the simplex algorithm to optimise the LP.
        // Returns a vector of -1s if unbounded, -2s if infeasible.
        pair<vector<double>, double> simplex();

        //reset function
        void reset();

    public:
        pair<vector<double>,double> input(matrix vorfaktoren, vector<double> max_gleichung,vector<double> rechte_seite);
        Simplex();

};