#include "lp.h"
#include <cmath>
#include <algorithm>
#include <limits>

using namespace arma;

/**
 * Two-Phase Simplex Method Implementation
 *
 * Solves: minimize c'x subject to Ax = b, x >= 0
 *
 * Phase 1: Find initial feasible basis using artificial variables
 * Phase 2: Optimize objective from feasible basis
 */

/**
 * Find entering variable using Bland's rule (smallest index with negative reduced cost)
 * Returns -1 if optimal
 */
static int find_entering_variable(const rowvec& reduced_costs, double tol) {
    int n = reduced_costs.n_elem;
    for (int j = 0; j < n; ++j) {
        if (reduced_costs(j) < -tol) {
            return j;
        }
    }
    return -1;
}

/**
 * Find leaving variable using minimum ratio test with Bland's rule for ties
 * Returns -1 if unbounded
 */
static int find_leaving_variable(const mat& tableau, int entering_col,
                                  const uvec& basis, double tol) {
    int m = tableau.n_rows - 1; // Exclude objective row
    int leaving_row = -1;
    double min_ratio = std::numeric_limits<double>::infinity();

    for (int i = 0; i < m; ++i) {
        double pivot_elem = tableau(i, entering_col);
        if (pivot_elem > tol) {
            double ratio = tableau(i, tableau.n_cols - 1) / pivot_elem;
            if (ratio < min_ratio - tol) {
                min_ratio = ratio;
                leaving_row = i;
            } else if (std::abs(ratio - min_ratio) < tol) {
                // Bland's rule: choose smallest basis index for ties
                if (leaving_row == -1 || basis(i) < basis(leaving_row)) {
                    leaving_row = i;
                }
            }
        }
    }

    return leaving_row;
}

/**
 * Perform pivot operation on tableau
 */
static void pivot(mat& tableau, int pivot_row, int pivot_col) {
    double pivot_elem = tableau(pivot_row, pivot_col);
    tableau.row(pivot_row) /= pivot_elem;

    int m = tableau.n_rows;
    for (int i = 0; i < m; ++i) {
        if (i != pivot_row) {
            double multiplier = tableau(i, pivot_col);
            tableau.row(i) -= multiplier * tableau.row(pivot_row);
        }
    }
}

/**
 * Phase 1: Find initial feasible basis using artificial variables
 *
 * Solves: minimize sum(artificial variables)
 * Returns: basis indices and whether feasible solution exists
 */
static bool phase1_simplex(const mat& A, const vec& b, uvec& basis,
                           mat& tableau, double tol, int max_iter = 10000) {
    uword m = A.n_rows;  // Number of constraints
    uword n = A.n_cols;  // Number of original variables

    // Build Phase 1 tableau: [A | I | b]
    // Variables: [original vars | artificial vars | RHS]
    // Objective: minimize sum of artificial variables
    tableau = mat(m + 1, n + m + 1, fill::zeros);

    // Constraint rows
    tableau(span(0, m-1), span(0, n-1)) = A;
    tableau(span(0, m-1), span(n, n+m-1)) = eye<mat>(m, m);
    tableau(span(0, m-1), n + m) = b;

    // Objective row: minimize sum of artificial variables
    // After adding constraints, reduced costs for artificial vars
    tableau(m, span(n, n+m-1)).fill(1.0);

    // Initial basis: artificial variables
    basis = uvec(m);
    for (uword i = 0; i < m; ++i) {
        basis(i) = n + i;
    }

    // Update objective row by subtracting constraint rows
    for (uword i = 0; i < m; ++i) {
        tableau.row(m) -= tableau.row(i);
    }

    // Run simplex on Phase 1 problem
    int iter = 0;
    while (iter < max_iter) {
        // Get reduced costs (objective row, excluding RHS)
        rowvec reduced_costs = tableau(m, span(0, n + m - 1));

        // Find entering variable
        int entering = find_entering_variable(reduced_costs, tol);
        if (entering == -1) {
            // Optimal solution found for Phase 1
            break;
        }

        // Find leaving variable
        int leaving = find_leaving_variable(tableau, entering, basis, tol);
        if (leaving == -1) {
            // Unbounded (shouldn't happen in Phase 1 with artificial vars)
            return false;
        }

        // Pivot
        pivot(tableau, leaving, entering);
        basis(leaving) = entering;

        ++iter;
    }

    if (iter >= max_iter) {
        return false;
    }

    // Check if feasible: objective value should be ~0
    double phase1_obj = -tableau(m, n + m);
    if (phase1_obj > tol) {
        // Infeasible
        return false;
    }

    // Remove artificial variables from basis if present
    for (uword i = 0; i < m; ++i) {
        if (basis(i) >= n) {
            // Artificial variable in basis - try to pivot it out
            for (uword j = 0; j < n; ++j) {
                if (std::abs(tableau(i, j)) > tol) {
                    pivot(tableau, i, j);
                    basis(i) = j;
                    break;
                }
            }
        }
    }

    return true;
}

/**
 * Phase 2: Optimize objective function from feasible basis
 */
static LPResult phase2_simplex(const vec& c, const mat& tableau_p1, const uvec& basis_p1,
                                uword n_orig, double tol, int max_iter = 10000) {
    LPResult result;
    uword m = basis_p1.n_elem;

    // Remove artificial variable columns from tableau
    // Keep only: original variables (0..n_orig-1) and RHS (last column)
    uword rhs_col = tableau_p1.n_cols - 1;
    mat tableau(m + 1, n_orig + 1);
    tableau(span(0, m-1), span(0, n_orig-1)) = tableau_p1(span(0, m-1), span(0, n_orig-1));
    tableau(span(0, m-1), n_orig) = tableau_p1(span(0, m-1), rhs_col);

    // Copy basis (but artificial vars should have been pivoted out already)
    uvec basis = basis_p1;

    // Build Phase 2 objective row
    tableau.row(m).zeros();
    tableau(m, span(0, n_orig - 1)) = c.t();

    // Update reduced costs for current basis
    for (uword i = 0; i < m; ++i) {
        if (basis(i) < n_orig) {
            tableau.row(m) -= c(basis(i)) * tableau.row(i);
        }
    }

    // Run simplex for Phase 2
    int iter = 0;
    while (iter < max_iter) {
        // Get reduced costs
        rowvec reduced_costs = tableau(m, span(0, n_orig - 1));

        // Find entering variable
        int entering = find_entering_variable(reduced_costs, tol);
        if (entering == -1) {
            // Optimal solution found
            result.status = LPStatus::OPTIMAL;
            result.objval = -tableau(m, n_orig);

            // Extract solution
            result.solution = vec(n_orig, fill::zeros);
            for (uword i = 0; i < m; ++i) {
                if (basis(i) < n_orig) {
                    result.solution(basis(i)) = tableau(i, n_orig);
                }
            }
            return result;
        }

        // Find leaving variable
        int leaving = find_leaving_variable(tableau, entering, basis, tol);
        if (leaving == -1) {
            // Unbounded
            result.status = LPStatus::UNBOUNDED;
            return result;
        }

        // Pivot
        pivot(tableau, leaving, entering);
        basis(leaving) = entering;

        ++iter;
    }

    // Max iterations reached
    result.status = LPStatus::ERROR;
    return result;
}

/**
 * Main LP solver using two-phase simplex
 */
LPResult solve_lp_simplex(const vec& c, const mat& A, const vec& b, double tol) {
    LPResult result;

    // Validate inputs
    if (A.n_rows != b.n_elem || A.n_cols != c.n_elem) {
        result.status = LPStatus::ERROR;
        return result;
    }

    uword m = A.n_rows;
    uword n = A.n_cols;

    // Check for negative RHS and flip constraints if needed
    mat A_work = A;
    vec b_work = b;
    for (uword i = 0; i < m; ++i) {
        if (b_work(i) < -tol) {
            A_work.row(i) = -A_work.row(i);
            b_work(i) = -b_work(i);
        }
    }

    uvec basis;
    mat tableau;

    // Phase 1: Find feasible basis
    if (!phase1_simplex(A_work, b_work, basis, tableau, tol)) {
        result.status = LPStatus::INFEASIBLE;
        return result;
    }

    // Phase 2: Optimize
    result = phase2_simplex(c, tableau, basis, n, tol);

    return result;
}

/**
 * Compute bounds on contrasts using custom LP solver (optimized version)
 */
std::tuple<mat, mat> bounds_lp_contrast_cpp(
    const mat& x,
    const mat& y,
    const mat& contr_m,
    double ub,
    double scale,
    double shift,
    bool sum_one,
    bool has_ub
) {
    uword n = x.n_rows;      // Number of observations
    uword n_x = x.n_cols;    // Number of predictors
    uword n_y = y.n_cols;    // Number of outcomes
    uword n_c = contr_m.n_cols; // Number of contrasts
    uword n_vars = n_y * n_x;   // Number of variables (entries of B)

    mat res_min(n, n_c);
    mat res_max(n, n_c);

    // Build constraint matrix structure (reused across observations)
    uword n_eq = n_y + (sum_one ? n_x : 0);
    uword n_ineq = has_ub ? n_vars : 0;
    uword n_constraints = n_eq + n_ineq;
    uword n_total_vars = n_vars + n_ineq; // Original + slack variables

    mat A(n_constraints, n_total_vars, fill::zeros);
    vec b(n_constraints);

    // Part 2: sum_one constraints (static across observations)
    if (sum_one) {
        for (uword k = 0; k < n_x; ++k) {
            for (uword j = 0; j < n_y; ++j) {
                A(n_y + k, j * n_x + k) = 1.0;
            }
            b(n_y + k) = scale * (1.0 - shift);
        }
    }

    // Part 3: Upper bound constraints (static across observations)
    if (has_ub) {
        for (uword idx = 0; idx < n_vars; ++idx) {
            A(n_eq + idx, idx) = 1.0;           // Original variable
            A(n_eq + idx, n_vars + idx) = 1.0;  // Slack variable
            b(n_eq + idx) = ub;
        }
    }

    // Solve LP for each observation and contrast
    for (uword j = 0; j < n_c; ++j) {
        vec contr = contr_m.col(j);
        double obj_offset = sum(contr) * shift;

        // Build objective for original variables (reused across observations)
        vec c(n_total_vars, fill::zeros);
        c(span(0, n_vars - 1)) = contr;

        for (uword i = 0; i < n; ++i) {
            // Part 1: Update observation-specific constraints B %*% x[i,] = y[i,]
            for (uword row = 0; row < n_y; ++row) {
                for (uword col = 0; col < n_x; ++col) {
                    A(row, row * n_x + col) = x(i, col);
                }
                b(row) = y(i, row);
            }

            // Solve for minimum and maximum using two-phase simplex
            LPResult sol_min = solve_lp_simplex(c, A, b);
            LPResult sol_max = solve_lp_simplex(-c, A, b);

            if (sol_min.status == LPStatus::OPTIMAL) {
                res_min(i, j) = sol_min.objval + obj_offset;
            } else {
                res_min(i, j) = datum::nan;
            }

            if (sol_max.status == LPStatus::OPTIMAL) {
                res_max(i, j) = -sol_max.objval + obj_offset;
            } else {
                res_max(i, j) = datum::nan;
            }
        }
    }

    // Apply inverse transformation
    res_min = res_min * scale + shift;
    res_max = res_max * scale + shift;

    return std::make_tuple(res_min, res_max);
}
