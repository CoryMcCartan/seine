#include "bounds.h"
#include <cmath>
#include <algorithm>
#include <tuple>

using namespace arma;

/**
 * Compute bounds on entries of matrix B satisfying B * x = y
 *
 * For each row of x and y, solves LP to minimize/maximize each entry of B
 * subject to: B %*% x = y, B in [bounds[0], bounds[1]]
 */
std::tuple<mat, mat> bounds_lp(const mat& x, const mat& y, const vec& bounds) {
    uword n = x.n_rows;
    uword p = x.n_cols;
    uword q = y.n_cols;
    
    // Initialize output matrices
    // B is q x p, so we have q*p entries in row-major order
    mat result_min(n, q * p);
    mat result_max(n, q * p);
    
    double lower_bound = bounds(0);
    double upper_bound = bounds(1);
    
    // Solve for each row
    for (uword i = 0; i < n; ++i) {
        rowvec xi = x.row(i);
        rowvec yi = y.row(i);
        
        // Solve for each entry of B (row-major order)
        // Row-major: B[0,0], B[0,1], ..., B[0,p-1], B[1,0], B[1,1], ..., B[1,p-1], ..., B[q-1,p-1]
        for (uword idx = 0; idx < q * p; ++idx) {
            // Convert linear index to (row, col) in row-major order
            uword b_row = idx / p;
            uword b_col = idx % p;
            
            // Compute min and max for B[b_row, b_col] in parallel
            // We want to optimize B[b_row, b_col] subject to:
            // 1. For each row j: sum_i B[j, i] * x[i] = y[j]
            // 2. bounds[0] <= B[j, i] <= bounds[1]
            
            // Start with box constraints
            double min_val = lower_bound;
            double max_val = upper_bound;
            
            if (std::abs(xi(b_col)) > 1e-8) {
                // From row constraint: sum_i B[b_row, i] * xi[i] = yi[b_row]
                // B[b_row, b_col] * xi[b_col] + sum_{i != b_col} B[b_row, i] * xi[i] = yi[b_row]
                // B[b_row, b_col] = (yi[b_row] - sum_{i != b_col} B[b_row, i] * xi[i]) / xi[b_col]
                // To minimize B[b_row, b_col]: maximize sum_{i != b_col} B[b_row, i] * xi[i]
                // To maximize B[b_row, b_col]: minimize sum_{i != b_col} B[b_row, i] * xi[i]
                
                double sum_max = 0.0;
                double sum_min = 0.0;
                
                for (uword i_col = 0; i_col < p; ++i_col) {
                    if (i_col == b_col) continue;
                    if (xi(i_col) > 0) {
                        sum_max += upper_bound * xi(i_col);
                        sum_min += lower_bound * xi(i_col);
                    }
                }
                
                double min_from_row = (yi(b_row) - sum_max) / xi(b_col);
                double max_from_row = (yi(b_row) - sum_min) / xi(b_col);
                
                min_val = std::max(min_val, min_from_row);
                max_val = std::min(max_val, max_from_row);
            }
            
            // Check feasibility
            if (!std::isfinite(min_val) || !std::isfinite(max_val) || min_val > max_val + 1e-10) {
                result_min(i, idx) = datum::nan;
                result_max(i, idx) = datum::nan;
            } else {
                result_min(i, idx) = min_val;
                result_max(i, idx) = max_val;
            }
        }
    }
    
    return std::make_tuple(result_min, result_max);
}
