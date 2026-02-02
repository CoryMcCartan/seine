// (c) 2025 Cory McCARTAN

#include "simplex.h"

namespace simplex {

// Create zero matrix
arma::mat zero_mat(int n, int m) {
    if (n > 0 && m > 0) {
        return arma::zeros<arma::mat>(n, m);
    }
    return arma::mat();
}

// Create identity matrix
arma::mat identity_mat(int n) {
    if (n > 0) {
        return arma::eye<arma::mat>(n, n);
    }
    return arma::mat();
}

// Perform tableau pivot operation
arma::mat pivot(const arma::mat& tab, int pivot_row, int pivot_col) {
    arma::mat result = tab;
    double pivot_val = tab(pivot_row, pivot_col);
    arma::vec pivot_col_vec = tab.col(pivot_col);
    
    // Update all rows except pivot row
    for (arma::uword i = 0; i < tab.n_rows; i++) {
        if (i != static_cast<arma::uword>(pivot_row)) {
            result.row(i) = tab.row(i) - (tab(i, pivot_col) / pivot_val) * tab.row(pivot_row);
        }
    }
    
    // Update pivot row
    result.row(pivot_row) = tab.row(pivot_row) / (-pivot_val);
    result(pivot_row, pivot_col) = 1.0 / pivot_val;
    
    // Update pivot column (except pivot row)
    for (arma::uword i = 0; i < tab.n_rows; i++) {
        if (i != static_cast<arma::uword>(pivot_row)) {
            result(i, pivot_col) = pivot_col_vec(i) / pivot_val;
        }
    }
    
    return result;
}

// Core simplex solver
SimplexResult simplex_core(
    const arma::vec& obj_coef,
    const arma::mat& constr_mat,
    const arma::vec& rhs,
    const arma::vec& init,
    arma::uvec basic,
    double obj_val,
    int stage,
    int n_orig,
    double eps,
    int max_iter
) {
    int n_vars = constr_mat.n_cols;
    int n_constr = constr_mat.n_rows;
    
    if (n_orig == -1) {
        n_orig = n_vars;
    }
    if (max_iter == -1) {
        max_iter = n_orig;
    }
    
    // Create nonbasic indices (complement of basic)
    arma::uvec nonbasic(n_vars - n_constr);
    int nb_idx = 0;
    for (int i = 0; i < n_vars; i++) {
        bool is_basic = false;
        for (arma::uword j = 0; j < basic.n_elem; j++) {
            if (basic(j) == static_cast<arma::uword>(i)) {
                is_basic = true;
                break;
            }
        }
        if (!is_basic) {
            nonbasic(nb_idx++) = i;
        }
    }
    
    // Build initial tableau
    arma::mat tab(n_constr, nonbasic.n_elem + 1);
    tab.col(0) = rhs;
    tab.cols(1, nonbasic.n_elem) = -constr_mat.cols(nonbasic);
    
    arma::vec obj_row;
    
    if (stage == 2) {
        // Stage 2: standard simplex
        tab = arma::join_cols(tab, arma::zeros<arma::rowvec>(nonbasic.n_elem + 1));
        tab(n_constr, 0) = obj_val;
        for (arma::uword i = 0; i < nonbasic.n_elem; i++) {
            tab(n_constr, i + 1) = obj_coef(nonbasic(i));
        }
        obj_row = tab.row(n_constr).cols(1, nonbasic.n_elem).t();
    } else {
        // Stage 1: auxiliary problem
        int n_aux_start = n_constr + n_orig - n_vars;
        arma::vec aux_obj_row = arma::sum(tab.rows(n_aux_start, n_constr - 1), 0).t();
        
        arma::mat top_row = arma::zeros<arma::rowvec>(nonbasic.n_elem + 1);
        top_row(0) = obj_val;
        for (arma::uword i = 0; i < nonbasic.n_elem; i++) {
            top_row(i + 1) = obj_coef(nonbasic(i));
        }
        
        tab = arma::join_cols(top_row, tab);
        tab = arma::join_cols(tab, aux_obj_row.t());
        
        obj_row = aux_obj_row.subvec(1, aux_obj_row.n_elem - 1);
    }
    
    // Main simplex iterations
    int iter = 1;
    while (!arma::all(obj_row > -eps) && iter <= max_iter) {
        // Find entering column (most negative objective coefficient)
        arma::uword min_idx = obj_row.index_min();
        int pivot_col = min_idx + 1;
        
        // Find leaving row (minimum ratio test)
        arma::uvec neg_indices;
        if (stage == 2) {
            neg_indices = arma::find(tab.submat(0, pivot_col, n_constr - 1, pivot_col) < -eps);
        } else {
            neg_indices = arma::find(tab.submat(1, pivot_col, n_constr, pivot_col) < -eps);
            if (neg_indices.n_elem > 0) {
                neg_indices = neg_indices + 1;
            }
        }
        
        if (neg_indices.n_elem == 0) {
            break; // unbounded
        }
        
        arma::vec ratios(neg_indices.n_elem);
        for (arma::uword i = 0; i < neg_indices.n_elem; i++) {
            ratios(i) = -tab(neg_indices(i), 0) / tab(neg_indices(i), pivot_col);
        }
        
        arma::uword min_ratio_idx = ratios.index_min();
        int pivot_row = neg_indices(min_ratio_idx);
        
        // Perform pivot
        tab = pivot(tab, pivot_row, pivot_col);
        
        // Update basic/nonbasic sets
        int temp;
        if (stage == 1) {
            temp = basic(pivot_row - 1);
            basic(pivot_row - 1) = nonbasic(pivot_col - 1);
            nonbasic(pivot_col - 1) = temp;
            obj_row = tab.row(n_constr + 1).cols(1, nonbasic.n_elem).t();
        } else {
            temp = basic(pivot_row);
            basic(pivot_row) = nonbasic(pivot_col - 1);
            nonbasic(pivot_col - 1) = temp;
            obj_row = tab.row(n_constr).cols(1, nonbasic.n_elem).t();
        }
        
        iter++;
    }
    
    SimplexResult result;
    
    if (stage == 1) {
        // Stage 1 results
        result.val_aux = tab(n_constr + 1, 0);
        
        // Check if auxiliary variables still in basis
        if (result.val_aux < eps) {
            for (arma::uword j = 0; j < basic.n_elem; j++) {
                if (basic(j) >= static_cast<arma::uword>(n_orig)) {
                    // Pivot out artificial variable
                    int pivot_row = j + 1;
                    arma::rowvec row_vals = arma::abs(tab.row(pivot_row).cols(1, nonbasic.n_elem));
                    arma::uvec valid = arma::find(row_vals > eps);
                    if (valid.n_elem > 0) {
                        int pivot_col = arma::index_min(nonbasic(valid)) + 1;
                        tab = pivot(tab, pivot_row, pivot_col);
                        
                        int temp = basic(pivot_row - 1);
                        basic(pivot_row - 1) = nonbasic(pivot_col - 1);
                        nonbasic(pivot_col - 1) = temp;
                    }
                }
            }
        }
        
        result.soln = arma::zeros<arma::vec>(n_vars);
        for (arma::uword i = 0; i < basic.n_elem; i++) {
            result.soln(basic(i)) = tab(i + 1, 0);
        }
        result.solved = -1;
        result.value = tab(0, 0);
        
        result.constr_mat = arma::zeros<arma::mat>(n_constr, n_vars);
        result.constr_mat.cols(basic) = identity_mat(n_constr);
        result.constr_mat.cols(nonbasic) = -tab.submat(1, 1, n_constr, nonbasic.n_elem);
        
        result.obj_coef = arma::zeros<arma::vec>(n_vars);
        result.obj_coef(nonbasic) = tab.row(0).cols(1, nonbasic.n_elem).t();
        
        result.obj_aux = arma::zeros<arma::vec>(n_vars);
        result.obj_aux(nonbasic) = tab.row(n_constr + 1).cols(1, nonbasic.n_elem).t();
        
        result.basic = basic;
    } else {
        // Stage 2 results
        result.soln = arma::zeros<arma::vec>(n_vars);
        for (arma::uword i = 0; i < basic.n_elem; i++) {
            result.soln(basic(i)) = tab(i, 0);
        }
        result.value = tab(n_constr, 0);
        
        result.constr_mat = arma::zeros<arma::mat>(n_constr, n_vars);
        result.constr_mat.cols(basic) = identity_mat(n_constr);
        result.constr_mat.cols(nonbasic) = tab.submat(0, 1, n_constr - 1, nonbasic.n_elem);
        
        result.obj_coef = arma::zeros<arma::vec>(n_vars);
        result.obj_coef(nonbasic) = tab.row(n_constr).cols(1, nonbasic.n_elem).t();
        
        result.solved = (iter <= max_iter) ? 1 : 0;
        result.basic = basic;
    }
    
    return result;
}

// Main simplex wrapper
SimplexResult simplex(
    const arma::vec& obj_coef,
    const arma::mat& A1,
    const arma::vec& b1,
    const arma::mat& A2,
    const arma::vec& b2,
    const arma::mat& A3,
    const arma::vec& b3,
    bool maximize,
    int max_iter,
    double eps
) {
    int m1 = (A1.n_elem > 0) ? A1.n_rows : 0;
    int m2 = (A2.n_elem > 0) ? A2.n_rows : 0;
    int m3 = (A3.n_elem > 0) ? A3.n_rows : 0;
    int m = m1 + m2 + m3;
    int n = obj_coef.n_elem;
    
    if (max_iter == -1) {
        max_iter = n + 2 * m;
    }
    
    arma::vec obj = obj_coef;
    if (maximize) {
        obj = -obj;
    }
    
    SimplexResult result;
    
    if (m2 + m3 == 0) {
        // Only <= constraints, use stage 2 directly
        arma::vec extended_obj = arma::join_cols(obj, arma::zeros<arma::vec>(m1));
        arma::mat extended_A = arma::join_rows(A1, identity_mat(m1));
        arma::vec init = arma::join_cols(arma::zeros<arma::vec>(n), b1);
        arma::uvec basic = arma::linspace<arma::uvec>(n, n + m1 - 1, m1);
        
        result = simplex_core(extended_obj, extended_A, b1, init, basic, 0.0, 2, n + m1, eps, max_iter);
    } else {
        // Need stage 1 (auxiliary problem)
        if (m2 > 0) {
            // Has >= constraints, need surplus and artificial variables
            arma::vec extended_obj = arma::join_cols(obj, arma::zeros<arma::vec>(m1 + 2*m2 + m3));
            
            // Build constraint matrix
            arma::mat A_all;
            if (m1 > 0 && m3 > 0) {
                A_all = arma::join_cols(arma::join_cols(A1, A2), A3);
            } else if (m1 > 0) {
                A_all = arma::join_cols(A1, A2);
            } else if (m3 > 0) {
                A_all = arma::join_cols(A2, A3);
            } else {
                A_all = A2;
            }
            
            // Build slack/surplus/artificial variable matrix
            arma::mat slack_part, surplus_part, artificial_part;
            
            // Slack variables column (for <= constraints)
            if (m1 > 0) {
                slack_part = arma::join_cols(identity_mat(m1), zero_mat(m2 + m3, m1));
            } else {
                slack_part = zero_mat(m, m1);
            }
            
            // Surplus variables column (for >= constraints, negative)
            if (m3 > 0) {
                surplus_part = arma::join_cols(zero_mat(m1, m2), 
                                              arma::join_cols(-identity_mat(m2), zero_mat(m3, m2)));
            } else {
                surplus_part = arma::join_cols(zero_mat(m1, m2), -identity_mat(m2));
            }
            
            // Artificial variables column (for >= and = constraints)
            artificial_part = arma::join_cols(zero_mat(m1, m2 + m3), identity_mat(m2 + m3));
            
            arma::mat slack_surplus_artificial;
            if (m1 > 0) {
                slack_surplus_artificial = arma::join_rows(arma::join_rows(slack_part, surplus_part), artificial_part);
            } else {
                slack_surplus_artificial = arma::join_rows(surplus_part, artificial_part);
            }
            
            arma::mat extended_A = arma::join_rows(A_all, slack_surplus_artificial);
            
            // Build RHS vector
            arma::vec rhs;
            if (m1 > 0 && m3 > 0) {
                rhs = arma::join_cols(arma::join_cols(b1, b2), b3);
            } else if (m1 > 0) {
                rhs = arma::join_cols(b1, b2);
            } else if (m3 > 0) {
                rhs = arma::join_cols(b2, b3);
            } else {
                rhs = b2;
            }
            
            // Build initial solution vector
            arma::vec init;
            if (m1 > 0 && m3 > 0) {
                init = arma::join_cols(arma::join_cols(arma::zeros<arma::vec>(n), b1),
                                      arma::join_cols(arma::zeros<arma::vec>(m2), 
                                                    arma::join_cols(b2, b3)));
            } else if (m1 > 0) {
                init = arma::join_cols(arma::zeros<arma::vec>(n), 
                                      arma::join_cols(b1, arma::join_cols(arma::zeros<arma::vec>(m2), b2)));
            } else if (m3 > 0) {
                init = arma::join_cols(arma::zeros<arma::vec>(n), 
                                      arma::join_cols(arma::zeros<arma::vec>(m2), arma::join_cols(b2, b3)));
            } else {
                init = arma::join_cols(arma::zeros<arma::vec>(n), arma::join_cols(arma::zeros<arma::vec>(m2), b2));
            }
            
            arma::uvec basic(m);
            for (int i = 0; i < m1; i++) {
                basic(i) = n + i;
            }
            for (int i = 0; i < m2 + m3; i++) {
                basic(m1 + i) = n + m1 + m2 + i;
            }
            
            SimplexResult stage1_result = simplex_core(extended_obj, extended_A, rhs, init, basic, 
                                                      0.0, 1, n + m1 + m2, eps, max_iter);
            
            if (stage1_result.val_aux > eps) {
                // Infeasible
                result = stage1_result;
            } else {
                // Continue to stage 2
                arma::vec stage2_obj = stage1_result.obj_coef.subvec(0, n + m1 + m2 - 1);
                arma::mat stage2_A = stage1_result.constr_mat.cols(0, n + m1 + m2 - 1);
                arma::vec stage2_init = stage1_result.soln.subvec(0, n + m1 + m2 - 1);
                arma::vec stage2_rhs = arma::zeros<arma::vec>(basic.n_elem);
                for (arma::uword i = 0; i < basic.n_elem; i++) {
                    stage2_rhs(i) = stage1_result.soln(stage1_result.basic(i));
                }
                
                result = simplex_core(stage2_obj, stage2_A, stage2_rhs, stage2_init, 
                                     stage1_result.basic, stage1_result.value, 2, n + m1 + m2, eps, max_iter);
            }
        } else {
            // Only = constraints (no >= constraints)
            arma::vec extended_obj = arma::join_cols(obj, arma::zeros<arma::vec>(m1 + m3));
            
            arma::mat A_all;
            arma::vec rhs, init;
            if (m1 > 0 && m3 > 0) {
                A_all = arma::join_cols(A1, A3);
                rhs = arma::join_cols(b1, b3);
                init = arma::join_cols(arma::zeros<arma::vec>(n), arma::join_cols(b1, b3));
            } else if (m1 > 0) {
                A_all = A1;
                rhs = b1;
                init = arma::join_cols(arma::zeros<arma::vec>(n), b1);
            } else {
                A_all = A3;
                rhs = b3;
                init = arma::join_cols(arma::zeros<arma::vec>(n), b3);
            }
            
            arma::mat extended_A = arma::join_rows(A_all, identity_mat(m1 + m3));
            arma::uvec basic = arma::linspace<arma::uvec>(n, n + m1 + m3 - 1, m1 + m3);
            
            SimplexResult stage1_result = simplex_core(extended_obj, extended_A, rhs, init, basic,
                                                      0.0, 1, n + m1, eps, max_iter);
            
            if (stage1_result.val_aux > eps) {
                // Infeasible
                result = stage1_result;
            } else {
                // Continue to stage 2
                arma::vec stage2_obj = stage1_result.obj_coef.subvec(0, n + m1 - 1);
                arma::mat stage2_A = stage1_result.constr_mat.cols(0, n + m1 - 1);
                arma::vec stage2_init = stage1_result.soln.subvec(0, n + m1 - 1);
                arma::vec stage2_rhs = arma::zeros<arma::vec>(basic.n_elem);
                for (arma::uword i = 0; i < basic.n_elem; i++) {
                    stage2_rhs(i) = stage1_result.soln(stage1_result.basic(i));
                }
                
                result = simplex_core(stage2_obj, stage2_A, stage2_rhs, stage2_init,
                                     stage1_result.basic, stage1_result.value, 2, n + m1, eps, max_iter);
            }
        }
    }
    
    if (maximize) {
        result.value = -result.value;
    }
    
    return result;
}

} // namespace simplex
