#' Simplex Algorithm for Linear Programming
#'
#' Solves linear programming problems using the simplex algorithm.
#'
#' @param a Numeric vector of objective function coefficients.
#' @param A1 Matrix of coefficients for <= constraints (or NULL).
#' @param b1 Numeric vector of right-hand sides for <= constraints.
#' @param A2 Matrix of coefficients for >= constraints (or NULL).
#' @param b2 Numeric vector of right-hand sides for >= constraints.
#' @param A3 Matrix of coefficients for = constraints (or NULL).
#' @param b3 Numeric vector of right-hand sides for = constraints.
#' @param maxi Logical; if TRUE, maximizes the objective (default FALSE).
#' @param n_iter Maximum number of iterations (default: n + 2*m).
#' @param eps Tolerance for numerical comparisons (default: 1e-10).
#'
#' @return A list with components:
#' \describe{
#'   \item{soln}{Solution vector}
#'   \item{solved}{Status: 1 = solved, 0 = max iterations, -1 = infeasible}
#'   \item{value}{Optimal objective value}
#'   \item{maxi}{Whether maximization was used}
#'   \item{slack}{Slack variables (if A1 provided)}
#'   \item{surplus}{Surplus variables (if A2 provided)}
#'   \item{artificial}{Artificial variables (if stage 1)}
#'   \item{obj}{Original objective coefficients}
#' }
#'
#' @export
simplex <- function(a,
                    A1 = NULL, b1 = NULL,
                    A2 = NULL, b2 = NULL,
                    A3 = NULL, b3 = NULL,
                    maxi = FALSE,
                    n_iter = NULL,
                    eps = 1e-10) {
    
    # Convert vectors to matrices if needed
    if (!is.null(A1) && !is.matrix(A1)) {
        A1 <- matrix(A1, nrow = 1)
    }
    if (!is.null(A2) && !is.matrix(A2)) {
        A2 <- matrix(A2, nrow = 1)
    }
    if (!is.null(A3) && !is.matrix(A3)) {
        A3 <- matrix(A3, nrow = 1)
    }
    
    # Calculate default n_iter if not provided
    n <- length(a)
    m <- 0
    if (!is.null(A1)) m <- m + nrow(A1)
    if (!is.null(A2)) m <- m + nrow(A2)
    if (!is.null(A3)) m <- m + nrow(A3)
    
    if (is.null(n_iter)) {
        n_iter <- n + 2 * m
    }
    
    # Call C++ implementation
    result <- simplex_cpp(
        a = a,
        A1 = A1,
        b1 = b1,
        A2 = A2,
        b2 = b2,
        A3 = A3,
        b3 = b3,
        maxi = maxi,
        n_iter = as.integer(n_iter),
        eps = eps
    )
    
    # Add names to solution
    names(result$soln) <- paste0("x", seq_len(n))
    names(result$obj) <- paste0("x", seq_len(n))
    
    # Set class
    class(result) <- "simplex"
    
    result
}
