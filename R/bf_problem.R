#' Count The r-by-n 0-1 Matrices Without Adjacent 1's in Any Row or Column
#'
#' @param r a positive integer
#' @param n a positive integer
#' @param vec logical, default value is TRUE
#'
#' @return a number (the n-th number) or a vector (\code{the_re[1:n]})
#' @importFrom utils combn
#' @export
#'
#' @examples
#' (count_grids(r = 1, n = 6))
#' (count_grids(r = 2, n = 6))
#' (count_grids(r = 3, n = 6))
#' (count_grids(r = 4, n = 6))
#' (count_grids(r = 5, n = 6))
#' (count_grids(r = 6, n = 6))
#' (count_grids(r = 7, n = 6))
#' (count_grids(r = 8, n = 6))
#' (count_grids(r = 9, n = 6))
count_grids <- function(r = 1, n = 1, vec = TRUE)
{the_re <- rep(0, n)
if(r == 1) {
  the_re <- vapply(1:n, function(i) fibo(i + 2), numeric(1)) } else {
    m <- fibo_matrix(r)
    d <- dim(m)[1]
    the_re[1] <- d
    the_m <- diag(d)
    for(i in 1:(n-1)) {
      the_m <- the_m %*% m
      the_re[i + 1] <- sum(the_m)
    }
  }
if(vec == TRUE) {
  return(the_re) } else {
    return(the_re[n])
  }
}

#' Calculate Coefficients of P(x)
#'
#' @param s_sq s-squared, a positive number
#' @param eigen_vec a vector
#'
#' @return a vector
#' @export
#'
calc_coef <- function(s_sq, eigen_vec)
{n <- length(eigen_vec)
 the_terms <- rep(0, n)
 for(i in 1:n) {
  the_combn_m <- utils::combn(eigen_vec, i)
  the_terms[i] <- sum(apply(the_combn_m, 2, prod)) * (-1)^i
 }
 the_re <- s_sq * c(1, the_terms)
}


#' Check of P(X) and Q(x)
#'
#' @param r row number
#' @param tol tolerance
#' @param P_of_x_coef Coefficients of P(x)
#' @param Q_of_x_coef Coefficients of Q(x)
#'
#' @return a list
#' @export
#'
#' @examples
#' check_P_and_Q_of_x(r = 5,
#'                    P_of_x_coef = c(13, 47, -37, -129, 68, 49, -23, -3, 1),
#'                    Q_of_x_coef = c(1, -4, -36, 0, 105, -15, -64, 20, 4, -1))
#'
#' check_P_and_Q_of_x(r = 6,
#'                    P_of_x_coef = c(21, 71, -215, -385, 668, 234, -400, 9, 49, -3, -1),
#'                    Q_of_x_coef = c(1, -8, -62, 78, 375, -300, -486, 385, 30, -52, 2, 1))
#'


check_P_and_Q_of_x <- function(r, tol = 1e-8,
                               P_of_x_coef,
                               Q_of_x_coef)
{the_m <- fibo_matrix(r)
 the_eigen <- eigen(the_m)
 the_s_sq <- Rmpfr::mpfr(colSums(the_eigen[[2]])^2, 80)
 which_to_keep <- which(the_s_sq > tol)
 the_top <- the_s_sq[which_to_keep]
 the_b <- the_eigen[[1]][which_to_keep]

 N <- length(the_b)
 the_re_P <- Rmpfr::mpfr(rep(0, N), 80)
 for(i in 1:N) {
   the_v <- the_b[-i]
   cal_v <- calc_coef(the_top[i], the_v)
   the_re_P <- the_re_P + cal_v
 }
 the_re_Q <- Rmpfr::mpfr(calc_coef(1, the_b), 80)
 diff_1 <- the_re_P - P_of_x_coef
 diff_2 <- the_re_Q - Q_of_x_coef
 chk_Px <- all(abs(diff_1) < tol)
 chk_Qx <- all(abs(diff_2) < tol)

 the_chk <- chk_Px && chk_Qx

 return(list("the_check" = the_chk,
             "top_diff" = max(abs(diff_1)),
             "bottom_diff" = max(abs(diff_2))))
}
