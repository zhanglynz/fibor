#' Fibonacci Numbers
#'
#' @param i a non-negative integer
#'
#' @return the ith Fibonacci number
#' @export
#'
#' @examples
#' sapply(0:5, fibo)
fibo <- function(i)
{if(i <= 2) {if(i > 0) 1 else 0
 } else {
  Recall(i - 1) + Recall(i - 2) }
}


#' Fibonacci Phi-Matrix Series
#'
#' @param n a positive integer
#'
#' @return the n-the Fibonacci Phi-matrix
#' @export
#'
#' @examples
#' (Phi_1 <- fibo_matrix(1))
#' (Phi_2 <- fibo_matrix(2))
#' (Phi_3 <- fibo_matrix(3))
fibo_matrix <- function(n = 1)
{d <- fibo(n + 2)
 fibo_m <- matrix(0, d, d)
 fibo_m[1, 1:2] <- c(1, 1)
 fibo_m[2, 1:2] <- c(1, 0)
 if(n >= 2) {
   d0 <- 2
   for(i in 2:n)
     {d1 <- fibo(i + 2)
      fibo_m[1:d0, (d0 + 1):d1] <- fibo_m[1:d0, 1:(fibo(i))]
      fibo_m[(d0 + 1):d1, 1:d0] <- fibo_m[1:(fibo(i)), 1:d0]
      d0 <- d1
     }
 }
 return(fibo_m)
}
