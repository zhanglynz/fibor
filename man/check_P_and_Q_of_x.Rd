% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bf_problem.R
\name{check_P_and_Q_of_x}
\alias{check_P_and_Q_of_x}
\title{Check of P(X) and Q(x)}
\usage{
check_P_and_Q_of_x(r, tol = 1e-08, P_of_x_coef, Q_of_x_coef)
}
\arguments{
\item{r}{row number}

\item{tol}{tolerance}

\item{P_of_x_coef}{Coefficients of P(x)}

\item{Q_of_x_coef}{Coefficients of Q(x)}
}
\value{
a list
}
\description{
Check of P(X) and Q(x)
}
\examples{
check_P_and_Q_of_x(r = 5,
                   P_of_x_coef = c(13, 47, -37, -129, 68, 49, -23, -3, 1),
                   Q_of_x_coef = c(1, -4, -36, 0, 105, -15, -64, 20, 4, -1))

check_P_and_Q_of_x(r = 6,
                   P_of_x_coef = c(21, 71, -215, -385, 668, 234, -400, 9, 49, -3, -1),
                   Q_of_x_coef = c(1, -8, -62, 78, 375, -300, -486, 385, 30, -52, 2, 1))

}
