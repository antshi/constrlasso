#' @title The function constrlasso
#'
#' @description This function fits a linearly constrained lasso regression,
#' using a predictor matrix X, a response y and a tuning parameter value lambda.
#' It results in a vector of coefficient estimates.
#' The function corresponds to lsq_constrsparsereg of the SparseReg MATLAB-toolbox by Zhou and Gaines,
#' see \href{http://hua-zhou.github.io/SparseReg/}{the project page}.
#' Only the Quadratic Programming Algorithm for ENET is implemented.
#' For more information, see
#' \href{http://hua-zhou.github.io/media/pdf/GainesKimZhou08CLasso.pdf}{\insertCite{gaines2018algorithms;textual}{constrlasso}}.
#'
#' @param X an nxp matrix of p regressors with n observations.
#' @param y an nx1 response vector with n observations.
#' @param lambda a tuning parameter value for the lasso penalty. Default value is lambda=0.
#' @param Aeq a cxp equality constraint matrix, containing c constraints for p regressors.
#' Default value is Aeq=NULL, no equality constraints.
#' @param beq a cx1 equality constraint vector. Default value is beq=NULL, no equality constraints.
#' @param A a cxp inequality constraint matrix, containing c constraints for p regressors.
#' Default value is A=NULL, no inequality constraints.
#' @param b a cx1 inequality constraint vector. Default value is b=NULL, no inequality constraints.
#' @param penidx a logical px1 vector, indicating which coefficients are to be penalized.
#' Default value is penidx=NULL and imposes penalty on all p coefficients.
#' @param method a character string, the method to be used.
#' Possible values are "QP" (default) for Quadratic Programming and "CD" for Coordinate Descent.
#' "CD" uses the glmnet package and only works without equality and inequality constraints.
#'
#' @section Details:
#' The Constrained Lasso as in \insertCite{gaines2018algorithms;textual}{constrlasso} minimizes
#' \deqn{0.5||y - X \beta ||^2_2 + \lambda||\beta||_1,}
#' subject to \eqn{Aeq \beta = beq} and \eqn{A \beta\le b}.
#'
#' @return betahat a px1 vector of estimated coefficients.
#' @return dual_eq equality duals. The function returns an empty vector for no equality constraints.
#' @return dual_neq inequality duals. The function returns an empty vector for no inequality constraints.
#'
#' @examples
#' library(constrlasso)
#' library(MASS)
#' set.seed(1234)
#' n <- 200
#' p <- 50
#' Xmat <- matrix(, n, p)
#' for (i in 1:p) {
#'   Xmat[, i] <- rnorm(n, runif(1, -3, 3), runif(1, 1, 2))
#' }
#' betas <- runif(p, -2, 2)
#' nonzeros <- sample(1:p, 20, replace = FALSE)
#' yvec <- Xmat[, nonzeros] %*% betas[nonzeros] + rnorm(n, 0, 2)
#' classoreg_results <- constrlasso(Xmat, yvec, lambda = 0)
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited
#'
#'
#' @export constrlasso
constrlasso <-
  function(X,
           y,
           lambda,
           Aeq = NULL,
           beq = NULL,
           A = NULL,
           b = NULL,
           penidx = NULL,
           method = "QP") {
    X <- as.matrix(X)
    n <- dim(X)[1]
    p <- dim(X)[2]
    y <- as.matrix(y)
    dim(y) <- c(n, 1)

    if (is.null(penidx)) {
      penidx <- matrix(TRUE, p, 1)
    }
    dim(penidx) <- c(p, 1)

    if (is.null(A)) {
      A <- matrix(0, 0, p)
      b <- rep(0, 0)
    }

    if (is.null(Aeq)) {
      Aeq <- matrix(0, 0, p)
      beq <- rep(0, 0)
    }

    m1 <- dim(Aeq)[1]
    m2 <- dim(A)[1]

    if (is.null(lambda)) {
      stop(cat("Please, enter a value for the penalty parameter lambda."))
    }

    ## no constraints
    if (m1 == 0 && m2 == 0) {
      # no penalization
      if (abs(lambda) < 1e-16) {
        wt <- matrix(1, n, 1)
        dim(wt) <- c(n, 1)
        Xwt <- X * as.numeric(sqrt(wt))
        ywt <- as.numeric(sqrt(wt)) * y

        betahat <- matrix(lm(ywt ~ Xwt - 1)$coef, p, 1)
        dual_eq <- rep(0, 0)
        dual_neq <- rep(0, 0)
      } else {
        # with penalization
        if (method == "CD") {
          glmnet_res <-
            glmnet(
              Xwt,
              ywt,
              alpha = 1,
              lambda = lambda,
              penalty.factor = penidx,
              maxit = 1000,
              intercept = FALSE,
              standardize = TRUE
            )

          betahat <- matrix(glmnet_res$beta, p, 1)
          dual_eq <- rep(0, 0)
          dual_neq <- rep(0, 0)
        } else if (method == "QP") {
          wt <- matrix(1, n, 1)
          dim(wt) <- c(n, 1)
          Xwt <- X * as.numeric(sqrt(wt))
          ywt <- as.numeric(sqrt(wt)) * y

          # quadratic coefficient
          H <- t(Xwt) %*% Xwt
          H <- rbind(cbind(H, -H), cbind(-H, H))

          # linear coefficient
          f <- -t(Xwt) %*% ywt
          f <- rbind(f, -f) + lambda * rbind(penidx, penidx)

          # optimizer
          x <- OP(Q_objective(H, L = t(f)))
          opt <- ROI_solve(x, solver = "qpoases")
          opt_sol <- opt$message$primal_solution

          # estimators
          betahat <-
            matrix(opt_sol[1:p] - opt_sol[(p + 1):length(opt_sol)], p, 1)
          dual_eq <- rep(0, 0)
          dual_neq <- rep(0, 0)
        }
      }
    } else {
      ## with constraints

      if (method == "CD") {
        warning("The CD method does not work with constraints. The solution is generated with QP.")
      }

      wt <- matrix(1, n, 1)
      dim(wt) <- c(n, 1)
      Xwt <- X * as.numeric(sqrt(wt))
      ywt <- as.numeric(sqrt(wt)) * y

      # quadratic coefficient
      H <- t(Xwt) %*% Xwt
      H <- rbind(cbind(H, -H), cbind(-H, H))

      # linear coefficient
      f <- -t(Xwt) %*% ywt
      f <- rbind(f, -f) + lambda * rbind(penidx, penidx)

      # constraints
      Amatrix <- rbind(cbind(Aeq, -Aeq), cbind(A, -A))
      bvector <- c(beq, b)

      # optimizer
      x <-
        OP(
          Q_objective(H, L = t(f)),
          L_constraint(
            L = Amatrix,
            dir = c(rep("==", m1), rep("<=", m2)),
            rhs = bvector
          )
        )
      opt <- ROI_solve(x, solver = "qpoases")
      opt_sol <- opt$message$primal_solution

      # estimators
      betahat <-
        matrix(opt_sol[1:p] - opt_sol[(p + 1):length(opt_sol)], p, 1)
      duals <- opt$message$dual_solution[-(1:(2 * p))]

      if (m1 != 0) {
        dual_eq <- -matrix(duals[1:m1], m1, 1)
      } else {
        dual_eq <- rep(0, 0)
      }
      if (m2 != 0) {
        dual_neq <- -matrix(duals[(m1 + 1):(m1 + m2)], m2, 1)
      } else {
        dual_neq <- rep(0, 0)
      }
    }

    return(list(
      "betahat" = betahat,
      "dual_eq" = dual_eq,
      "dual_neq" = dual_neq
    ))
  }


#' @title The function constrlasso_path
#'
#' @description This function performs a Constrained Lasso Solution Path
#' as in \insertCite{gaines2018algorithms;textual}{constrlasso}.
#' It computes the solution path for the constrained lasso problem, using a predictor matrix X and a response y.
#' The constrained lasso solves the standard lasso \insertCite{tibshirani1996regression;textual}{constrlasso}
#' subject to the linear equality constraints \eqn{Aeq \beta = beq}
#' and linear inequality constraints \eqn{A \beta \le b}.
#' The result lambda_path contains the values of the tuning parameter along the solution path and beta_path -
#' the estimated regression coefficients for each value of lambda_path.
#' The function corresponds to lsq_classopath of the SparseReg MATLAB-toolbox of by Zhou and Gaines,
#' see \href{http://hua-zhou.github.io/SparseReg/}{the project page}.
#' For more information, see \href{http://hua-zhou.github.io/media/pdf/GainesKimZhou08CLasso.pdf}{\insertCite{gaines2018algorithms;textual}{constrlasso}}.
#'
#' @param X an nxp matrix with p regressors with n observations.
#' @param y an nx1 response vector with n observations.
#' @param Aeq a cxp equality constraint matrix, containing c constraints for p regressors.
#' Default value is Aeq=NULL, no equality constraints.
#' @param beq a cx1 equality constraint vector. Default value is beq=NULL, no equality constraints.
#' @param A a cxp inequality constraint matrix, containing c constraints for p regressors.
#' Default value is A=NULL, no inequality constraints.
#' @param b a cx1 inequality constraint vector. Default value is b=NULL, no inequality constraints.
#' @param penidx a logical px1 vector, indicating which coefficients are to be penalized.
#' Default value is penidx=NULL and allows all p coefficients to be penalized.
#' @param init_method a character string, the initializing method to be used.
#' Possible values are "QP" (default) for Quadratic Programming and "LP" for Linear Programming.
#' "LP" is recommended only when it's reasonable to assume that all coefficient estimates initialize at zero.
#' @param epsilon a tuning parameter for ridge penalty in case of high-dimensional (n>p) regressors matrix X.
#' Default value is 1e-4.
#' @param stop_lambda_tol a tolerance value for the tuning lasso parameter.
#' The algorithm stops when hitting this tolerance. Default value is 1e-7.
#' @param ceiling_tol a tolerance value for the change in subgradients. Default value is 1e-10.
#' @param zeros_tol a tolerance value for the zero equality of coefficients. Default value is 1e-20.
#' @param verbose a logical parameter. TRUE prints along the constraint lasso solution path. Default value is FALSE.
#'
#' @section Details:
#' The Constrained Lasso as in \insertCite{gaines2018algorithms;textual}{constrlasso} minimizes
#' \deqn{0.5||y - X \beta ||^2_2 + \lambda||\beta||_1,}
#' subject to \eqn{Aeq \beta = beq} and \eqn{A \beta\le b}.
#'
#' @return lambda_path a vector of the tuning parameter values along the solution path.
#' @return beta_path  a matrix with estimated regression coefficients for each value of lambda_path
#' @return df_path a vector with degrees of freedom along the solution path
#' @return objval_path a vector with values of the objective function for each value of lambda_path
#' @examples
#'
#' library(constrlasso)
#' set.seed(1234)
#' n <- 200
#' p <- 50
#' Xmat <- matrix(, n, p)
#' for (i in 1:p) {
#'   Xmat[, i] <- rnorm(n, runif(1, -3, 3), runif(1, 1, 2))
#' }
#' betas <- runif(p, -2, 2)
#' nonzeros <- sample(1:p, 20, replace = FALSE)
#' yvec <- Xmat[, nonzeros] %*% betas[nonzeros] + rnorm(n, 0, 2)
#' classopath_results <- constrlasso_path(Xmat, yvec)
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited
#'
#' @export constrlasso_path
constrlasso_path <-
  function(X,
           y,
           Aeq = NULL,
           beq = NULL,
           A = NULL,
           b = NULL,
           penidx = NULL,
           init_method = "QP",
           epsilon = 1e-4,
           stop_lambda_tol = 1e-7,
           ceiling_tol = 1e-10,
           zeros_tol = 1e-20,
           verbose = FALSE) {
    if (verbose) {
      printer <- print
    } else {
      printer <- function(x) {
      }
    }

    X <- as.matrix(X)
    n <- dim(X)[1]
    n_orig <- n
    p <- dim(X)[2]

    y <- as.matrix(y)
    dim(y) <- c(n, 1)

    if (is.null(penidx)) {
      penidx <- matrix(TRUE, p, 1)
    }
    dim(penidx) <- c(p, 1)

    if (is.null(A)) {
      A <- matrix(0, 0, p)
      b <- rep(0, 0)
    }

    if (is.null(Aeq)) {
      Aeq <- matrix(0, 0, p)
      beq <- rep(0, 0)
    }

    m1 <- dim(Aeq)[1]
    m2 <- dim(A)[1]

    # check if a ridge penalty has to be included
    if (n < p) {
      warning("Adding a small ridge penalty (default is 1e-4) since n < p.")

      if (epsilon <= 0) {
        warning(
          "The ridge tuning parameter epsilon must be positive,
          switching to default value (1e-4)."
        )
        epsilon <- 1e-4
      }

      # create augmented data for n<p with ridge penalty
      y <- rbind(y, matrix(rep.int(0, p)))
      X <- rbind(X, sqrt(epsilon) * diag(1, p))
      n_orig <- n
    } else {
      # make sure X has a full column rank
      R <- qr(X)$qr[1:p, ]
      X_fullrank <-
        sum(matrix(abs(diag(R))) > abs(R[1]) * max(n, p) * .Machine$double.eps)
      X_fullrank <- qr(X)$rank

      if (X_fullrank != p) {
        warning("Adding a small ridge penalty (default is 1e-4) since X is rank deficient")
        if (epsilon <= 0) {
          warning(
            "The ridge tuning parameter epsilon must be positive,
            switching to default value (1e-4)."
          )
          epsilon <- 1e-4
        }
        # create augmented data for n<p with ridge penalty
        y <- rbind(y, matrix(rep.int(0, p)))
        X <- rbind(X, sqrt(epsilon) * diag(1, p))
        n_orig <- n
      }
    }

    #### define iterations and path help parameters

    max_iters <- 5 * (p + m2)
    beta_path <- matrix(0, p, max_iters)
    dualpath_eq <- matrix(0, m1, max_iters)
    dualpath_ineq <- matrix(0, m2, max_iters)
    lambda_path <- matrix(0, 1, max_iters)
    df_path <- matrix(Inf, 1, max_iters)
    objval_path <- matrix(0, 1, max_iters)
    violations_path <- matrix(Inf, 1, max_iters)

    #### initialization
    H <- t(X) %*% X

    if (init_method == "LP") {
      warning("LP is used for initialization, assumes initial solution is unique.")

      obj_fun <- matrix(1, 2 * p, 1)
      lb <- 0
      lb_mat <- diag(1, 2 * p)

      constr_mat <- rbind(cbind(Aeq, -Aeq), cbind(A, -A), lb_mat)
      constr_rhs <- c(beq, b, matrix(rep.int(lb, 2 * p), 2 * p, 1))
      constr_dir <- c(rep("=", m1), rep("<=", m2), rep(">=", 2 * p))

      lpsol <-
        lp(
          direction = "min",
          obj_fun,
          constr_mat,
          constr_dir,
          constr_rhs,
          compute.sens = TRUE
        )
      lpres <- matrix(lpsol$solution, 2 * p, 1)

      beta_path[, 1] <- lpres[1:p] - lpres[(p + 1):(2 * p)]

      dual_eqlin <- matrix(lpsol$duals[1:m1], m1, 1)
      if (m2 != 0) {
        dual_ineqlin <- matrix(lpsol$duals[(m1 + 1):(m1 + m2)], m2, 1)
      } else {
        dual_ineqlin <- rep(0, 0)
      }

      dualpath_eq[, 1] <- -dual_eqlin
      dualpath_ineq[, 1] <- dual_ineqlin
    } else if (init_method == "QP") {
      obj_fun <- matrix(1, 2 * p, 1)
      lb <- 0
      lb_mat <- diag(1, 2 * p)

      constr_mat <- rbind(cbind(Aeq, -Aeq), cbind(A, -A), lb_mat)
      constr_rhs <- c(beq, b, matrix(rep.int(lb, 2 * p), 2 * p, 1))
      constr_dir <- c(rep("=", m1), rep("<=", m2), rep(">=", 2 * p))

      lpsol <-
        lp(
          direction = "min",
          obj_fun,
          constr_mat,
          constr_dir,
          constr_rhs,
          compute.sens = TRUE
        )
      lpres <- matrix(lpsol$solution, 2 * p, 1)

      beta_path[, 1] <- lpres[1:p] - lpres[(p + 1):(2 * p)]

      dual_eqlin <- matrix(lpsol$duals[1:m1], m1, 1)
      if (m2 != 0) {
        dual_ineqlin <- matrix(lpsol$duals[(m1 + 1):(m1 + m2)], m2, 1)
      } else {
        dual_ineqlin <- rep(0, 0)
      }

      dualpath_eq[, 1] <- -dual_eqlin
      dualpath_ineq[, 1] <- dual_ineqlin

      # initialize sets
      dualpath_ineq[which(dualpath_ineq[, 1] < 0), 1] <- 0
      sets_active <- abs(beta_path[, 1]) > 1e-4 | !penidx
      beta_path[!sets_active, 1] <- 0

      # initialize subgradient vector and find lambda_max
      resid <- y - X %*% beta_path[, 1]
      subgrad <-
        t(X) %*% resid - t(Aeq) %*% dualpath_eq[, 1] - t(A) %*% dualpath_ineq[, 1]
      lambda_max <- max(abs(subgrad))

      # Use QP at lambda_max to initialize
      constrlasso_sol <-
        constrlasso(
          X,
          y,
          lambda = lambda_max,
          Aeq = Aeq,
          beq = beq,
          A = A,
          b = b,
          penidx = penidx,
          method = "QP"
        )

      beta_path[, 1] <- constrlasso_sol$betahat
      dualpath_eq[, 1] <- constrlasso_sol$dual_eq
      dualpath_ineq[, 1] <- constrlasso_sol$dual_neq
    }

    # initialize sets
    dualpath_ineq[which(dualpath_ineq[, 1] < 0), 1] <- 0
    sets_active <- abs(beta_path[, 1]) > 1e-4 | !penidx
    beta_path[!sets_active, 1] <- 0

    resid_ineq <- A %*% beta_path[, 1] - b
    set_ineq_border <- resid_ineq == 0
    n_ineq_border <- length(which(set_ineq_border))

    # initialize subgradient vector and find lambda_path
    resid <- y - X %*% beta_path[, 1]
    subgrad <-
      t(X) %*% resid - t(Aeq) %*% dualpath_eq[, 1] - t(A) %*% dualpath_ineq[, 1]
    lambda_path[, 1] <- max(abs(subgrad))
    idx <- which.max(abs(subgrad))
    subgrad[sets_active] <- sign(beta_path[sets_active, 1])
    subgrad[!sets_active] <- subgrad[!sets_active] / lambda_path[, 1]
    sets_active[idx] <- TRUE
    n_active <- length(which(sets_active))
    n_inactive <- length(which(!sets_active))

    # calculate value for the objective function
    objval_path[, 1] <-
      0.5 * (sum((y - X %*% beta_path[, 1])^2)) + lambda_path[, 1] * sum(abs(beta_path[, 1]))

    # calculate degrees of freedom
    Aeq_rank <- qr(Aeq)$rank
    df_path[, 1] <- n_active - Aeq_rank - n_ineq_border

    # set initial violations counter to 0
    violations_path[1] <- 0

    # direction of the algorithm (sign)
    dir_sign <- -1


    #### MAIN LOOP FOR PATH FOLLOWING

    for (k in 2:max_iters) {
      printer(k)
      printer(lambda_path[, k - 1])

      # threshold near-zero lambdas to zero and stop algorithm
      if (lambda_path[, k - 1] <= (0 + stop_lambda_tol)) {
        lambda_path[, k - 1] <- 0
        printer(paste("BREAK. Previous Lambda < ", stop_lambda_tol, ".",
          sep =
            ""
        ))
        break
      }

      M <- cbind(
        H[sets_active, sets_active], t(matrix(Aeq[, sets_active], ncol = n_active)),
        t(matrix(A[set_ineq_border, sets_active], ncol = n_active))
      )

      M <-
        rbind(M, cbind(
          rbind(matrix(Aeq[, sets_active], ncol = n_active), matrix(A[set_ineq_border, sets_active],
            ncol =
              n_active
          )),
          matrix(0, m1 + n_ineq_border, m1 + n_ineq_border)
        ))

      ## calculate derivative
      # try using a regular inverse first, otherwise the Moore-Penrose Inverse
      inv_calc_error <- function(e) {
        dir <-
          -(ginv(M) %*% rbind(
            matrix(subgrad[sets_active], n_active, 1),
            matrix(0, m1 + n_ineq_border, 1)
          ))
        printer("Moore-Penrose-Inverse used.")
        dir
      }
      dir <-
        tryCatch(
          dir_sign * (solve(M, rbind(
            matrix(subgrad[sets_active], n_active, 1), matrix(0, m1 + n_ineq_border, 1)
          ))),
          error = inv_calc_error
        )

      if (n_inactive != 0) {
        dir_subgrad <-
          -cbind(
            matrix(H[!sets_active, sets_active], ncol = n_active),
            t(matrix(Aeq[, !sets_active], ncol = n_inactive)),
            t(matrix(A[set_ineq_border, !sets_active], ncol = n_inactive))
          ) %*% dir
      } else {
        dir_subgrad <- matrix(0, 0, m1 + n_ineq_border)
      }

      ### check additional events related to potential subgradient violations

      ## inactive coefficients moving too slowly

      # negative subgradient
      inact_slow_neg_idx <- which((1 * dir_sign - ceiling_tol) <= subgrad[!sets_active] &
        subgrad[!sets_active] <= (1 * dir_sign + ceiling_tol) &
        1 * dir_sign < dir_subgrad)
      # positive subgradient
      inact_slow_pos_idx <-
        which((-1 * dir_sign - ceiling_tol) <= subgrad[!sets_active] &
          subgrad[!sets_active] <= (-1 * dir_sign + ceiling_tol) &
          dir_subgrad < -1 * dir_sign)

      ## "Active" coeficients estimated as 0 with potential sign mismatch
      # Positive subgradient but negative derivative
      sign_mismatch_pos_idx <- which((0 - ceiling_tol) <= subgrad[sets_active] &
        subgrad[sets_active] <= (1 + ceiling_tol) &
        dir_sign * dir[1:n_active] <= (0 - ceiling_tol) &
        beta_path[sets_active, k - 1] == 0)
      # Negative subgradient but positive derivative
      sign_mismatch_neg_idx <-
        which((-1 - ceiling_tol) <= subgrad[sets_active] &
          subgrad[sets_active] <= (0 + ceiling_tol) &
          (0 + ceiling_tol) <= dir_sign * dir[1:n_active] &
          beta_path[sets_active, k - 1] == 0)

      # reset violation counter (to avoid infinite loops)
      violation_counter <- 0

      ### Outer while loop for checking all conditions together

      while (length(inact_slow_neg_idx) != 0 ||
        length(inact_slow_pos_idx) != 0 ||
        length(sign_mismatch_pos_idx) != 0 || length(sign_mismatch_neg_idx) != 0) {
        printer("VIOLATIONS DUE TO SLOW ALGO MOVEMENT OR POS/NEG MISMATCH")

        ## Monitor & fix condition 1 violations
        while (length(inact_slow_neg_idx) != 0) {
          printer("Violation inact_slow_neg_idx")

          # Identify & move problem coefficient

          inactive_coefs <-
            which(!sets_active) # indices corresponding to inactive coefficients
          viol_coeff <-
            inactive_coefs[inact_slow_neg_idx] # identify problem coefficient
          sets_active[viol_coeff] <-
            TRUE # put problem coefficient back into active set
          n_active <-
            length(which(sets_active)) # determine new number of active/inactive coefficients
          n_inactive <- length(which(!sets_active))
          n_ineq_border <-
            length(which(set_ineq_border)) # determine number of active/binding inequality constraints

          # Recalculate derivative for coefficients & multipliers

          M <-
            cbind(
              H[sets_active, sets_active],
              t(matrix(Aeq[, sets_active], ncol = n_active)),
              t(matrix(A[set_ineq_border, sets_active], ncol = n_active))
            )
          M <- rbind(M, cbind(
            rbind(
              matrix(Aeq[, sets_active], ncol = n_active),
              matrix(A[set_ineq_border, sets_active], ncol = n_active)
            ), matrix(0, m1 + n_ineq_border, m1 + n_ineq_border)
          ))
          inv_calc_error <- function(e) {
            dir <-
              -(ginv(M) %*% rbind(
                matrix(subgrad[sets_active], n_active, 1),
                matrix(0, m1 + n_ineq_border, 1)
              ))
            dir
          }
          dir <-
            tryCatch(
              dir_sign * (solve(M, rbind(
                matrix(subgrad[sets_active], n_active, 1),
                matrix(0, m1 + n_ineq_border, 1)
              ))),
              error = inv_calc_error
            )

          # calculate derivative for lambda*subgradient
          if (n_inactive != 0) {
            dir_subgrad <-
              -cbind(matrix(H[!sets_active, sets_active], ncol = n_active), t(matrix(Aeq[, !sets_active],
                ncol =
                  n_inactive
              )), t(matrix(A[set_ineq_border, !sets_active], ncol = n_inactive))) %*%
              dir
          } else {
            dir_subgrad <- matrix(0, 0, m1 + n_ineq_border)
          }

          # check for violations again

          # Negative subgradient
          inact_slow_neg_idx <-
            which((1 * dir_sign - ceiling_tol) <= subgrad[!sets_active] &
              subgrad[!sets_active] <= (1 * dir_sign + ceiling_tol) &
              1 * dir_sign < dir_subgrad)
          # Positive subgradient
          inact_slow_pos_idx <-
            which((-1 * dir_sign - ceiling_tol) <= subgrad[!sets_active] &
              subgrad[!sets_active] <= (-1 * dir_sign + ceiling_tol) &
              dir_subgrad < -1 * dir_sign)
          # Positive subgrad but negative derivative
          sign_mismatch_pos_idx <- which((0 - ceiling_tol) <= subgrad[sets_active] &
            subgrad[sets_active] <= (1 + ceiling_tol) &
            dir_sign * dir[1:n_active] <= (0 - ceiling_tol) &
            beta_path[sets_active, k - 1] == 0)
          # Negative subgradient but positive derivative
          sign_mismatch_neg_idx <-
            which((-1 - ceiling_tol) <= subgrad[sets_active] &
              subgrad[sets_active] <= (0 + ceiling_tol) &
              (0 + ceiling_tol) <= dir_sign * dir[1:n_active] &
              beta_path[sets_active, k - 1] == 0)

          # update violation counter
          violation_counter <- violation_counter + 1
          if (violation_counter >= max_iters) {
            printer("Too many violations.")
            break
          }
        }

        ## Monitor & fix subgradient condition 2 violations

        while (length(inact_slow_pos_idx) != 0) {
          printer("violation inact_slow_pos_idx")

          # Identify & move problem coefficient

          inactive_coefs <-
            which(!sets_active) # indices corresponding to inactive coefficients
          viol_coeff <-
            inactive_coefs[inact_slow_pos_idx] # identify problem coefficient
          sets_active[viol_coeff] <-
            TRUE # put problem coefficient back into active set;
          n_active <-
            length(which(sets_active)) # determine new number of active/inactive coefficients
          n_inactive <- length(which(!sets_active))
          n_ineq_border <-
            length(which(set_ineq_border)) # determine number of active/binding inequality constraints

          # Recalculate derivative for coefficients & multipliers

          M <-
            cbind(
              H[sets_active, sets_active],
              t(matrix(Aeq[, sets_active], ncol = n_active)),
              t(matrix(A[set_ineq_border, sets_active], ncol = n_active))
            )
          M <-
            rbind(M, cbind(
              rbind(matrix(Aeq[, sets_active], ncol = n_active), matrix(A[set_ineq_border, sets_active],
                ncol =
                  n_active
              )),
              matrix(0, m1 + n_ineq_border, m1 + n_ineq_border)
            ))

          inv_calc_error <- function(e) {
            dir <-
              -(ginv(M) %*% rbind(
                matrix(subgrad[sets_active], n_active, 1),
                matrix(0, m1 + n_ineq_border, 1)
              ))
            printer("Moore-Penrose-Inverse used.")
            dir
          }
          dir <-
            tryCatch(
              dir_sign * (solve(M, rbind(
                matrix(subgrad[sets_active], n_active, 1),
                matrix(0, m1 + n_ineq_border, 1)
              ))),
              error = inv_calc_error
            )

          if (n_inactive != 0) {
            dir_subgrad <-
              -cbind(matrix(H[!sets_active, sets_active], ncol = n_active), t(matrix(Aeq[, !sets_active],
                ncol =
                  n_inactive
              )), t(matrix(A[set_ineq_border, !sets_active], ncol = n_inactive))) %*%
              dir
          } else {
            dir_subgrad <- matrix(0, 0, m1 + n_ineq_border)
          }

          # check for violations again

          # Positive subgradient
          inact_slow_pos_idx <-
            which((-1 * dir_sign - ceiling_tol) <= subgrad[!sets_active] &
              subgrad[!sets_active] <= (-1 * dir_sign + ceiling_tol) &
              dir_subgrad < -1 * dir_sign)
          # Positive subgrad but negative derivative
          sign_mismatch_pos_idx <- which((0 - ceiling_tol) <= subgrad[sets_active] &
            subgrad[sets_active] <= (1 + ceiling_tol) &
            dir_sign * dir[1:n_active] <= (0 - ceiling_tol) &
            beta_path[sets_active, k - 1] == 0)
          # Negative subgradient but positive derivative
          sign_mismatch_neg_idx <-
            which((-1 - ceiling_tol) <= subgrad[sets_active] &
              subgrad[sets_active] <= (0 + ceiling_tol) &
              (0 + ceiling_tol) <= dir_sign * dir[1:n_active] &
              beta_path[sets_active, k - 1] == 0)

          # update violation counter
          violation_counter <- violation_counter + 1

          if (violation_counter >= max_iters) {
            printer("Too many violations.")
            break
          }
        }

        ## Monitor & fix subgradient condition 3 violations

        while (length(sign_mismatch_pos_idx) != 0) {
          printer("violation sign_mismatch_pos_idx")

          # Identify & move problem coefficient

          active_coefs <-
            which(sets_active) # indices corresponding to inactive coefficients
          viol_coeff <-
            active_coefs[sign_mismatch_pos_idx] # identify problem coefficient
          sets_active[viol_coeff] <-
            FALSE # put problem coefficient back into active set;
          n_active <-
            length(which(sets_active)) # determine new number of active/inactive coefficients
          n_inactive <- length(which(!sets_active))
          n_ineq_border <-
            length(which(set_ineq_border)) # determine number of active/binding inequality constraints


          # Recalculate derivative for coefficients & multipliers

          M <-
            cbind(H[sets_active, sets_active], t(matrix(Aeq[, sets_active], ncol = n_active)), t(matrix(A[set_ineq_border, sets_active],
              ncol =
                n_active
            )))
          M <-
            rbind(M, cbind(
              rbind(matrix(Aeq[, sets_active], ncol = n_active), matrix(A[set_ineq_border, sets_active],
                ncol =
                  n_active
              )),
              matrix(0, m1 + n_ineq_border, m1 + n_ineq_border)
            ))

          inv_calc_error <- function(e) {
            dir <-
              -(ginv(M) %*% rbind(
                matrix(subgrad[sets_active], n_active, 1),
                matrix(0, m1 + n_ineq_border, 1)
              ))
            printer("Moore-Penrose-Inverse used.")
            dir
          }
          dir <-
            tryCatch(
              dir_sign * (solve(M, rbind(
                matrix(subgrad[sets_active], n_active, 1),
                matrix(0, m1 + n_ineq_border, 1)
              ))),
              error = inv_calc_error
            )

          if (n_inactive != 0) {
            dir_subgrad <-
              -cbind(matrix(H[!sets_active, sets_active], ncol = n_active), t(matrix(Aeq[, !sets_active],
                ncol =
                  n_inactive
              )), t(matrix(A[set_ineq_border, !sets_active], ncol = n_inactive))) %*%
              dir
          } else {
            dir_subgrad <- matrix(0, 0, m1 + n_ineq_border)
          }

          # check for violations again

          # Positive subgrad but negative derivative
          sign_mismatch_pos_idx <- which((0 - ceiling_tol) <= subgrad[sets_active] &
            subgrad[sets_active] <= (1 + ceiling_tol) &
            dir_sign * dir[1:n_active] <= (0 - ceiling_tol) &
            beta_path[sets_active, k - 1] == 0)

          # Negative subgradient but positive derivative
          sign_mismatch_neg_idx <-
            which((-1 - ceiling_tol) <= subgrad[sets_active] &
              subgrad[sets_active] <= (0 + ceiling_tol) &
              (0 + ceiling_tol) <= dir_sign * dir[1:n_active] &
              beta_path[sets_active, k - 1] == 0)

          # update violation counter
          violation_counter <- violation_counter + 1

          if (violation_counter >= max_iters) {
            printer("Too many violations.")
            break
          }
        }

        ## Monitor & fix subgradient condition 4 violations

        while (length(sign_mismatch_neg_idx) != 0) {
          printer("violation sign_mismatch_neg_idx")

          # identify & move problem coefficient

          active_coefs <-
            which(sets_active) # indices corresponding to inactive coefficients
          viol_coeff <-
            active_coefs[sign_mismatch_neg_idx] # identify problem coefficient
          sets_active[viol_coeff] <-
            FALSE # put problem coefficient back into active set
          n_active <-
            length(which(sets_active)) # determine new number of active/inactive coefficients
          n_inactive <- length(which(!sets_active))
          n_ineq_border <-
            length(which(set_ineq_border)) # determine number of active/binding inequality constraints


          # recalculate derivative for coefficients & multipliers

          M <-
            cbind(H[sets_active, sets_active], t(matrix(Aeq[, sets_active], ncol = n_active)), t(matrix(A[set_ineq_border, sets_active],
              ncol =
                n_active
            )))
          M <-
            rbind(M, cbind(
              rbind(matrix(Aeq[, sets_active], ncol = n_active), matrix(A[set_ineq_border, sets_active],
                ncol =
                  n_active
              )),
              matrix(0, m1 + n_ineq_border, m1 + n_ineq_border)
            ))

          inv_calc_error <- function(e) {
            dir <-
              -(ginv(M) %*% rbind(
                matrix(subgrad[sets_active], n_active, 1),
                matrix(0, m1 + n_ineq_border, 1)
              ))
            printer("Moore-Penrose-Inverse used.")
            dir
          }
          dir <-
            tryCatch(
              dir_sign * (solve(M, rbind(
                matrix(subgrad[sets_active], n_active, 1),
                matrix(0, m1 + n_ineq_border, 1)
              ))),
              error = inv_calc_error
            )

          if (n_inactive != 0) {
            dir_subgrad <-
              -cbind(matrix(H[!sets_active, sets_active], ncol = n_active), t(matrix(Aeq[, !sets_active],
                ncol =
                  n_inactive
              )), t(matrix(A[set_ineq_border, !sets_active], ncol = n_inactive))) %*%
              dir
          } else {
            dir_subgrad <- matrix(0, 0, m1 + n_ineq_border)
          }

          # check for violations again

          # Negative subgradient but positive derivative
          sign_mismatch_neg_idx <-
            which((-1 - ceiling_tol) <= subgrad[sets_active] &
              subgrad[sets_active] <= (0 + ceiling_tol) &
              (0 + ceiling_tol) <= dir_sign * dir[1:n_active] &
              beta_path[sets_active, k - 1] == 0)

          # update violation counter
          violation_counter <- violation_counter + 1

          if (violation_counter >= max_iters) {
            printer("Too many violations.")
            break
          }
        }

        # update violation trackers to see if any issues persist
        # Negative subgradient
        inact_slow_neg_idx <-
          which((1 * dir_sign - ceiling_tol) <= subgrad[!sets_active] &
            subgrad[!sets_active] <= (1 * dir_sign + ceiling_tol) &
            1 * dir_sign < dir_subgrad)
        # Positive subgradient
        inact_slow_pos_idx <-
          which((-1 * dir_sign - ceiling_tol) <= subgrad[!sets_active] &
            subgrad[!sets_active] <= (-1 * dir_sign + ceiling_tol) &
            dir_subgrad < -1 * dir_sign)
        # Positive subgrad but negative derivative
        sign_mismatch_pos_idx <- which((0 - ceiling_tol) <= subgrad[sets_active] &
          subgrad[sets_active] <= (1 + ceiling_tol) &
          dir_sign * dir[1:n_active] <= (0 - ceiling_tol) &
          beta_path[sets_active, k - 1] == 0)
        # Negative subgradient but positive derivative
        sign_mismatch_neg_idx <-
          which((-1 - ceiling_tol) <= subgrad[sets_active] &
            subgrad[sets_active] <= (0 + ceiling_tol) &
            (0 + ceiling_tol) <= dir_sign * dir[1:n_active] &
            beta_path[sets_active, k - 1] == 0)

        if (violation_counter >= max_iters) {
          printer("Too many violations.")
          break
        }
      } ### end of violation check while-loop

      # store number of violations
      violations_path[k] <- violation_counter

      # calculate derivative for residual inequality
      dirresid_ineq <-
        matrix(A[!set_ineq_border, sets_active],
          nrow = length(which(!set_ineq_border)), ncol =
            n_active
        ) %*% dir[1:n_active]

      ### Determine lambda for next event (via delta lambda)

      next_lambda_beta <- matrix(Inf, p, 1)

      ## Events based on changes in coefficient status

      # Active coefficient going inactive
      next_lambda_beta[sets_active, ] <- -dir_sign * beta_path[sets_active, k - 1] / dir[1:n_active]

      # inactive coefficient becoming positive
      t1 <-
        dir_sign * lambda_path[, k - 1] * (1 - subgrad[!sets_active]) / (dir_subgrad - 1)
      # threshold values hitting ceiling
      t1[t1 <= (0 + ceiling_tol)] <- Inf

      # inactive coefficient becoming negative
      t2 <- -dir_sign * lambda_path[, k - 1] * (1 + subgrad[!sets_active]) / (dir_subgrad + 1)
      # threshold values hitting ceiling
      t2[t2 <= (0 + ceiling_tol)] <- Inf

      # choose smaller delta lambda out of t1 and t2
      next_lambda_beta[!sets_active, ] <- pmin(t1, t2)

      # ignore delta lambdas numerically equal to zero
      next_lambda_beta[next_lambda_beta <= ceiling_tol | !penidx, ] <- Inf

      ## Events based inequality constraints

      # clear previous values
      next_lambda_ineq <- matrix(Inf, m2, 1)

      # inactive inequality constraint becoming active
      next_lambda_ineq[!set_ineq_border, ] <-
        as.matrix(-dir_sign * resid_ineq[!set_ineq_border], length(which(!set_ineq_border)), 1) /
          (as.matrix(dirresid_ineq, length(which(!set_ineq_border)), 1))

      # active inequality constraint becoming inactive
      next_lambda_ineq[set_ineq_border, ] <-
        as.matrix(-dir_sign * dualpath_ineq[set_ineq_border, k - 1]) / as.matrix(dir[-(1:(n_active +
          m1))], n_ineq_border, 1)

      # ignore delta lambdas equal to zero
      next_lambda_ineq[next_lambda_ineq <= ceiling_tol, ] <- Inf

      # find smallest lambda
      chg_lambda <-
        min(rbind(next_lambda_beta, next_lambda_ineq), na.rm = TRUE)

      # find all indices corresponding to this chg_lambda
      idx <-
        which((rbind(next_lambda_beta, next_lambda_ineq) - chg_lambda) <= ceiling_tol)

      # terminate path following if no new event found
      if (is.infinite(chg_lambda)) {
        chg_lambda <- lambda_path[, k - 1]
      }

      # update values at new lambda and move to next lambda, make sure is not negative

      if ((lambda_path[, k - 1] + dir_sign * chg_lambda) < 0) {
        chg_lambda <- lambda_path[, k - 1]
      }

      printer(chg_lambda)
      # calculate new value of lambda
      lambda_path[, k] <- lambda_path[, k - 1] + dir_sign * chg_lambda

      ## update parameter and subgradient values

      # new coefficient estimates
      beta_path[sets_active, k] <-
        beta_path[sets_active, k - 1] + dir_sign * chg_lambda * dir[1:n_active]

      # force near-zero coefficients to be zero (helps with numerical issues)
      beta_path[abs(beta_path[, k]) < zeros_tol, k] <- 0

      # new subgradient estimates
      subgrad[!sets_active] <-
        as.matrix(lambda_path[, k - 1] * subgrad[!sets_active, ] + dir_sign * chg_lambda *
          dir_subgrad) / lambda_path[, k]

      ## update dual variables

      # update duals1 (lagrange multipliers for equality constraints)
      dualpath_eq[, k] <-
        dualpath_eq[, k - 1] + dir_sign * chg_lambda * as.matrix(dir[(n_active + 1):(n_active +
          m1)], m1, 1)

      # update duals2 (lagrange multipliers for inequality constraints)
      dualpath_ineq[set_ineq_border, k] <-
        dualpath_ineq[set_ineq_border, k - 1] + dir_sign * chg_lambda * as.matrix(dir[(n_active +
          m1 + 1):nrow(dir)], n_ineq_border, 1)

      # update residual inequality
      resid_ineq <- A %*% beta_path[, k] - b

      # update sets

      for (j in seq_along(idx)) {
        curidx <- idx[j]
        if (curidx <= p && sets_active[curidx]) {
          # an active coefficient hits 0, or
          sets_active[curidx] <- FALSE
        } else if (curidx <= p && !sets_active[curidx]) {
          # a zero coefficient becomes nonzero
          sets_active[curidx] <- TRUE
        } else if (curidx > p) {
          # an inequality on boundary becomes strict, or a strict inequality hits boundary
          set_ineq_border[curidx - p] <- !set_ineq_border[curidx - p]
        }
      }

      # determine new number of active coefficients
      n_active <- length(which(sets_active))
      n_inactive <- length(which(!sets_active))

      # determine number of active/binding inequality constraints
      n_ineq_border <- length(which(set_ineq_border))

      # calculate value of the objective function
      objval_path[k] <- 0.5 * (sum((y - X %*% beta_path[, k])^2)) + lambda_path[k] *
        sum(abs(beta_path[, k]))

      # calculate degrees of freedom (dfs)
      df_path[k] <- n_active - Aeq_rank - n_ineq_border

      # break algorithm when dfs are exhausted
      if (df_path[k - 1] >= n_orig) {
        printer("BREAK. No more degrees of freedom.")
        break
      }

      if (length(which(sets_active)) == p) {
        printer("BREAK. All coefficients active. No further Sparsity.")
        break
      }
    }

    beta_path <- beta_path[, 1:(k - 1)]
    lambda_path <- lambda_path[1:(k - 1)]
    objval_path <- objval_path[1:(k - 1)]
    df_path <- df_path[1:(k - 1)]
    df_path[df_path < 0] <- 0

    return(
      list(
        "beta_path" = beta_path,
        "lambda_path" = lambda_path,
        "objval_path" = objval_path,
        "df_path" = df_path
      )
    )
  }
