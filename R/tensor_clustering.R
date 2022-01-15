#' Weighted higher-order initialization
#'
#' Weighted higher-order initialization for multiway spherical clustering under degree-corrected tensor block model.
#' This function takes the tensor/matrix observation, the number of clusters, and a logic variable indicating the symmetry
#' as input. Output is the estimated clustering assignment.
#'
#'
#' @param Y     array/matrix, order-3 tensor/matrix observation
#' @param r     vector, the number of clusters on each mode; see "details"
#' @param asymm logic variable, if "TRUE", assume the clustering assignment differs in different modes; if "FALSE", assume all the modes share the same clustering assignment
#' @return a list containing the following:
#'
#' \code{z0} { a list of vectors recording the estimated clustering assignment }
#'
#' \code{s0} { a list of vectors recording the index of degenerate entities with random clustering assignment}
#'
#' @details   \code{r} should be a length 2 vector for matrix and length 3 vector for tensor observation;
#'
#' all the elements in \code{r} should be integer larger than 1;
#'
#' matrix case and symmetric case only allow \code{r} with the same number of clusters on each mode;
#'
#' observations with non-identical dimension on each mode are only applicable with \code{asymm = T}.
#'
#'
#'
#' @export
#' @examples
#'
#' test_data = sim_dTBM(seed = 1, imat = FALSE, asymm = FALSE, p = c(50,50,50), r = c(3,3,3),
#'                     core_control = "control", s_min = 0.05, s_max = 1,
#'                     dist = "normal", sigma = 0.5,
#'                     theta_dist = "pareto", alpha = 4, beta = 3/4)
#'
#' initialization <- wkmeans(test_data$Y, r = c(3,3,3), asymm = FALSE)
#'


wkmeans <- function(Y, r, asymm) {
  imat <- F

  if (length(dim(Y)) == 2) {
    cat("matrix case \n")
    dim(Y) <- c(dim(Y), 1)
    imat <- T
  }

  if (imat == T & length(r) != 2) {
    warning("need to input a length 2 vector for the number of clusters")
    return()
  }

  if (imat == F & length(r) != 3) {
    warning("need to input a length 3 vector for the number of clusters")
    return()
  }

  r1 <- r[1]
  r2 <- r[2]

  if (imat == T) {
    r3 <- 1

    if (r1 != r2) {
      warning("matrix case requires the same number of clusters on two modes")
      return()
    }

    if (r1 <= 1 | r2 <= 1) {
      warning("all the numbers of clusters should be larger than 1")
      return()
    }

    if (sum(dim(Y)[1:2]) / dim(Y)[1] != 2 & asymm == F) {
      warning("use asymmetric algorithm for observation with non-identical dimension on each mode")
      return()
    }
  } else if (imat == F) {
    r3 <- r[3]

    if (asymm == F) {
      if (r1 != r2 | r2 != r3 | r1 != r3) {
        warning("symmetric case requires the same number of clusters on every mode")
        # r3 = r2 = r1 = min(c(r1,r2,r3))
        return()
      }

      if (sum(dim(Y)) / dim(Y)[1] != 3) {
        warning("use asymmetric algorithm for observation with non-identical dimension on each mode")
        return()
      }
    }

    if (r1 <= 1 | r2 <= 1 | r3 <= 1) {
      warning("all the numbers of clusters should be larger than 1")
      return()
    }
  }

  ### two-step SVD
  Y <- as.tensor(Y)
  # first SVD
  u1 <- svd(unfold(Y, 1, c(3, 2))@data)$u[, 1:r1]
  u2 <- svd(unfold(Y, 2, c(1, 3))@data)$u[, 1:r2]
  if (imat == F) {
    u3 <- svd(unfold(Y, 3, c(1, 2))@data)$u[, 1:r3]
  } else if (imat == T) {
    u3 <- as.matrix(1)
  }

  # second SVD

  ##### for matrix case, second SVD asks r1 = r2
  hu1 <- svd(unfold(ttl(Y, list(t(u2), t(u3)), ms = c(2, 3)), 1, c(3, 2))@data)$u[, 1:r1]
  hu2 <- svd(unfold(ttl(Y, list(t(u1), t(u3)), ms = c(1, 3)), 2, c(1, 3))@data)$u[, 1:r2]
  if (imat == F) {
    hu3 <- svd(unfold(ttl(Y, list(t(u1), t(u2)), ms = c(1, 2)), 3, c(1, 2))@data)$u[, 1:r3]
  } else if (imat == T) {
    hu3 <- as.matrix(1)
  }

  ### estimated X
  X1 <- hu1 %*% t(hu1) %*% unfold(ttl(Y, list(hu2 %*% t(hu2), hu3 %*% t(hu3)), ms = c(2, 3)), 1, c(3, 2))@data
  X2 <- hu2 %*% t(hu2) %*% unfold(ttl(Y, list(hu1 %*% t(hu1), hu3 %*% t(hu3)), ms = c(1, 3)), 2, c(1, 3))@data
  X3 <- hu3 %*% t(hu3) %*% unfold(ttl(Y, list(hu1 %*% t(hu1), hu2 %*% t(hu2)), ms = c(1, 2)), 3, c(1, 2))@data

  if (asymm == T) {
    res1 <- single_wkmeans(X1, r1)
    res2 <- single_wkmeans(X2, r2)
    if (imat == F) {
      res3 <- single_wkmeans(X3, r3)

      z0 <- list(as.numeric(res1$z), as.numeric(res2$z), as.numeric(res3$z))
      s0 <- list(s01 = res1$s0, s02 = res2$s0, s03 = res3$s0)
    } else if (imat == T) {
      z0 <- list(as.numeric(res1$z), as.numeric(res2$z))
      s0 <- list(s01 = res1$s0, s02 = res2$s0)
    }

    return(list(z0 = z0, s0 = s0))
  } else if (asymm == F) {
    z <- rep(0, dim(Y)[1])
    sc <- 1:length(z)
    ### normalization
    l2 <- apply(X1, 1, function(x) sqrt(sum(x^2))) # l2 norm

    if (length(which(l2 == 0) > 0)) {
      sc <- which(l2 != 0)
      X1 <- X1[sc, ]
      l2 <- l2[sc]
    }

    Xs <- diag(l2^(-1)) %*% X1

    ## weighted kmeans
    diss <- dist(Xs, method = "euclidean", p = 2, upper = T, diag = T)

    z[sc] <- wcKMedoids(diss^2, k = r1, weights = l2^2, method = "PAMonce", cluster.only = T)
    z[-sc] <- sample(unique(z[sc]), length(z[-sc]), replace = T)


    s0 <- setdiff(1:length(z), sc)
    z <- as.factor(z)
    levels(z) <- 1:r1

    z0 <- as.numeric(z)

    if (imat == F) {
      z0 <- list(z0, z0, z0)
      s0 <- list(s0, s0, s0)
    } else if (imat == T) {
      z0 <- list(z0, z0)
      s0 <- list(s0, s0)
    }

    return(list(z0 = z0, s0 = s0))
  }
}





#' Angle-based iteration
#'
#' Angle-based iteration for multiway spherical clustering under degree-corrected tensor block model.
#' This function takes the tensor/matrix observation, initial clustering assignment, and a logic variable indicating the symmetry
#' as input. Output is the refined clustering assignment.
#'
#'
#' @param Y         array/matrix, order-3 tensor/matrix observation
#' @param z0        a list of vectors, initial clustering assignment; see "details"
#' @param max_iter  integer, max number of iterations if update does not converge
#' @param alpha1    number, substitution of degenerate core tensor; see "details"
#' @param asymm     logic variable, if "TRUE", assume the clustering assignment differs in different modes; if "FALSE", assume all the modes share the same clustering assignment
#'
#' @return a list containing the following:
#'
#' \code{z} {a list of vectors recording the estimated clustering assignment}
#'
#' \code{s_deg} {logic variable, if "TRUE", degenerate estimated core tensor/matrix occurs during the iteration; if "FALSE", otherwise}
#'
#' @details \code{z0} should be a length 2 list for matrix and length 3 list for tensor observation;
#' observations with non-identical dimension on each mode are only applicable with \code{asymm = T};
#'
#' When the estimated core tensor has a degenerate slice, i.e., a slice with all zero elements, randomly pick an entry in the degenerate slice with value \code{alpha1}.
#'
#' @export
#' @examples
#' test_data = sim_dTBM(seed = 1, imat = FALSE, asymm = FALSE, p = c(50,50,50), r = c(3,3,3),
#'                     core_control = "control", s_min = 0.05, s_max = 1,
#'                     dist = "normal", sigma = 0.5,
#'                     theta_dist = "pareto", alpha = 4, beta = 3/4)
#'
#' initialization <- wkmeans(test_data$Y, r = c(3,3,3), asymm = FALSE)
#'
#' iteration <- angle_iteration(test_data$Y, initialization$z0, max_iter = 20, asymm = FALSE)


angle_iteration = function(Y, z0, max_iter, alpha1 = 0.01, asymm){
  #Y = test_data$Y; z0 = test_ini$z0
  imat <- F
  s_deg <- F

  # z0 should be a list
  if (length(dim(Y)) == 2) {
    cat("matrix case \n")
    dim(Y) <- c(dim(Y), 1)
    imat <- T

    if(sum(dim(Y)[1:2])/dim(Y)[1] != 2 & asymm == F){
      warning("use asymmetric algorithm for observation with non-identical dimension on each mode")
      return()
    }
  }

  if(sum(dim(Y))/dim(Y)[1] != 3 & asymm == F & imat == F){
    warning("use asymmetric algorithm for observation with non-identical dimension on each mode")
    return()
  }

  z <- lapply(z0, renumber)
  if (imat == T) {
    z[[3]] <- 1
  }

  for (iter in 1:max_iter) {
    cat("iter = ", iter, "\n")

    # estimate S
    est_S <- updateS(Y, z, imat)

    # update z1
    # reducde Y
    Y1 <- Cal_Y1(Y, z, imat)
    re1 <- single_Aiteration(unfold(as.tensor(est_S), 1, c(3, 2))@data, unfold(as.tensor(Y1), 1, c(3, 2))@data, alpha1)
    z1_new = re1$z
    z1_new <- renumber(z1_new)

    if(asymm == F){

      if(imat == T){
        z_new = list(z1_new, z1_new,as.vector(1))
      }else if(imat == F){
        z_new = list(z1_new, z1_new, z1_new)
      }

      if(re1$s_deg == T){
        s_deg = T
      }

      if (identical(z_new, z)) {
        break
      }

      z = z_new

    }else if(asymm == T){

      # update z2
      # reducde Y
      Y2 <- Cal_Y2(Y, z,imat)
      re2 <- single_Aiteration(unfold(as.tensor(est_S), 2, c(1, 3))@data, unfold(as.tensor(Y2), 2, c(1, 3))@data, alpha1)
      z2_new = re2$z
      z2_new <- renumber(z2_new)

      #update z3
      # reducde Y
      if (imat == T) {
        z3_new <- 1

        if(re1$s_deg == T | re2$s_deg == T){
          s_deg = T
        }

      } else if (imat == F) {
        Y3 <- Cal_Y3(Y, z)
        re3 <- single_Aiteration(unfold(as.tensor(est_S), 3, c(1, 2))@data, unfold(as.tensor(Y3), 3, c(1, 2))@data, alpha1)
        z3_new = re3$z
        z3_new <- renumber(z3_new)

        if(re1$s_deg == T | re2$s_deg == T|re3$s_deg == T){
          s_deg = T
        }

      }

      z_new_list = list(z1_new,z2_new, z3_new)

      if (identical(z_new_list, z)) {
        break
      }

      z <- z_new_list

    }

  }

  if(imat == T){
    z[[3]] = NULL
  }

  return(list(z = z, s_deg = s_deg))
}


#' Simulation of degree-corrected tensor block models
#'
#' Generate order-3 tensor/matrix observations with degree heterogeneity under degree-corrected tensor block models.
#'
#' @param seed          number, random seed for generating data
#' @param imat          logic variable, if "TRUE", generate matrix data; if "FALSE", generate order-3 tensor data
#' @param asymm         logic variable, if "TRUE", clustering assignment differs in different modes; if "FALSE", all the modes share the same clustering assignment
#' @param p             vector, dimension of the tensor/matrix observation
#' @param r             vector, number of clusters on each mode
#' @param core_control  character, the way to control the generation of core tensor/matrix; see "details"
#' @param delta         number, Frobenius norm of the slices in core tensor if \code{core_control = "control"}
#' @param s_min         number, value of off-diagonal elements in original core tensor/matrix if \code{core_control = "control"}
#' @param s_max         number, value of diagonal elements in original core tensor/matrix if \code{core_control = "control"}
#' @param dist          character, distribution of tensor/matrix observation; see "details"
#' @param sigma         number, standard deviation of Gaussian noise if \code{dist = "normal"}
#' @param theta_dist    character, distribution of degree heterogeneity; see "details"
#' @param alpha         number, shape parameter in pareto distribution if \code{theta_dist = "pareto"}
#' @param beta          number, scale parameter in pareto distribution if \code{theta_dist = "pareto"}
#'
#' @return  a list containing the following:
#'
#' \code{Y} {array ( if \code{imat = F} )/matrix ( if \code{imat = T} ), simulated tensor/matrix observations with dimension \code{p}  }
#'
#' \code{X} {array ( if \code{imat = F} )/matrix ( if \code{imat = T} ), mean tensor/matrix of the observation, i.e., the expectation of \code{Y}}
#'
#' \code{S} {array ( if \code{imat = F} )/matrix ( if \code{imat = T} ), core tensor/matrix recording the block effects with dimension \code{r}}
#'
#' \code{theta} {a list of vectors, degree heterogeneity on each mode}
#'
#' \code{z} {a list of vectors, clustering assignment on each mode}
#'
#'
#' @details  The general tensor observation is generated as
#'
#' \code{Y = S x1 Theta1 M1 x2 Theta2 M2 x3 Theta3 M3 + E,}
#'
#' where \code{S} is the core tensor, \code{Thetak} is a diagonal matrix with elements in the \code{k}-th vector of \code{theta},
#' \code{Mk} is the membership matrix based on the clustering assignment in the \code{k}-th vector of \code{z} with \code{r[k]} clusters,
#' \code{E} is the mean-zero noise tensor, and \code{xk} refers to the matrix-by-tensor product on the \code{k}-th mode, for \code{k = 1,2,3}.
#'
#' If \code{imat = T}, \code{Y,S,E} degenerate to matrix and \code{Y = Theta1 M1 S M2^T Theta2^T + E}.
#'
#' If \code{asymm = F}, \code{Thetak = Theta} and \code{Mk = M} for all \code{k = 1,2,3}.
#'
#' \code{core_control} specifies the way to generate \code{S}:
#'
#' If \code{core_control = "control"}, first generate \code{S} as a diagonal tensor/matrix with diagonal elements \code{s_max} and off-diagonal elements \code{s_min};
#' then scale the original core such that Frobenius norm of the slices equal to \code{delta}, i.e, \code{delta = sqrt(sum(S[1,,]^2))} or \code{delta = sqrt(sum(S[1,]^2))};
#' ignore the scaling if \code{delta = NULL}; option \code{"control"} is only applicable for symmetric case \code{asymm = F}.
#'
#' If \code{core_control = "random"}, generate \code{S} with random entries following uniform distribution U(0,1).
#'
#' \code{dist} specifies the distribution of \code{E}: \code{"normal"} for Gaussian and \code{"binary"} for Bernoulli distribution; \code{sigma} specifices the standard deviation if \code{dist = "normal"}.
#'
#' \code{theta_dist} firstly specifies the distribution of \code{theta}: \code{"non"} for constant 1, \code{"abs_normal"} for absoulte normal distribution, \code{"pareto"} for pareto distribution; \code{alpha, beta} specify the shape and scale parameter if \code{theta_dist = "pareto"};
#' then scale \code{theta} to have mean equal to one in each cluster.
#'
#'
#' @export
#' @examples
#'
#' test_data = sim_dTBM(seed = 1, imat = FALSE, asymm = FALSE, p = c(50,50,50), r = c(3,3,3),
#'                     core_control = "control", s_min = 0.05, s_max = 1,
#'                     dist = "normal", sigma = 0.5,
#'                     theta_dist = "pareto", alpha = 4, beta = 3/4)

sim_dTBM = function(seed = NA,imat = F,asymm = F, p, r,
                    core_control = c("random", "control"), delta = NULL, s_min = NULL, s_max =NULL,
                    dist = c("normal", "binary"), sigma = 1,
                    theta_dist = c("abs_normal", "pareto", "non"), alpha = NULL, beta = NULL){

  # p, r are vectors

  if (is.na(seed) == FALSE) set.seed(seed)

  if(imat == T){
    cat("generate matrix data \n")
    r = r[1:2]
    p = p[1:2]
  }

  if(asymm == F){

    if(sum(p)/p[1] != 3 & imat == F){
      warning("all the modes share the same dimension in symmetric case")
      return()
    }

    if(sum(p)/p[1] != 2 & imat == T){
      warning("all the modes share the same dimension in symmetric case")
      return()
    }

    if(sum(r)/r[1] != 3 & imat == F){
      warning("all the modes share the same number of clusters in symmetric case")
      return()
    }

    if(sum(r)/r[1] != 2 & imat == T){
      warning("all the modes share the same number of clusters in symmetric case")
      return()
    }
  }

  if(core_control == "control"){
    if(asymm == T){
      warning("core control is only applicable for symmetric case")
      return()
    }

    S <- sim_S(r[1], s_min, s_max, delta, imat)
  }else if(core_control == "random"){
    S = array(runif(prod(r)), dim = r)
  }

  if(asymm == F){
    # generate assignment
    z <- sample(1:r[1], p[1], replace = T)
    z <- renumber(z)

    theta = generate_theta(p[1],theta_dist,z,alpha, beta)

    # generate mean tensor
    if(imat == F){
      X <- ttl(as.tensor(S[z, z, z]), list(diag(theta), diag(theta), diag(theta)), ms = c(1, 2, 3))@data

      z = list(z,z,z)
      theta = list(theta,theta,theta)
    }else if(imat == T){
      X <- ttl(as.tensor(S[z, z]), list(diag(theta), diag(theta)), ms = c(1, 2))@data

      z = list(z,z)
      theta = list(theta,theta)

    }


  }else if(asymm == T){

    z1 = renumber( sample(1:r[1], p[1], replace = T))
    z2 = renumber( sample(1:r[2], p[2], replace = T))

    theta1 = generate_theta(p[1], theta_dist,z1,alpha, beta)
    theta2 = generate_theta(p[2], theta_dist,z2,alpha, beta)

    if(imat == F){
      z3 = renumber( sample(1:r[3], p[3], replace = T))
      theta3 = generate_theta(p[3], theta_dist,z3,alpha, beta)

      theta = list(theta1, theta2, theta3)
      z = list(z1,z2,z3)

      X <- ttl(as.tensor(S[z1, z2, z3]), list(diag(theta1), diag(theta2), diag(theta3)), ms = c(1, 2, 3))@data
    }else if(imat == T){
      X <- ttl(as.tensor(S[z1, z2]), list(diag(theta1), diag(theta2)), ms = c(1, 2))@data

      theta = list(theta1, theta2)
      z = list(z1,z2)

    }

  }

  # generate data
  if (dist == "normal") {
    Y <- X + array(rnorm(prod(dim(X)), 0, sigma), dim = dim(X))
  } else if (dist == "binary") {
    X[which(X > 1)] <- 1
    Y <- array(rbinom(prod(dim(X)), 1, as.vector(X)), dim = dim(X))
  }

  return(list(Y = Y, X = X, S = S, theta = theta, z = z))

}



#' Number of clusters selection
#'
#' Estimate the number of clusters in the degree-corrected tensor block model based on BIC criterion. The choice of BIC
#' aims to balance between the goodness-of-fit for the data and the degree of freedom in the population model.
#' This function is restricted for the Gaussian observation.
#'
#' @param Y         array/matrix, order-3 Gaussian tensor/matrix observation
#' @param r_range   matrix, candidates for the number of clusters on each row; see "details"
#' @param asymm         logic variable, if "TRUE", clustering assignment differs in different modes; if "FALSE", all the modes share the same clustering assignment
#'
#' @return a list containing the following:
#'
#' \code{r} {vector, the number of clusters among the candidates with minimal BIC value}
#'
#' \code{bic} {vector, the BIC value for each candidiate}
#'
#' @details \code{r_range} should be a two-column matrix for matrix and three-column matrix for tensor observation;
#'
#' all the elements in \code{r_range} should be integer larger than 1;
#'
#' matrix case and symmetric case only allow candidates with the same number of clusters on each mode;
#'
#' observations with non-identical dimension on each mode are only applicable with \code{asymm = T}.
#'
#' @export
#' @examples
#'
#' test_data = sim_dTBM(seed = 1, imat = FALSE, asymm = FALSE, p = c(50,50,50), r = c(3,3,3),
#'                     core_control = "control", s_min = 0.05, s_max = 1,
#'                     dist = "normal", sigma = 0.5,
#'                     theta_dist = "pareto", alpha = 4, beta = 3/4)
#'
#' r_range = rbind(c(2,2,2), c(3,3,3),c(4,4,4),c(5,5,5))
#' selection <- select_r(test_data$Y, r_range, asymm = FALSE)


select_r = function(Y,r_range,asymm = F){

  imat = F

  if (length(dim(Y)) == 2) {
    cat("matrix case \n")
    imat <- T


    r_range = r_range[,1:2]
  }



  if(asymm == F){

    if(sum(rowSums(r_range)/r_range[,1]) != 3*dim(r_range)[1] & imat == F){
      warning("all the modes share the same number of clusters in symmetric case")
      return()
    }
  }

  if(sum(rowSums(r_range)/r_range[,1]) != 2*dim(r_range)[1] & imat == T){
    warning("all the modes share the same number of clusters in matrix case")
    return()
  }


  bic = rep(0,dim(r_range)[1])
  p = dim(Y)

  for (i in 1:dim(r_range)[1]) {

    r =  r_range[i,]
    cat("given r = ",r, "\n")

    # obtain z_hat
    initial = wkmeans(Y, r, asymm = asymm)
    z_hat = angle_iteration(Y,initial$z0, max_iter = 20, asymm = asymm)$z

    # obtain S_hat
    if(imat == T){
      dim(Y) <- c(dim(Y), 1)
    }
    S_hat = updateS(Y,z_hat,imat)

    # obtain theta_hat
    theta_hat = theta_estimate(Y,z_hat, imat)

    # obtain bic
    if(imat == T){

      X_hat = ttl(as.tensor(S_hat[z_hat[[1]], z_hat[[2]],]), list(diag(theta_hat[[1]]),diag(theta_hat[[2]])), ms = c(1,2))@data
      dim(X_hat) = c(dim(X_hat),1)
      if(asymm == F){
        bic[i] = p[1]^2*log(sum((X_hat - Y)^2)) +(r[1]^2 + p[1]*log(r[1]) + p[1] - r[1])*log(p[1]^2)
      }else if(asymm == T){
        bic[i] = prod(p)*log(sum((X_hat - Y)^2)) + (prod(r) + sum(p*log(r) + p)- sum(r))*log(prod(p))
      }
    }else if(imat == F){
      X_hat = ttl(as.tensor(S_hat[z_hat[[1]], z_hat[[2]], z_hat[[3]]]), list(diag(theta_hat[[1]]),diag(theta_hat[[2]]), diag(theta_hat[[3]])), ms = c(1,2,3))@data

      if(asymm == F){
        bic[i] = p[1]^3*log(sum((X_hat - Y)^2)) +(r[1]^3 + p[1]*log(r[1]) + p[1] - r[1])*log(p[1]^3)
      }else if(asymm == T){
        bic[i] = prod(p)*log(sum((X_hat - Y)^2)) + (prod(r) + sum(p*log(r) + p)- sum(r))*log(prod(p))
      }

    }

    if(imat == T){
      dim(Y) = dim(Y)[1:2]
    }

  }# for r
  return(list(r = r_range[which.min(bic),], BIC = bic))

}
