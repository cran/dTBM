########## Bricks for dTBM ############
library(tensorregress)
library(WeightedCluster) # wcKMedoids
library(EnvStats) # patero



renumber <- function(z) {
  z_re <- rep(0, length(z))
  uniq <- unique(z)

  for (i in 1:length(uniq)) {
    z_re[z == uniq[i]] <- i
  }
  return(z_re)
}

single_wkmeans = function(X, r){
  z <- rep(0, dim(X)[1])
  sc <- 1:length(z)

  l2 <- apply(X, 1, function(x) sqrt(sum(x^2))) # l2 norm

  if (length(which(l2 == 0) > 0)) {
    sc <- which(l2 != 0)
    X <- X[sc, ]
    l2 <- l2[sc]
  }

  Xs <- diag(l2^(-1)) %*% X

  ## weighted kmeans
  diss <- dist(Xs, method = "euclidean", p = 2, upper = T, diag = T)

  z[sc] <- wcKMedoids(diss^2, k = r, weights = l2^2, method = "PAMonce", cluster.only = T)
  z[-sc] <- sample(unique(z[sc]), length(z[-sc]), replace = T)


  s0 <- setdiff(1:length(z), sc)
  z <- as.factor(z)
  levels(z) <- 1:r

  return(list(z = as.numeric(z), s0 = s0))
}


updateS <- function(Y, z, imat) {
  # z is a list

  r1 = length(unique(z[[1]]))
  r2 = length(unique(z[[2]]))


  #r <- length(unique(z))

  if(imat == F){
    r3 = length(unique(z[[3]]))
    S <- array(0, dim = c(r1, r2, r3))

    for (a in 1:r1) {
      for (b in 1:r2) {
        for (c in 1:r3) {
          S[a, b, c] <- mean(Y[z[[1]] == a, z[[2]] == b, z[[3]] == c], na.rm = T)
        }
      }
    }
  }else if(imat == T){
    S <- array(0, dim = c(r1, r2, 1))

    for (a in 1:r1) {
      for (b in 1:r2) {
        S[a, b, 1] <- mean(as.matrix(Y[z[[1]] == a, z[[2]] == b,1]), na.rm = T)
      }
    }
  }

  return(S)
}


Cal_Y1 <- function(Y, z, imat) {
  r2 = length(unique(z[[2]]))


  # imat = F
  # if(length(dim(Y)) == 2){
  #   imat = T
  # }

  p <- dim(Y)[1]

  if(imat == F){
    r3 = length(unique(z[[3]]))
    Y1 <- array(0, dim = c(p, r2, r3))
    for (b in 1:r2) {
      for (c in 1:r3) {
        Y_cut = Y[, z[[2]] == b, z[[3]] == c]
        dim(Y_cut) = c(p, sum(z[[2]] == b), sum(z[[3]] == c))
        Y1[, b, c] <- apply(Y_cut, 1, mean)
      }
    }
  }else if(imat == T){
    Y1 <- array(0, dim = c(p, r2, 1))
    for (b in 1:r2) {
      Y1[, b, 1] <- apply(as.matrix(Y[, z[[2]] == b,1]), 1, mean)
      #Y1[, b, 1] <- mean(Y[, z == b,1])
    }
  }

  return(Y1)
}




Cal_Y2 <- function(Y, z, imat) {
  r1 = length(unique(z[[1]]))

  # imat = F
  # if(length(dim(Y)) == 2){
  #   imat = T
  # }

  p <- dim(Y)[2]

  if(imat == F){
    r3 = length(unique(z[[3]]))
    Y2 <- array(0, dim = c(r1, p, r3))
    for (a in 1:r1) {
      for (c in 1:r3) {
        Y_cut = Y[z[[1]] == a, , z[[3]] == c]
        dim(Y_cut) = c( sum(z[[1]] == a),p, sum(z[[3]] == c))
        Y2[a, , c] <- apply(Y_cut, 2, mean)
      }
    }
  }else if(imat == T){
    Y2 <- array(0, dim = c( r1, p, 1))
    for (a in 1:r1) {
      Y2[a, , 1] <- apply(as.matrix(Y[z[[1]] == a,,1]), 2, mean)
      #Y1[, b, 1] <- mean(Y[, z == b,1])
    }
  }

  return(Y2)
}

Cal_Y3 <- function(Y, z) {
  r1 = length(unique(z[[1]]))
  r2 = length(unique(z[[2]]))

  p <- dim(Y)[3]

  Y3 <- array(0, dim = c(r1, r2, p))
  for (a in 1:r1) {
    for (b in 1:r2) {
      Y_cut = Y[z[[1]] == a, z[[2]] == b , ]
      dim(Y_cut) = c( sum(z[[1]] == a),sum(z[[2]] == b),p)
      Y3[a, b, ] <- apply(Y_cut, 3, mean)
    }
  }

  return(Y3)
}


theta_estimate = function(Y, z, imat){

  # imat = F
  # if(length(dim(Y)) == 2){
  #   imat = T
  # }

  r1 = length(unique(z[[1]]))
  r2 = length(unique(z[[2]]))

  # theta 1
  Y1_unfold <- unfold(as.tensor(Cal_Y1(Y, z, imat)), 1, c(3, 2))@data

  mtheta_hat1 <- apply(Y1_unfold, 1, function(x) sqrt(sum(x^2)))
  for (a in 1:r1) {
    ind = z[[1]] == a
    mtheta_hat1[ind] = mtheta_hat1[ind]*sum(ind)/sum(mtheta_hat1[ind])
  }


  # theta 2
  Y2_unfold <- unfold(as.tensor(Cal_Y2(Y, z, imat)), 2, c(1, 3))@data

  mtheta_hat2 <- apply(Y2_unfold, 1, function(x) sqrt(sum(x^2)))
  for (a in 1:r2) {
    ind = z[[2]] == a
    mtheta_hat2[ind] = mtheta_hat2[ind]*sum(ind)/sum(mtheta_hat2[ind])
  }

  if(imat == T){
    theta_hat = list(mtheta_hat1, mtheta_hat2)
  }else if(imat == F){

    r3 = length(unique(z[[3]]))
    # theta 3
    Y3_unfold <- unfold(as.tensor(Cal_Y3(Y, z)), 3, c(1, 2))@data

    mtheta_hat3 <- apply(Y3_unfold, 1, function(x) sqrt(sum(x^2)))
    for (a in 1:r3) {
      ind = z[[3]] == a
      mtheta_hat3[ind] = mtheta_hat3[ind]*sum(ind)/sum(mtheta_hat3[ind])
    }

    theta_hat = list(mtheta_hat1, mtheta_hat2, mtheta_hat3)
  }

  return(theta_hat)
}


cosin <- function(v1, v2) {
  v1_norm <- sqrt(sum(v1^2))
  v2_norm <- sqrt(sum(v2^2))
  if (v1_norm == 0 | v2_norm == 0) {
    return(cosin = 1)
  } else {
    return(cosin = t(v1) %*% v2 / (v1_norm * v2_norm))
  }
}

updatez <- function(Y1_unfold, S_unfold) {
  z_new <- rep(0, dim(Y1_unfold)[1])

  for (j in 1:length(z_new)) {
    dist <- c()
    for (a in 1:dim(S_unfold)[1]) {
      dist <- c(dist, cosin(Y1_unfold[j, ], S_unfold[a, ]))
    }
    z_new[j] <- which.max(dist)
  }

  return(z_new)
}

single_Aiteration = function(S_unfold,Y_unfold,alpha1){

  z <- rep(0, dim(Y_unfold)[1])
  sc <- 1:length(z)
  s_deg <- F

  lS <- apply(S_unfold, 1, function(x) sum(x^2))
  index <- which(lS == 0)
  if (length(index) > 0) {
    for (i in index) {
      S_unfold[i, sample(1:dim(S_unfold)[2],1)] <- alpha1
      s_deg <- T
    }
  }

  l2 <- apply(Y_unfold, 1, function(x) sum(x^2))
  if (length(which(l2 == 0) > 0)) {
    sc <- which(l2 != 0)
    Y_unfold <- Y_unfold[sc, ]
  }

  # update cluster via anlge distance
  z[sc] <- updatez(Y_unfold, S_unfold)
  z[-sc] <- sample(unique(z[sc]), length(z[-sc]), replace = T)

  return(list(z = z, s_deg = s_deg))
}




sim_S = function(r, s_min, s_max, delta = NULL, imat = F){

  if(imat == F){
    S <- array(s_min, dim = rep(r,3))
    for (i in 1:r) {
      S[i, i, i] <- s_max
    }

    if(!is.null(delta)){
      delta_min = sqrt( sum( S[1,,]^2 ))
      S = S * delta / delta_min
    }
  }else if(imat == T){
    S <- array(s_min, dim = rep(r,2))
    for (i in 1:r) {
      S[i, i] <- s_max
    }

    if(!is.null(delta)){
      delta_min = sqrt( sum(S[1,]^2 ))
      S = S * delta / delta_min
    }

  }

  return(S)
}



generate_theta = function(p, theta_dist, z, alpha = NULL, beta = NULL){
  r = length(unique(z))

  if (theta_dist == "abs_normal") {
    theta <- abs(rnorm(p, 0, 0.5)) + 1 - sqrt(1 / (2 * pi))
  } else if (theta_dist == "pareto") {
    theta <- rpareto(p, location = beta, shape = alpha) # choose alpha*beta = alpha - 1
  } else if (theta_dist == "non"){
    theta = rep(1, p)
  }

  # in group normalization
  for (a in 1:r) {
    index_a = z == a
    theta[index_a] = theta[index_a]*sum(index_a)/sum(theta[index_a])
  }

  return(theta)

}

