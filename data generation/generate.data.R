generate.data <-
  function(n1,
           n2,
           p,
           rate.connect = 0.1,
           rate.drop = 0.8,
           m.pert = 10,
           umin = 0.3,
           umax = 0.8,
           T = 10,
           diffmode) {
    source('./data generation/SFNG.R')
    source('./data generation/theta2partial.R')
    library(MASS)


    net.structure <- SFNG(p, 2, 1)
    diag(net.structure) <- 0



    dense = matrix(runif(p * p), p)
    dense = (dense - 0.5) / 0.5 * (umax - umin) + umin * sign(dense - 0.5)
    dense[upper.tri(dense)] <- 0
    dense <- dense + t(dense)

    theta1 <- net.structure * dense
    theta2 <- theta1
    theta_arr <- list()
    delta_arr <- list()
    sigma_arr <- list()
    X_arr <- list()
    index.pert_arr <- list()
    theta_arr[[1]] <- theta1

    for(theta_idx in 2:T){
        theta_init = theta_arr[[theta_idx - 1]]
        theta_now = theta_init

        all.index <- 1:p

        #group genes
        score <- rep(0, p)
        for (i in 1:p) {
            score[i] <- sum(net.structure[, i])
        }
        index.sort <- order(score, decreasing = TRUE)

        index.vital <- index.sort[1:ceiling(0.2 * p)]
        index.moderate <-
        index.sort[(ceiling(0.2 * p) + 1):ceiling(0.4 * p)]
        index.original <- index.sort[(ceiling(0.4 * p) + 1):p]

        # select diffgene
        if (diffmode == 1) {
            index.pert <- sample(index.vital, m.pert, replace = FALSE)
        } 
        else if (diffmode == 2) {
            index.pert <- sample(index.moderate, m.pert, replace = FALSE)
        } 
        else if (diffmode == 3) {
            index.pert <- sample(index.original, m.pert, replace = FALSE)
        } 
        else if (diffmode == 4) {
            index.pert <- sample(index.vital, m.pert / 2, replace = FALSE)
            index.pert <-c(index.pert, sample(index.original, m.pert / 2, replace = FALSE))
        } 
        else if (diffmode == 5) {
            index.pert <- sample(index.vital, m.pert / 10, replace = FALSE)
            index.pert <- c(index.pert, sample(index.moderate, m.pert / 10 * 4, replace = FALSE))
            index.pert <- c(index.pert, sample(index.original, m.pert / 10 * 5, replace = FALSE))
        }


        #  pertubation
        for (pert in index.pert) {
            connect.now <- which(net.structure[pert, ] != 0)
            disconnect.now <- all.index[-c(connect.now, pert)]

            disconnect.index <- which(runif(p) < rate.drop)
            disconnect.index <- intersect(connect.now, disconnect.index)

            theta_now[pert, disconnect.index] <- 0
            theta_now[disconnect.index, pert] <- 0

            connect.index <- which(runif(p) < rate.connect)
            connect.index <- intersect(disconnect.now, connect.index)

            l3 <- length(connect.index)
            connect.dense <- runif(l3)
            connect.dense <-
                (connect.dense - 0.5) / 0.5 * umin + (umax - umin) * sign(connect.dense -
                                                                            0.5)
            theta_now[pert, connect.index] <- connect.dense
            theta_now[connect.index, pert] <- connect.dense
            }
        if (sum((theta_init - theta_now) != 0) == 0) {
            for (pert in index.pert) {
                connect.now <- which(net.structure[pert, ] != 0)
                disconnect.now <- all.index[-c(connect.now, pert)]

                disconnect.index <- which(runif(p) < rate.drop)
                disconnect.index <-
                intersect(connect.now, disconnect.index)

                theta_now[pert, disconnect.index] <- 0
                theta_now[disconnect.index, pert] <- 0

                connect.index <- which(runif(p) < rate.connect)
                connect.index <- intersect(disconnect.now, connect.index)

                l3 <- length(connect.index)
                connect.dense <- runif(l3)
                connect.dense <-
                (connect.dense - 0.5) / 0.5 * umin + (umax - umin) * sign(connect.dense -
                                                                            0.5)
                theta_now[pert, connect.index] <- connect.dense
                theta_now[connect.index, pert] <- connect.dense
            }
        }

        eig1 <- eigen(theta_init)$values
        eig2 <- eigen(theta_now)$values
        eig.min <- min(c(eig1, eig2))
        eig.diag <- (abs(eig.min) + 0.1) * diag(p)
        theta_init <- theta_init + eig.diag
        theta_now <- theta_now + eig.diag

        X1 <- mvrnorm(n1, rep(0, p), solve(theta_init))
        X2 <- mvrnorm(n2, rep(0, p), solve(theta_now))

        sigma1 <- cov(X1)
        sigma2 <- cov(X2)

        thetap <- lapply(list(theta_init, theta_now), theta2partial)
        delta.true <- thetap[[1]] - thetap[[2]]
        if(theta_idx == 2){
            X_arr[[1]] = X1
            X_arr[[2]] = X2
            theta_arr[[1]] = theta_init
            theta_arr[[2]] = theta_now
            delta_arr[[1]] = delta.true
            sigma_arr[[1]] = sigma1
            sigma_arr[[2]] = sigma2
            index.pert_arr[[1]] = index.pert
        }
        else{
            X_arr[[theta_idx]] = X2
            theta_arr[[theta_idx]] = theta_now
            delta_arr[[theta_idx]] = delta.true
            sigma_arr[[theta_idx]] = sigma2
            index.pert_arr[[theta_idx]] = index.pert
        }
    }
    return(list(
        theta_arr = theta_arr,
        sigma_arr = sigma_arr,
        delta_arr = delta_arr,
        X_arr = X_arr,
        index.pert_arr = index.pert_arr
    ))
  }



# n1=1000
# n2=1000
# p=50
# rate.connect = 0.1
# rate.drop = 0.8
# m.pert = 10
# umin = 0.3
# umax = 0.8
# diffmode=3
# T = 10
# source('./data generation/SFNG.R')
# source('./data generation/theta2partial.R')
# library(MASS)


# net.structure <- SFNG(p, 2, 1)
# diag(net.structure) <- 0



# dense = matrix(runif(p * p), p)
# dense = (dense - 0.5) / 0.5 * (umax - umin) + umin * sign(dense - 0.5)
# dense[upper.tri(dense)] <- 0
# dense <- dense + t(dense)

# theta1 <- net.structure * dense
# theta2 <- theta1
# theta_arr <- list()
# delta_arr <- list()
# sigma_arr <- list()
# X_arr <- list()
# index.pert_arr <- list()
# theta_arr[[1]] <- theta1

# for(theta_idx in 2:T){
#     theta_init = theta_arr[[theta_idx - 1]]
#     theta_now = theta_init

#     all.index <- 1:p

#     #group genes
#     score <- rep(0, p)
#     for (i in 1:p) {
#         score[i] <- sum(net.structure[, i])
#     }
#     index.sort <- order(score, decreasing = TRUE)

#     index.vital <- index.sort[1:ceiling(0.2 * p)]
#     index.moderate <-
#     index.sort[(ceiling(0.2 * p) + 1):ceiling(0.4 * p)]
#     index.original <- index.sort[(ceiling(0.4 * p) + 1):p]

#     # select diffgene
#     if (diffmode == 1) {
#         index.pert <- sample(index.vital, m.pert, replace = FALSE)
#     } 
#     else if (diffmode == 2) {
#         index.pert <- sample(index.moderate, m.pert, replace = FALSE)
#     } 
#     else if (diffmode == 3) {
#         index.pert <- sample(index.original, m.pert, replace = FALSE)
#     } 
#     else if (diffmode == 4) {
#         index.pert <- sample(index.vital, m.pert / 2, replace = FALSE)
#         index.pert <-c(index.pert, sample(index.original, m.pert / 2, replace = FALSE))
#     } 
#     else if (diffmode == 5) {
#         index.pert <- sample(index.vital, m.pert / 10, replace = FALSE)
#         index.pert <- c(index.pert, sample(index.moderate, m.pert / 10 * 4, replace = FALSE))
#         index.pert <- c(index.pert, sample(index.original, m.pert / 10 * 5, replace = FALSE))
#     }


#     #  pertubation
#     for (pert in index.pert) {
#         connect.now <- which(net.structure[pert, ] != 0)
#         disconnect.now <- all.index[-c(connect.now, pert)]

#         disconnect.index <- which(runif(p) < rate.drop)
#         disconnect.index <- intersect(connect.now, disconnect.index)

#         theta_now[pert, disconnect.index] <- 0
#         theta_now[disconnect.index, pert] <- 0

#         connect.index <- which(runif(p) < rate.connect)
#         connect.index <- intersect(disconnect.now, connect.index)

#         l3 <- length(connect.index)
#         connect.dense <- runif(l3)
#         connect.dense <-
#             (connect.dense - 0.5) / 0.5 * umin + (umax - umin) * sign(connect.dense -
#                                                                         0.5)
#         theta_now[pert, connect.index] <- connect.dense
#         theta_now[connect.index, pert] <- connect.dense
#         }
#     if (sum((theta_init - theta_now) != 0) == 0) {
#         for (pert in index.pert) {
#             connect.now <- which(net.structure[pert, ] != 0)
#             disconnect.now <- all.index[-c(connect.now, pert)]

#             disconnect.index <- which(runif(p) < rate.drop)
#             disconnect.index <-
#             intersect(connect.now, disconnect.index)

#             theta_now[pert, disconnect.index] <- 0
#             theta_now[disconnect.index, pert] <- 0

#             connect.index <- which(runif(p) < rate.connect)
#             connect.index <- intersect(disconnect.now, connect.index)

#             l3 <- length(connect.index)
#             connect.dense <- runif(l3)
#             connect.dense <-
#             (connect.dense - 0.5) / 0.5 * umin + (umax - umin) * sign(connect.dense -
#                                                                         0.5)
#             theta_now[pert, connect.index] <- connect.dense
#             theta_now[connect.index, pert] <- connect.dense
#         }
#     }

#     eig1 <- eigen(theta_init)$values
#     eig2 <- eigen(theta_now)$values
#     eig.min <- min(c(eig1, eig2))
#     eig.diag <- (abs(eig.min) + 0.1) * diag(p)
#     theta_init <- theta_init + eig.diag
#     theta_now <- theta_now + eig.diag

#     X1 <- mvrnorm(n1, rep(0, p), solve(theta_init))
#     X2 <- mvrnorm(n2, rep(0, p), solve(theta_now))

#     sigma1 <- cov(X1)
#     sigma2 <- cov(X2)

#     thetap <- lapply(list(theta_init, theta_now), theta2partial)
#     delta.true <- thetap[[1]] - thetap[[2]]
#     if(theta_idx == 2){
#         X_arr[[1]] = X1
#         X_arr[[2]] = X2
#         theta_arr[[1]] = theta_init
#         theta_arr[[2]] = theta_now
#         delta_arr[[1]] = delta.true
#         sigma_arr[[1]] = sigma1
#         sigma_arr[[2]] = sigma2
#         index.pert_arr[[1]] = index.pert
        
#     }
#     else{
#         X_arr[[theta_idx]] = X2
#         theta_arr[[theta_idx]] = theta_now
#         delta_arr[[theta_idx]] = delta.true
#         sigma_arr[[theta_idx]] = sigma2
#         index.pert_arr[[theta_idx]] = index.pert
#     }
# }


