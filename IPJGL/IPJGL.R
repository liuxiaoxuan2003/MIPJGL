IPJGL <- function(
    data,
    lambda1 = 4.6,
    lambda2 = 5.2,
    err.threshold = 1e-4,
    step.max = 1e2,
    truncate = 1e-5,
    normtype = '2'){
    source('./IPJGL/soft.tau1.R')
    source('./IPJGL/soft.tauV.R')
    source('./IPJGL/solve.theta.R')
    source('./IPJGL/Dnorm.R')
    library('progress')
    # 初始化存储每次迭代结果的列表
    theta_arr <- list()
    Z_arr <- list()
    V_arr <- list()
    W_arr <- list()
    O_arr <- list()
    F_arr <- list()
    G_arr <- list()
    Q_arr <- list()
    sigma_arr <- list()
    theta_arr_old <- list()
    for (i in 1:T) {  
        # 从X_arr中获取数据生成sigma_arr
        sigma_arr[[i]] <- cov(data$X_arr[[i]])

        # 生成初始值
        n <- nrow(data$X_arr[[i]])
        p <- ncol(data$X_arr[[i]])

        I <- diag(p)

        theta <- I
        Z <- I
        V <- I
        W <- I
        O <- matrix(0, p, p)
        F <- matrix(0, p, p)
        G <- matrix(0, p, p)
        Q <- matrix(0, p, p)

        # 存储初始值到相应的列表中
        theta_arr[[i]] <- theta
        theta_arr_old[[i]] <- I
        Z_arr[[i]] <- Z
        V_arr[[i]] <- V
        W_arr[[i]] <- W
        O_arr[[i]] <- O
        F_arr[[i]] <- F
        G_arr[[i]] <- G
        Q_arr[[i]] <- Q
    }

    rho0 <- 50

    miu <- 5

    rho.max <- 1e10

    step.max <- 1e2

    rho <- rho0


    rho.steps <- floor(log10(rho.max) / log10(miu))

    pb <- progress_bar$new(total = rho.steps * step.max)

    # ADMM算法
    for (j in 1:rho.steps) {
        for (step in 1:step.max) {
            pb$tick()
            for (i in 1:T) {
                if (i == 1){
                    theta_arr[[i]] = solve.theta((n * sigma_arr[[i]] + F_arr[[i]] + Q_arr[[i]]) - rho * (theta_arr[[i+1]] + (V_arr[[i]] + W_arr[[i]]) + Z_arr[[i]]),
                            V_arr[[i]],
                            rho,
                            n,
                            lambda2,
                            normtype)
                }
                else if (i == T){
                    theta_arr[[i]] = solve.theta((n * sigma_arr[[i]] - F_arr[[i - 1]] + Q_arr[[i - 1]]) - rho * (theta_arr[[i - 1]] - (V_arr[[i - 1]] + W_arr[[i - 1]]) + Z_arr[[i]]),
                        V_arr[[i - 1]],
                        rho,
                        n,
                        lambda2,
                        normtype)
                }
                else{
                    theta_arr[[i]] = solve.theta1((n * sigma_arr[[i]] + F_arr[[i]] - F_arr[[i - 1]] + Q_arr[[i]] - Q_arr[[i - 1]]) - rho * (theta_arr[[i - 1]] + theta_arr[[i + 1]] + (V_arr[[i]] + W_arr[[i]]) - (V_arr[[i - 1]] + W_arr[[i - 1]]) + Z_arr[[i]]),
                        V_arr[[i - 1]],
                        V_arr[[i]],
                        rho,
                        n,
                        lambda2,
                        normtype)
                }
            }
            for (i in 1: (T - 1)){
                V_arr[[i]] <- soft.tauV(F_arr[[i]] - G_arr[[i]] + rho * (theta_arr[[i]] - theta_arr[[i+1]] - W_arr[[i]] + t(W_arr[[i]])),
                                    theta_arr[[i]],
                                    theta_arr[[i+1]],
                                    lambda2,
                                    rho,
                                    normtype)
            }

            for (i in 1: T){
                Z_arr[[i]] <- soft.tau1(theta_arr[[i]] + Q_arr[[i]] / rho, lambda1 / rho)
            }
            for (i in 1: (T - 1)){
                W_arr[[i]] <- 0.5 * (t(V_arr[[i]]) - V_arr[[i]] + theta_arr[[i]] - theta_arr[[i+1]]) + 0.5 / rho * (F_arr[[i]] + t(G_arr[[i]]))
            }
            for (i in 1: (T - 1)){
                F_arr[[i]] <- F_arr[[i]] + rho * (theta_arr[[i]] - theta_arr[[i+1]] - (V_arr[[i]] + W_arr[[i]]))
            }
            for (i in 1: (T - 1)){
                G_arr[[i]] <- G_arr[[i]] + rho * (V_arr[[i]] - t(W_arr[[i]]))
            }

            for (i in 1: T){
                Q_arr[[i]] <- Q_arr[[i]] + rho * (theta_arr[[i]] - Z_arr[[i]])
            }
            err <- norm(as.matrix(theta_arr[[1]] - theta_arr_old[[1]]), 'f') / norm(theta_arr_old[[1]], 'f')
            for (i in 2: T){
                err = max(err, norm(as.matrix(theta_arr[[i]] - theta_arr_old[[i]]), 'f') / norm(theta_arr_old[[i]], 'f'))
            }
            
            if (err <= err.threshold) {
                break
            }
            for (i in 1: T){
                theta_arr_old[[i]] = theta_arr[[i]]
            }
        }
        rho <- rho * miu
    }

    # 打印算法结束信息
    cat(paste0('\n algorithm is done! relative error=', round(err, 8)))

    # 对结果进行处理
    for (i in 1: T){
        theta_arr[[i]] <- (theta_arr[[i]] + t(theta_arr[[i]])) / 2
    }
    for (i in 1: T){
        Z_arr[[i]] <- (Z_arr[[i]] + t(Z_arr[[i]])) / 2
    }
    for (i in 1: T){
        theta_arr[[i]][abs(theta_arr[[i]]) < truncate] <- 0
    }
    for (i in 1: T){
        Z_arr[[i]][abs(Z_arr[[i]]) < truncate] <- 0
    }
    for (i in 1: (T - 1)){
        V_arr[[i]][abs(V_arr[[i]]) < truncate] <- 0
    }
    return(
        list(
            theta_arr = theta_arr,
            Z_arr = Z_arr,
            V_arr = V_arr
        )
    )
}
# IPJGL <- function(
#     data,
#     lambda1 = 4.6,
#     lambda2 = 5.2,
#     err.threshold = 1e-4,
#     step.max = 1e2,
#     truncate = 1e-5,
#     normtype = '2'){
#     source('./IPJGL/soft.tau1.R')
#     source('./IPJGL/soft.tauV.R')
#     source('./IPJGL/solve.theta.R')
#     source('./IPJGL/Dnorm.R')
#     library('progress')
#     # 初始化存储每次迭代结果的列表
#     theta_arr <- list()
#     Z_arr <- list()
#     V_arr <- list()
#     W_arr <- list()
#     O_arr <- list()
#     F_arr <- list()
#     G_arr <- list()
#     Q_arr <- list()
#     sigma_arr <- list()
#     theta_arr_old <- list()
#     for (i in 1:T) {  
#         # 从X_arr中获取数据生成sigma_arr
#         sigma_arr[[i]] <- cov(data$X_arr[[i]])

#         # 生成初始值
#         n <- nrow(data$X_arr[[i]])
#         p <- ncol(data$X_arr[[i]])

#         I <- diag(p)

#         theta <- I
#         Z <- I
#         V <- I
#         W <- I
#         O <- matrix(0, p, p)
#         F <- matrix(0, p, p)
#         G <- matrix(0, p, p)
#         Q <- matrix(0, p, p)

#         # 存储初始值到相应的列表中
#         theta_arr[[i]] <- theta
#         theta_arr_old[[i]] <- I
#         Z_arr[[i]] <- Z
#         V_arr[[i]] <- V
#         W_arr[[i]] <- W
#         O_arr[[i]] <- O
#         F_arr[[i]] <- F
#         G_arr[[i]] <- G
#         Q_arr[[i]] <- Q
#     }

#     rho0 <- 50

#     miu <- 5

#     rho.max <- 1e10

#     step.max <- 1e2

#     rho <- rho0


#     rho.steps <- floor(log10(rho.max) / log10(miu))

#     pb <- progress_bar$new(total = rho.steps * step.max)

#     # ADMM算法
#     for (j in 1:rho.steps) {
#         for (step in 1:step.max) {
#             pb$tick()
#             for (i in 1:T) {  # 对每一组数据进行更新
#                 if (i == 1){
#                     theta_arr[[i]] = solve.theta((n * sigma_arr[[i]] + F_arr[[i]] + Q_arr[[i]]) - rho * (theta_arr[[i+1]] + (V_arr[[i]] + W_arr[[i]]) + Z_arr[[i]]),
#                             V_arr[[i]],
#                             rho,
#                             n,
#                             lambda2,
#                             normtype)
#                 }
#                 else if (i == T){
#                     theta_arr[[i]] = solve.theta((n * sigma_arr[[i]] - F_arr[[i - 1]] + Q_arr[[i - 1]]) - rho * (theta_arr[[i - 1]] - (V_arr[[i - 1]] + W_arr[[i - 1]]) + Z_arr[[i]]),
#                         V_arr[[i - 1]],
#                         rho,
#                         n,
#                         lambda2,
#                         normtype)
#                 }
#                 else{
#                     theta_arr[[i]] = solve.theta((n * sigma_arr[[i]] + F_arr[[i]] - F_arr[[i - 1]] + Q_arr[[i]]) - rho * (theta_arr[[i - 1]] + theta_arr[[i + 1]] + (V_arr[[i]] + W_arr[[i]]) - (V_arr[[i - 1]] + W_arr[[i - 1]]) + Z_arr[[i]]),
#                         V_arr[[i]] + V_arr[[i - 1]],
#                         rho,
#                         n,
#                         lambda2,
#                         normtype)
#                 }
#             }
#             for (i in 1: (T - 1)){
#                 V_arr[[i]] <- soft.tauV(F_arr[[i]] - G_arr[[i]] + rho * (theta_arr[[i]] - theta_arr[[i+1]] - W_arr[[i]] + t(W_arr[[i]])),
#                                     theta_arr[[i]],
#                                     theta_arr[[i+1]],
#                                     lambda2,
#                                     rho,
#                                     normtype)
#             }

#             for (i in 1: T){
#                 Z_arr[[i]] <- soft.tau1(theta_arr[[i]] + Q_arr[[i]] / rho, lambda1 / rho)
#             }
#             for (i in 1: (T - 1)){
#                 W_arr[[i]] <- 0.5 * (t(V_arr[[i]]) - V_arr[[i]] + theta_arr[[i]] - theta_arr[[i+1]]) + 0.5 / rho * (F_arr[[i]] + t(G_arr[[i]]))
#                 F_arr[[i]] <- F_arr[[i]] + rho * (theta_arr[[i]] - theta_arr[[i+1]] - (V_arr[[i]] + W_arr[[i]]))
#                 G_arr[[i]] <- G_arr[[i]] + rho * (V_arr[[i]] - t(W_arr[[i]]))
#             }

#             for (i in 1: T){
#                 Q_arr[[i]] <- Q_arr[[i]] + rho * (theta_arr[[i]] - Z_arr[[i]])
#             }
#             err <- norm(as.matrix(theta_arr[[1]] - theta_arr_old[[1]]), 'f') / norm(theta_arr_old[[1]], 'f')
#             for (i in 2: T){
#                 err = max(err, norm(as.matrix(theta_arr[[i]] - theta_arr_old[[i]]), 'f') / norm(theta_arr_old[[i]], 'f'))
#             }
            
#             if (err <= err.threshold) {
#                 break
#             }
#             for (i in 1: T){
#                 theta_arr_old[[i]] = theta_arr[[i]]
#             }
#         }
#         rho <- rho * miu
#     }

#     # 打印算法结束信息
#     cat(paste0('\n algorithm is done! relative error=', round(err, 8)))

#     # 对结果进行处理
#     for (i in 1: T){
#         theta_arr[[i]] <- (theta_arr[[i]] + t(theta_arr[[i]])) / 2
#     }
#     for (i in 1: T){
#         Z_arr[[i]] <- (Z_arr[[i]] + t(Z_arr[[i]])) / 2
#     }
#     for (i in 1: T){
#         theta_arr[[i]][abs(theta_arr[[i]]) < truncate] <- 0
#     }
#     for (i in 1: T){
#         Z_arr[[i]][abs(Z_arr[[i]]) < truncate] <- 0
#     }
#     for (i in 1: (T - 1)){
#         V_arr[[i]][abs(V_arr[[i]]) < truncate] <- 0
#     }
#     return(
#         list(
#             theta_arr = theta_arr,
#             Z_arr = Z_arr,
#             V_arr = V_arr
#         )
#     )
# }