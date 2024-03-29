rm(list = ls(all = TRUE))
source('source.scripts.R')

data <- generate.data(
  n1,
  n2,
  p,
  rates.connect,
  rates.drop,
  m.pert,
  umin = 0.3,
  umax = 0.8,
  diffmode = 5,
)
delta_true_arr <- data$delta_arr
print(length(data$delta_arr))
##### IPJGL #####
results.IPJGL <-
  IPJGL(
    data = data
  )
theta_pre_arr = results.IPJGL$theta_arr
Z_pre_arr = results.IPJGL$Z_arr
V_pre_arr = results.IPJGL$V_arr

keep.edges <- p

delta_pre_arr <- list()
for (i in 2: T){
    delta_pre_arr[[i]] = theta2partial(Z_pre_arr[[i - 1]]) - theta2partial(Z_pre_arr[[i]])
}

p_true_arr <- list()
p_pred_arr <- list()
for (i in 2: T){
    p_true_arr[[i]] = showmatrix(keep.largest.N(delta_true_arr[[i]], keep.edges), main = 'truth')
    p_pred_arr[[i]] = showmatrix(keep.largest.N(delta_pre_arr[[i]], keep.edges), main = 'MIPJGL')
}



# delta.IPJGL <- list()
# print(length(results.IPJGL$Z_arr))
# # 循环遍历相邻的两个部分
# for (i in 2: T ) {
#   Z1 <- results.IPJGL$Z_arr[[i-1]]  
#   Z2 <- results.IPJGL$Z_arr[[i]]
  
#   # 计算theta2partial
#   delta <- theta2partial(Z1) - theta2partial(Z2)
  
#   # 将计算结果附加到delta.IPJGL中
#   delta.IPJGL[[i]] <- delta
# }




# rm(list = ls(all = TRUE))
# source('source.scripts.R')

# data <- generate.data(
#   n1,
#   n2,
#   p,
#   rates.connect,
#   rates.drop,
#   m.pert,
#   umin = 0.3,
#   umax = 0.8,
#   diffmode = 3,
# )
# delta.true <- data$delta.true

# ##### IPJGL #####
# results.IPJGL <-
#   IPJGL(
#     data = data
#   )
# delta.IPJGL <- theta2partial(results.IPJGL$Z1)-theta2partial(results.IPJGL$Z2)

# keep.edges <- p
# p1 <- showmatrix(keep.largest.N(delta.true,keep.edges),main = 'truth')
# p2 <- showmatrix(keep.largest.N(delta.IPJGL,keep.edges),main = 'IPJGL')

