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
##### TVIPJGL #####
results.TVIPJGL <-
  TVIPJGL(
    data = data
  )
theta_pre_arr = results.TVIPJGL$theta_arr
Z_pre_arr = results.TVIPJGL$Z_arr
V_pre_arr = results.TVIPJGL$V_arr

keep.edges <- p

delta_pre_arr <- list()
for (i in 2: T){
    delta_pre_arr[[i]] = theta2partial(Z_pre_arr[[i - 1]]) - theta2partial(Z_pre_arr[[i]])
}

p_true_arr <- list()
p_pred_arr <- list()
for (i in 2: T){
    save_true_path = sprintf("./fig/delta_true%d-%d.pdf", i-1, i)
    save_pred_path = sprintf("./fig/delta_pred%d-%d.pdf", i-1, i)
    p_true_arr[[i]] = showmatrix(keep.largest.N(delta_true_arr[[i]], keep.edges), main = 'truth', save_path=save_true_path)
    p_pred_arr[[i]] = showmatrix(keep.largest.N(delta_pre_arr[[i]], keep.edges), main = 'MTVIPJGL', save_path=save_pred_path)
}

# 
# delta.TVIPJGL <- list()
# print(length(results.TVIPJGL$Z_arr))
# # 循环遍历相邻的两个部分
# for (i in 2: T ) {
#   Z1 <- results.TVIPJGL$Z_arr[[i-1]]  
#   Z2 <- results.TVIPJGL$Z_arr[[i]]
  
#   # 计算theta2partial
#   delta <- theta2partial(Z1) - theta2partial(Z2)
  
#   # 将计算结果附加到delta.TVIPJGL中
#   delta.TVIPJGL[[i]] <- delta
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

# ##### TVIPJGL #####
# results.TVIPJGL <-
#   TVIPJGL(
#     data = data
#   )
# delta.TVIPJGL <- theta2partial(results.TVIPJGL$Z1)-theta2partial(results.TVIPJGL$Z2)

# keep.edges <- p
# p1 <- showmatrix(keep.largest.N(delta.true,keep.edges),main = 'truth')
# p2 <- showmatrix(keep.largest.N(delta.TVIPJGL,keep.edges),main = 'TVIPJGL')

