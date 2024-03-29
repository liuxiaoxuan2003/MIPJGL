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
delta.true <- data$delta_arr
print(length(data$delta_arr))
##### IPJGL #####
results.IPJGL <-
  IPJGL(
    data = data
  )
#delta.IPJGL <- theta2partial(results.IPJGL$Z1)-theta2partial(results.IPJGL$Z2)

delta.IPJGL <- list()
print(length(results.IPJGL$Z_arr))
# 循环遍历相邻的两个部分
for (i in 2: T ) {
  Z1 <- results.IPJGL$Z_arr[[i-1]]  
  Z2 <- results.IPJGL$Z_arr[[i]]
  
  # 计算theta2partial
  delta <- theta2partial(Z1) - theta2partial(Z2)
  
  # 将计算结果附加到delta.IPJGL中
  delta.IPJGL[[i]] <- delta
}

p<-50
keep.edges <- p

is.numeric(keep.edges)

# 循环绘制图形
for (i in 1:T) {
  # 从 delta_true 中获取第 i 组数据
  delta_true_first <- delta.true[[i]]
  
  # 检查 delta_true_first 和 keep.edges 是否符合要求
  if (!is.numeric(delta_true_first) || !is.numeric(keep.edges)) {
    cat("Skipping plot for Delta True and Delta IPJGL", i, "due to invalid parameters.\n")
    next  # 跳过当前循环，继续下一个
  }
  
  # 从 delta_IPJGL 中获取第 i 组数据
  delta_IPJGL_first <- delta.IPJGL[[i]]
  
  # 检查 delta_IPJGL_first 是否符合要求
  if (!is.numeric(delta_IPJGL_first)) {
    cat("Skipping plot for Delta IPJGL", i, "due to invalid parameters.\n")
    next  # 跳过当前循环，继续下一个
  }
  
  # 绘制 Delta True 图形
  p1 <- showmatrix(keep.largest.N(delta_true_first, keep.edges), main = paste("Delta True", i))
  print(p1)

}


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

