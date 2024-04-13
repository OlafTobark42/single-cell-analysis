install.packages("parallel")
library("parallel")


## 一个简单的例子
system.time(for(i in 1:4){Sys.sleep(2)})
## 或者用lapply改写成：
system.time(lapply(1:4, function(i) Sys.sleep(2)))
## 设置并行环境
library(parallel)
## 检测系统可用的核数
detectCores()
## 默认返回的结构逻辑的核数，需修改logical=FALSE，返回物理核数
detectCores(logical=FALSE)
## 建立2核的集群
cl <- makeCluster(20)

## 不使用并行计算
system.time(lapply(1:4, function(i) Sys.sleep(2)))

## 使用parallel包，运行时间减半
## 在非windows系统下，使用mclapply函数
system.time(
  mclapply(1:4, function(i) Sys.sleep(2),mc.cores=20)
)
## 在windows系统下，使用parlapply函数
system.time(
  parlapply(cl, 1:4,function(i) Sys.sleep(2))
)

##关闭集群
stopCluster(cl)


library("future")
plan()
plan("multiprocess",workers=30)
plan("default")
