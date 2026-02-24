rm(list = ls())

library(dplyr)
library(tidyr)

data  = read.csv("Windows_Ave_example_1.csv", check.names = FALSE, header = FALSE)


val_indices <- grep("val", data[, 1], ignore.case = TRUE)
num_vals <- length(val_indices)
result <- data.frame(matrix(0, ncol=100,dimnames=list(NULL,paste0('Time',1:100))))
for (i in 1:num_vals) {
  if (i<100){
    a <- val_indices[i]+1
    b <- val_indices[i+1]-1
    for (j in 1:ncol(data)){
      print(i)
      print(j)
      c <- val_indices[j]+1
      d <- val_indices[j+1]-1
      result[c:d, i] <- data[a:b, j] 
    }
  }
  else {
    a <- val_indices[i]+1
    e <- nrow(data)
    for (j in 1:ncol(data)){
      print(i)
      print(j)
      c <- val_indices[j]+1
      d <- val_indices[j+1]-1
      result[c:d, i] <- data[a:e,j]
  }
  }
}

result <- result[-1,]
write.csv(result,file = "T-Windows_Ave_Example_1.csv")
save(result,file = "output/output.Rdata")

