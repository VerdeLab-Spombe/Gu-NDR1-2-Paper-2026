
rm(list=ls())

# Read document for the vector x,y coordinate 
doc <- read.csv('Example_image.csv', header = T)

#Input the cell migration vector
x_start <- as.numeric(doc[1,2])
y_start <- as.numeric(doc[1,3])

x_end <- as.numeric(doc[1,4])
y_end <- as.numeric(doc[1,5])

#Input the nuclear and Golgi apparatus vector 

output <- data.frame(number = doc[,"Trailing"])
for (i in 1: nrow(doc)) {
 
   print(i)
  
  x_nuclear <- as.numeric(doc[i,7])
  y_nuclear <- as.numeric(doc[i,8])

  x_golgi <- as.numeric(doc[i,9])
  y_golgi <- as.numeric(doc[i,10])
  
  x1 <- x_end - x_start #(x1,y1): the caculated migration vector
  y1 <- y_end - y_start
  
  x2 <- x_golgi - x_nuclear #(x2,y2): the calculated the centroid of nuclear to Golgi apparatus
  y2 <- y_golgi - y_nuclear
  dot = x1*x2 + y1*y2      # dot product
  det = x1*y2 - y1*x2      # determinant
  angle = atan2(det, dot)*180/pi  # atan2(y, x) or atan2(sin, cos)
  output[i,2] <- angle
  
  if (angle < 0) {
    output[i,3] <- 360-abs(angle)
  }
  if (angle >= 0) {
    output[i,3] <- angle
  }
  if (abs(angle) <= 60) {
    output[i,4] <- 'Polarity'
  }
  if (abs(angle) > 60) {
    output[i,4] <- 'Unpolarity'
  }
    
  end
}

for (i in 1: nrow(doc)) {
  
  print(i)
  
  x_nuclear <- as.numeric(doc[i,12])
  y_nuclear <- as.numeric(doc[i,13])
  
  x_golgi <- as.numeric(doc[i,14])
  y_golgi <- as.numeric(doc[i,15])
  
  x1 <- x_end - x_start #(x1,y1): the caculated migration vector
  y1 <- y_end - y_start
  
  x2 <- x_golgi - x_nuclear #(x2,y2): the calculated the centroid of nuclear to Golgi apparatus
  y2 <- y_golgi - y_nuclear
  dot = x1*x2 + y1*y2      # dot product
  det = x1*y2 - y1*x2      # determinant
  angle = atan2(det, dot)*180/pi  # atan2(y, x) or atan2(sin, cos)
  output[i,5] <- angle
  
  if (angle < 0) {
    output[i,6] <- 360-abs(angle)
  }
  if (angle >= 0) {
    output[i,6] <- angle
  }
  if (abs(angle) <= 60) {
    output[i,7] <- 'Polarity'
  }
  if (abs(angle) > 60) {
    output[i,7] <- 'Unpolarity'
  }
  end
}

colnames(output)[2:7] <- c('Leading_1','Leading_2','L_Polarity?','Trailing_1','Trailing_2','T_Polarity?')

write.csv(output,file = "angels_total_results.csv")
save(output,file = "output.Rdata")




