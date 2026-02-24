
rm(list=ls())

# Read document for the vector x,y coordinate ,
doc <- data.table::fread('Testing-stackreg in µm per min.csv ', header = T)

#Input the cell migration vector
x_start <- 0
y_start <- 0

x_end <- 2
y_end <- 0

#Input the tn and tn-1 vector 

output <- data.frame(number = doc[,"Slice"])
for (i in 1: nrow(doc)) {
 
  print(i)
  
  x_tn_1 <- as.numeric(doc[i,3])
  y_tn_1 <- as.numeric(doc[i,4])

  x_tn <- as.numeric(doc[i+1,3])
  y_tn <- as.numeric(doc[i+1,4])
  
  x1 <- x_end - x_start #(x1,y1): the caculated migration vector
  y1 <- y_end - y_start
  
  x2 <- x_tn - x_tn_1 #(x2,y2): the calculated the tn-1 to tn
  y2 <- y_tn - y_tn_1
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
  if (abs(angle)<=90){
    output[i,5] = 90-abs(angle)
  }
  else
    output[i,5] = -(abs(angle)-90)
  end
}


colnames(output)[2:5] <- c('Original_angle','Modified_angel','T_Polarity?','Adjusted_angel')

write.csv(output,file = "Testing-angels_total_results_90 as polarity.csv")
save(output,file = "output_90aspolarity.Rdata")





