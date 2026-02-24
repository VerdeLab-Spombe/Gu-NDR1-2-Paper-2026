rm(list = ls())
library(sojourner)

# if(!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("sojourner")
# Specify folder with data
folder1=system.file("extdata","Test_1",package="sojourner")
folder2=system.file("extdata","Test_2",package="sojourner")

# Create track list
trackll1<-createTrackll(folder=folder1, interact=FALSE, input=3, ab.track=FALSE, cores=1, frameRecord=TRUE)
trackll2<-createTrackll(folder=folder2, interact=FALSE, input=3, ab.track=FALSE, cores=1, frameRecord=TRUE)


# Take a look at the list by 
str(trackll1,1)
# Filter/choose tracks 10 frames or longer for all analysis
trackll.fi1<-filterTrack(trackll=trackll1, filter=c(min=10,max=Inf))
trackll.fi2<-filterTrack(trackll=trackll2, filter=c(min=10,max=Inf))

# Merge tracks from different image files in the folder
trackll.fi.me1=c(mergeTracks(folder1, trackll.fi1))
trackll.fi.me2=c(mergeTracks(folder2, trackll.fi2))


str(trackll.fi.me1,1)
# calculate MSD for all tracks longer than 10 frames
#resolution is the  ratio of pixel to uM, our is 1.3043um/pixel



data_individual1 <- msd(trackll.fi.me1,dt=45,resolution=1.3043,summarize=F,cores=4,plot=TRUE,output=TRUE)
data_summarized1 <- msd(trackll.fi.me1,dt=45,resolution=1.3043,summarize=TRUE,cores=4,plot=TRUE,output=TRUE)

data_individual2 <- msd(trackll.fi.me2,dt=45,resolution=1.3043,summarize=F,cores=4,plot=TRUE,output=TRUE)
data_summarized2 <- msd(trackll.fi.me2,dt=45,resolution=1.3043,summarize=TRUE,cores=4,plot=TRUE,output=TRUE)



head(data_individual1[[1]])

individual1 <- data_individual1[[1]]
summarized1 <- data_summarized1[[1]]

individual2 <- data_individual2[[1]]
summarized2 <- data_summarized2[[1]]



write.csv(individual1,file = "Test_1_MSD_individual.csv")
write.csv(summarized1,file = "Test_1_MSD_summarized.csv")

write.csv(individual2,file = "Test_2_MSD_individual.csv")
write.csv(summarized2,file = "Test_2_MSD_summarized.csv")

