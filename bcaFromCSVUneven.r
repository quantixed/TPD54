library(bootstrap)
library(beeswarm)
# import data
myrawdata <- read.csv(file='data.txt', header=TRUE, stringsAsFactors=FALSE)
# find the mean of the control values (assume it is in 1st column)
ctrlmean <- mean(myrawdata[,1], na.rm=TRUE)
# subtract mean from the data
mydiffdata <- myrawdata - ctrlmean
# uneven number of rows mean NaNs, deal with them
mydata <- lapply(mydiffdata, function(col)col[!is.na(col)])
# don't need to work on the first column
cols <- length(mydata) - 1
# make empty vector to take the mean
mymeans <- rep(NA, cols)
mybcalo <- rep(NA, cols)
mybcahi <- rep(NA, cols)
# use a for loop to get for each mean diff, low and high BCa CI 
for (i in 1:cols ){
  mymeans[i] <- mean(mydata[[i+1]])
  bca = bcanon(mydata[[i+1]],10000,mean,alpha=0.025)
  mybcalo[i] <- bca$conf[1,2]
  bca = bcanon(mydata[[i+1]],10000,mean,alpha=0.975)
  mybcahi[i] <- bca$conf[1,2]
}
#beeswarm(myrawdata,pch=16,method="swarm")
#write.csv(mydata, file = "mydata.csv",row.names=FALSE)
# make dataframe of the result
result <- data.frame(mymeans,mybcalo,mybcahi)
# save to wd
write.csv(result, file = "result.csv",row.names=FALSE)
