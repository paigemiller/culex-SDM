#go to extracted file on local directory
#setwd("~/Desktop/ghcn-data/ghcnm.v3.3.0.20160216")

#read table, weird format, blank spaces give it problems, fill=TRUE to just let blank spaces be a thing
ghcn.vars <- c() 
ghcn <- read.table("ghcnm.tavg.v3.3.0.20160216.qca.dat", na.strings = "-9999", fill=TRUE, stringsAsFactors=FALSE, header=FALSE)

# separate v1
col1 <- ghcn$V1 
ghcn.col1 <- data.frame(matrix(NA, nrow = length(ghcn[,1]), ncol = 3))
colnames(ghcn.col1) <- c("ID", "year", "var")

for (i in 1:length(col1)){
  split.i <- unlist(strsplit(col1[i], split=""))
  ghcn.col1$ID[i] <- paste(split.i[1:11], sep="", collapse="")
  ghcn.col1$year[i] <- paste(split.i[12:15], sep="", collapse="")
  ghcn.col1$var[i] <- paste(split.i[16:19], sep="", collapse="")
}

#function that's faster --working on it
split.col1 <- function(ID){
  ID <- unlist(strsplit(ID, split=""))
  
}

# remove col1
ghcn <- ghcn[,2:25]

# join ghcn.col with ghcn
require(dplyr)
GHCN <- bind_cols(ghcn.col1, ghcn)

# see README for vars



