
# Loading dataset and Dictionary 
data <- read.csv("~/Desktop/changepointdetection/inputs/processed_o3_data.csv")
Fmatrix <- read.csv("~/Desktop/changepointdetection/inputs/Fmatrix.csv")
FmatrixBA <- read.csv("~/Desktop/changepointdetection/inputs/FmatrixO3_BA.csv")
FmatrixComp <- read.csv("~/Desktop/changepointdetection/inputs/FmatrixO3_Comp.csv")
FmatrixERA <- read.csv("~/Desktop/changepointdetection/inputs/FmatrixO3_ERA.csv")
FmatrixFlora <- read.csv("~/Desktop/changepointdetection/inputs/FmatrixO3_Flora.csv")
FmatrixPance <- read.csv("~/Desktop/changepointdetection/inputs/FmatrixO3_Pance.csv")

dim(Fmatrix)

vec<-processed_o3_data[,"O3_BA"]

t<-c(1:length(vec))
t<-t[-which(is.na(vec))]

Fmatrix <- Fmatrix[-which(is.na(vec)),]


plot(Fmatrix[,"Haar0.71"],type="l")
lines(Fmatrix[,"Haar0.773"],col="black")
lines(Fmatrix[,"cos(2pi6posx)"],col="black")
lines(FmatrixBA[,"Haar0.71"],col="blue")
lines(FmatrixBA[,"Haar0.773"],col="blue")
lines(FmatrixBA[,"cos(2pi6posx)"],col="red")

plot(Fmatrix[,"cos.2pi6posx."],type="l")

lines(FmatrixComp[,"cos.2pi6posx."],col="blue")
lines(Fmatrix0_BA[,"cos.(2pi6posx)"],col="red")
colnames(FmatrixBA)
colnames(Fmatrix)[90]


colnames(Fmatrix0_BA)

colnames(Fmatrix)
colnames(FmatrixBA)



write.csv(FmatrixBA,"FmatrixO3_BA.csv",row.names = FALSE)


