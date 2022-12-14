# Emilie Lebarbier (Paris) and Meili Baragatti (Montpellier), 27/06/2017

library(wavethresh)
if (!require('fda')) install.packages('fda'); library('fda')

#
debug = FALSE


# Creation of a dictionary with the Haar functions
# arguments J and R: we want the Haar functions with resolutions from 2^R to 2^J
Dico <- function(J,R){
  N=2^(J+1)
  HaarMat=GenW(N,filter.number=1,family="DaubExPhase")
  seq=seq(J,R,-1)
  Nb.bases.j=c(2^seq)
  sum=sum(Nb.bases.j)
  Fh=HaarMat[,2:(sum[length(sum)]+1)]
  return(Dico.Funct=Fh)
}

# Evaluation of the functions in Dic.Funct at positions posx
Eval.posx<- function(Dic.Funct,posx,N){
  position.Funct=seq(1/N,1,1/N)  
  M=dim(Dic.Funct)[2]
  FhNew=matrix(0,ncol=M,nrow=length(posx))
  for (i in 1:length(posx)){
    FhNew[i,]=Dic.Funct[which(position.Funct>=posx[i])[1],]
  }
  return(Fhx=FhNew)
}

# Remove from the dictionary the bases which are null at all evaluated positions
Bases.NonNulles<- function(F.bases){
  nbases=dim(F.bases)[2]
  PresentBases=c(1:nbases)
  sumcolonnes = apply(F.bases,2,sum)
  toremove = which(sumcolonnes == 0)
  
  if (length(toremove) !=0){
    F.bases.Present=F.bases[,-toremove]
    PresentBases = PresentBases[-toremove]
  } else {
    F.bases.Present=F.bases
    PresentBases = PresentBases
  }
  list(Fh.P=F.bases.Present, PresentBases=PresentBases)
}

########################################
# Dico with cste, Fourier, Haar, x^{1/3}, x^{1/2}, 
# x et x^2, log(x), splines ordres 3: example truffes
########################################

CreationBases<-function(t,J,R,invperiode){
  n=length(t) 
  posx=t/max(t)
  
  ## Fonction cste
  cste <- rep(1,n)
  
  ## Fonctions fourier
  namesfourier <- c()
  Ff      = c()
  for (i in 1:invperiode){
    Ff = cbind(Ff,sin(2*pi*i*posx),cos(2*pi*i*posx))
    namesfourier <- c(namesfourier,paste("sin(2pi*",i,"posx)",sep=""),paste("cos(2pi",i,"posx)",sep=""))
  }
  
  ## Fonctions Haar
  Dic.Funct=Dico(J,R)
  Dic.Funct=abs(Dic.Funct)
  Fh=Eval.posx(Dic.Funct,posx,2^(J+1))
  namesHaar <- c()
  for (i in 1:dim(Fh)[2]){
    namesHaar <- c(namesHaar,paste("Haar",floor(i/2^J*1000)/1000,sep="")) 
  }
  
  ## x^{1/3}, x^{1/2}, x et x^2
  x13 <- posx^{1/3}
  x12 <- sqrt(posx)
  x <- posx
  xx <- posx*posx
  ## log(x)
  logx <- log(posx+1)
  
  ## Splines
  numKnots = 13
  knots <- quantile(unique(t),seq(0,1,length=(numKnots+2))[-c(1,(numKnots+2))])
  B2  = create.bspline.basis(c(min(t),max(t)),norder=2,breaks=knots)
  Fs2 = eval.basis(t,B2)
  namesFs2 <- B2$names
  B3  = create.bspline.basis(c(min(t),max(t)),norder=3,breaks=knots)
  Fs3 = eval.basis(t,B3)
  namesFs3 <- B3$names
  B4  = create.bspline.basis(c(min(t),max(t)),norder=4,breaks=knots)
  Fs4 = eval.basis(t,B4)
  namesFs4 <- B4$names
  
  ## Creation matrice F
  # F.Bases=cbind(cste,Fh,Ff,Fs3,x13,x12,x,xx,logx)
  F.Bases=cbind(cste,Fh,Ff,Fs3,x,xx)
  colnames(F.Bases) = c("cste",namesHaar,namesfourier,namesFs3,"x","xx")
  #F.Bases=cbind(cste,Ff,Fh,Fs2,Fs3,Fs4,x13,x12,x,xx,logx)
  #colnames(F.Bases) = c("cste",namesfourier,namesHaar,namesFs2,namesFs3,namesFs4,"x13","x12","x","xx","logx")
  
  ## On enleve les colonnes correspondant aux fonctions de base nulles partout
  temp <- Bases.NonNulles(F.Bases)
  F.Bases <- temp[[1]]
  FunctionsNonNulles <- temp[[2]]
  
  list(F.Bases=F.Bases,FunctionsNonNulles=FunctionsNonNulles)
}

main <- function() {
  
  filename <- paste0(getwd(), "/inputs/processed_o3_data.csv")
  if(debug)  print(filename)
  
  data <- read.csv(filename)
  vec <- c()
  cnames <- colnames(data)
  for (i in 1:3) {
    vec<-data[,i]
    
    t<-c(1:length(vec))
    t<-t[-which(is.na(vec))]
    
    #Fmatrix_0 <- CreationBases(t,7,7,10)[[1]]
    Fmatrix <- CreationBases(t,7,7,10)[[1]]
    #Fmatrix <- Fmatrix_0[-which(is.na(vec)),]
    
    output <- paste0(getwd(),"/inputs/Fmatrix",cnames[i],".csv")
    if(debug) print("Output file:")
    if(debug) print(output)
    write.csv(Fmatrix, output) ## 
  }
  print(" Success!! ")
}

main()
