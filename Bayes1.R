#authors: Bruno Jerkovic, Josip Koprcina, Christos Mavrikis

#If you are running this script for the first time
#please uncomment and install the packages below

#install.packages(c("remotes","pROC","naivebayes"))
#remotes::install_github("jtextor/bayesianNetworks")
#install.packages("lavaan", repos="http://cran.at.r-project.org/")
#install.packages("dummies")
#install.packages("bnlearn")
#install.packages("dagitty")
#install.packages("Rgraphviz",dependenacies = TRUE)
#install.packages("pcalg",dependenacies = TRUE)
#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install(c("Rgraphviz", "RBGL"))
#include used libs for assignment
#BiocManager::valid()   
#library(dummies)
#library(tibble)
#library(naivebayes)
#install.packages("corrplot")
library(corrplot)
library(lavaan)
library(bnlearn)
library(dagitty)
library(bayesianNetworks)
library(pROC)
library(BiocManager)
library(pcalg)

#Error when trying to use bayesian network ci.test library so i imported it manually
ci.test2 <- function( x, y, z=NULL, data ) {
  require(ranger)
  
  if( !is.data.frame(data) ){
    stop("Please supply your data in a data frame!")
  }
  
  if( length(z) == 0 ){
    z <- NULL
  }	
  
  n <- nrow(data)
  
  if( is.ordered( data[,x] ) ){
    data[,x] <- as.integer(data[,x])
    tx <- matrix( data[,x],ncol=1 )
  } else if( is.numeric( data[,x] ) ){
    data[,x] <- rank( data[,x] )
    tx <- matrix( data[,x],ncol=1 )
  } else if( is.factor( data[,x] ) ){ 
    tx <- table( seq_len(n), data[,x] )
  } else {
    stop( "unsupported data type: x ")
  }
  
  if( is.ordered( data[,y] ) ){
    data[,y] <- as.integer(data[,y])
    ty <- matrix( data[,y], ncol=1 )
  } else if( is.numeric( data[,y] ) ){
    data[,y] <- rank( data[,y] )
    ty <- matrix( data[,y], ncol=1 )
  } else if( is.factor( data[,y] ) ){ 
    ty <- table( seq_len(n), data[,y] )
  } else {
    stop( "unsupported data type: y ")
  }
  
  # residualize
  if( is.null(z) ){
    tx <- sweep( tx, 2, colMeans(tx) )
    ty <- sweep( ty, 2, colMeans(ty) )
  } else {
    forest.x <- ranger::ranger( x=data[,z,drop=FALSE], y=data[,x],
                                probability=ncol(tx)>1 )
    tx <- tx - predict( forest.x, data=data )$predictions
    
    forest.y <- ranger::ranger( x=data[,z,drop=FALSE], y=data[,y],
                                probability=ncol(ty)>1 )
    ty <- ty - predict( forest.y, data=data )$predictions
  }
  
  
  if( is.factor(data[,x]) ){
    tx <- tx[,-ncol(tx),drop=FALSE]
  }
  if( is.factor(data[,y]) ){
    ty <- ty[,-ncol(ty),drop=FALSE]
  }
  
  k <- ncol(tx)
  r <- ncol(ty)	
  
  m <- 1
  R <- matrix( 0, nrow=n, ncol=(k)*(r) )
  
  # this takes long and can probably be done better
  keep.col <- rep( TRUE, ncol(R) )
  for( i in 1:(k) ){
    for( j in 1:(r) ){
      keep.col[m] <- (k==1 || i<k) && (r==1 || j<r)
      R[,m] <- (tx[,i]*ty[,j]) / sqrt(mean(tx[,i]^2)) / sqrt(mean(ty[,j]^2))
      m <- m+1
    }
  }
  
  cm <- colMeans(R)
  effect <- abs(cm)[which.max( cm^2 )]
  
  R <- R[,keep.col,drop=FALSE]
  p <- ncol(R)
  x2c <- (n-p) / p / (n-1) * n * t(colMeans(R)) %*% MASS::ginv( cov(R) ) %*% colMeans(R)
  list( statistic=c(x2c), effect=effect, parameter=max(k-1,1)*max(r-1,1),
        p.value=pchisq(c(x2c), max(k-1,1)*max(r-1,1),lower.tail=FALSE) )
}
#pass our dataset into a variable 
obese <- read.csv("ObesityDataSet_raw_and_data_sinthetic.csv")
#print some instances to be sure everything is OK
head(obese[1:10,])

obese<-obese[sample(nrow(obese)),]


#DAG - plot BN - Initial assumption

g <- dagitty('
dag {
bb="0,0,1,1"
Age [pos="0.050,0.037"]
CAEC [pos="0.258,0.368"]
CALC [pos="0.948,0.314"]
CH2O [pos="0.250,0.277"]
FAF [pos="0.791,0.234"]
FAVC [pos="0.394,0.364"]
FCVC [pos="0.502,0.348"]
Gender [pos="0.249,0.037"]
Height [pos="0.140,0.685"]
MTRANS   [pos="0.653,0.475"]
NCP [pos="0.381,0.523"]
NObeyesdad [pos="0.258,0.935"]
SCC [pos="0.817,0.119"]
SMOKE [pos="0.027,0.326"]
TUE [pos="0.121,0.277"]
Weight [pos="0.511,0.726"]
family_history_with_overweight [pos="0.547,0.038"]
Age -> CALC
Age -> CAEC
Age -> MTRANS
Age -> FAF
Age -> FAVC
Age -> FCVC
Age -> Height
Age -> SCC
Age -> SMOKE
Age -> TUE
Age -> Weight
Age -> CH2O
CAEC -> Weight
CALC-> Weight
CH2O -> Weight
FAF -> CH2O
FAF -> FCVC
FAF -> NCP
FAF -> Weight
FAVC -> Weight
FCVC -> Weight
Gender -> FAF
Gender -> SCC
Gender -> SMOKE
Gender -> TUE
Gender -> Weight
Gender -> Height
Height -> NObeyesdad
MTRANS -> Weight
NCP -> FAVC
NCP -> Weight
SCC -> CH2O
SCC -> FAF
SCC -> FAVC
SCC -> FCVC
SCC -> NCP
SCC -> Weight
SCC -> CALC
SCC -> SMOKE
SMOKE -> Weight
TUE -> Weight
Weight -> NObeyesdad
family_history_with_overweight -> CAEC
family_history_with_overweight -> FAVC
family_history_with_overweight -> FCVC
family_history_with_overweight -> NCP
family_history_with_overweight -> SCC
family_history_with_overweight -> Weight
}
')

plot(((g)))

#If you want a clearer or larger picture of the plot
#you can export the plot using the Export button

#Start of preprocessing 
obese$Age    <- round(obese$Age)
obese$Height <- round(obese$Height,digits=2)
obese$Weight <- round(obese$Weight)
obese$FCVC   <- round(obese$FCVC)
obese$NCP    <- round(obese$NCP)
obese$FAF    <- round(obese$FAF)
obese$TUE    <- round(obese$TUE)
obese$CH2O   <- round(obese$CH2O)
table(obese$Weight)
#Dividing into Age groups
s <-rep("<20",nrow(obese))
s[obese$Age>=20 &obese$Age<35] <- "20-34"
s[obese$Age>=35 &obese$Age<50] <- "35-49"
s[obese$Age>=50] <- "50>="
obese$Age <-as.numeric(ordered(s)) # <20 = 1 , 20-34 = 2 , 35-49 = 3 , 50>= = 4
#print column $Age

table(obese$CALC)
#Dividing into Height groups
s <-rep("<1.60",nrow(obese))
s[obese$Height>=1.60 &obese$Height<1.75] <- "1.60-1.74"
s[obese$Height>=1.75 &obese$Height<1.85] <- "1.75-1.84"
s[obese$Height>=1.85] <- "1.85>="
obese$Height <-as.numeric(ordered(s)) # <1.60 = 1 , 1.60-1.74 = 2 , 1.75-1.84 = 3 , 1.85>= = 4
#print column $Height
#obese$Height

#Dividing into Weight groups
s <-rep("<50",nrow(obese))
s[obese$Weight>=50 &obese$Weight<70] <- "50-69"
s[obese$Weight>=70 &obese$Weight<90] <- "70-89"
s[obese$Weight>=90] <- "90>="
obese$Weight <-as.numeric(ordered(s))  # <50 = 1 , 50-69 = 2 , 70-89 = 3 ,90>= = 4
#print column $Weight
#obese$Weight

#Breakdown Weight levels as well
obese$NObeyesdad <- as.character(obese$NObeyesdad)
obese$NObeyesdad[obese$NObeyesdad %in% c("Insufficient_Weight","Normal_Weight")] <- "Normal_Weight"
obese$NObeyesdad[obese$NObeyesdad %in% c("Overweight_Level_I","Overweight_Level_II")] <- "Overweight"
obese$NObeyesdad[obese$NObeyesdad %in% c("Obesity_Type_I","Obesity_Type_II","Obesity_Type_III")] <- "Obese"
obese$NObeyesdad <- as.numeric(ordered(obese$NObeyesdad))#Normal_Weight = 1 Overweight = 2 Obese= 3
#print column $NObeyesdad
#obese$NObeyesdad


#as.numeric ordered
obese$family_history_with_overweight <- as.numeric(ordered(obese$family_history_with_overweight,c("no","yes"))) #no = 1 &  yes = 2
obese$FAVC                           <- as.numeric(ordered(obese$FAVC,c("no","yes"))) #no = 1 &  yes = 2
obese$SCC                            <- as.numeric(ordered(obese$SCC,c("no","yes"))) #no = 1 &  yes = 2
obese$SMOKE                          <- as.numeric(ordered(obese$SMOKE,c("no","yes"))) #no = 1 &  yes = 2
obese$Gender                         <- as.numeric(ordered(obese$Gender,c("Female","Male"))) #Female = 1 &  Male = 2
obese$CAEC                           <- as.numeric(ordered(obese$CAEC, levels=c("no","Sometimes","Frequently","Always"))) #no = 1 Sometimes = 2 Frequently = 3 Always = 4
obese$CALC                           <- as.numeric(ordered(obese$CALC, levels=c("no","Sometimes","Frequently","Always"))) #no = 1 Sometimes = 2 Frequently = 3 Always = 4
obese$MTRANS                         <- as.numeric(ordered(obese$MTRANS, c("Walking","Bike","Public_Transportation","Motorbike","Automobile")))

#End of preprocessing

#Fitting the Model With Obese dataset and checking the Coefficients

net <- model2network(toString(g,"bnlearn"))
net

#--------- manually check at coeffiecents this has been analyzed in the cover letter----
#fit1 <- bn.fit(net,obese,method="mle") #change
#summary(fit1)
#table(obese$NCP)
#coefficients(fit1)
##lm( Height ~ Age, data=obese )

#Any coefficient with a score ranging between -0.1 to +0.1 have been manually removed
g_cutoff <- dagitty('
dag {
bb="0,0,1,1"
Age [pos="0.050,0.037"]
CAEC [pos="0.258,0.368"]
CALC [pos="0.948,0.314"]
CH2O [pos="0.250,0.277"]
FAF [pos="0.791,0.234"]
FAVC [pos="0.394,0.364"]
FCVC [pos="0.502,0.348"]
Gender [pos="0.249,0.037"]
Height [pos="0.140,0.685"]
MTRANS   [pos="0.653,0.475"]
NCP [pos="0.381,0.523"]
NObeyesdad [pos="0.258,0.935"]
SCC [pos="0.817,0.119"]
SMOKE [pos="0.027,0.326"]
TUE [pos="0.121,0.277"]
Weight [pos="0.511,0.726"]
family_history_with_overweight [pos="0.547,0.038"]
Age -> MTRANS
Age -> TUE
Age -> Weight
CAEC -> Weight
CALC-> Weight
FAF -> NCP
FAVC -> Weight
FCVC -> Weight
Gender -> FAF
Gender -> Weight
Gender -> Height
Height -> NObeyesdad
SCC -> FAVC
Weight -> NObeyesdad
family_history_with_overweight -> CAEC
family_history_with_overweight -> FAVC
family_history_with_overweight -> Weight
}
')

plot(g_cutoff)

#Check for effect
#for(x in names(g_cutoff)){
 # px <- parents(g_cutoff,x)
  #for(y in px){
   # tst <- ci.test2(x,y,setdiff(px,y),data=obese)
    #if(tst$effect > 0.1)
     # print(paste(y,'->',x,'effect:',tst$effect,'p.value',tst$p.value))
  #}
#}


#PREDICTIONS - Perform 10 fold cross validation - random sampling

#Create 10 equally size folds
folds <- cut(seq(1,nrow(obese)),breaks=10,labels=FALSE)
net <- model2network(toString(g_cutoff,"bnlearn"))
net

#initialize matrix for each fold auc score
mat<-matrix(list(), nrow=10, ncol=2)
colnames(mat) <-c("Train","Test")
for (k in 1:10) {
mat[k,1] <- 0
mat[k,2] <- 0
}
#print(mat)
localTests(g, obese, type ="cis.chisq", max.conditioning.variables = 2)[,4,drop=FALSE]
localTests(g_cutoff, obese, type ="cis.chisq", max.conditioning.variables = 2)[,4,drop=FALSE]
for (k in 1:10) {
  #folds and fitting using BNLEARN
  testIndexes <- which(folds==k,arr.ind=TRUE)
  train.set <- obese[-testIndexes, ]
  test.set <- obese[testIndexes, ]
  fit1 <- bn.fit(net,train.set,method="mle") #change
  fit1
  #---------- Train -------------------
  p <-  predict(fit1, node="NObeyesdad", data=train.set[c("Height","Weight")], method="bayes-lw")
  mat[k,1] <- auc(train.set$NObeyesdad,p)
  roc1<-roc(train.set$NObeyesdad,p) 

  #---------- Test -------------------
  p1 <-  predict(fit1, node="NObeyesdad",data=test.set[c("Height","Weight")], method="bayes-lw")
  mat[k,2] <- auc(test.set$NObeyesdad,p1)
  roc2<-roc(test.set$NObeyesdad,p1)

  title_k <- paste(c("Plotting fold", k))

  plot(roc1,col='red',main=title_k) #train
  legend("topleft",c("test","train"),fill=c("blue","red"))
  plot(roc2,add=TRUE, col='blue') #test
}

print("Results per fold")
print(mat)
average_of_test <- matrix(unlist(mat), ncol=length(mat) )
sprintf("Average of Test: %f" ,mean(average_of_test[11:20]))

