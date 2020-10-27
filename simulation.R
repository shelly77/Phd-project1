
  #load TMB template 
  require('TMB')
  require("mvtnorm")
  compile("great_tit_Final.cpp")
  dyn.load(dynlib("great_tit_Final"))

  source("functions.R")
  source("great_tit_data_pre.R")

####################################
  ##compare two models' AIC valuesm provided that both of the models are converged
  powtest <- function(listA,listB){
    if (listA[[1]] !=0 || listB[[1]] !=0)
      return(NA)
    else ifelse(listB[[3]]+2 < listA[[3]],1,0) #the third element of the list is AIC
  }


  ## create empty matrix and list to store the simulation restults
  powFlu <- matrix(NA,nrow=nrep,ncol=length(autolist))
  powZip <- matrix(NA,nrow=nrep,ncol=length(autolist))
  powPhi <- matrix(NA,nrow=nrep,ncol=length(autolist))
  powAllPhi <- matrix(NA,nrow=nrep,ncol=length(autolist))
  powFem <- matrix(NA,nrow=nrep,ncol=length(autolist))

  bestfit <- list()
  
  nrep=100 # number of simulations
  autolist <- c(0,0.1,0.25,0.5,0.75,0.9) # vector of auto-correlation values

for (iauto in autolist ){
  for (irep in 1:nrep){

    auto <- autolist[iauto]
    Phi=matrix(c(auto,0,0,0,auto,0,0,0,0), nrow=3)
    
    #simulate data
    n <- 100
    tmax <- 50
    mu_alpha <- 1.8  #mean of eta_alpha
    mu_theta <- 18 #mean of eta_theta
    mu_omega <- log(45) #mean of eta_omega (the standardized strength of stabilizing selection is 0.267)
    sigma1 <- 0.2 # standard deviaiton of errors of alpha 
    sigma2 <- 18  # standard deviaiton of errors of theta
    sigma3 <- 0.2 # standard deviaiton of errors of omega
    sigma <- c(sigma1,sigma2,sigma3)
    
    #to simulate the VAR1 process at different time points
    #we first calculate the covariance matrix of alpha, theta and omega (Gamma0)
    # when the entry in sigma is 0, it indicates that there is no random effects of alpha, theta or omega
    Sigma <- diag(c(sigma[1]^2,sigma[2]^2,sigma[3]^2))# Sigma is covariance matrix of the white noises for VAR1 alpha, theta and omega
    
    d=3
    vecSigma <- c(rep(NA,d*d))
    for (i in 1:d)
      for (j in 1:d){
        vecSigma[i+(j-1)*d] = Sigma[i,j]}
    
    A <- matrix(NA,d*d,d*d)
    for (i in 1:d)
      for (j in 1:d)
        for (k in 1:d)
          for (l in 1:d)
            A[(i-1)*d+k, (j-1)*d+l] = -Phi[i,j]*Phi[k,l]
    
    A <- diag(1,d*d)+A
    Ainv = solve(A)
    #Gamma0 is the stationary covariance matrix for VAR1 alpha, theta, omega
    vecGamma0 <- Ainv%*%vecSigma
    Gamma0 <- matrix(NA,d,d)
    for (i in 1:d)
      for (j in 1:d){
        Gamma0[i,j] = vecGamma0[i+(j-1)*d]}
    ##simulate laying dates, number of fledglings and year t
    #set.seed(irep)
    t <- rep(1:tmax,rpois(tmax,n))   # simulate on average mean_n individuals every gen.
    tfactor <- factor(t)
    ntotal <- length(t)   # total number of individual observations 
    #set.seed(irep)
    components <- sample(1:2,prob=c(0.7,0.3),size=ntotal,replace=TRUE)
    mus <- c(23,62)
    sds <-c(7.5,10.5)
    z <- floor(rnorm(n=ntotal,mean=mus[components],sd=sds[components]))

    varPro <- matrix(NA,tmax,3)
    # simulate VAR1 process alpha, theta and omega
    #set.seed(irep)
    varPro[1,] <- rmvnorm(1,mean=c(0,0,0),sigma=Gamma0)
    for (j in 2:tmax) varPro[j,] <- t(Phi%*%varPro[j-1,])+rmvnorm(1,mean=c(0,0,0),sigma=Sigma)
    
    eta_alpha = mu_alpha +  varPro[,1]
    eta_theta = mu_theta +  varPro[,2]
    eta_omega = mu_omega +  varPro[,3]
    
    #simulate random mother effects, take true mothers ID
    motherID <- mytitbroodFinal$MothersID[1:ntotal]
    epsi <- rnorm(length(unique(motherID)),0,sd=1)
    idNum <- as.numeric(factor(motherID))

    logw <- vector()
    motherEff <- vector()
    for(i in 1:ntotal){
      #add random mother effects to the log expected fitness
      motherEff[i]<- epsi[idNum[i]]
      #the sd of mother effects is set to 0.05
      logw[i] <- eta_alpha[t[i]]-(z[i]-eta_theta[t[i]])^2/(2*(exp(eta_omega[t[i]]))^2)+0.05*motherEff[i]
      }
    w <- exp(logw) # expected fitness 
    
    #zero-inflated probability 
    beta0_p <- 2
    prob_fun <- beta0_p 
    prob <- 1/(1+exp(prob_fun))
    
    #set.seed(irep)
    y<-vector()
    for (i in 1:ntotal){
      y[i] <- rzipois(1,w[i],prob)}
  
    alldata <- data.frame(y=y,t=t,z=z,tfactor=tfactor,motherID=motherID)
    
    dataTMB = 
      list(
        LayDateApril  = alldata$z, # lay date for each brood
        NumFle = alldata$y, # number of fledglings for each brood
        t = as.factor(alldata$t),# index of year from 1955 to 2015
        j = as.factor(alldata$motherID)
      )
    
  
    #############################################
    # Initialize parameters
      parameters <- list(
      alpha        = rep(0,nlevels(dataTMB$t)),
      theta        = rep(0,nlevels(dataTMB$t)),
      omega        = rep(0,nlevels(dataTMB$t)),
      muAlpha      = rep(3,2),
      muTheta      = rep(10,2),
      muOmega      = rep(3,2),
      logSigma     = c(-2,1,-1),
      parCor       = rep(0,3),
      Phi          = matrix(rep(0,9),nrow=3),
      beta0_P      = 0,
      beta1_P      = 0,
      u0_P         = rep(1,nlevels(dataTMB$t)),
      u1_P         = rep(0.1,nlevels(dataTMB$t)),
      beta0_W      = 0,
      beta1_W      = 0,
      u0_W         = rep(1,nlevels(dataTMB$t)),
      u1_W         = rep(1,nlevels(dataTMB$t)),
      logSdRanP    = c(0,-3),
      thetaRanP    = -0.5,
      logSdRanW    = rep(-2,2),
      thetaRanW    = 0,
      logTauMother = rep(-2,2),
      epsi         = rep(0.1,nlevels(dataTMB$j))
    )

    
    #############################################
    map0 = list(
      alpha          = factor(rep(NA,nlevels(dataTMB$t))),
      theta          = factor(rep(NA,nlevels(dataTMB$t))),
      omega          = factor(rep(NA,nlevels(dataTMB$t))),
      muAlpha        = factor(rep(NA,2)),
      muTheta        = factor(rep(NA,2)),
      muOmega        = factor(rep(NA,2)),
      logSigma       = factor(rep(NA,3)),
      parCor         = factor(rep(NA,3)),
      Phi            = factor(matrix(rep(NA,9),nrow=3)), 
      beta0_P        = factor(NA),
      beta1_P        = factor(NA),
      u0_P           = factor(rep(NA,nlevels(dataTMB$t))),
      u1_P           = factor(rep(NA,nlevels(dataTMB$t))),
      beta0_W        = factor(NA),
      beta1_W        = factor(NA),
      u0_W           = factor(rep(NA,nlevels(dataTMB$t))),
      u1_W           = factor(rep(NA,nlevels(dataTMB$t))),
      logSdRanP      = factor(rep(NA,2)),
      thetaRanP      = factor(NA),
      logSdRanW      = factor(rep(NA,2)),
      thetaRanW      = factor(NA),
      logTauMother   = factor(rep(NA,2)),
      epsi           = factor(rep(NA,nlevels(dataTMB$j)))
    )
    

    ########################################################
    #model selection
    #######################################################

    model0 <- list(parameters=parameters, map=map0, 
                    args.MakeADFun=list(data = dataTMB, DLL = "great_tit_Final", 
                                        random = c("alpha","theta","omega","epsi",
                                                   "u0_P","u1_P","u0_W", "u1_W"),
                                        silent = FALSE))
    class(model0) <- "tmbmodel" # null model
    
    
    model1 <- update(model0,optimizer="nlminb",list(muAlpha = factor(c(NA,1)),
                                                    muTheta = factor(c(NA,1)),
                                                    muOmega = factor(c(NA,1)),
                                                    alpha=NULL,
                                                    theta=NULL,
                                                    omega=NULL,
                                                    logSigma=factor(c(1,NA,NA)),
                                                    beta0_P = NULL,
                                                    epsi=NULL,
                                                    logTauMother=factor(c(NA,1))),
                     cv=0.2,seed=123)

    list1 <- listmodel(model1)
    
    model2 <- update(model1,optimizer="nlminb",list(logSigma=factor(c(1,2,3))))

    list2 <- listmodel(model2)
    
    model3 <- update(model2,optimizer="nlminb",list(Phi=factor(matrix(c(1,NA,NA,NA,2,NA,NA,NA,NA),nrow=3))))
    
    list3 <- listmodel(model3)
    
    model4 <- update(model3,optimizer="nlminb",list( logSigma=factor(c(1,2,NA))))

    list4 <- listmodel(model4)

    model5 <- update(model3,optimizer="nlminb",list(beta1_P=NULL,u0_p=NULL,logSdRanP=factor(c(1,NA))))

    list5 <- listmodel(model5)

    model6 <- update(model3,optimizer="nlminb",list(Phi=factor(matrix(c(1,2,NA,3,4,NA,NA,NA,NA),nrow=3))))

    list6 <- listmodel(model6)

    model7 <- update(model3,optimizer="nlminb",list(epsi=factor(rep(NA,nlevels(dataTMB$j))),
                                                    logTauMother=factor(c(NA,NA))))

    list7 <- listmodel(model7)

    #distinguish properly models with and without fluctuating omega
    powFlu[irep,iauto] <- powtest(list4,list3)

    #identify zero-inflated probability as a parameter over a selection episode
    powZip[irep,iauto] <- powtest(list5,list3)

    #identify the transition parameters in $\Phi$.
    powPhi[irep,iauto] <- powtest(list2,list3)
    
    # identify the correlation of alpha and theta
    powAllPhi[irep,iauto] <- powtest(list6,list3)
    
    #identify the random mother effects
    powFem[irep,iauto] <- powtest(list7,list3)

    if (model3$convergence==0){ dataset <- summary(model3$sdreport,"report")
    bestfit[[irep+(iauto-1)*nrep]] <- list(subset(dataset, rownames(dataset)
    == "PhiB"))} else bestfit[[irep+(iauto-1)*nrep]] <- NA
  }
}

 simPowerTest <- list(powFlu,powZip,powPhi,powAllPhi,powFem)
 save(simPowerTest,file="simPowerTest.Rdata")
 save(bestfit, file="simBestfit.Rdata")


# Figure S1
################################
#plot percent of the true model selected over alternative model
 par(mfrow=c(1,2))
 plot(NULL,NULL,xlim=c(0,0.95),ylim=c(0,1),xlab="Auto-correlation",ylab="Percent of true model selected over alternative model",
     cex.lab=1.2,xaxt='n')
 for (i in 1:length(autolist)){
  lines(autolist,apply(na.omit(powZip),2,mean),lty=1,xaxt='n',pch=14,type="b")
  lines(autolist,apply(na.omit(powFlu),2,mean),lty=2,xaxt='n',pch=15,type="b")
  lines(autolist,apply(na.omit(powPhi),2,mean),lty=3,xaxt='n',pch=16,type="b")
  lines(autolist,apply(na.omit(powAllPhi),2,mean),lty=4,xaxt='n',pch=17,type="b")
  lines(autolist,apply(na.omit(powFem),2,mean),lty=5,xaxt='n',pch=18,type="b")
 }
 axis(1,at=autolist)
 legend(0.18,0.15,
       title="Percent of true model selected over:", title.adj=0,
       legend = c("model with fluctuating zip", expression(paste("model with fixed ", omega)),
                  expression(paste("model without autocorrelation in ",alpha, " and ", theta)),
                  expression(paste("model with correlated ",alpha," and ", theta)),
                  "model without random mother effects"),  
       cex = 0.9, lty = c(1,2,3,4,5),pch=c(14,15,16,17,18),box.lty=0)

################################
 #plot estimated auto-correlation
 autoAlpha <- matrix(NA,nrow=nrep,ncol=length(autolist))
 autoTheta<- matrix(NA,nrow=nrep,ncol=length(autolist))
 for (iauto in 1:6){    
   for (rep in 1:nrep){ 
     autoAlpha[rep,iauto] <- unlist(bestfit[[rep+(iauto-1)*nrep]])[1]
     autoTheta[rep,iauto] <- unlist(bestfit[[rep+(iauto-1)*nrep]])[5]
   }
 }
 
 autoStandAlpha <- apply(autoAlpha,2,sd)
 autoStandTheta <- apply(autoTheta,2,sd)
 
 require("plotrix")
 plotCI(autolist,colMeans(autoAlpha), uiw=autoStandAlpha,xlab="", cex.lab=1.2,
        ylab="Estimated auto-correlation",ylim =c(-0.2,1),gap=0,col="red",xaxt='n',pch=24,pt.bg="red")
 lines(0:1,0:1,lty=2,lwd=2,col="grey")
 par(new=TRUE) 
 plotCI(autolist,colMeans(autoTheta), uiw=autoStandTheta,xlab="Auto-correlation",
        ylab="",ylim =c(-0.2,1),gap=0,col="blue",xaxt='n',pch=21,pt.bg="blue",cex.lab=1.2)
 axis(1,at=autolist)
 legend(
   x=0.5, # x coordinate of the top left of the legend
   y=0, # y coordinate of the top left of the legend
   box.lty=0, # line type to surround the legend box (0 for none)
   legend=c(expression(paste("auto-correlation of ", alpha)),
            expression(paste("auto-correlation of ", theta))), # sequence of text for the legend
   pch=c(24,21), # sequence of point types for the legend; -1 is a nonexistent point
   pt.bg=c("red","blue"), # sequence of fill colours for the points
   col=c("red","blue"),
   pt.cex=c(1.2,1.2),
   lty=c(1,1),
   cex=1.2# sequence of line types for the legend
 )




