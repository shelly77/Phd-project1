
 #fit the model defined by C++ template with our great tit data

 ################################################################
 #load TMB template 
 library(TMB)
 compile("great_tit_Final.cpp")
 dyn.load(dynlib("great_tit_Final"))

 #################################################################
 # load data
 source("great_tit_data_pre.R")

 #################################################################
 # Creat data list
 data = 
   list(
     LayDateApril  = as.numeric(as.character(mytitbroodFinal$LayDateApril)), # lay date for each brood
     NumFle = as.numeric(as.character(mytitbroodFinal$NumberFledged)), # number of fledglings for each brood
     t = as.factor(mytitbroodFinal$YearOfBreeding), # index of year from 1955 to 2015
     j = as.factor(mytitbroodFinal$MothersID), # mother identity
     makePred = 0, # 0 represent not making predictions
     predDates = rep((seq(0,100,by=10)),61), #only used for making predictiongs
     tPred = as.factor(rep(c(1:61),each=11)) #only used for making predictiongs
   )

 ##################################################################
 # Initialize parameters
  parameters <- list(
   alpha        = rep(1,nlevels(data$t)),
   theta        = rep(5,nlevels(data$t)),
   omega        = rep(1,nlevels(data$t)),
   muAlpha      = rep(3,2),
   muTheta      = rep(10,2),
   muOmega      = rep(3,2),
   logSigma     = c(-2,1,-1),
   parCor       = rep(0,3),
   Phi          = matrix(rep(0,9),nrow=3),
   beta0_P      = 0,
   beta1_P      = 0,
   u0_P         = rep(1,nlevels(data$t)),
   u1_P         = rep(0.1,nlevels(data$t)),
   beta0_W      = 0,
   beta1_W      = 0,
   u0_W         = rep(1,nlevels(data$t)),
   u1_W         = rep(1,nlevels(data$t)),
   logSdRanP    = rep(-2,2),
   thetaRanP    = 0,
   logSdRanW    = rep(-2,2),
   thetaRanW    = 0,
   logTauMother = rep(-2,2),
   epsi         = rep(0.1,nlevels(data$j))
 )


#######################################################
 map00 = list(
   alpha          = factor(rep(NA,nlevels(data$t))),
   theta          = factor(rep(NA,nlevels(data$t))),
   omega          = factor(rep(NA,nlevels(data$t))),
   muAlpha        = factor(rep(NA,2)),
   muTheta        = factor(rep(NA,2)),
   muOmega        = factor(rep(NA,2)),
   logSigma       = factor(rep(NA,3)),
   parCor         = factor(rep(NA,3)),
   Phi            = factor(matrix(rep(NA,9),nrow=3)), 
   beta0_P        = factor(NA),
   beta1_P        = factor(NA),
   u0_P           = factor(rep(NA,nlevels(data$t))),
   u1_P           = factor(rep(NA,nlevels(data$t))),
   beta0_W        = factor(NA),
   beta1_W        = factor(NA),
   u0_W           = factor(rep(NA,nlevels(data$t))),
   u1_W           = factor(rep(NA,nlevels(data$t))),
   logSdRanP      = factor(rep(NA,2)),
   thetaRanP      = factor(NA),
   logSdRanW      = factor(rep(NA,2)),
   thetaRanW      = factor(NA),
   logTauMother   = factor(rep(NA,2)),
   epsi           = factor(rep(NA,nlevels(data$j)))
 )
 
#########################################################
 #model selection
 library("numDeriv")
 source("functions.R")

 model00 <- list(parameters=parameters, map=map00, 
                args.MakeADFun=list(data = data, DLL = "great_tit_Final", 
                                    random = c("alpha","theta","omega","epsi",
                                               "u0_P","u1_P","u0_W",
                                               "u1_W"),
                                    silent = FALSE))
  class(model00) <- "tmbmodel" # null model


#eta_alpha, eta_theta and eta_omega are fixed across years
 model1 <- update(model00,optimizer="nlminb",list(muAlpha=factor(c(1,2)),
                                                  muTheta= factor(c(1,2)),
                                                  muOmega=factor(c(1,2))),
                  cv=0.2,seed=123)
 model1
 #aic = 28174.88 
 #p = 6

 #######################
 # directional selection for P
 model2 <- update(model1,optimizer="nlminb",list(muAlpha=factor(c(NA,1)),
                                                 muTheta= factor(c(NA,1)),
                                                 muOmega=factor(c(NA,1)),
                                                 beta0_P=NULL,
                                                 beta1_P=NULL),
                  cv=0.1,seed=123)
 model2
 #aic = 28173.9 
 #p = 5
 
 #directional selection for P and P has random intercept for each year
 model2a <- update(model2,optimizer="nlminb",list( u0_P=NULL,
                                                   u1_P=NULL,
                                                   logSdRanP =factor(c(1,NA))),
                   cv=0.2,seed=123)
 model2a
 #aic = 27939.3 
 #p = 6
 
 #directional selection for P and P also has random slope for each year
 model2b <- update(model2a,optimizer="nlminb",list(logSdRanP =factor(c(1,2))),
                   cv=0.2,seed=123)
 model2b
 #aic = 27847.1 
 #p = 7
 
 #directional selection for P and  correlated annual random intercept and slope
 model2c <- update(model2b,optimizer="nlminb",list(thetaRanP=NULL),
                   cv=0.01,seed=123)
 
 model2c
 #aic = 2785.26
 #p = 8
 
 #########################
 # directional selection for W
 model3 <- update(model1,optimizer="nlminb",list(muAlpha=factor(c(1,NA)),
                                                 muTheta= factor(c(1,NA)),
                                                 muOmega=factor(c(1,NA)),
                                                 beta0_W=NULL,
                                                 beta1_W=NULL),
                  cv=0.2,seed=123)
 model3
 #aic = 28175.43  
 #p = 5
 
 #directional selection for W and W has random intercept for each year
 model3a<- update(model3,optimizer="nlminb",list( u0_W=NULL,
                                                  u1_W=NULL,
                                                  logSdRanW =factor(c(1,NA))),
                  cv=0.2,seed=123)
 model3a
 #aic = 27636.37 
 #p = 6
 
 #directional selection for W, W also has random slope for each year
 model3b <- update(model3a,optimizer="nlminb",list(logSdRanW =factor(c(1,2))),
                   cv=0.2,seed=123)
 model3b
 #aic = 27437.87 
 #p = 7
 
 #directional selection for W and  correlated annual random intercept and slope
 model3c <- update(model3b,optimizer="nlminb",list(thetaRanW=NULL),
                   cv=0.1,seed=123)
 model3c
 #aic =  27365.22 
 #p = 8
 
 #########################
 # random alpha theta and omega based on model 1, 2c and 3c
 
 # model 4a hard to converge, so update model4a based on model 4
 model4 <- update(model1, optimizer="nlminb",list(alpha=NULL,
                                                  theta = NULL,
                                                  omega = NULL,
                                                  logSigma = factor(c(1,NA,2))),
                  cv=0.2,seed=123)
 model4
 #aic = 27169.13 
 #p = 8
 
 # stablized selection for P and W, random alpha, theta and omega
 model4a <- update(model4, optimizer="nlminb",list(logSigma = NULL),
                   cv=0.1,seed=123)                              
 model4a
 #aic = 27042.34  
 #p = 9
 
 # model 4c hard to converge, so update model4c based on model 4b
 model4b <- update(model2c,optimizer="nlminb",list(alpha=NULL,
                                                   theta = NULL,
                                                   omega = NULL,
                                                   logSigma = factor(c(1,NA,3))),
                   cv=0.02,seed=123)
 model4b 
 #aic = 27041.97
 #p = 10
 
 # directional selection for P, stabilized selection for W, random alpha, theta and omega
 model4c <- update(model4b,optimizer="nlminb",list(logSigma = factor(c(1,2,3))),
                   cv=0.1,seed=123)
 model4c 
 #aic = 27038 
 #p = 11
 
 # model 4e hard to converge, so update model4e based on model 4d
 model4d <- update(model3c,optimizer="nlminb",list(alpha=NULL,
                                                   theta = NULL,
                                                   omega = NULL,
                                                   logSigma = factor(c(1,NA,2))),
                   cv=0.01,seed=123)
 model4d
 #aic =  27057.22 
 #p = 10
 
 # directional selection for W, stabilized selection for P, random alpha, theta and omega
 model4e <- update(model4d,optimizer="nlminb",list(logSigma = NULL),
                   cv=0.01,seed=123)
 model4e
 #aic = 27024.37 
 #p = 11
 
 #########################
 # AR1 alpha theta and omega based on model 4a, 4c and 4e
 
 #stabilising selection for P and W, AR1 for alpha, theta and omega
 model5 <- update(model4a,optimizer="nlminb",list(Phi=factor(matrix(c(1,NA,NA,NA,2,NA,NA,NA,3),nrow=3))),
                  cv=0.2,seed=123)
 model5
 #aic = 27009.63  
 #p = 12
 
 # directional selection for P, AR1 for alpha, theta and omega
 model5a <- update(model4c,optimizer="nlminb",list(Phi=factor(matrix(c(1,NA,NA,NA,2,NA,NA,NA,3),nrow=3))),
                   cv=0.55,seed=123)
 model5a
 #aic = 26968.73 
 #p = 14
 
 # directional selection for W, AR1 for alpha, theta and omega
 model5b <- update(model4e,optimizer="nlminb",list(Phi=factor(matrix(c(1,NA,NA,NA,2,NA,NA,NA,3),nrow=3))),
                   cv=0.1,seed=123)
 
 model5b
 #aic = 27002.46
 #p = 14
 
 #########################
 # VAR1 alpha theta and omega based on model 4a, 4c and 4e
 
 #stabilising selection for P and W, VAR1 for alpha, theta and omega
 model6 <- update(model4a,optimizer="nlminb",list(Phi=NULL),
                  cv=0.2,seed=123)
 model6
 #aic = 27016.91 
 #p = 18
 #only entries 1,5 are siginicant
 
 # directional selection for P, VAR1 for alpha, theta and omega
 model6a <- update(model4c,optimizer="nlminb",list(Phi=NULL),
                   cv=0.2,seed=123)
 model6a
 #aic = 27016.72 
 #p = 20
 #no significant entries
 
 # directional selection for W, VAR1 for alpha, theta and omega, cannot get optimization 
 #converged with full Phi, exclude one entry  of Phi and later add it back
 model6b <- update(model4e,optimizer="nlminb",list(Phi=factor(matrix(c(1,2,3,4,5,NA,6,7,8),nrow=3))),
                   cv=0.2,seed=123)
 model6b
 #aic = 27008.48 
 #p = 19
 #only first entry in Phi is significant 
 
 
 #########################
 # based on model 6, 5a and 6b, keep only significant (at statistics 0.05) entries in Phi,
 # add the entry not estimated by model 6b back
 
 #stabilising selection for P and W, AR1 for alpha, theta 
 model7 <- update(model6,optimizer="nlminb",list(Phi=factor(matrix(c(1,NA,NA,NA,2,NA,NA,NA,NA),nrow=3))), 
                  cv=0.01,seed=123)
 model7
 #aic = 27010.06 
 #p = 11
 
 # directional selection for P, AR1 for alpha, theta 
 model7a <- update(model5a,optimizer="nlminb",list(Phi=factor(matrix(c(1,NA,NA,NA,2,NA,NA,NA,NA),nrow=3))), 
                   cv=0.1,seed=123)
 model7a
 #aic = 26967.02 
 #p = 13
 
 # directional selection for W, AR1 for alpha, theta 
 model7b <- update(model6b,optimizer="nlminb",list(Phi=factor(matrix(c(1,NA,NA,NA,NA,2,NA,NA,NA),nrow=3))), 
                   cv=0.1,seed=123)
 model7b
 #aic =  27004.37 
 #p = 13
 
 # directional selection for W, AR1 for alpha, theta  (based on mode 7b, only keep the first significant entry)
 model7c <- update(model6b,optimizer="nlminb",list(Phi=factor(matrix(c(1,NA,NA,NA,NA,NA,NA,NA,NA),nrow=3))), 
                   cv=0.1,seed=123)
 model7c
 #aic =  27003.99 
 #p = 12
 
 ################################
 # so far for stabilising selection for P and W, model 5 is the best
 # for directional selection for P, model 7a is the best
 # for directional selection for W, model 5b is the best
 
 # next we add correlations to the errors of alpha, theta and omega processes based on model 5, 7 and 7a
 
 ##stabilising selection for P and W, VAR1 for alpha, theta and omega with correlated errors
 model8 <- update(model5,optimizer="nlminb",list(parCor=factor(c(1,2,3))),
                  cv=0.1,seed=123)
 model8
 #aic = 26985.69 
 #p = 15
 #does improve model fitting
 
 #directional selection for P, VAR1 for alpha, theta and omega with correlated errors
 model8a <- update(model7a,optimizer="nlminb",list(parCor=factor(c(1,2,3))),
                   cv=0.1,seed=123)
 model8a
 #aic = 26946.99
 #p = 16
 #does improve model fitting
 
 #directional selection for W, VAR1 for alpha, theta and omega with correlated errors
 model8b <- update(model5b,optimizer="nlminb",list(parCor=factor(c(1,2,3))),
                   cv=0.01,seed=123)
 model8b
 #aic = 27002.71 
 #p = 17
 # does not improve model fitting
 
 ################################
 #so far model 8a is the best, add mother effect
 
 #stabilising selection for P and W, AR1 for alpha, theta and omega, mother effects included
 model9 <- update(model8a,optimizer="nlminb",list(epsi=NULL,logTauMother=NULL),
                  cv=0.01,seed=123)
 
 model9
 #aic =  26927.89 
 #p = 18
 
 #########################
 # directional selection for P and W
 model10 <- update(model4c,optimizer="nlminb",list(muAlpha=factor(c(NA,NA)),
                                                   muTheta= factor(c(NA,NA)),
                                                   muOmega=factor(c(NA,NA)),
                                                   beta0_W=NULL,
                                                   beta1_W=NULL,
                                                   u0_W=NULL,
                                                   u1_W=NULL,
                                                   logSdRanW =factor(c(1,2)),
                                                   thetaRanW = NULL),
                   cv=0.01,seed=123)
 model10
 #aic = 27032
 #p = 13
 
 ################################
 #so far model 9 is the best,  play around it next
 
 bestmodel <- model9
 #aic = 26927.89 
 
 # add all the entries back to Phi again to see if the model fit improves
 model11 <- update(bestmodel,optimizer="nlminb",list(Phi=NULL),
                   cv=0.1,seed=123)
 model11
 #aic = 26931.965
 #p = 25
 
 # only keep the significant entry of Phi and parCor
 model12 <- update(bestmodel,optimizer="nlminb",list(Phi=factor(matrix(c(1,NA,NA,NA,NA,NA,NA,NA,NA),nrow=3)), 
                                                     parCor=factor(c(1,NA,NA))),
                   cv=0.1,seed=123)
 model12
 #aic = 26942.21 
 #p = 15
 
 
 #only keep entry 1 and 5 of Phi and exclude unsignificant entries in parCor
 model13 <- update(bestmodel,optimizer="nlminb",list(parCor=factor(c(1,NA,NA))),
                   cv=0.01,seed=123)
 model13
 #aic = 26928.43  
 #p = 16
 
 
 ################################
 #now model 9a is still the best,  play around it next
 
 #change random alpha into fixed, theta and omega are random
 model14 <- update(bestmodel,optimizer="nlminb",list(logSigma=factor(c(NA,1,2)),
                                                     Phi=factor(matrix(rep(NA,9),nrow=3)),
                                                     parCor=factor(c(NA,NA,NA))),
                   cv=0.1,seed=123)
 model14
 #aic = 27138.36 
 #p = 12
 
 #change random alpha into fixed, AR1 theta and omega 
 model14a <- update(bestmodel,optimizer="nlminb",list(logSigma=factor(c(NA,1,2)),
                                                      Phi= factor(matrix(c(NA,NA,NA,NA,1,NA,NA,NA,2),nrow=3)),
                                                      parCor=factor(c(NA,NA,NA))),
                    cv=0.1,seed=123)
 model14a
 #aic =  27131.87
 #p = 14
 
 #change random theta into fixed, alpha and omega are random
 model14b <- update(bestmodel,optimizer="nlminb",list(logSigma=factor(c(1,NA,2)),
                                                      Phi=factor(matrix(rep(NA,9),nrow=3)),
                                                      parCor=factor(c(NA,NA,NA))),
                    cv=0.01,seed=123)
 model14b
 #aic =  27024.88 
 #p = 12
 
 #change random theta into fixed, AR1 alpha and omega 
 model14c <- update(bestmodel,optimizer="nlminb",list(logSigma=factor(c(1,NA,2)),
                                                      Phi=factor(matrix(c(1,NA,NA,NA,NA,NA,NA,NA,2),nrow=3)),
                                                      parCor=factor(c(NA,NA,NA))),
                    cv=0.2,seed=123)
 model14c
 #aic = 27013.36 
 #p = 14
 
 #change random omega into fixed, AR1 alpha and theta 
 model14d <- update(bestmodel,optimizer="nlminb",list(logSigma=factor(c(1,2,NA)),
                                                      Phi=factor(matrix(c(1,NA,NA,NA,2,NA,NA,NA,NA),nrow=3)),
                                                      parCor=factor(c(NA,NA,NA))),
                    cv=0.1,seed=123)
 model14d
 #aic = 26952.73 
 #p = 14
 
 #change random omega into fixed, VAR1 alpha and theta 
 model14e <- update(bestmodel,optimizer="nlminb",list(logSigma=factor(c(1,2,NA)),
                                                      Phi=factor(matrix(c(1,2,NA,3,4,NA,NA,NA,NA),nrow=3)),
                                                      parCor=factor(c(NA,NA,NA))),
                    cv=0.1,seed=123)
 model14e
 #aic = 26956.53 
 #p = 16
 
 #change random omega into fixed, VAR1 alpha and theta, add correlation to errors of alpha and theta
 model14f <- update(bestmodel,optimizer="nlminb",list(logSigma=factor(c(1,2,NA)),
                                                      Phi=factor(matrix(c(1,2,NA,3,4,NA,NA,NA,NA),nrow=3)),
                                                      parCor=factor(c(1,NA,NA))),
                    cv=0.1,seed=123)
 model14f
 #aic = 26934.99 
 #p = 17
 
 #change random omega into fixed, VAR1 alpha and theta with significant entries in Phi, add correlation to errors of alpha and theta
 model14g <- update(bestmodel,optimizer="nlminb",list(logSigma=factor(c(1,2,NA)),
                                                      Phi=factor(matrix(c(1,NA,NA,NA,2,NA,NA,NA,NA),nrow=3)),
                                                      parCor=factor(c(1,NA,NA))),
                    cv=0.01,seed=123)
 model14g
 #aic =  26931.16 
 #p = 15
 #model 14g with fixed omega is very close to the best model!
 
 
 #remove the mother effect from P
 model14h <- update(bestmodel,optimizer="nlminb",list(logTauMother=factor(c(NA,1))),
                    cv=0.01,seed=123)
 model14h
 #aic = 26948.99 
 #p = 17
 
 #remove the mother effect from W
 model14i <- update(bestmodel,optimizer="nlminb",list(logTauMother=factor(c(1,NA))),
                    cv=0.1,seed=123)
 model14i
 #aic = 26934.33 
 #p = 17
 
 
 #now we are sure model 9  is indeed the best model, we can directly arrive at it based on the null model
 #################################
 #arrive at the best model based on null model 
 BestModel <- update(model00,optimizer="nlminb",list(alpha = NULL,
                                                     theta = NULL,
                                                     omega = NULL,
                                                     logSigma = NULL,
                                                     muAlpha = factor(c(NA,1)),
                                                     muTheta = factor(c(NA,1)),
                                                     muOmega = factor(c(NA,1)),
                                                     Phi = factor(matrix(c(1,NA,NA,NA,2,NA,NA,NA,NA),nrow=3)),
                                                     parCor = factor(c(1,2,3)),
                                                     beta0_P= NULL,
                                                     beta1_P= NULL,
                                                     u0_P= NULL,
                                                     u1_P = NULL,
                                                     logSdRanP = factor(c(1,2)),
                                                     thetaRanP = NULL,
                                                     epsi = NULL,
                                                     logTauMother = NULL),
                     cv=0.02,seed=123)
 BestModel
 #aic = 26927.89 
 #p = 18
 #save(BestModel, file="BestModel.Rdata")
 #summary(BestModel$sdreport,"fixed",p.value = TRUE)
 #summary(BestModel$sdreport,t,p.value = TRUE)
 
 
 
 ##############################################################################
 # the code below is for fitting the selected model (model9) with data from
 # 1973-2015.
 #############################################################################
 
 #load TMB template 
 library(TMB)
 compile("great_tit_Final.cpp")
 dyn.load(dynlib("great_tit_Final"))
 
 #################################################################
 # load data
 source("great_tit_data_pre.R")
 
 # Creat data list
 dataPart = 
   list(
     LayDateApril  = as.numeric(as.character(mytitbroodPart$LayDateApril)),
     NumFle = as.numeric(as.character(mytitbroodPart$NumberFledged)),
     t = as.factor(mytitbroodPart$YearOfBreeding),
     j = as.factor(mytitbroodPart$MothersID)
   )
 
 
 # Initialize parameters
 parametersPart <- list(
   alpha        = rep(1,nlevels(dataPart$t)),
   theta        = rep(5,nlevels(dataPart$t)),
   omega        = rep(1,nlevels(dataPart$t)),
   muAlpha      = rep(2,2),
   muTheta      = rep(20,2),
   muOmega      = rep(4,2),
   logSigma     = c(-2,3,-2),
   parCor       = rep(0,3),
   Phi          = matrix(rep(0,9),nrow=3),
   beta0_P   = 3,
   beta1_P   = 0,
   u0_P   = rep(1,nlevels(dataPart$t)),
   u1_P    = rep(1,nlevels(dataPart$t)),
   beta0_W = 0,
   beta1_W = 0,
   u0_W = rep(1,nlevels(dataPart$t)),
   u1_W = rep(1,nlevels(dataPart$t)),
   logSdRanP    = rep(-2,2),
   thetaRanP    = 0,
   logSdRanW = rep(-2,2),
   thetaRanW = 0,
   logTauMother= rep(-2,2),
   epsi         = rep(0.1,nlevels(dataPart$j))
 )
 
 map00Part = list(
   alpha          = factor(rep(NA,nlevels(dataPart$t))),
   theta          = factor(rep(NA,nlevels(dataPart$t))),
   omega          = factor(rep(NA,nlevels(dataPart$t))),
   muAlpha        = factor(rep(NA,2)),
   muTheta        = factor(rep(NA,2)),
   muOmega        = factor(rep(NA,2)),
   logSigma       = factor(rep(NA,3)),
   parCor         = factor(rep(NA,3)),
   Phi            = factor(matrix(rep(NA,9),nrow=3)), 
   beta0_P     = factor(NA),
   beta1_P     = factor(NA),
   u0_P     = factor(rep(NA,nlevels(dataPart$t))),
   u1_P      = factor(rep(NA,nlevels(dataPart$t))),
   beta0_W   = factor(NA),
   beta1_W   = factor(NA),
   u0_W   = factor(rep(NA,nlevels(dataPart$t))),
   u1_W   = factor(rep(NA,nlevels(dataPart$t))),
   logSdRanP      = factor(rep(NA,2)),
   thetaRanP      = factor(NA),
   logSdRanW   = factor(rep(NA,2)),
   thetaRanW   = factor(NA),
   logTauMother  = factor(rep(NA,2)),
   epsi           = factor(rep(NA,nlevels(dataPart$j)))
 )
 
 #######################################################
 library("numDeriv")
 source("functions.R")
 
 model00Part <- list(parameters=parametersPart, map=map00Part, 
                     args.MakeADFun=list(data = dataPart, DLL = "great_tit_Final", 
                                         random = c("alpha","theta","omega","epsi",
                                                    "u0_P","u1_P","u0_W",
                                                    "u1_W"),
                                         silent = FALSE))
 class(model00Part) <- "tmbmodel" # null model
 
 
 BestModelPart <- update(model00Part,optimizer="nlminb",list(alpha = NULL,
                                                             theta = NULL,
                                                             omega = NULL,
                                                             logSigma = NULL,
                                                             muAlpha = factor(c(NA,1)),
                                                             muTheta = factor(c(NA,1)),
                                                             muOmega = factor(c(NA,1)),
                                                             Phi = factor(matrix(c(1,NA,NA,NA,2,NA,NA,NA,NA),nrow=3)),
                                                             parCor = factor(c(1,2,3)),
                                                             beta0_P= NULL,
                                                             beta1_P= NULL,
                                                             u0_P= NULL,
                                                             u1_P = NULL,
                                                             logSdRanP = factor(c(1,2)),
                                                             thetaRanP = NULL,
                                                             epsi = NULL,
                                                             logTauMother = NULL),
                         cv=0.1,seed=123)
 BestModelPart
 #aic= 20190.21 
 #p= 18
 #save(BestModelPart, file="BestModelPart.Rdata")
 #summary(BestModelPart$sdreport,"fixed",p.value = TRUE)
 #summary(BestModelPart$sdreport,"report")
 
 
 # For comparison with estimates from Chevin et.al 2015, we fit the best model with fixed omega with partial data set
 # the estimates are reported in the last column in Table S1 
 modeltest <- update(BestModelPart,optimizer="nlminb",list(logSigma=factor(c(1,2,NA)),
                                                           Phi=factor(matrix(c(1,NA,NA,NA,2,NA,NA,NA,NA),nrow=3)),
                                                           parCor=factor(c(1,NA,NA))),
                     cv=0.01,seed=123)
 modeltest
 #summary(modeltest$sdreport,"fixed",p.value = TRUE)
 #summary(modeltest$sdreport,"report")
 
 
 ##make predictions inside the template
 ##############################
 
 data = 
   list(
     LayDateApril  = as.numeric(as.character(mytitbroodFinal$LayDateApril)), # lay date for each brood
     NumFle = as.numeric(as.character(mytitbroodFinal$NumberFledged)), # number of fledglings for each brood
     t = as.factor(mytitbroodFinal$YearOfBreeding), # index of year from 1955 to 2015
     j = as.factor(mytitbroodFinal$MothersID), # mother identity
     makePred = 1,
     predDates = rep((seq(0,100,by=10)),61),
     tPred = as.factor(rep(c(1:61),each=11))
   )
 
 model00 <- list(parameters=parameters, map=map00, 
                 args.MakeADFun=list(data = data, DLL = "great_tit_Final", 
                                     random = c("alpha","theta","omega","epsi",
                                                "u0_P","u1_P","u0_W",
                                                "u1_W"),
                                     silent = FALSE))
 class(model00) <- "tmbmodel" # null model
 
 
 BestModel <- update(model00,optimizer="nlminb",list(alpha = NULL,
                                                     theta = NULL,
                                                     omega = NULL,
                                                     logSigma = NULL,
                                                     muAlpha = factor(c(NA,1)),
                                                     muTheta = factor(c(NA,1)),
                                                     muOmega = factor(c(NA,1)),
                                                     Phi = factor(matrix(c(1,NA,NA,NA,2,NA,NA,NA,NA),nrow=3)),
                                                     parCor = factor(c(1,2,3)),
                                                     beta0_P= NULL,
                                                     beta1_P= NULL,
                                                     u0_P= NULL,
                                                     u1_P = NULL,
                                                     logSdRanP = factor(c(1,2)),
                                                     thetaRanP = NULL,
                                                     epsi = NULL,
                                                     logTauMother = NULL),
                     cv=0.05,seed=123)
 BestModel
 #aic = 26927.89 
 #p = 18
 #save(BestModel, file="BestModel.Rdata")
 #summary(BestModel$sdreport,"fixed",p.value = TRUE)
 #summary(BestModel$sdreport,"report",p.value = TRUE)
 
 




