
 #Tables and figures in the paper#
 #####################################################################################
 #to run the following codes, you have to run the code for model selection, get the
 #candidate models fitted, and the best model (BestModel) object was saved as "BestModel.Rdata".
 load("BestModel.Rdata")
 source("great_tit_data_pre.R")
 rep_fixed <- summary(BestModel$sdreport,"fixed",p.value = TRUE)
 rep_random <- summary(BestModel$sdreport,"random")
 rep_value <- summary(BestModel$sdreport,"report")

 #estimates of fixed parameters and reported estimates for directional selection, calculate the confidence interval
 t <- as.factor(mytitbroodFinal$YearOfBreeding)
 CI_fixed_left <- rep_fixed[,1]-1.96*rep_fixed[,2]
 CI_fixed_right <- rep_fixed[,1]+1.96*rep_fixed[,2]
 rep_fixed <- round(cbind(rep_fixed,CI_fixed_left,CI_fixed_right),3)
 CI_value_left <- rep_value[,1]-1.96*rep_value[,2]
 CI_value_right <- rep_value[,1]+1.96*rep_value[,2]
 rep_value <- round(cbind(rep_value,CI_value_left,CI_value_right),3)

 attach(mytitbroodFinal)
 #calculation of average within-year phenotypic standard deviation#######
 sd_all <- aggregate(LayDateApril~YearOfBreeding, data=mytitbroodFinal, FUN=function(x) 
  c(sd <-sd(x)))
 sd_mean <- mean(sd_all[,2])
 #17.54047
 #calculation of annual selection strenghth#######
 estOmega <- unname(subset(rep_random,row.names(rep_random)=="omega"))
 dataOmega <- estOmega[,1]+unname(subset(rep_fixed,row.names(rep_fixed)=="muOmega")[,1])
 selStren <- sd_all[,2]^2/(exp(dataOmega)^2+sd_all[,2]^2)
 #min(selStren) #0.006 2013
 #max(selStren) #0.305  1962
 #mean(selStren) #0.135
 detach(mytitbroodFinal)

 #calculation of average within-year phenotypic standard deviation#######
 # with partial data from 1973 to 2015
 attach(mytitbroodPart)
 sd_all <- aggregate(LayDateApril~YearOfBreeding, data=mytitbroodPart, FUN=function(x) 
   c(sd <- sd(x)))
 sd_mean <- mean(sd_all[,2])
 #16.39423
 detach(mytitbroodPart)

 #########plot directional selection gradients for P (Figure 1)#############
 beta1_P <- subset(rep_fixed,row.names(rep_fixed)=="beta1_P")[,1]
 u1_P <- subset(rep_random,row.names(rep_random)=="u1_P")
 beta1_P_vector=rep(beta1_P,100)
 sdVecRanP <-subset(rep_value,row.names(rep_value)=="sdVecRanP")[2,1]

 layout(matrix(1:2,1,2),widths=c(5.2,2),heights=c(5.2))
 par(mar=c(5,5,4,2))
 require(plotrix)

 plotCI(c(1:61),unname(u1_P[,1]+beta1_P), uiw=u1_P[,2],xaxt="n",yaxt="n",xlab="Year",
       ylab="Selection gradient (the fixed slope + random slope)",ylim =c(-0.1,0.05),barcol="grey",gap=0,cex.lab=1.2)
 ticks=c(-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04)
 axis(2,at=ticks,labels=ticks)
 axis(1,at=c(1,6,11,16,21,26, 31,36,41,46,51,56,61),
     labels=c("1955","1960","1965","1970","1975","1980","1985","1990","1995","2000","2005","2010","2015"))
 abline(a=0,b=0,lty=2,col="red")

 par(mar=c(5,1,4,2))
 ylim <- c(-0.1,0.05)
 plot(NA,ylim=ylim,xlim=c(0,1),xlab="Distribution of selection gradients",ylab="",yaxt="n")
 ticks=c(-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04)
 axis(2,at=ticks,labels=ticks)
 box()
 zz <- seq(ylim[1],ylim[2],length=100)
 lines(exp(-(zz-beta1_P_vector)^2/(2*sdVecRanP^2)),zz)
 abline(a=beta1_P,b=0,lty=2,col="red")
 arrows(0.2,-0.025,0.2,-0.012, length = 0.1, angle = 30, lty=1,code=1, lwd=2, col="grey")
 legend(-0.11,0.005, legend = "the estimated fixed slope ", bty = "n", col=c("grey"), cex=1)
 legend(0,0, legend = expression(hat(beta)[p]^(1) == -0.025), bty = "n", col=c("grey"), cex=1)

 # number of years with positive directional selection gradients.
 # length(which(unname(u1_P[,1]+beta1_P)>0))

########### plot optimum movement (Figure 2)#########
 t <- as.factor(mytitbroodFinal$YearOfBreeding)
 estFixed <- summary(BestModel$sdreport,"fixed")
 muTheta <- subset(estFixed,rownames(estFixed)=="muTheta")[1,1]
 estTheta <- muTheta+unname(summary(BestModel$sdreport,"random")[,1])[62:(62+60)]
 #mean(estTheta)
 year <- c(1:range(as.numeric(t))[2])
 plot(estTheta~year,type="l",ylab="Optimal laying date, mean laying date",xlab="Year",ylim=range(-37:80),
      xaxt="n",yaxt="n",col="deepskyblue3",lwd=1.5,cex.lab=1.5)
 axis(1,at=seq(1,61,by=5),
      labels=c("1955","1960","1965","1970","1975","1980","1985","1990","1995","2000","2005","2010","2015"))
 axis(1,at=year,labels=FALSE,tcl=-.25)
 axis(2,at=c(1,1+30,1+30+31),labels=c("April 1", "May 1", "June 1"),line=-.8)
 #calculate mean laying date by year
 meanz <- tapply(as.numeric(mytitbroodFinal$LayDateApril),
                as.factor(mytitbroodFinal$YearOfBreeding),FUN=mean)
 #mean(meanz)
 points(year,unname(meanz),pch=20)
 #calculate 95% confidence interval
 left <- estTheta-1.96*unname(summary(BestModel$sdreport,"random")[,2][62:(62+60)])
 right <- estTheta+1.96*unname(summary(BestModel$sdreport,"random")[,2][62:(62+60)])
 matlines(year,cbind(left,right),lty=2,col="deepskyblue4")

##############################################
##plot the observed number of fledglings and estimated number of fledglings conditional on
#zero random mother effects (Figure 3)
 expect_y = subset(sdreportAll, rownames(sdreportAll) == "expect_y")
 predDatay <- data.frame(data$tPred, data$predDates, unname(expect_y[,1:2]))

 layout(matrix(c(1:66), 11, 6,byrow=F))
 par(mar=c(0.1,0.1,0.1,0.1), oma=c(5, 7, 1, 1),mgp=c(3,0.4,0))
 for (i in 1:61){
   
   estYearDat <- mytitbroodFinal[which(mytitbroodFinal$YearOfBreeding==1954+i),]
   predYearDat <- predDatay[which(predDatay$data.tPred==i),]
  
   myx <- c(0:100) 
   colType <- estYearDat$BroodType
   plot(estYearDat$LayDateApril,estYearDat$NumberFledged,pch=20,xlab='',ylab='',cex = 0.5,
       cex.axis=0.7,xaxt="n", yaxt="n",ylim=c(-0.5,15),xlim=c(0,100),col=colType+10 )
   lwMean <- loess(estYearDat$NumberFledged~estYearDat$LayDateApril)
   conIntMean <-predict(lwMean, newdata=myx, se=T)
   text(5,14,1954+i,cex=1)
  
   lines(myx, conIntMean$fit,col="darkgrey",lwd=1.5)
   lines(myx, conIntMean$fit - qt(0.975,conIntMean$df)*conIntMean$se, lty=2,col="darkgrey")
   lines(myx, conIntMean$fit + qt(0.975,conIntMean$df)*conIntMean$se, lty=2,col="darkgrey")
  
   lines(predYearDat$data.predDates,predYearDat$X1,lwd=1.5)
   lines(predYearDat$data.predDates,predYearDat$X1+1.96*predYearDat$X2, lty=2)
   lines(predYearDat$data.predDates,predYearDat$X1-1.96*predYearDat$X2, lty=2)
  
   if (estYearDat$YearOfBreeding[1]  %in% c(1965,1976,1987,1998,2009,2015)){
     axis(1,at=c(1,1+30,1+30+31,1+30+31+30),labels=c("April 1", "May 1", "June 1","July 1"),cex.axis=0.8)}
   if (estYearDat$YearOfBreeding[1]  %in% c(1955:1965)){
     axis(2,cex.axis=0.6)}
 }
 legend(1,13,legend = c("First brood", "Replacement brood (first brood failed)", "Second brood after successful first brood"),  
        cex = 0.67, col = c("red","green",  "blue"),pch = c(16,16,16),box.lty=0)

 mtext("Observed and predicted number of fledglings",side=2,line=1.8,outer = T,cex=1)
 mtext("Laying date",side=1,line=1.8,outer = T,cex=1)


##############################################
#########Figures in supporting document######

#######plot the index of successful broods and predicted
##probability of successful broods for each year (Fgiure S2)
 sdreportAll <- summary(BestModel$sdreport,"report",p.value = TRUE)
 p_nonzeroy <- subset(sdreportAll, rownames(sdreportAll) == "p_nonzeroy")
 predDatap <- data.frame(data$tPred, data$predDates, unname(p_nonzeroy[,1:2]))

 layout(matrix(c(1:66), 11, 6,byrow=F))
  par(mar=c(0.1,0.1,0.1,0.1), oma=c(5, 7, 1, 1),mgp=c(3,0.4,0))
for (i in 1:61){
  
  estYearDat <- mytitbroodFinal[which(mytitbroodFinal$YearOfBreeding==1954+i),]
  predYearDat <- predDatap[which(predDatap$data.tPred==i),]
  myx <- c(0:100)
  
  colType <- estYearDat$BroodType
  with(estYearDat,plot(LayDateApril, sucBroodIndex, pch=20,xlab='', ylab='',cex = 0.8,cex.axis=0.7,xaxt="n", yaxt="n",
                       ylim=c(0,1.1),xlim=c(0,100), col=colType+10))
  lwP <- loess(estYearDat$sucBroodIndex~estYearDat$LayDateApril)
  conInt<-predict(lwP,newdata = myx,se=T)
  text(5,0.1,1954+i,side=3,cex=1)
  lines(myx,conInt$fit, col="darkgrey",lwd=1.5)
  lines(myx,conInt$fit - qt(0.975,conInt$df)*conInt$se, lty=2,col="darkgrey")
  lines(myx,conInt$fit + qt(0.975,conInt$df)*conInt$se, lty=2,col="darkgrey")
  
  lines(predYearDat$data.predDates,predYearDat$X1,lwd=1.5)
  lines(predYearDat$data.predDates,predYearDat$X1+1.96*predYearDat$X2, lty=2)
  lines(predYearDat$data.predDates,predYearDat$X1-1.96*predYearDat$X2, lty=2)
  
  if (estYearDat$YearOfBreeding[1]  %in% c(1965,1976,1987,1998,2009,2015)){
    axis(1,at=c(1,1+30,1+30+31,1+30+31+30),labels=c("April 1", "May 1", "June 1","July 1"),cex.axis=0.8)}
  if (estYearDat$YearOfBreeding[1]  %in% c(1955:1965)){
    axis(2,cex.axis=0.6)}
  if (estYearDat$YearOfBreeding[1] == 1955){
    legend(1,0.8,legend = c("First brood", "Replacement brood (first brood failed)", "Second brood after successful first brood"),  
           cex = 0.67, col = c("red","green",  "blue"),pch = c(16,16,16),box.lty=0)}
}
 mtext("Observed index and predicted probability of successful-brooding",side=2,line=1.8,outer = T,cex=1)
 mtext("Laying date",side=1,line=1.8,outer = T,cex=1)
 
 #######plot the nonzero number of fledglings and predicted
 ##number of nonzero fledglings (Fgiure S3)
 sdreportAll <- summary(BestModel$sdreport,"report",p.value = TRUE)
 w_nonzeroy <- subset(sdreportAll, rownames(sdreportAll) == "expect_nonzeroy")
 predDataw <- data.frame(data$tPred, data$predDates, unname(w_nonzeroy[,1:2]))
 
 layout(matrix(c(1:66), 11, 6,byrow=F))
 par(mar=c(0.1,0.1,0.1,0.1), oma=c(5, 7, 1, 1),mgp=c(3,0.4,0))
 for (i in 1:61){
   
   estYearDat <- mytitbroodFinal[which(mytitbroodFinal$YearOfBreeding==1954+i & mytitbroodFinal$NumberFledged>0),]
   predYearDat <- predDataw[which(predDataw$data.tPred==i),]
   myx <- c(0:100)
   
   colType <- estYearDat$BroodType
   with(estYearDat,plot(LayDateApril, NumberFledged, pch=20,xlab='', ylab='',cex = 0.8,cex.axis=0.7,xaxt="n", yaxt="n",
                        ylim=c(1,15),xlim=c(0,100), col=colType+10))
   lwP <- loess(estYearDat$NumberFledged~estYearDat$LayDateApril)
   conInt<-predict(lwP,newdata = myx,se=T)
   text(5,14,1954+i,cex=1)
   lines(myx,conInt$fit, col="darkgrey",lwd=1.5)
   lines(myx,conInt$fit - qt(0.975,conInt$df)*conInt$se, lty=2,col="darkgrey")
   lines(myx,conInt$fit + qt(0.975,conInt$df)*conInt$se, lty=2,col="darkgrey")
   
   lines(predYearDat$data.predDates,predYearDat$X1,lwd=1.5)
   lines(predYearDat$data.predDates,predYearDat$X1+1.96*predYearDat$X2, lty=2)
   lines(predYearDat$data.predDates,predYearDat$X1-1.96*predYearDat$X2, lty=2)
   
   if (estYearDat$YearOfBreeding[1]  %in% c(1965,1976,1987,1998,2009,2015)){
     axis(1,at=c(1,1+30,1+30+31,1+30+31+30),labels=c("April 1", "May 1", "June 1","July 1"),cex.axis=0.8)}
   if (estYearDat$YearOfBreeding[1]  %in% c(1955:1965)){
     axis(2,cex.axis=0.6)}
   if (estYearDat$YearOfBreeding[1] == 1955){
     legend(1,0.8,legend = c("First brood", "Replacement brood (first brood failed)", "Second brood after successful first brood"),  
            cex = 0.67, col = c("red","green",  "blue"),pch = c(16,16,16),box.lty=0)}
 }
 mtext("Observed index and predicted probability of successful-brooding",side=2,line=1.8,outer = T,cex=1)
 mtext("Laying date",side=1,line=1.8,outer = T,cex=1)
 
 
 ############# plot optimum movement estimated with partial data (Figure S4)
 rep_randomPart <- summary(BestModelPart$sdreport,"random")
 t <- as.factor(mytitbroodPart$YearOfBreeding)
 estFixedPart <- summary(BestModelPart$sdreport,"fixed")
 muThetaPart <- subset(estFixedPart,rownames(estFixedPart)=="muTheta")[1,1]
 estThetaPart <- muThetaPart+unname(subset(rep_randomPart,row.names(rep_randomPart)=="theta"))[,1]
 year <- c(1:range(as.numeric(t))[2])
 plot(estThetaPart~year,type="l",ylab="Optimal laying date, mean laying date with partial data",xlab="Year",ylim=range(-37:80),
      xaxt="n",yaxt="n",col="deepskyblue3",lwd=1.5)
 axis(1,at=c(1,8,13,18,23,28,33,38,43),
      labels=c("1973","1980","1985","1990","1995","2000","2005","2010","2015"))
 axis(1,at=year,labels=FALSE,tcl=-.25)
 axis(2,at=c(1,1+30,1+30+31),labels=c("April 1", "May 1", "June 1"),line=-.8)
 #calculate mean laying date by year
 meanzPart<-tapply(as.numeric(mytitbroodPart$LayDateApril),
                   as.factor(mytitbroodPart$YearOfBreeding),FUN=mean)
 points(year,unname(meanzPart),pch=20)
 #calculate 95% confidence interval
 left <- estThetaPart-1.96*unname(subset(rep_randomPart,row.names(rep_randomPart)=="theta"))[,2]
 right <- estThetaPart+1.96*unname(subset(rep_randomPart,row.names(rep_randomPart)=="theta"))[,2]
 matlines(year,cbind(left,right),lty=2,col="deepskyblue4")

 ########## plot movement of omega, alpha, probability of successful-brooding 
 ##and mean number of fledglings (Figure S5)
 par(mfrow=c(2,2))
 #plot width movement
 estOmega <- unname(subset(rep_random,row.names(rep_random)=="omega"))
 year <- c(1:61)
 dataOmega <- estOmega[,1]+unname(subset(rep_fixed,row.names(rep_fixed)=="muOmega")[,1])
 lwOmega <- loess(exp(dataOmega)~c(1:61))
 plot(exp(dataOmega)~year,type="l",xlab="Year",ylab="Movement of width of fitness function",xaxt="n")
 axis(1,at=c(1,6,11,16,21,26, 31,36,41,46,51,56,61),
      labels=c("1955","1960","1965","1970","1975","1980","1985","1990","1995","2000","2005","2010","2015"))
 axis(1,at=year,labels=FALSE,tcl=-.25)
 lines(year,lwOmega$fitted,col="darkgreen",lwd=2)
 
 #plot alpha against year
 estAlpha <- unname(subset(rep_random,row.names(rep_random)=="alpha"))
 year <- levels(as.factor(mytitbroodFinal$YearOfBreeding))
 dataAlpha <- estAlpha[,1]+unname(subset(rep_fixed,row.names(rep_fixed)=="muAlpha")[,1])
 lwAlpha <- loess(exp(dataAlpha)~c(1:61))
 plot(exp(dataAlpha)~year,type="l",xlab="Year",ylab="Estimated maximum number of fledglings")
 lines(year,lwAlpha$fitted,col="red",lwd=2)
 
 #plot mean 1-p against year
 estSuc <- 1-unname(subset(rep_value,row.names(rep_value)=="prob"))[,1]
 year <- as.factor(mytitbroodFinal$YearOfBreeding)
 dataP <- cbind(estSuc,year)
 meanP <- aggregate(estSuc~year, data=dataP, FUN=function(x) 
   c(mean=mean(x), count=length(x)))[,2][,1]
 dataP <- as.data.frame(dataP)
 lwP <- loess(meanP~unique(year),data=dataP)
 plot(meanP~levels(year),type="l",xlab="Year",ylab="Mean probability of successful-brooding")
 lines(levels(year),lwP$fitted,col="blue",lwd=2)
 
 #plot the mean number of fledglings against year
 estMean <- unname(subset(rep_value,row.names(rep_value)=="poi_lamb"))[,1]
 year <- as.factor(mytitbroodFinal$YearOfBreeding)
 dataMean <- cbind(estMean,year)
 meanW <- aggregate(estMean~year, data=dataMean, FUN=function(x)
   c(mean=mean(x), count=length(x)))[,2][,1]
 dataMean <- as.data.frame(dataMean)
 lwMean <- loess(meanW~unique(year),data=dataMean)
 plot(meanW~levels(year),type="l",xlab="Year",ylab="Mean number of fledglings")
 lines(levels(year),lwMean$fitted,col="yellow",lwd=2)

 ####### plot autocorrelation of the random slopes for directional selection via P (Figure S6)
 u1_P <- subset(rep_random,row.names(rep_random)=="u1_P")
 AutoCorrelation <- acf(unname(u1_P[,1]), plot = FALSE)
 plot(AutoCorrelation, main = "ACF for random slopes")

 
 
 

