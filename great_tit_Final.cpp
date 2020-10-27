
// This file is C++ code defining the model as a template
#include <TMB.hpp>

template<class Type>

Type objective_function<Type>::operator() ()
{
  using namespace density;
  
  DATA_VECTOR(LayDateApril); // laying dates
  DATA_VECTOR(NumFle); // observed number of fledglings 
  DATA_FACTOR(t); // years 
  DATA_FACTOR(j); // mothers
  
  DATA_INTEGER(makePred); // an integer telling if we want to make predictions
  DATA_VECTOR(predDates); // on what layings dates we want to make predicitons
  DATA_FACTOR(tPred); // for which years we want to make predictions
  
  PARAMETER_VECTOR(alpha); 
  PARAMETER_VECTOR(theta);
  PARAMETER_VECTOR(omega);
  
  PARAMETER_VECTOR(muAlpha);
  PARAMETER_VECTOR(muTheta);
  PARAMETER_VECTOR(muOmega); vector<Type> expMuOmega=exp(muOmega); ADREPORT(expMuOmega);
  //sigma is standard deviation vector of errors for AR1/VAR1 alpha, theta and omega
  PARAMETER_VECTOR(logSigma); vector<Type> sigma = exp(logSigma); ADREPORT(sigma);
  //vector<Type> varVec = sigma*sigma; ADREPORT(varVec);
  
  //vector of entries of cholesky decomposition of correlation matrix
  PARAMETER_VECTOR(parCor);
  //to guarantee phi falls between -1 and 1, transform Phi into PhiB, the reported
  //estimates of Phi in the paper is actually PhiB here.
  PARAMETER_MATRIX(Phi);  matrix<Type>  PhiB =2/(1+exp(-2*Phi.array()))-1; ADREPORT(PhiB);
  
  PARAMETER(beta0_P); //fixed intercept of episode P
  PARAMETER(beta1_P); //fixed slope of episode P
  PARAMETER_VECTOR(u0_P); //random intercepts
  PARAMETER_VECTOR(u1_P); //random slopes
  
  PARAMETER(beta0_W); //fixed intercept of episode W
  PARAMETER(beta1_W); //fixed slope of episode W
  PARAMETER_VECTOR(u0_W); //random intercepts
  PARAMETER_VECTOR(u1_W); //random slopes
  
  //sdVecRanP is vector of standard deviations of random intercepts and slopes for component P
  PARAMETER_VECTOR(logSdRanP);    vector<Type> sdVecRanP = exp(logSdRanP); ADREPORT(sdVecRanP);
  //thetaRanP is correlation between random intercepts and slopes for component P
  PARAMETER(thetaRanP);  
  //sdVecRanLamb is vector of standard deviations of random intercepts and slopes for component W
  PARAMETER_VECTOR(logSdRanW); vector<Type> sdVecRanLamb = exp(logSdRanW); ADREPORT(sdVecRanLamb);
  //thetaRanW is unstructured correlation between random intercepts and slopes for component W
  PARAMETER(thetaRanW); 
  
  //tauMother is the vector of standard deviations of mother effects for the two components
  PARAMETER_VECTOR(logTauMother); vector<Type> tauMother = exp(logTauMother); ADREPORT(tauMother);
  PARAMETER_VECTOR(epsi); //random mother effect 
  
  int nObs  = NumFle.size();
  int nYear = NLEVELS(t);
  int nMother = NLEVELS(j);
  
  int nPred  = predDates.size();
  
  Type nll=Type(0);
  //through if statement, we turn on or off alpha, theta or omega
  if (CppAD::Variable(logSigma(0)) || CppAD::Variable(logSigma(1)) || CppAD::Variable(logSigma(2)) )
  {
    int d = 3;
    
    UNSTRUCTURED_CORR_t<Type> neg_log_density(parCor);
    // SigmaCor is the correlation matrix of the white noises for VAR1 alpha, theta and omega
    matrix<Type> SigmaCor(d,d);
    SigmaCor=neg_log_density.cov(); 
    
    // Sigma is covariance matrix of the white noises for VAR1 alpha, theta and omega
    matrix<Type> Sigma(d,d); 
    Sigma(0,0)=sigma(0)*sigma(0);
    Sigma(1,1)=sigma(1)*sigma(1);
    Sigma(2,2)=sigma(2)*sigma(2);
    Sigma(0,1)=sigma(0)*sigma(1)*SigmaCor(0,1);
    Sigma(0,2)=sigma(0)*sigma(2)*SigmaCor(0,2);
    Sigma(1,2)=sigma(1)*sigma(2)*SigmaCor(1,2);
    Sigma(1,0)=Sigma(0,1);
    Sigma(2,0)=Sigma(0,2);
    Sigma(2,1)=Sigma(1,2);
    
    vector<Type> vecSigma(d*d);
    
    for (int i=0; i<d; i++)
      for (int j=0; j<d; j++) {
        
        vecSigma(i+j*d) = Sigma(i,j);}
      
    ADREPORT(SigmaCor);
    ADREPORT(vecSigma);
    ADREPORT(Sigma);
    
    matrix<Type> A(d*d,d*d);
    for (int i=0; i<d; i++)
      for (int j=0; j<d; j++)
        for (int k=0; k<d; k++)
          for (int l=0; l<d; l++)
            A(i*d+k, j*d+l) = -PhiB(i,j)*PhiB(k,l);
    for (int i=0; i<d*d; i++)
      A(i,i) += 1;
    matrix<Type> Ainv = A.inverse();
    
    //Gamma0 is the stationary covariance matrix for VAR1 alpha, theta, omega
    vector<Type> vecGamma0 = Ainv*vecSigma;
    matrix<Type> Gamma0(d,d);
    for (int i=0; i<d; i++)
      for (int j=0; j<d; j++){
        Gamma0(i,j) = vecGamma0(i+j*d);}
      ADREPORT(Gamma0);
    matrix<Type> stdGamma0=sqrt(Gamma0.array()); // report standard deviation instead of variance for alpha, theta and omega
    ADREPORT(stdGamma0);
    
    matrix<Type> errorMax(nYear,d);
    for (int t=0; t<nYear; t++){
      errorMax(t,0)= alpha(t);
      errorMax(t,1)= theta(t);
      errorMax(t,2)= omega(t);
    }
    
    nll += MVNORM(Gamma0)(errorMax.row(0)); 
    
    for (int t=1; t<errorMax.rows(); t++) {
      vector<Type> errorMaxt = errorMax.row(t);
      vector<Type> PhierrorMax = errorMax.row(t-1)*PhiB;
      vector<Type> resid = errorMaxt - PhierrorMax;
      
      nll += MVNORM(Sigma)(resid);
    }
  }
  
  //mother effect contribution to likelihood
  for (int j=0; j<nMother; ++j){
    if (CppAD::Variable(epsi(j)))
      nll-=dnorm(epsi(j),Type(0),Type(1),1); //iid N(0,1) female id effects
  }
  
  vector<Type> betaRanP(2);
  vector<Type> betaRanLamb(2);
  
  matrix<Type> CorRanP(2,2);  
  CorRanP(0,0) = 1; CorRanP(0,1) = thetaRanP;
  CorRanP(1,0) = thetaRanP; CorRanP(1,1) = 1;
  
  matrix<Type> CorRanLamb(2,2);  
  CorRanLamb(0,0) = 1; CorRanLamb(0,1) = thetaRanW;
  CorRanLamb(1,0) = thetaRanW; CorRanLamb(1,1) = 1;
  
  
  //contribution of random intercepts and slopes of directional selection to likelihood
  for (int t=0; t<nYear; ++t){
    
    if (CppAD::Variable(logSdRanP(0))){
      
      betaRanP(0)=u0_P(t);
      betaRanP(1)=u1_P(t);
      nll += VECSCALE(MVNORM(CorRanP),sdVecRanP)(betaRanP); }
    
    if (CppAD::Variable(logSdRanW(0))){
      betaRanLamb(0)=u0_W(t);
      betaRanLamb(1)=u1_W(t);
      nll += VECSCALE(MVNORM(CorRanLamb),sdVecRanLamb)(betaRanLamb);}
  }
  
  matrix<Type> etaAlpha(3,nYear);
  matrix<Type> etaTheta(3,nYear);
  matrix<Type> etaOmega(3,nYear);
  
  for (int s=0; s<2; s++){
    for (int t=0; t<nYear; t++){
      
      etaAlpha(s,t) = muAlpha(s);
      etaTheta(s,t) = muTheta(s);
      etaOmega(s,t) = muOmega(s);
      
      if (CppAD::Variable(logSigma(0)))
        etaAlpha(s,t) += alpha(t);
      
      if (CppAD::Variable(logSigma(1)))
        etaTheta(s,t) += theta(t);
      
      if (CppAD::Variable(logSigma(2)))
        etaOmega(s,t) += omega(t);
    }
  }
  
  vector<Type> prob_fun(nObs);
  vector<Type> poi_lamb_fun(nObs);
  
  vector<Type> prob(nObs);
  vector<Type> poi_lamb(nObs);
  
  //either directional selection or stabilizing selection on each episode
  for (int i=0; i<nObs; ++i){
    
    if (CppAD::Variable(muAlpha(0))) //turn stabilizing selection on
      prob_fun(i) = etaAlpha(0,t(i))
      -pow((LayDateApril(i)-etaTheta(0,t(i))),2)/(2*pow(exp(etaOmega(0,t(i))),2));
    if (CppAD::Variable(beta0_P)) //turn directional selection on with only fixed effect
      prob_fun(i) = beta0_P;
      if (CppAD::Variable(beta1_P)) //turn directional selection on with only fixed effect
        prob_fun(i) +=  beta1_P*LayDateApril(i);
    if (CppAD::Variable(logSdRanP(0)))
      prob_fun(i) +=  u0_P(t(i)); //add random intercepts to the directional selection
    if (CppAD::Variable(logSdRanP(1)))
      prob_fun(i) += u1_P(t(i))*LayDateApril(i); //add random slopes to the directional selection
    if (CppAD::Variable(logTauMother(0))) //turn the mother effect on
      prob_fun(i) += tauMother(0)*epsi(j(i));
    
    
    if (CppAD::Variable(muAlpha(1)))
      poi_lamb_fun(i) = etaAlpha(1,t(i))
      -pow((LayDateApril(i)-etaTheta(1,t(i))),2)/(2*pow(exp(etaOmega(1,t(i))),2));
    if (CppAD::Variable(beta0_W))
      poi_lamb_fun(i) = beta0_W + beta1_W*LayDateApril(i);
    if (CppAD::Variable(logSdRanW(0)))
      poi_lamb_fun(i) +=  u0_W(t(i)) ;
    if (CppAD::Variable(logSdRanW(1)))
      poi_lamb_fun(i) += u1_W(t(i))*LayDateApril(i);
    if (CppAD::Variable(logTauMother(1)))
      poi_lamb_fun(i) += tauMother(1)*epsi(j(i));
    
 
    prob(i) = 1/(1+exp(prob_fun(i)));
    
    poi_lamb(i) = exp(poi_lamb_fun(i));

    nll-= dzipois(NumFle(i), poi_lamb(i), prob(i), true);

  }
  
  
  if (makePred == 1){
    
    vector<Type> prob_pred(nPred);
    vector<Type> p_pred(nPred);
    
    vector<Type> poi_lamb_pred(nPred);
    vector<Type> lamb_pred(nPred);
    
    vector<Type> p_nonzeroy(nPred);
    vector<Type> expect_nonzeroy(nPred);
    
    vector<Type> expect_y(nPred);
    
    for (int i=0; i< nPred; ++i){
    //ignore the mother effects and make predictions
    prob_pred(i) = beta0_P + beta1_P*predDates(i) +  u0_P(tPred(i)) + u1_P(tPred(i))*predDates(i);
    p_pred(i) = 1/(1+exp(prob_pred(i)));
    
    poi_lamb_pred(i) = etaAlpha(1,tPred(i))
      -pow((predDates(i)-etaTheta(1,tPred(i))),2)/(2*pow(exp(etaOmega(1,tPred(i))),2));
    
    lamb_pred(i) = exp(poi_lamb_pred(i));
    
    //calculate p(Y|>0)
    p_nonzeroy(i) =  (1-p_pred(i))*(1-exp(-lamb_pred(i)));
    
    //calculate E(Y|Y>0)
    expect_nonzeroy(i) =  lamb_pred(i)/(1-exp(-lamb_pred(i)));
    
    //calculate E(Y|w,p)=(1-p)*w
    expect_y(i) =  (1-p_pred(i))*lamb_pred(i);
    
      }
    ADREPORT(p_nonzeroy);
    ADREPORT(expect_nonzeroy);
    ADREPORT(expect_y);
    }

    ADREPORT(prob);
    ADREPORT(poi_lamb);
    ADREPORT(exp(muOmega));
  return nll;
}








