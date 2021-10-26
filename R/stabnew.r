#' Check-based yield stability analysis
#' @param y	A response varible vector used for stability analysis
#' @param Gen	A vector of genotypes.
#' @param Env	A vector of environments.
#' @param times	Times of resampling used for stability analysis.
#' @param check	One or more checks used for stability analysis.
#' @param Rep	An argument with replication: Rep=TRUE or with replication: Rep=FALSE
#' @param X	A vector or matrix of other predictable variables. Default is NULL.
#' @param alpha	A nominal probability values used for statistical tests. Default is NULL, 0.05
#' @return A list of yield stability results
#' @references
#' Finlay, K.W., G.N. Wilkinson 1963. The analysis of adaptation in a plant breeding programme. Australian Journal of Agricultural Research 14: 742-754.
#'
#' Wu, J., K. Glover, W. Berzonsky, 2012. Statistical tests for stability analysis with resampling techniques. 25th Conference of Applied Statistics in Agriculture. p88-108. April 29-May 01, 2012. Manhattan, KS
#'
#' Wu, J., K. Glover, and N. Mueller 2014. Check based stability analysis method and its application to winter wheat variety trials," Conference on Applied Statistics in Agriculture. P102-114. https://doi.org/10.4148/2475-7772.1006
#' @examples
#' data(maize)
#' #names(maize)
#' Geno=as.vector(maize$Cultivar)
#' Env=paste(maize$Location,maize$Year,sep=":")
#' y=maize$Yld
#' res=stab.fw.check(y,Gen=Geno,Env=Env,times=10,check=c("Hai He"),Rep=FALSE)
#' res
##end
#' @export

stab.fw.check=function(y,Gen,Env,times,check,Rep,X=NULL,alpha=NULL){
  if(is.null(alpha))alpha=0.05
  if(is.null(X))X=NULL
  if(Rep==TRUE){
    dat=GetGEMean(y,Gen,Env,X)
    y1=dat$y
    if(is.null(X))X1=NULL
    else X1=dat[,-c(1:3)]
    Gen1=dat$Gen
    Env1=dat$Env
    result=fw.check(y1,Gen1,Env1,check,times,X1,alpha)
  }
  else result=fw.check(y,Gen,Env,times,check,X,alpha)
  return(result)
}


## F-W based stability analysis
#' F-W Regression Based Yield Stability Analysis
#' @param y A vector of yield data
#' @param Gen A vector of Genotypes
#' @param Env A vector of Environments
#' @param times Replication number for resampling
#' @param Rep Replication included or not included
#' @param X Independent variables matrix or vector
#' @param alpha Preset alpha value
#' @return A list of yield stability results
#'
#' @references
#' Finlay, K.W., G.N. Wilkinson 1963. The analysis of adaptation in a plant breeding programme. Australian Journal of Agricultural Research 14: 742-754.
#'
#' Wu, J., K. Glover, W. Berzonsky, 2012. Statistical tests for stability analysis with resampling techniques. 25th Conference of Applied Statistics in Agriculture. p88-108. April 29-May 01, 2012. Manhattan, KS
#'
#' Wu, J., K. Glover, and N. Mueller 2014. Check based stability analysis method and its application to winter wheat variety trials," Conference on Applied Statistics in Agriculture. https://doi.org/10.4148/2475-7772.1006

#' @examples
#' require(genstab)
#' data(maize)
#' #names(maize)
#' Geno=as.vector(maize$Cultivar)
#' Env=paste(maize$Location,maize$Year,sep=":")
#' y=maize$Yld
#' res=stab.fw(y,Gen=Geno,Env=Env,times=10,Rep=TRUE)
#' res
#' ##end
#' @export

stab.fw=function(y,Gen,Env,times,Rep,X=NULL,alpha=NULL){
  if(is.null(alpha))alpha=0.05
  if(is.null(X))X=NULL
  if(Rep==TRUE){
    dat=GetGEMean(y,Gen,Env,X)
    y1=dat$y
    if(is.null(X))X1=NULL
    else X1=dat[,-c(1:3)]
    Gen1=dat$Gen
    Env1=dat$Env
    result=fw_reg(y1,Gen1,Env1,times,X1,alpha)
  }
  else result=fw_reg(y,Gen,Env,times,X,alpha)
  return(result)
}


#' Group variances with resampling
#' @description Group variance calculation with two resampling techniques:permuation and bootstraping
#' @param Y	A matrix including One or more traits
#' @param class	A vector of the first factor for calculating variance. For example, a vector of genotypes.
#' @param cls2	A vector of the second factor used within-group bootstraping for variance. It can be default
#' @param resample	Resampling technique option. resample="Boot" is for bootstrapping. resample="Perm" is for permutation.
#' @param times	Number of resampling used. The default number is 1000.
#' @param alpha	A nomimal probability used for statistical test. The default value is 0.05.
#' @return A list of variances and confidence intervals for genotypes or environments
#' @author Jixiang Wu <jixiang.wu@sdstate.edu>
#' @references
#' Finlay, K.W., G.N. Wilkinson 1963. The analysis of adaptation in a plant breeding programme. Australian Journal of Agricultural Research 14: 742-754.
#'
#' Wu, J., K. Glover, W. Berzonsky, 2012. Statistical tests for stability analysis with resampling techniques. 25th Conference of Applied Statistics in Agriculture. p88-108. April 29- May 01, 2012. Manhattan, KS
#' @examples
#' data(maize)
#' #names(maize)
#' Geno=as.vector(maize$Cultivar)
#' Env=paste(maize$Location,maize$Year,sep=":")
#' y=maize$Yld
#' res=stab.var(y,class=Geno,cls2=Env,resample="Boot",times=100)
#' res
#' res=stab.var(y,class=Geno,resample="Perm",times=100)
#' res
#' @export

stab.var=function(Y,class,cls2=NULL,resample,times=NULL,alpha=NULL){
  Y=as.matrix(Y)
  tn=ncol(Y)
  if(is.null(alpha))alpha=0.05
  if(is.null(times))times=1000
  if(is.null(cls2))cls2=NULL
  if(tn==1)return(Var2(Y[,1],class,cls2,resample,times,alpha))
  else{
    RES=list()
    for(i in 1:tn)RES[[i]]=Var2(Y[,i],class,cls2,resample,times,alpha)
    names(RES)=colnames(Y)
    return(RES)
  }
}
#' Group means and ranks with resampling
#' @description Group mean and rank calculation with two resampling techniques:permuation and bootstraping
#' @param Y	A matrix including One or more traits
#' @param class	A vector of the first factor for calculating variance. For example, a vector of genotypes.
#' @param cls2	A vector of the second factor used within-group bootstraping for variance. It can be default
#' @param resample	Resampling technique option. resample="Boot" is for bootstrapping. resample="Perm" is for permutation.
#' @param times	Number of resampling used. The default number is 1000.
#' @param alpha	A nomimal probability used for statistical test. The default value is 0.05.
#' @return A list of variances and confidence intervals for genotypes or environments
#' @author Jixiang Wu <jixiang.wu@sdstate.edu>
#' @references
#' Finlay, K.W., G.N. Wilkinson 1963. The analysis of adaptation in a plant breeding programme. Australian Journal of Agricultural Research 14: 742-754.
#'
#' Wu, J., K. Glover, W. Berzonsky, 2012. Statistical tests for stability analysis with resampling techniques. 25th Conference of Applied Statistics in Agriculture. p88-108. April 29- May 01, 2012. Manhattan, KS
#' @examples
#' data(maize)
#' #names(maize)
#' Geno=as.vector(maize$Cultivar)
#' Env=paste(maize$Location,maize$Year,sep=":")
#' y=maize$Yld
#' res=stab.mean(y,class=Geno,cls2=Env,resample="Boot",times=100)
#' res
#' res=stab.mean(y,class=Geno,resample="Perm",times=100)
#' res
#' @export

stab.mean=function(Y,class,cls2=NULL,resample,times=NULL,alpha=NULL){
  Y=as.matrix(Y)
  tn=ncol(Y)
  if(is.null(alpha))alpha=0.05
  if(is.null(times))times=1000
  if(is.null(cls2))cls2=NULL
  if(tn==1)return(Rank2(Y[,1],class,cls2,resample,times,alpha))
  else{
    RES=list()
    for(i in 1:tn)RES[[i]]=Rank2(Y[,i],class,cls2,resample,times,alpha)
    names(RES)=colnames(Y)
    return(RES)
  }
}


reg.boot=function(y,X,times=NULL,boot=NULL,ALPHA=NULL){
  if(is.null(boot))boot="original"
  if(is.null(ALPHA))ALPHA=0.05
  if(is.null(times))times=1000
  p1=ALPHA/2
  p2=1-p1

  dat=data.frame(y,X)
  dat=na.omit(dat)
  #ok=complete.cases(dat)
  #dat=dat[ok,]
  n=length(dat$y)
  mod0=lm(y~.,data=dat)
  #attributes(mod0)
  r0=summary(mod0)$r.squared
  res0=mod0$residual
  yhat0=mod0$fitted.values
  b0=mod0$coefficients
  nc=length(b0)
  B=matrix(0,times,nc)
  r1=numeric(times)
  if(boot=="original"){
    bootmethod="boot_original"
    for(i in 1:times){
      id=sample(1:n,replace=TRUE)
      dat1=dat[id,]
      mod1=lm(y~.,data=dat1)
      B[i,]=mod1$coefficients
      r1[i]=summary(mod1)$r.squared
    }
  }

  else if(boot=="residual"){
    bootmethod="boot_residual"
    X1=dat[,-1]
    for(i in 1:times){
      id=sample(1:n,replace=T)
      y1=yhat0+res0[id]
      dat1=data.frame(y1,X1)
      #dat1=dat[id,]
      mod1=lm(y1~.,data=dat1)
      B[i,]=mod1$coefficients
      r1[i]=summary(mod1)$r.squared
    }
  }
  B=data.frame(B,r1)

  p=numeric(nc)
  classnames=c(names(b0),"r.squared")
  CI=matrix(0,nc+1,2)
  b1=numeric(nc+1)
  b0=c(b0,r0)
  for(i in 1:(nc+1)){
    #i=1
    v=B[,i]
    v=na.omit(v)
    #
    #ok=complete.cases(v)
    #v=v[ok]
    tm=length(v)
    CI[i,]=quantile(v,p=c(p1,p2),na.rm=TRUE)
    b1[i]=mean(v)
    if(b0[i]==0)p[i]=1.000
    else if(b0[i]>0)p[i]=length(which(v<0))/tm
    else if(b0[i]<0)p[i]=length(which(v>0))/tm
  }
  p=p*2
  mean=data.frame(b0,b1,p,CI)
  colnames(mean)=c("Orig","boot","Prob","LL","UL")
  rownames(mean)=classnames
  res=list(parameters=mean,alpha=ALPHA,bootmethod=bootmethod)
  return(res)
}

GetGEMean=function(y,Gen,Env,X=NULL){
  GE=tapply(y,list(Env,Gen),mean)
  GE=as.data.frame.table(GE)
  colnames(GE)=c("Env","Gen","y")
  if(is.null(X)==FALSE){
    X=as.matrix(X)
    nx=ncol(X)
    for(i in 1:nx){
      GE1=tapply(X[,i],list(Env,Gen),mean)
      GE1=as.data.frame.table(GE1)
      GE=cbind(GE,GE1[,3])
    }
    colnames(GE)=c("Env","Gen","y",colnames(X))
  }
  dat=data.frame(GE)
  dat=na.omit(dat)
  #ok=complete.cases(dat)
  #dat=dat[ok,]
  return(dat)
}

fw_reg=function(y,Gen,Env,times,X=NULL,alpha=NULL){
  if(is.null(alpha))alpha=0.05
  p1=alpha/2
  p2=1-p1
  if(is.null(X)){
    dat1=data.frame(Gen,Env,y)
    colnames(dat1)=c("Gen","Env","y")
  }
  else {
    dat1=data.frame(Gen,Env,y,X)
    colnames(dat1)=c("Gen","Env","y",colnames(X))
  }
  dat1=dat1[order(dat1$Gen,dat1$Env),]
  #return(dat1)
  #
  EI=tapply(dat1$y,dat1$Env,mean)
  name0=names(EI)
  gnames=unique(Gen)
  gn=length(unique(Gen))
  r.square=numeric(gn)
  b=numeric(gn)

  DAT=list()
  for(i in 1:gn){
    #i=1
    cat("Interation=",i,"\n")
    id=which(dat1$Gen==gnames[i])
    len=length(EI)
    #if(is.null(X)==F){
    if(length(id)==len)dat2=data.frame(dat1[id,],EI)
    else if(length(id)<length(EI)){
      ei=dat1$y[id]
      name1=dat1$Env[id]
      #name1=names(ei)
      intername=intersect(name1,name0)
      id1=numeric(length(intername))
      id2=numeric(length(intername))
      for(j in 1:length(id1)){
        id0=which(name1==intername[j])
        id1[j]=id0
        id0=which(name0==intername[j])
        id2[j]=id0
      }
      ei1=ei[id1]
      ei2=EI[id2]
      dat2=data.frame(dat1[id,],ei2)
    }
    #}
    colnames(dat2)=c(colnames(dat1),"EI")
    DAT[[i]]=data.frame(dat2)
  }

  RES=list()
  if(is.null(X)){
    for(i in 1:gn){
      dat2=DAT[[i]]
      y1=dat2[,3]
      X1=dat2[,-c(1:3)]
      RES[[i]]=reg.boot(y1,X1,times,boot="residual")
    }
  }
  else{
    for(i in 1:gn){

      dat2=DAT[[i]]
      nc=ncol(dat2)
      #dat2=dat2[,-nc]
      #y1=dat2[,3]
      #X1=dat2[,-c(1:3)]
      y1=dat2[,nc]
      X1=dat2[,c(4:(nc-1))]

      RES[[i]]=reg.boot(y1,X1,times,boot="residual")
    }
  }
  names(RES)=gnames
  return(RES)
}

Rank2=function(y,class,cls2=NULL,resample,times=NULL,alpha=NULL){
  if(is.null(alpha))alpha=0.05
  if(is.null(times))times=1000
  p1=alpha/2
  p2=1-p1
  m0=tapply(y,class,mean,na.rm=TRUE)
  nc=length(m0)
  Mean=matrix(0,times,nc)
  n=length(y)
  if(is.null(cls2))cls2=rep(1,n)
  cnames=unique(cls2)
  cn=length(cnames)
  for(i in 1:times){
    #i=1
    if(resample=="Boot"){
      y1=NULL
      class1=NULL
      for(j in 1:cn){
        id=which(cls2==cnames[j])
        y0=y[id]
        c0=class[id]
        n0=length(y0)
        id0=sample(n0,replace=TRUE)
        y1=c(y1,y0[id0])
        class1=c(class1,c0[id0])
      }
      #
      #
      #id=sample(n,replace=T)
      #y1=y[id]
      #class1=class[id]
    }
    else{
      id=sample(n,replace=FALSE)
      y1=y[id]
      class1=class
    }
    m1=tapply(y1,class1,mean,na.rm=TRUE)
    if(length(m1)<nc)m1=vectcomp(m0,m1)
    Mean[i,]=m1  ##tapply(y1,class1,mean,na.rm=T)
  }
  classnames=names(m0)
  #cn=length(m)
  CI=matrix(0,nc,2)
  m1=numeric(nc)
  #dim(Mean)
  for(i in 1:nc){
    #i=1
    v=Mean[,i]
    CI[i,]=quantile(v,p=c(p1,p2),na.rm=TRUE)
    m1[i]=mean(v)
  }
  mean=data.frame(m0,m1,CI)
  colnames(mean)=c("Orig",resample,"LL","UL")
  rownames(mean)=classnames

  order0=order(m0)
  rank0=numeric(nc)
  rank0[order0]=1:nc
  Rank1=matrix(0,times,nc)
  v=numeric(nc)
  for(i in 1:times){
    id=order(Mean[i,])
    v[id]=1:nc
    Rank1[i,]=v
  }
  CI=matrix(0,nc,2)
  order1=numeric(nc)
  for(i in 1:nc){
    v=Rank1[,i]
    CI[i,]=quantile(v,p=c(p1,p2),na.rm=TRUE)
    order1[i]=mean(v)
  }

  rnk=data.frame(rank0,order1,CI)
  colnames(rnk)=c("Orig",resample,"LL","UL")
  rownames(rnk)=classnames
  result=list(mean=mean,rank=rnk,alpha=alpha)
  return(result)
}


Var2=function(y,class,cls2=NULL,resample,times=NULL,alpha=NULL){
  if(is.null(alpha))alpha=0.05
  if(is.null(times))times=1000
  p1=alpha/2
  p2=1-p1
  m0=tapply(y,class,var,na.rm=TRUE)
  nc=length(m0)
  Mean=matrix(0,times,nc)
  n=length(y)
  if(is.null(cls2))cls2=rep(1,n)
  cnames=unique(cls2)
  cn=length(cnames)
  for(i in 1:times){
    if(resample=="Boot"){
      y1=NULL
      class1=NULL
      for(j in 1:cn){
        id=which(cls2==cnames[j])
        y0=y[id]
        c0=class[id]
        n0=length(y0)
        id0=sample(n0,replace=TRUE)
        y1=c(y1,y0[id0])
        class1=c(class1,c0[id0])
      }
      #
      #
      #id=sample(n,replace=T)
      #y1=y[id]
      #class1=class[id]
    }
    else{
      id=sample(n,replace=FALSE)
      y1=y[id]
      class1=class
    }
    m1=tapply(y1,class1,var,na.rm=TRUE)
    if(length(m1)<nc)m1=vectcomp(m0,m1)
    Mean[i,]=m1  ##tapply(y1,class1,mean,na.rm=T)
  }
  classnames=names(m0)
  #cn=length(m)
  CI=matrix(0,nc,2)
  m1=numeric(nc)
  #dim(Mean)
  for(i in 1:nc){
    #i=1
    v=Mean[,i]
    CI[i,]=quantile(v,p=c(p1,p2),na.rm=TRUE)
    m1[i]=mean(v)
  }
  mean=data.frame(m0,m1,CI)
  colnames(mean)=c("Orig",resample,"LL","UL")
  rownames(mean)=classnames
  result=list(Var=mean,alpha=alpha)
  return(result)
}



fw.check=function(y,Gen,Env,times,check,X=NULL,alpha=NULL){
  if(is.null(alpha))alpha=0.05
  p1=alpha/2
  p2=1-p1
  y=as.vector(y)
  Env=as.vector(Env)
  Gen=as.vector(Gen)
  if(is.null(X)){
    dat1=data.frame(Gen,Env,y)
    colnames(dat1)=c("Gen","Env","y")
  }
  else {
    dat1=data.frame(Gen,Env,y,X)
    colnames(dat1)=c("Gen","Env","y",colnames(X))
  }
  dat1=dat1[order(dat1$Gen,dat1$Env),]

  id=get.check.data(dat1$Gen,check)
  EI=tapply(dat1$y[id],dat1$Env[id],mean) ##get check based EI
  name0=names(EI)
  Gen1=as.vector(dat1$Gen[-id])

  gnames=unique(Gen1)
  gn=length(unique(Gen1))
  r.square=numeric(gn)
  b=numeric(gn)
  #name0=names(EI)
  DAT=list()
  for(i in 1:gn){
    #i=2
    cat("Interation=",i,"\n")
    id=which(dat1$Gen==gnames[i])
    len=length(EI)
    #if(is.null(X)==F){
    if(length(id)==len)dat2=data.frame(dat1[id,],EI)
    else if(length(id)<length(EI)){
      ei=dat1$y[id]
      name1=dat1$Env[id]
      #length(name1)
      #name1=names(ei)
      intername=intersect(name1,name0)
      id1=numeric(length(intername))
      id2=numeric(length(intername))
      for(j in 1:length(id1)){
        id0=which(name1==intername[j])
        id1[j]=id0
        id0=which(name0==intername[j])
        id2[j]=id0
      }
      ei1=ei[id1]
      ei2=EI[id2]
      dat2=data.frame(gnames[i],name1,ei1,ei2)
    }
    #}
    colnames(dat2)=c(colnames(dat1),"EI")
    DAT[[i]]=data.frame(dat2)
  }
  ##DAT[[1]]


  #DAT[[1]]
  #ww03[2,]
  RES=list()
  for(i in 1:gn){
    #i=1
    dat2=DAT[[i]]
    y1=dat2[,3]
    X1=as.matrix(dat2[,-c(1:3)])
    colnames(X1)=colnames(dat2)[-c(1:3)]
    RES[[i]]=reg.boot(y1,X1,times,boot="residual")
  }
  names(RES)=gnames
  return(RES)
}


get.check.data=function(Genotype,GenComNames,...){
  gcn=length(GenComNames)
  for(i in 1:gcn){
    if(i==1)id=which(Genotype==GenComNames[i])
    else id=c(id,which(Genotype==GenComNames[i]))
  }
  return(id)
}

vectcomp=function(v0,v1){
  name0=names(v0)
  name1=names(v1)
  vn=length(v0)
  v2=numeric(vn)
  for(i in 1:vn){
    id=which(name0[i]==name1)
    if(length(id)>0)v2[i]=v1[id]
  }
  names(v2)=name0
  return(v2)
}


#' Maize yield trial data
#' @references
#' Fan X.M., Kang M.S., Chen H.M., Zhang Y.D., Tan J., Xu C.X. (2007) Yield stability of maize hybrids evaluated in multi-environment trials in Yunnan, China. Agronomy Journal.99:220-228
#'
#'@examples
#'str(maize)
"maize"
