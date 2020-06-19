#RScript to pre-process data to run the Bayesian Hierarchical Mixture model, scripted by Brian Leung, Jun 19, 2020 (brian.leung2@mcgill.ca)

#CLEAN VERSION
#args <- commandArgs(trailingOnly = TRUE)
#
print("For first run, need to compile STAN - make sure you set initiate_stan: in terminal: init_stan=T, then turn off afterwards: init_stan=F. Do not need to do it again, unless change STAN code")
library(rstan)
library(foreach)
library(doParallel)
if(!exists("init_stan"))
	init_stan=F

registerDoParallel()
options(cores=30) #set number of cores to use, depends on how many computer has
options(mc.cores=1)
rstan_options(auto_write = TRUE)

#use info and growth_rates, separate by system. Select all ones where ts>2, and where sd>0. apply stan. For ts=2, fit r & sd and use as priors.
#what does the mean growth rate indicate though?
iter=10000
warmup=3000
chains=3
nr=chains*(iter-warmup)
info_nm="info.rds"
gr1<-readRDS("growth_rates_example.rds") 
info<-readRDS(info_nm)
info=info[which(info$id %in% names(gr1)),] #ensure data sets match

info$combo=as.character(info$combo)
u<-unique(info$combo)

nrep=length(u)
tr=1:4 #testing all combos of the clusters
#**********FUNCTIONS

get_dat<-function(rws)
{
	pos=which(info$combo==rws)
	id=as.character(info$id[pos])
	gr=gr1[id]
	
	pop=unlist(lapply(gr,length),use.names=F)
	v<-unlist(lapply(gr,var),use.names=F)
	s1=(sum(v*(pop-1),na.rm=T)/sum((pop-1)))^.5 #pooled SD

	sdgr=s1/(pop)^.5 #get SE - i.e., distr. of sample means
	y=unlist(lapply(gr,mean),use.names=F)
	om<-which(is.na(y))
	if(length(om)>0)
	{
		y=y[-om]
		sdgr=sdgr[-om]
		pop=pop[-om]
	}
	J=length(y)
	l=list(y=y,J=J,pop=pop, s=s1,sigma1=sdgr)
	return(l)
}
DIC<-function(theta,tau,frac,sigma, tK,l1)
{
	Dm=matrix(0,nrow=nr,ncol=l1$J)
	for(j in 1:l1$J)
	{
		if(option==1)
		{
			s=sigma/(l1$pop[j]^.5)
		}else{
			s=l1$sigma1[j]
		}

		sd1=(tau^2+s^2)^.5 #tau is zero if K is zero
		if(tK<2)
		{
			Dm[,j]=dnorm(l1$y[j]-theta,sd=sd1)
		}else{
			Dm[,j]=apply(frac*dnorm(l1$y[j]-theta,sd=sd1),1,sum)
		}
	}

	Dm1=-2*apply(log(Dm),1,sum)
	DEp1=var(Dm1)/2
	return(c(mean(Dm1),DEp1,mean(Dm1)+DEp1))

}
do_fit<-function(fit,l1,K,rep)
{
	#hyper mean
	theta <- matrix(unlist(extract(fit, 'theta'), use.names=F), ncol=K,nrow=nr)
	tau <- matrix(unlist(extract(fit, 'tau'),use.names=F),ncol=K,nrow=nr)
	
	if(K>1)
	{
		frac<-matrix(unlist(extract(fit, 'frac'),use.names=F),ncol=K,nrow=nr)
	}else{
		frac=1
	}
	ss<-summary(fit)$summary
	rn<-names(ss[,1])
	rpos=grep("theta\\[",rn)
	rpos=c(rpos,grep("tau\\[",rn))
	rpos=c(rpos,grep("frac\\[",rn))
	rn=rn[rpos]
	summary<-ss[rpos,c(1,3,4,8,10)] #include columns we're interested in from summary table
	cn=colnames(summary)
	snm=paste(rep(cn,each=length(rn)),rep(rn, length(cn)),sep="_")
	snm=gsub("\\[","_",snm)
	snm=gsub("\\]","",snm)
	snm=gsub("\\%","",snm)

	snm=paste(snm,K,sep="_")
	ret<-as.vector(summary)
	
	D=DIC(theta,tau,frac,sigma, K,l1)
	ret=c(ret,D)
	snm=c(snm,paste(c("lk","df","dic"),K,sep="_"))
	names(ret)=snm
	return(ret)
}

ll<-list()
for(rep in 1:nrep)
{
	ll[[rep]]<-get_dat(u[rep])
}
if(init_stan==T)
{ #running once, to make sure stan program is compiled, before parellization
	print("initiating")
	l=ll[[1]]
	l$K=1
	l$bound1=-10
	l$bound2=mean(l$y)-sd(l$y)
	l$bound3=mean(l$y)+sd(l$y)
	fit=stan(file='BHM.stan', data=l, iter=iter, warmup=warmup, chains=chains,control=list(adapt_delta=.9))
	end() #then need to re-run again
}
print("running now")
tmp<-foreach(rep = 1:nrep) %dopar%
{
	l<-ll[[rep]]
	if(length(l$y)==0)
		return(NULL)
	ret=c(system=u[rep], pop=length(l$y))
	for(k in tr) #iterate through all K
	{
		l$K=k
		l$bound3=mean(l$y)+sd(l$y)
		l$bound1=-10
		l$bound2=mean(l$y)-sd(l$y)
		if(k==4)#change this later - this is to check the upper bound only
		{
			l$K=2
			l$bound1=l$bound3
			l$bound2=10
		}		
		sink("tmp_stanoutput") #we don't want the output
		fit=stan(file='BHM.stan', data=l, iter=iter, warmup=warmup, chains=chains,control=list(adapt_delta=.9))
		sink()
		ret=c(ret,do_fit(fit,l,l$K,rep))
	}
	print(rep)
	return(ret)	
} #return complex output file for later analysis
system('rm tmp_stanoutput')

dat=t(matrix(unlist(tmp),nrow=length(tmp[[1]]),ncol=length(tmp)))
dat<-as.data.frame(dat)
names(dat)=names(tmp[[1]])
#because K=4 gets converted to K=2, must change back
pos=(which(names(dat)=="dic_3")+1):ncol(dat)
names(dat)[pos]=gsub("_2","_4",names(dat)[pos])
names(dat)[pos]=gsub("4_4","2_4",names(dat)[pos])
for(i in 2:ncol(dat))
{
	dat[,i]=as.numeric(as.character(dat[,i]))
}

pos=grep("dic",names(dat))
if(length(pos)>1) 
{
	dat$which_min<-apply(dat[,pos],1,which.min) #gives you a value order 1 to 4, for which scenario is lowest
	ret_vals<-function()
	{
		col=grep("mean_theta_1",names(dat))
		v<-rep(NA,nrow(dat))
		for(i in 1:nrow(dat))
			v[i]=dat[i,col[dat$which_min[i]]]

		return(v)
	}
	dat$min_main=ret_vals()
}

saveRDS(dat,"dat.rds",version=2) 
print("NOTE: Check Rhat values. if Rhat deviates substantially from one, convergence is suspect. Please re-run analysis for those systems")

