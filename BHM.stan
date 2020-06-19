//STAN code to run the Bayesian Hierarchical Mixture model, scripted by Brian Leung, Jun 19, 2020
data {
	int J; //the number of groups
	int K;//the number of clusters 1-3
	real y[J]; // Observations - group means
	real<lower=0> sigma1[J]; //within group SE 
	real bound1; 
	real bound2;
	real bound3;	
}
parameters {

	real<lower=0, upper=10> tau[K]; //variance of all pops: hyper-parameter

	real<lower=-10,upper=10> theta1; //mean of all pops: hyper-parameter
	real<lower=bound1,upper=bound2> theta2[(K>1) ? 1:0]; 
	real<lower=bound3,upper=10> theta3[(K==3) ? 1:0]; 
	real<lower=0,upper=.49> frac2[(K>1) ? 1:0]; //fraction for each element. must sum to 1 frac, betw 0,1, frac between 0, 1-k
	real<lower=0,upper=.49> frac3[(K==3) ? 1:0]; 

}

transformed parameters {
	real theta[K];
	real frac[K];
	if(K>0)
	{
		theta[1]=theta1;
		frac[1]=1;
	}
	if(K>1)
	{
		theta[2]=theta2[1];
		frac[2]=frac2[1];
		frac[1]=1-frac[2];
	}
	if(K==3)
	{
		theta[3]=theta3[1];
		frac[3]=frac3[1];
		frac[1]=frac[1]-frac3[1];
	}	


}


model {

	if(K>1)
	{
		real log_frac[K];
		log_frac=log(frac);
		for(j in 1:J)
		{
			vector[K] ls;
			for(k in 1:K)
			{
				real s;
				s=(tau[k]^2+sigma1[j]^2)^.5;
				ls[k]=log_frac[k]+normal_lpdf(y[j] | theta[k],s);
			}
			target+=log_sum_exp(ls);
	
		}
	}	
	if(K==1)
	{
		vector[J] s;
		for(j in 1:J)
		{
			s[j]=(tau[1]^2+sigma1[j]^2)^.5;
		}
		y~normal(theta1,s);
	}
	for(k in 1:K) //priors
	{
		target+= -0.5*log(tau[k]);
		target += -0.5*log(frac[k]);
	}
	
}

