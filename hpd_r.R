n=3
x=c(2,6,1)
sum_x=sum(x)
alpha=2
beta=5
alpha_post=n+alpha #posterior alpha
beta_post=sum_x-n+beta #posterior beta

highest_posterior_density_beta<-function(a,b,target){
  max=optimize(f=function(x) dbeta(x,a,b),lower=0,upper=1,maximum=TRUE)
  max_x=max$maximum
  max_func=max$objective
  
  search_grid=rev(seq(0,max_func,by=0.00001))
  
  results=matrix(NA,nrow=length(search_grid),ncol=4,
                 dimnames=list(c(search_grid),c('Lower','Upper','Integrate','Diff')))
  
  func_for_roots<-function(x,y,a,b){dbeta(x,a,b)-y}
  
  for(i in 1:length(search_grid)){
    results[i,1]=uniroot(func_for_roots,y=search_grid[i],a=a,b=b,lower=0,upper=max_x)$root
    results[i,2]=uniroot(func_for_roots,y=search_grid[i],a=a,b=b,lower=max_x,upper=1)$root
    results[i,3]=integrate(f=function(x,a,b)dbeta(x,a,b),a=a,b=b,
                           lower=results[i,1],upper=results[i,2])$value
    results[i,4]=abs(results[i,3]-target)
  }
  
  results_new=cbind(results,search_grid)
  
  output=results_new[results_new[,4]==min(results_new[,4])]
  list(lower=output[1],upper=output[2],int=output[3],diff=output[4],density=output[5])
}
hpd=highest_posterior_density_beta(alpha_post,beta_post,0.95)
hpd_int=c(hpd$lower,hpd$upper)
hpd_int

ggplot(data=data.frame(x=seq(0,1,length.out=1000)),aes(x))+
  stat_function(fun=function(x) dbeta(x,n+alpha,sum_x-n+beta),size=1.5)+
  geom_vline(xintercept=hpd$lower,col='red',size=1.5,linetype=2)+
  geom_vline(xintercept=hpd$upper,col='red',size=1.5,linetype=2)+
  geom_hline(yintercept=hpd$density,col='blue',size=1.5)+
  labs(title='HPD',x='X',y='Density')+
  annotate('text',x=0.85,y=0.5,label='HPD Interval Density Value',size=5)+
  annotate('text',x=0.55,y=3,label='95% HPD Interval',size=5)