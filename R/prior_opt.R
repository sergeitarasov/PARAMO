pgamma(2, shape=1.8, scale=1.5)
qgamma(c(.1, .9), shape=1.8, scale=1.5)

rexp(1000, rgamma(1000,shape=.383, scale=274)) %>% hist(., breaks=150)
rexp(1000, rgamma(1000,shape=.383, scale=274)) %>% round(3) %>% mean()
rexp(1000, rgamma(1000,shape=0.3763443, scale=1.8618148)) %>% mean()
rgamma(1000,shape=1, scale=.85) %>% hist(., breaks=150)
rgamma(1000,shape=0.3231829, scale=5.2959589) %>% round(2)

rexp(10000, 50) %>% hist(., breaks=150)
rexp(1000, .2) %>% round(2) %>% min()
rinvgamma(10000, shape = 2.2, rate=0.024) %>% hist(., breaks=150)
rinvgamma(10000, shape = 22, rate=0.4) %>% mean()

rinvgamma(10000, shape = 2.2, rate=0.024) %>% hist(., breaks=150)
qinvgamma(c(.025, .975), shape = 2.2, rate=0.024)

rinvgamma(10000, shape = 2.01, rate=0.02) %>% hist(., breaks=150)
qinvgamma(c(.025, .975), shape = 2.01, rate=0.02)
rinvgamma(1000, shape = 3, rate=0.0001) %>% hist(., breaks=150)
rinvgamma(1000, shape = 2, rate=0.0001) %>% var()

1/.01=10
rexp(1000, 300) %>% hist(., breaks=150)
1/100=0.01
1/14

qexp(c(.025, .975), 500) %>% as.character()
qexp(c(.025, .975), 50) %>% as.character()
dexp(c(.02), 50)

dgamma(.003, shape=.5, scale=20)
dgamma(.1, shape=.5, scale=20)
dgamma(20, shape=.5, scale=20)
dgamma(5, shape=.9, scale=10)
.9/100

data=c(0.6466898, 5.3839229)
gamma_opt_prior<-function(param, quant=c(.1, .9), data){
  #qgamma(c(.1, .9), shape=1.8, scale=1.5)-data
    qgamma(quant, shape=param[1], scale=param[2])-data->tmp
  sum(tmp^2)
}

gamma_opt_prior<-function(param, quant=c(.1, .9), data){
#qgamma(c(.1, .9), shape=1.8, scale=1.5)-data
  if (param[1]>1){
  qgamma(quant, shape=param[1], scale=param[2])-data->tmp
  }
  else tmp<-100
  
  sum(tmp^2)
}

nlm(gamma_opt_prior, p=c(5,5), data=c(0.6466898, 5.3839229))

nlm(gamma_opt_prior, p=c(4,3), data=c(0.5, 300))

qgamma(c(.1, .9), shape=1.009, scale=3)
dgamma(300, shape=1.009, scale=30)

beta_opt_prior<-function(param, quant=c(.1, .9), data){
  #qgamma(c(.1, .9), shape=1.8, scale=1.5)-data
  qbeta(quant, shape1=param[1], shape2=param[2])-data->tmp
  sum(tmp^2)
}

nlm(beta_opt_prior, p=c(3,3), data=c(0.1, .9))
qbeta(c(.1, .9), 1, 1)
dexp(.01, .1)

rb=rbeta(10000, 1.0058, 1)
delta=(300-.5)/(0.9-.1)
hist(rb*delta+.5, breaks=150)
