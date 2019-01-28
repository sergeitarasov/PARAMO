part_scheme=list(c(1, 2), c(3), c(4))
part_scheme=list(c(1), c(2), c(3,4))
part_scheme=list(c(1,2), c(3,4))



rows2rate_matrix(M.temp, row.vec, dt_rates[,7])

n.lump=apply(dt_rates, 2, function(x)
  rows2rate_matrix(M.temp, row.vec, x) %>% is_strg_lumpable(., part_scheme)
  )
which(n.lump==T) %>% length

apply(dt_rates[,which(n.lump==T)], 2, function(x) length(unique(x))) %>% table


n.lump1=apply(dt_rates, 2, function(x)
  rows2rate_matrix(M.temp, row.vec, x) %>% is_strg_lumpable(., part_scheme)
)

setdiff(which(n.lump==T), which(n.lump1==T))

ncol=ncol(dt_rates)
hist(dt_rates[2,which(n.lump==T)], breaks=10)
hist(dt_rates[2, which(n.lump==F)], breaks=10)
hist(dt_rates[2,], breaks=10)

rate_coll=1
plot(density(dt_rates[rate_coll,which(n.lump==T)]),type="l",col="darkgreen",lwd=4)
lines(density(dt_rates[rate_coll,which(n.lump==F)]),type="l",col="red",lwd=4)
lines(density(dt_rates[rate_coll,]),type="l",col="blue",lwd=4)

rate_coll1=4
rate_coll2=7
plot(density(dt_rates[rate_coll1,which(n.lump==T)]+dt_rates[rate_coll2,which(n.lump==T)]),type="l",col="darkgreen",lwd=4)
lines(density(dt_rates[rate_coll1,which(n.lump==F)]+dt_rates[rate_coll2,which(n.lump==F)]),type="l",col="red",lwd=4)
lines(density(dt_rates[rate_coll,]),type="l",col="blue",lwd=4)


#####
max_err=apply(dt_rates, 2, function(x)
  rows2rate_matrix(M.temp, row.vec, x) %>% lump_max_error(., part_scheme, 1)
)
round(max_err, 5)
hist(max_err,breaks=20)
?hist

which(round(max_err, 4)==0) %>% length


MCMCpack::rdirichlet(1, c(1,1,1) )

