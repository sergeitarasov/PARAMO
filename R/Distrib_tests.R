# Kolmogorov-Smirnov Tests

x <- rnorm(50)
y <- runif(30)
y <- rnorm(50)

ks.test(x, y)

wilcox.test(rnorm(10), rnorm(10, 2), paired=T, conf.int = TRUE)
wilcox.test(rnorm(10), rnorm(20), paired=F)

wilcox.test(rnorm(200, mean=2, sd=10), rnorm(200, mean=2, sd=10), paired=F)
wilcox.test(rnorm(2000, mean=2, sd=10), rnorm(2000, mean=2, sd=1), paired=F)
wilcox.test(rnorm(2000, mean=2, sd=1), rexp(2000, 0.5), paired=F)
rexp(2000, 0.5) %>% mean()

hist(rnorm(2000, mean=2, sd=1)-rexp(2000, 0.5), breaks = 50)
hist(rexp(2000, 1)-rexp(2000, 0.5), breaks = 50)

t.test(rnorm(2000, mean=2.5, sd=1), rexp(2000, 0.5))
t.test(rnorm(2000, mean=3, sd=1), rnorm(2000, mean=2, sd=1))

library("wBoot", lib.loc="~/.local/R/site-library")

# Driving distances, in yards, for independent samples of drives off a
# 2-3/4" wooden tee and off a 3" Stinger Competition golf tee.
data("tees")
str(tees)
attach(tees)
# Note that the data are unstacked.

# 99% confidence interval for the difference between the mean driving
# distances of the two types of tees. Name variable DISTANCE.
boot.two.bca(REGULAR, STINGER, mean, stacked = FALSE, variable = "DISTANCE",
             conf.level = 0.99)

boot.two.bca(rnorm(100, mean=2.5, sd=1), rnorm(100, mean=2, sd=1), mean, stacked = FALSE, null.hyp=0,
             conf.level = 0.99)

boot.two.bca(rnorm(100, mean=2, sd=1), rnorm(100, mean=2, sd=1), mean, stacked = FALSE, null.hyp=0,
             conf.level = 0.99)

boot.two.bca(rnorm(100, mean=3, sd=1), rexp(2000, 0.5), mean, stacked = FALSE, null.hyp=0,
             conf.level = 0.99)

boot.two.bca(sample(ST$cranium, 5000), sample(unlist(ST, use.names = FALSE), 10000), mean, stacked = FALSE, null.hyp=0,
             conf.level = 0.99, R = 1000)

boot.two.bca(sample(ST$metanotum, 5000), sample(unlist(ST, use.names = FALSE), 10000), mean, stacked = FALSE, null.hyp=0,
             conf.level = 0.99, R = 1000)

boot.two.bca(sample(ST$cranium, 5000), sample(ST$mesonotum, 5000), mean, stacked = FALSE, null.hyp=0,
             conf.level = 0.99, R = 1000)

boot.two.bca(ST$cranium, ST$mesonotum, mean, stacked = FALSE, null.hyp=0,
             conf.level = 0.99, R = 5000)

boot.two.bca(ST$`fore wing`, ST$`hind wing`, mean, stacked = FALSE, null.hyp=0,
             conf.level = 0.99, R = 1000)

F
