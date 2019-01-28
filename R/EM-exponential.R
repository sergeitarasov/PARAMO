# EM for exponential

# Intitialization
data<-c(rexp(100, 0.1), rexp(100, 10))
data=v
data.kmeans <- kmeans(data, 2)
data.kmeans.cluster <- data.kmeans$cluster
data.df <- data_frame(x = data, cluster = data.kmeans.cluster)

data.df %>%
  mutate(num = row_number()) %>%
  ggplot(aes(y = num, x = x, color = factor(cluster))) +
  geom_point() +
  ylab("Values") +
  ylab("Data Point Number") +
  scale_color_discrete(name = "Cluster") +
  ggtitle("K-means Clustering")

# Stattistics
data.summary.df <- data.df %>%
  group_by(cluster) %>%
  summarize(mu = mean(x), size = n())

data.summary.df <- data.summary.df %>%
  mutate(lambda=1/mu ,alpha = size / sum(size))
data.summary.df


1/mean(data[1:100])
1/mean(data[101:200])

dexp(data[1:100], 0.089) %>% log() %>% sum
dexp(data[1:100], 0.1) %>% log() %>% sum

dexp(data[101:200], 0.089) %>% log() %>% sum

dexp(data[1:100], 10) %>% log() %>% sum
dexp(data[101:200], 10) %>% log() %>% sum

c(dexp(data[1:100], .1), dexp(data[101:200], 10)) %>% log %>% sum

#######
# log likelihood of exp
data=rexp(1000, 1)
lambda<-1
L_exp(1, data)
dexp(data, 1) %>% log%>% sum

start_time <- Sys.time()
for (i in 1:10^6) {L_exp(1, data)}
end_time <- Sys.time()
end_time - start_time


L_exp<-function(lambda, data){
  log(lambda^length(data))-(lambda*length(data)*mean(data))
}

###
# likelihood of every cateory (exponential distribution) over branch
br=z$maps[[172]]
lambda<-1

L_exp_br(1, br)

L_exp_br<-function(lambda, br){
  c(dexp(br[1:length(br)-1]), 1-exp(-lambda*br[length(br)]) )
}
###

#' Expectation Step of the EM Algorithm
#'
#' Calculate the posterior probabilities (soft labels) that each component
#' has to each data point.
#'
#' @param sd.vector Vector containing the standard deviations of each component
#' @param sd.vector Vector containing the mean of each component
#' @param alpha.vector Vector containing the mixing weights  of each component
#' @return Named list containing the loglik and posterior.df
e_step <- function(x, lambda.vector, alpha.vector) {
  comp1.prod <- dexp(x, lambda.vector[1]) * alpha.vector[1]
  comp2.prod <- dexp(x, lambda.vector[2]) * alpha.vector[2]

  sum.of.comps <- comp1.prod + comp2.prod

  comp1.post <- comp1.prod / sum.of.comps

  
  comp2.post <- comp2.prod / sum.of.comps

  
  sum.of.comps.ln <- log(sum.of.comps, base = exp(1))
  sum.of.comps.ln.sum <- sum(sum.of.comps.ln)
  
  list("loglik" = sum.of.comps.ln.sum,
       "posterior.df" = cbind(comp1.post, comp2.post))
}

e_step(x=data, lambda.vector=m.step[["lambda"]], alpha.vector=m.step[["alpha"]])
e.step=e_step(data, data.summary.df[["lambda"]], data.summary.df[["alpha"]])

#' Maximization Step of the EM Algorithm
#'
#' Update the Component Parameters
#'
#' @param x Input data.
#' @param posterior.df Posterior probability data.frame.
#' @return Named list containing the mean (mu), variance (var), and mixing
#'   weights (alpha) for each component.
m_step <- function(x, posterior.df) {
  comp1.n <- sum(posterior.df[, 1])
  comp2.n <- sum(posterior.df[, 2])
  
  comp1.lambda <- comp1.n / sum(posterior.df[, 1] * x)
  comp2.lambda <- comp2.n / sum(posterior.df[, 2] * x)
  
  comp1.alpha <- comp1.n / length(x)
  comp2.alpha <- comp2.n / length(x)
  
  list("lambda" = c(comp1.lambda, comp2.lambda),
             "alpha" = c(comp1.alpha, comp2.alpha))
}

m_step(data, e.step[["posterior.df"]])

i=1
for (i in 1:50) {
  if (i == 1) {
    # Initialization
    
    e.step <- e_step(data, data.summary.df[["lambda"]], data.summary.df[["alpha"]])
    
    m.step <- m_step(data, e.step[["posterior.df"]])
    
    cur.loglik <- e.step[["loglik"]]
    loglik.vector <- e.step[["loglik"]]
  } else {
    # Repeat E and M steps till convergence
  
    e.step <- e_step(data, m.step[["lambda"]], m.step[["alpha"]])
    
    m.step <- m_step(data, e.step[["posterior.df"]])
    
    loglik.vector <- c(loglik.vector, e.step[["loglik"]])
    
    loglik.diff <- abs((cur.loglik - e.step[["loglik"]]))
    if (loglik.diff < 1e-6) {
      break
    } else {
      cur.loglik <- e.step[["loglik"]]
    }
  }
}
loglik.vector
hist(data)

############
#' Plot a Mixture Component
#' 
#' @param x Input ata.
#' @param mu Mean of component.
#' @param sigma Standard of component.
#' @param lam Mixture weight of component.
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

data.frame(x = wait) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(m.step$mu[1], sqrt(m.step$var[1]), 
                            lam = m.step$alpha[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(m.step$mu[2], sqrt(m.step$var[2]), 
                            lam = m.step$alpha[2]),
                colour = "blue", lwd = 1.5) +
  ylab("Density") +
  xlab("Values") +
  ggtitle("Final GMM Fit")