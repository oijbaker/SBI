n <- 100000
eps <- 20
p <- 2
accepted_samples <- matrix(nrow=n, ncol=p) # stores accepted parameter samples
count <- 0
# Initialise
x_sim <- list()
dist <- list()
theta_sim <- list()


#################
# ABC algorithm #
#################

for(i in 1:n){
  theta_sim[i] <- runif(p, 0, 1) # sample from prior
  x_sim[i] <- ... # simulate data from SIR model
  dist[i] <- calc_dist(x_sim[i], x) # calculate distance
  # Accept-reject step
  if(dist[i] <= eps){
    accepted_samples[i,] <- theta_sim[i]
    count <- count + 1
  }
  
  cat("Acceptance rate: ", count/n, "\n")
}

