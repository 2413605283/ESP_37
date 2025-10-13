n=1000
hmax=5
##################1
h <- sample(rep(1:n, times = sample(1:hmax, n, replace = TRUE))[1:n])


beta <- runif(n)  # everyone's beta value，U(0,1)

######################2
get.net <- function(beta, nc = 15) {
  
  # arguments：
  #   beta : everyone's beta value，U(0,1)
  #   nc   : 平均社交连接数numbers of average network population
  # return：
  #   adj  : everyone's network index
  # ----------------------------------------
  
  n <- length(beta)
  beta_bar <- mean(beta)
  
  # pij :prob of i and j have network
  pij <- outer(beta, beta, "*") * nc / (beta_bar^2 * (n - 1))
  pij[pij > 1] <- 1   # 概率不能超过 1
  
  # random matrix 
  rand <- matrix(runif(n^2), n, n)
  
  # compare Pij and random matrix, get yes/no matrix
  link <- rand < pij
  
  # everyone cannot hae network with themselves
  diag(link) <- FALSE
  
  # make it symmetric, eg: i has connection with j, 
  # implies j also has connection with i
  link[lower.tri(link)] <- t(link)[lower.tri(link)]
  
  # transfer matrix(above we get) into adj(for each one, "true" means "has connection")
  # summary all true and get a subset which record the index of other people
  adj <- apply(link, 1, function(x) which(x))
  
  return(adj)
}

alink <- get.net(beta, nc = 15)
