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

## ---- Q3: SEIR with households + contacts + random mixing ----
nseir <- function(beta, h, alink,alpha = c(.1, .01, .01),delta=.2,gamma=.4,
                  nc = 15, nt = 100, pinf = .005){
  n <- length(beta); b_mean <- mean(beta)
  ah <- alpha[1]; ac <- alpha[2]; ar <- alpha[3]
  x <- rep(0, n); x[1:max(1, round(pinf*n))] <- 2  # 0=S,1=E,2=I,3=R
  S <- E <- I <- R <- rep(0, nt); S[1] <- sum(x==0); I[1] <- sum(x==2)
  hh <- split(1:n, h)
  for(t in 2:nt){
    u <- runif(n)
    x[x==2 & u<delta] <- 3            # I -> R
    x[x==1 & u<gamma] <- 2            # E -> I
    p_h <- p_n <- p_r <- numeric(n)
    
    ## household infection
    if (ah > 0)
      for (members in hh) {
        hi <- sum(x[members] == 2)
        if (hi > 0) p_h[members] <- 1 - (1 - ah)^hi
      }
    ## network infection
    if (ac > 0) {
      n_infectious <- which(x == 2)                # all infectious people
      inf_neighbour_count <- numeric(n)            # m_I(j) for each person j
      
      # count how many infectious neighbours each person has
      for (i in n_infectious) {
        neighbours <- alink[[i]]
        if (length(neighbours)) 
          inf_neighbour_count[neighbours] <- inf_neighbour_count[neighbours] + 1L
      }
      
      # convert neighbour counts into infection probabilities
      p_n <- 1 - (1 - ac)^inf_neighbour_count
    } else {
      p_n <- numeric(n)
    }
    
    ## random mixing (simple linear add of tiny pairwise risks)
    p_r <- numeric(n)
    if (ar > 0) {
      sumB <- sum(beta[x == 2])            # sum of beta over infectives
      C <- ar * nc / (b_mean^2 * (n - 1))
      p_r <- C * beta * sumB               # p_rnd_j = C * beta_j * sumB
      p_r[p_r > 1] <- 1                   # clamp
    }
    
    
    ## total S->E transition
    p <- 1 - (1 - p_h) * (1 - p_n) * (1 - p_r)
    s <- which(x==0)
    if(length(s)) x[s[runif(length(s)) < p[s]]] <- 1
    S[t] <- sum(x==0); E[t] <- sum(x==1); I[t] <- sum(x==2); R[t] <- sum(x==3)
  }
  list(S=S,E=E,I=I,R=R,beta=beta)
}