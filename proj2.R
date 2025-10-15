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

## ---------------------------------------------------------------
## nseir()
## Purpose
##   Simulate a stochastic SEIR epidemic with three infection ways:
##   (1) infected by household members, (2) infected by regular network contacts
##   incfection, and (3)random mixing infection.
##
## Inputs:
##   beta   : length-n numeric; individual variability in the transmission rate
##   h      : length-n integer; household IDs 
##   alink  : list of integer vectors; alink[[i]] are regular (non-household) contacts of person i
##   alpha  : c(ah, ac, ar). daily infection probs for household, network, and residual random mixing
##   delta  : daily probability I -> R
##   gamma  : daily probability E -> I
##   nc     : target mean degree (used to scale random-mixing term)
##   nt     : number of simulated days
##   pinf   : the proportion of the initial population to randomly start in the I state
##
## Output:
##   a list with elements S, E, I, R and t 
##   giving the total population in each class each day and the day, respectively.
##
## How it works:
##   Each individual j has a disease state x_j ∈ {S=0, E=1, I=2, R=3}.
##   The simulation proceeds in daily time steps. On each day:
##      1) Disease progression:
##         • Infectious individuals (I) recover with probability delta (I→R);
##         • Exposed individuals (E) become infectious with probability gamma (E→I).
##      2) Infection probabilities:
##          For each susceptible person j, compute their probability 
##          of becoming exposed through three independent transmission routes:
##            • Household: increases with the number of infectious household members;
##            • Network: increases with the number of infectious social contacts;
##            • Random mixing: background exposure proportional to beta_j and the
##              total infectious “mass” (Σ beta_i over infectives).
##
##          Then, combine the three independent infection risks as
##          p_j = 1 - (1 - p_hh)(1 - p_net)(1 - p_rnd). Finally, for each susceptible j,
##          Generate a random number U_j ~ Uniform(0,1). Compare it to p_j.
##          If U_j < p_j, person j becomes newly exposed (state changes S→E), otherwise,
##          j is still susceptible (state S).
##
##    After all transitions, record population totals of S, E, I, and R for that day.

## ---------------------------------------------------------------
nseir <- function(beta, h, alink,
                  alpha = c(.1, .01, .01), delta = .2, gamma = .4,
                  nc = 15, nt = 100, pinf = .005) {
  
  n <- length(beta)
  b_mean <- mean(beta)
  ah <- alpha[1]; ac <- alpha[2]; ar <- alpha[3]
  
  # Initialise all individuals as susceptible (state 0),
  # then randomly select a small fraction (pinf of the population) to be infectious (state 2).
  x <- rep(0, n)
  x[sample.int(n, max(1, round(pinf * n)))] <- 2
  
  
  
  # Create empty vectors to store daily totals of S, E, I, and R
  S <- E <- I <- R <- rep(0, nt)
  S[1] <- sum(x == 0); I[1] <- sum(x == 2)
  
  # Create a list where each element contains the indices of people in the same household.
  hh <- split(1:n, h)
  
  # ---- Daily simulation loop ---------------------------------------------
  for (t in 2:nt) {
    
    # Generate one random number per person to decide all their possible state changes
    u <- runif(n)
    
    # Disease progression:
    # Infectious recover with probability delta and exposed become infectious with gamma.
    x[x == 2 & u < delta] <- 3               # I -> R
    x[x == 1 & u < gamma] <- 2               # E -> I
    
    # The infection probabilities of each way for each individual j.
    # Naming: p_h (household), p_n (network), p_r (random mixing).
    p_h <- p_n <- p_r <- numeric(n)
    
    # Household infection:
    # If a household has hi infectious members, every member's daily risk is
    # p_h = 1 - (1 - ah)^hi, assuming independent transmission events.
    if (ah > 0) {
      # Loop through each household (list of member indices)
      for (members in hh) {                 
        hi <- sum(x[members] == 2)   # Count how many people in this household are infectious      
        if (hi > 0) p_h[members] <- 1 - (1 - ah)^hi   # assign infection probability to all household members
      }
    }
    
    # Contact-network infection:
    # Each person j faces infection risk from infectious contacts in their network.
    # If they have n_i infectious neighbours, their daily probability is p_n = 1 - (1 - ac)^n_i.
    if (ac > 0) {
      # Find indices of all currently infectious individuals
      n_infectious <- which(x == 2)
      
      # Initialise vector: n_i[j] counts infectious neighbours of j
      n_i <- numeric(n)
      
      # Loop through every infectious person
      for (i in n_infectious) {
        neighbours <- alink[[i]]
        
        # Check if they have any contacts,
        if (length(neighbours))
          
          # increase neighbour counts for those contacts
          n_i[neighbours] <- n_i[neighbours] + 1L
      }
      # Infection probabilities:
      p_n <- 1 - (1 - ac)^n_i
    }
    
    # Residual random mixing:
    # Calculate the probability of infection not captured by households or edges.
    if (ar > 0) {
      sumB <- sum(beta[x == 2])                    # sum of beta for all infectious people
      
      # daily probability of random mixing 
      p_r <- (ar * nc / (b_mean^2 * (n - 1)))* beta * sumB
      p_r[p_r > 1] <- 1                            # numerical safety (rare)
    }
    
    # Combine the 3 independent infection ways:
    p <- 1 - (1 - p_h) * (1 - p_n) * (1 - p_r)
    
    # Apply S->E only to susceptibles: Bernoulli(U < p_j) per person.
    s <- which(x == 0)
    if (length(s)) x[s[runif(length(s)) < p[s]]] <- 1
    
    # Record the number of individuals in each disease state (S, E, I, R) at day t.
    S[t] <- sum(x == 0); E[t] <- sum(x == 1)
    I[t] <- sum(x == 2); R[t] <- sum(x == 3)
  }
  
  # Return daily counts and beta.
  list(S = S, E = E, I = I, R = R, beta = beta, t=seq_len(nt))
}
