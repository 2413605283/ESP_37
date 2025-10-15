# Practical 2 — Social Structure in SEIR Models 
#
#Yunhan Zhang: s2176155. Xiyu Wu: s2799746. Tianyu Wang: s2794991
#Yunhan did question 3. Xiyu did question 4 and 5. Tianyu did question 1 and 2.
#Each member did roughly the same amount of work.
#
#
# Overview
#   This script simulates the spread of an infectious disease in a 
#   population where people are infected in three ways:
#     • within households,
#     • through regular social network of contacts ,
#     • through random mixing across the whole population.
#
#   The model follows the stochastic SEIR model, where each individual
#   is in one of four states:
#       S – Susceptible
#       E – Exposed (infected but not yet infectious)
#       I – Infectious
#       R – Recovered (immune)
#
#   Each day, individuals may progress (E→I, I→R) or become newly exposed
#   (S→E) based on their infectious contacts in these three ways.
# ======================================================================



# Q1: Generate the household vector "h".
#
# Purpose:
#   Create a length-n vector, showing which household each person belongs to.
#
#   e.g. if h[i] = h[j], it indicates i and j are in the same household.
#
# Inputs:
#   n: total population.
#   hmax: maximum household size.
#
# Output:
#   A vector "h", whose elements are numbers, recording each person's household index.
#  
# How it works:
#   (1) Rdomly assign each household a size between 1 and hmax.
#   (2) Rpeat (1) for n times to get n , each household size is U(0,hmax). 
#   (3) According to each household size, generate the same number of household index.
#   (4) Keep the first n entries so that the total population is n.
#   (5) Randomize the order.

n = 10000
hmax = 5

h <- sample(rep(1:n, times = sample(1:hmax, n, replace = TRUE))[1:n])



# Q2: generate the regular contacts vector "alink".
# 
# Arguments:
#   beta: a vector, whose element beta_i showing i's beta value.
#   beta_bar: the average value of all beta_i.
#   h: the household vector "h", the result of Q1.
#   nc: the average number of contacts per person.
#
# Function: get.net()
#
# Purpose:
#   Build the (non-household) social network for the population.
#   Each pair of people (i, j) is connected with probability p_ij,
#   p_ij = nc * beta_i * beta_j / (beta_bar^2 * (n - 1)).
#
#   e.g. if i's social network includes j, then j's social network must includes i,
#        (i and j are not in the same household).
#
# Inputs:
#   beta: a vector, whose element beta_i showing i's beta value.
#   nc: the average number of contacts per person.
#
# Output:
#   adj: a list of vectors; adj[[i]] contains the index of person i’s regular contacts
#
# How it works:
#   (1) For each person i, exclude i and his household members from possible contacts,
#       the remaining people are denoted as "candidate".
#   (2) Consider only people j > i to avoid making duplicate links.
#   (3) Compute link probabilities p_ij and generate random numbers from 0 to 1, 
#       if random numbers < p_ij, record a connection.
#   (5) After the loop, make the network symmetric, that is, pij = pji.
#   (6) Convert the logical matrix into an adjacency list.

beta <- runif(n)  # each person's beta value，follows U(0,1)
beta_bar <- mean(beta) #average value


# Note:
# For each person i, we only consider j where j>i when computing pij,
# then we set pji = pij, because social links are undirected.

# This indicates that we only need to compute the upper triangle of the matrix "link".

get.net <- function(beta, nc = 15) {

  # create an n*n logical matrix to record possible links
  link <- matrix(FALSE, n, n)
  
  # do not need to consider i=n, because no index is bigger than n
  for (i in 1:(n-1)) {
    
    # exclude i's household members and i himself
    family_i <- c(h[[i]], i)
    candidates <- setdiff((i+1):n, family_i)
    
    if (length(candidates) > 0) {
      # compute probability pij(cannot exceed 1) 
      pij <- beta[i] * beta[candidates] * nc / (beta_bar^2 * (n - 1))
      pij[pij > 1] <- 1
      
      # rand is a vector, containing random numbers from 0 to 1
      # generate random numbers and decide which pairs form a connection
      rand <- runif(length(candidates))
      connected <- candidates[rand < pij]
      
      # fill in "TRUE/FALSE" in the matrix "link"
      if (length(connected) > 0)
        link[i, connected] <- TRUE
    }
  }
  
  # make the matrix symmetric, as pij = pji
  link <- link | t(link)
  
  # convert the matrix to a list
  # if pij = TRUE, then adj[i] contains j 
  adj <- apply(link, 1, function(x) which(x))
  
  return(adj)
}

# calculate the input for Q3
alink <- get.net(beta, nc = 15)


#Q3:Runs the SEIR simulation.
# ---------------------------------------------------------------
# nseir()
# Purpose
#   Simulate a stochastic SEIR epidemic with three infection ways:
#   (1) infected by household members, (2) infected by regular network contacts
#   incfection, and (3)random mixing infection.
#
# Inputs:
#   beta   : length-n numeric; individual variability in the transmission rate
#   h      : length-n integer; household IDs 
#   alink  : list of integer vectors; alink[[i]] are regular (non-household) contacts of person i
#   alpha  : c(ah, ac, ar). daily infection probs for household, network, and residual random mixing
#   delta  : daily probability I -> R
#   gamma  : daily probability E -> I
#   nc     : target mean degree (used to scale random-mixing term)
#   nt     : number of simulated days
#  pinf   : the proportion of the initial population to randomly start in the I state
#
# Output:
#   a list with elements S, E, I, R and t 
#   giving the total population in each class each day and the day, respectively.
#
# How it works:
#   Each individual j has a disease state x_j ∈ {S=0, E=1, I=2, R=3}.
#   The simulation proceeds in daily time steps. On each day:
#      1) Disease progression:
#         • Infectious individuals (I) recover with probability delta (I→R);
#         • Exposed individuals (E) become infectious with probability gamma (E→I).
#      2) Infection probabilities:
#          For each susceptible person j, compute their probability 
#          of becoming exposed through three independent transmission routes:
#            • Household: increases with the number of infectious household members;
#            • Network: increases with the number of infectious social contacts;
#            • Random mixing: background exposure proportional to beta_j and the
#              total infectious “mass” (Σ beta_i over infectives).
#
#          Then, combine the three independent infection risks as
#          p_j = 1 - (1 - p_hh)(1 - p_net)(1 - p_rnd). Finally, for each susceptible j,
#          Generate a random number U_j ~ Uniform(0,1). Compare it to p_j.
#          If U_j < p_j, person j becomes newly exposed (state changes S→E), otherwise,
#          j is still susceptible (state S).
#
#    After all transitions, record population totals of S, E, I, and R for that day.

# ---------------------------------------------------------------
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

# Q4  Plotting function: visualize SEIR time series clearly
#
# Purpose
#   Given the SEIR simulation results from nseir() (S, E, I, R),
#   draw a clear, well-formatted line plot:
#      consistent colours and axes across scenarios
#      separate title() call to control font size and offset

# Inputs
#   res       : list containing res$S, res$E, res$I, res$R
#   main      : string, plot title
#   ylim      : numeric vector [ymin, ymax]; if NULL, set automatically
#   col       : colour vector for S/E/I/R (black/blue/red/darkgreen)
#   cex_main  : title text scaling factor (<1 slightly smaller)
#   line_main : vertical offset of the title (>1 moves it down slightly)

# Outputs
#   None (produces a plot on the current graphics device)

# Outline of logic
#   1) Use plot() once to draw S and set up axes (no title yet)
#   2) Add E, I, R with lines() to overlay in same coordinate frame
#   3) Add title() and legend() for clarity and consistency
# ------------------------------------------------------------
plot_epi <- function(res, main = "SEIR",
                     ylim = NULL,
                     col = c("black", "blue", "red", "darkgreen"),
                     cex_main = 0.95, line_main = 1.2,
                     legend_cex = 0.9,
                     legend_yintersp = 0.9) {
  
  # Auto-select y-axis range
  if (is.null(ylim))
    ylim <- c(0, max(res$S, res$E, res$I, res$R))
  
  # Plot S curve first to set up axes
  plot(seq_along(res$S), res$S, type = "l", lwd = 2, col = col[1],
       xlab = "day", ylab = "population", ylim = ylim, main = "")
  
  # Add main title
  title(main = main, cex.main = cex_main, line = line_main)
  
  # Overlay the E, I, R curves
  lines(seq_along(res$E), res$E, lwd = 2, col = col[2])
  lines(seq_along(res$I), res$I, lwd = 2, col = col[3])
  lines(seq_along(res$R), res$R, lwd = 2, col = col[4])
  
  # ---- Legend vertically placed beside the title ----
  op <- par(xpd = NA)  # allow drawing outside the plot region
  on.exit(par(op), add = TRUE)
  
  # Determine coordinates for legend placement
  usr <- par("usr")  # c(xmin, xmax, ymin, ymax)
  x_pos <- usr[2] + 0.08 * diff(usr[1:2])   # a bit right of the plot
  y_pos <- usr[4] - 0.02 * diff(usr[3:4])   # slightly above the top
  
  legend(x = x_pos, y = y_pos,
         legend = c("S", "E", "I", "R"),
         col = col, lwd = 2,
         bty = "n", horiz = FALSE,
         cex = legend_cex,
         y.intersp = legend_yintersp,
         x.intersp = 0.6,
         seg.len = 1.1)
}
# Q5  Run four SEIR scenarios and plot results

# Purpose
#   This function runs the four epidemic scenarios and plots them in a 2×2 grid for direct visual comparison:
#     (1) Full model: household + fixed contacts + random mixing
#     (2) Random mixing only (α_h=α_c=0, α_r=0.04)
#     (3) Full model with constant beta
#     (4) Random mixing only with constant beta

# Inputs
#   beta, h, alink : individual-level parameters:
#                    - beta: individual infectiousness/contact propensities
#                    - h: household membership vector (same length as beta)
#                    - alink: fixed contact network (list of neighbours)

#   nt  : number of time steps to simulate (default 100)
#         Each iteration corresponds to one “day” in the model
#         A value of 100 is sufficient to observe a full epidemic wave (rise–peak–decay) while keeping computation time reasonable

#   nc  : expected number of fixed contacts per person (15)

#   pinf: initial infection proportion (default 0.005)
#         Sets the fraction of individuals initially infectious
#         For n=10000, pinf=0.005 then 50 initial infectives
#         This value is large enough to seed an outbreak but small enough to observe the early exponential growth phase

#   alpha_full      : (α_h, α_c, α_r) used for the full model scenario
#   alpha_rand_only : (α_h, α_c, α_r) used for random-mixing-only scenario

# Outputs
#   Invisibly returns a list with results for each of the four scenarios:
#     $full, $rand_only, $full_const_beta, $rand_const_beta
#   Each element is a list of time series vectors (S, E, I, R, beta)
#   The function produces a 2×2 panel plot as side effect

# Implementation outline
#   - For each scenario, call nseir() with appropriate alpha and beta settings
#   - For constant-beta scenarios, rebuild the contact network because get.net() depends on beta heterogeneity
#   - Use par(mfrow=c(2,2)) to draw all four plots in one figure.
#   - Use larger top margin to prevent titles being clipped in RStudio.

#Comment
#  When the model includes family and social network structures (Full model), the epidemic rises faster, peaks higher, 
#  and comes earlier because the infection spreads rapidly within closely connected small groups. 
#  However, the cumulative number of infections at the end of the entire epidemic is often slightly lower than that of the random mixing model, 
#  as the clustered structure causes the epidemic to "burn out" quickly within local groups, reducing the chances of cross-group transmission.
#  In contrast, the curve of the random mixing model is smoother, slower, with a lower peak but a longer duration. 
#  The scenario using a constant β also indicates that these differences are mainly caused by structural effects rather than individual differences.
# ------------------------------------------------------------
run_four_scenarios <- function(beta, h, alink,
                               nt = 100, nc = 15, pinf = 0.005,
                               alpha_full = c(0.1, 0.01, 0.01),
                               alpha_rand_only = c(0, 0, 0.04)) {
  n <- length(beta)
  
  # (1) Full model – all three transmission channels active
  epi1 <- nseir(beta = beta, h = h, alink = alink,
                alpha = alpha_full, delta = .2, gamma = .4,
                nc = nc, nt = nt, pinf = pinf)
  
  # (2) Random mixing only – disable household and network links
  epi2 <- nseir(beta = beta, h = h, alink = alink,
                alpha = alpha_rand_only, delta = .2, gamma = .4,
                nc = nc, nt = nt, pinf = pinf)
  
  # (3) Full model with constant beta:
  #     Replace individual betas with their mean and rebuild network
  beta_bar <- mean(beta)
  beta_cst <- rep(beta_bar, n)
  alink_cst <- get.net(beta_cst, nc = nc)  # rebuild network
  
  epi3 <- nseir(beta = beta_cst, h = h, alink = alink_cst,
                alpha = alpha_full, delta = .2, gamma = .4,
                nc = nc, nt = nt, pinf = pinf)
  
  # (4) Constant beta + random mixing only:
  #     Use the same constant-beta network but disable household and network links
  epi4 <- nseir(beta = beta_cst, h = h, alink = alink_cst,
                alpha = alpha_rand_only, delta = .2, gamma = .4,
                nc = nc, nt = nt, pinf = pinf)
  
  # ---- Plotting section ----
  # Create a 2×2 grid of subplots; enlarge top margin to keep titles visible
  op <- par(mfrow = c(2,2), mar = c(4,4,5.4,4))
  on.exit(par(op), add = TRUE)  # restore settings afterwards
  
  # Optional enhancement: compute a common y-axis limit (Y) for fair visual comparison
  # Y <- max(epi1$S, epi1$E, epi1$I, epi1$R,
  #          epi2$S, epi2$E, epi2$I, epi2$R,
  #          epi3$S, epi3$E, epi3$I, epi3$R,
  #          epi4$S, epi4$E, epi4$I, epi4$R)
  
  # Plot each scenario with descriptive two-line titles
  plot_epi(epi1, main = "Full model\n(Household + Network + Random)")
  plot_epi(epi2, main = "Random mixing only\n(α_h=α_c=0, α_r=0.04)")
  plot_epi(epi3, main = "Full model with constant beta")
  plot_epi(epi4, main = "Random mixing only with constant beta")
  
  # Return all results invisibly
  invisible(list(full = epi1,
                 rand_only = epi2,
                 full_const_beta = epi3,
                 rand_const_beta = epi4))
}
dev.new(width = 13, height = 9)
res4 <- run_four_scenarios(beta = beta, h = h, alink = alink,
                           nt = 100, nc = 15, pinf = .005)
