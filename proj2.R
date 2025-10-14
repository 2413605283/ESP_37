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
                     col = c("black","blue","red","darkgreen"),
                     cex_main = 0.95, line_main = 1.2) {
  
  # Automatically choose y-axis range so that all curves are visible
  if (is.null(ylim))
    ylim <- c(0, max(res$S, res$E, res$I, res$R))
  
  # Draw the S curve first to establish axes.
  # Leave 'main' blank here and control title formatting later
  plot(seq_along(res$S), res$S, type = "l", lwd = 2, col = col[1],
       xlab = "day", ylab = "population", ylim = ylim, main = "")
  
  # Add the title separately to adjust size and offset precisely
  title(main = main, cex.main = cex_main, line = line_main)
  
  # Overlay the remaining three compartments on the same axes
  lines(seq_along(res$E), res$E, lwd = 2, col = col[2])
  lines(seq_along(res$I), res$I, lwd = 2, col = col[3])
  lines(seq_along(res$R), res$R, lwd = 2, col = col[4])
  
  # Add a legend (top right, no box) matching colours and labels
  legend("topright", bty = "n", lwd = 2, col = col,
         legend = c("S","E","I","R"))
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
  op <- par(mfrow = c(2,2), mar = c(4,4,5.4,1))
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
res4 <- run_four_scenarios(beta = beta, h = h, alink = alink,
                           nt = 100, nc = 15, pinf = .005)