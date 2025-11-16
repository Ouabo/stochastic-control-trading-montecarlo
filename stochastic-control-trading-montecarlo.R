############################################################
# Stochastic Optimal Control with Monte Carlo Approximation
# R implementation matching the Article description
############################################################
#Load the following packages: for Visualization##
library(ggplot2)
library(dplyr)
library(gganimate)
library(gifski)   
library(png) 
library(tidyr)
##############################

#set.seed(123)

##############################
# 1. Load or define price data
##############################

# Example: you would replace this with real data:
# Exple: S_hist <- your vector of historical daily prices 
# Here we create a dummy price series just for illustration:
N_hist <- 1000
S_hist <- cumsum(c(100, rnorm(N_hist, mean = 0, sd = 1)))  # random walk around 100

##########################################
# 2. Calibration of sigma from price data
##########################################

Delta_S <- diff(S_hist)            # increments
Delta_S_bar <- mean(Delta_S)
sigma_hat <- sqrt(sum((Delta_S - Delta_S_bar)^2) / (length(Delta_S) - 1))

cat("Calibrated sigma:", sigma_hat, "\n")

###################################################
# 3. Model parameters and discretization of states
###################################################

T_horizon <- 100        # number of decision steps (t = 0,...,T-1)
x_max     <- 30       # maximum inventory
a_max     <- 10        # maximum action per step
c_cost    <- 0.02     # trading cost coefficient
M_mc      <- 2000     # number of Monte Carlo samples

# Inventory grid: x_i in {0,1,...,x_max}
x_grid <- 0:x_max
I <- length(x_grid)

# Price grid: choose min/max around historical range
S_min <- quantile(S_hist, 0.05)
S_max <- quantile(S_hist, 0.95)
J     <- 21                          # number of price grid points
S_grid <- seq(S_min, S_max, length.out = J)

cat("Price grid from", S_min, "to", S_max, "with", J, "points.\n")

###########################################
# 4. Containers for V_t and optimal policy
###########################################

# We store V_t for t = 0,...,T.
# V_list[[t+1]] is the matrix of V_t, dimension I x J
V_list <- vector("list", T_horizon + 1)

# Also store optimal action a_t^* at each grid point (same dimension)
A_list <- vector("list", T_horizon)

#############################################
# 5. Terminal condition: V_T(x,S) = 0 for all
#############################################

V_T <- matrix(0, nrow = I, ncol = J)
V_list[[T_horizon + 1]] <- V_T

#####################################
# 6. Helper: linear interpolation in S
#####################################

# Given V_next (I x J) = V_{t+1}(x_i, S^{(j)}),
# x_index in 1:I, and a price S_val,
# return an approximate V_{t+1}(x_i, S_val) by linear interpolation.
interp_V_next <- function(V_next, x_index, S_val, S_grid) {
  # If S_val is outside the grid, clamp to edge
  if (S_val <= S_grid[1]) {
    return(V_next[x_index, 1])
  }
  if (S_val >= S_grid[length(S_grid)]) {
    return(V_next[x_index, length(S_grid)])
  }
  # Find j such that S_grid[j] <= S_val <= S_grid[j+1]
  j <- max(which(S_grid <= S_val))
  if (j == length(S_grid)) {
    return(V_next[x_index, j])
  }
  S_low <- S_grid[j]
  S_high <- S_grid[j + 1]
  lambda <- (S_high - S_val) / (S_high - S_low)
  # Linear interpolation
  V_low  <- V_next[x_index, j]
  V_high <- V_next[x_index, j + 1]
  V_interp <- lambda * V_low + (1 - lambda) * V_high
  return(V_interp)
}

############################################################
# 7. Backward induction: compute V_t and optimal policy A_t
############################################################

for (t in (T_horizon-1):0) {
  cat("Computing V_", t, "...\n", sep = "")
  
  V_next <- V_list[[t + 2]]   # V_{t+1}, matrix I x J
  V_t    <- matrix(0, nrow = I, ncol = J)
  A_t    <- matrix(0, nrow = I, ncol = J)  # optimal actions
  
  # Loop over all grid states (x_i, S^{(j)})
  for (i in 1:I) {
    x_val <- x_grid[i]
    
    for (j in 1:J) {
      S_val <- S_grid[j]
      
      # Action set A(x_i)
      a_max_i <- min(a_max, x_val)
      actions <- 0:a_max_i
      
      # If no inventory, only action 0
      if (a_max_i == 0) {
        V_t[i, j] <- 0
        A_t[i, j] <- 0
        next
      }
      
      # For each action a, compute Q-value
      Q_values <- numeric(length(actions))
      
      # Pre-simulate shocks for this state if you want re-use
      z_vec <- rnorm(M_mc, mean = 0, sd = 1)
      
      for (k in seq_along(actions)) {
        a <- actions[k]
        x_next <- x_val - a
        
        
        immediate <- a * S_val - c_cost * a^2
        
       
        V_future_sum <- 0
        
       
        x_next_index <- which(x_grid == x_next)
        
        for (m in 1:M_mc) {
          S_next <- S_val + sigma_hat * z_vec[m]
          V_future_sum <- V_future_sum +
            interp_V_next(V_next, x_next_index, S_next, S_grid)
        }
        
        V_future_mc <- V_future_sum / M_mc
        
        Q_values[k] <- immediate + V_future_mc
      }
      
      # Optimal value and action
      best_idx <- which.max(Q_values)
      V_t[i, j] <- Q_values[best_idx]
      A_t[i, j] <- actions[best_idx]
    }
  }
  
  V_list[[t + 1]] <- V_t
  A_list[[t + 1]] <- A_t
}

cat("Backward induction complete.\n")

############################################
#############Visualization##############################

get_opt_action <- function(t, x, S, A_list, x_grid, S_grid) {
  # t: time index (0,...,T_horizon-1)
  # x: current inventory (must be in x_grid)
  # S: current price
  # returns optimal action a_t^*(x,S) using nearest price grid
  
  A_t <- A_list[[t + 1]]          # matrix I x J
  
  
  i <- which(x_grid == x)
  if (length(i) == 0) stop("x not on grid")
  
  
  j <- which.min(abs(S_grid - S))
  
  a_opt <- A_t[i, j]
  return(a_opt)
}
##########################

# Choose initial state
x0 <- 20                    # start with 20 shares
S0 <- S_grid[ceiling(J/2)]  # start at mid price grid (or any price you like)

T_steps <- T_horizon        # number of decision steps

# Storage
time_vec <- 0:T_steps
x_path   <- numeric(T_steps + 1)
S_path   <- numeric(T_steps + 1)
a_path   <- numeric(T_steps)    # actions for t=0,...,T-1
r_path   <- numeric(T_steps)    # rewards

# init
x_path[1] <- x0
S_path[1] <- S0

for (t in 0:(T_steps - 1)) {
  i <- t + 1
  
  # get optimal action for current (t, x, S)
  a_t <- get_opt_action(
    t = t,
    x = x_path[i],
    S = S_path[i],
    A_list = A_list,
    x_grid = x_grid,
    S_grid = S_grid
  )
  
  # enforce not to sell more than inventory (just in case)
  a_t <- min(a_t, x_path[i])
  
  # immediate reward: a_t * S_t - c * a_t^2
  r_t <- a_t * S_path[i] - c_cost * a_t^2
  
  # simulate next price
  eps  <- rnorm(1)
  S_next <- S_path[i] + sigma_hat * eps
  
  # update inventory
  x_next <- x_path[i] - a_t
  
  # store
  a_path[i] <- a_t
  r_path[i] <- r_t
  S_path[i + 1] <- S_next
  x_path[i + 1] <- x_next
}

traj_df <- data.frame(
  t  = time_vec,
  x  = x_path,
  S  = S_path
)

traj_actions <- data.frame(
  t  = 0:(T_steps - 1),
  a  = a_path,
  r  = r_path
)
#(a) Price and inventory over time
# Price path
p1 <- ggplot(traj_df, aes(x = t, y = S)) +
  geom_line(color = "steelblue", size = 1.2) +
  labs(title = "Price Path ",
       x = "Time step", y = "Price S_t") +
  theme_minimal()

# Inventory path
p2 <- ggplot(traj_df, aes(x = t, y = x)) +
  geom_step(color = "darkred", size = 1.2) +
  labs(title = "Inventory Path under Optimal Policy",
       x = "Time step", y = "Inventory x_t") +
  theme_minimal()

p1
p2

#(b) Actions over time (how much you sell each step)
p3 <- ggplot(traj_actions, aes(x = t, y = a)) +
  geom_col(fill = "darkgreen") +
  labs(title = "Optimal Actions over Time",
       x = "Time step", y = "Action a_t (shares sold)") +
  theme_minimal()

p3

traj_actions$cum_r <- cumsum(traj_actions$r)

# (d) Reward over time (instantaneous + cumulative)
p4 <- ggplot(traj_actions, aes(x = t)) +
  geom_col(aes(y = r), fill = "steelblue", alpha = 0.6) +
  geom_line(aes(y = cum_r), color = "darkred", size = 1.1) +
  labs(title = "Reward Over Time",
       x = "Time step",
       y = "Reward r_t (bars) / Cumulative (line)") +
  theme_minimal()

#(c) Combined view: price + inventory + actions
#library(patchwork)

(p1 | p2) / (p3 | p4)

#######################
####Animation ##########
# Add cumulative reward
# 1) Add cumulative reward
traj_actions$cum_r <- cumsum(traj_actions$r)

# 2) Long data for line panels: Price, Inventory, Cumulative reward
lines_df <- bind_rows(
  traj_df %>%
    transmute(t, panel = "Price",             value = S),
  traj_df %>%
    transmute(t, panel = "Inventory",         value = x),
  traj_actions %>%
    transmute(t, panel = "Cumulative reward", value = cum_r)
)

# 3) Long data for actions (bars) in "Action" panel
actions_df <- traj_actions %>%
  transmute(t, panel = "Action (Shares Sold)", value = a)

# 4) Build frame-by-frame data so past values stay visible
max_frame <- max(traj_df$t)   # typically T_steps

lines_anim_df <- bind_rows(
  lapply(0:max_frame, function(f) {
    lines_df %>%
      filter(t <= f) %>%
      mutate(frame = f)
  })
)

actions_anim_df <- bind_rows(
  lapply(0:max_frame, function(f) {
    actions_df %>%
      filter(t <= f) %>%
      mutate(frame = f)
  })
)

# 5) Animation: 4 panels, bars + lines, all past values visible
p_anim <- ggplot() +
  # Actions (bars) in "Action" panel
  geom_col(
    data = actions_anim_df,
    aes(x = t, y = value, group = t),
    fill = "darkgreen",
    alpha = 0.7,
    width = 0.7
  ) +
  
  # Price, Inventory, Cumulative reward as lines in their panels
  geom_line(
    data = lines_anim_df,
    aes(x = t, y = value, color = panel, group = panel),
    size = 1.2
  ) +
  
  facet_wrap(~ panel, scales = "free_y", ncol = 2) +
  
  scale_color_manual(values = c(
    "Price"             = "steelblue",
    "Inventory"         = "darkred",
    "Cumulative reward" = "black"
  )) +
  
  labs(
    title = "Evolution under Optimal Policy — Time t = {current_frame}",
    x = "Time step",
    y = "",
    color = ""
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    strip.text      = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  
  # Each 'frame' value is a full snapshot up to that time
  transition_manual(frame)

# 6) Render GIF
anim <- animate(
  p_anim,
  renderer = gifski_renderer(),
  fps      = 10,
  nframes  = max_frame + 1,  # one frame per time step
  width    = 800,
  height   = 600
)

anim_save("optimal_execution_animation.gif", animation = anim)


