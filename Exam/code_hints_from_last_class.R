##################################################
## Project:         FIRS sampling exam hints
## Script purpose:  These are hints for the exam based on the code we 
##                  developed on the projector during the last day of class.
## Date:            2026-09-30
## Author:          Alexander Massey (alexander.massey@ign.fr)
##
## Notes:           I cleaned it up to be a bit more readable and generalized 
##                  the systematic grid sampling function we worked on. 
##                  Try different n_sim and n (per the suggested ones).
##################################################


#####################
# Sampling Exercise #
#####################

library(forestinventory)
data(zberg)

zberg_forest <- zberg |>
  filter(!is.na(basal)) |>
  select(cluster, x_terr, y_terr, stade, melange, couver, stem, basal) |>
  mutate(plot_id = row_number())

# Proposal sample 1
sample_SRS_plots <- sample(zberg_forest$plot_id, 50)
zberg_forest$in_sample_SRS_plots <- zberg_forest$plot_id %in% sample_SRS_plots

plot(
  zberg_forest |> select(x_terr, y_terr),
  col = "lightgray",
  pch = 16
)

points(
  zberg_forest[zberg_forest$in_sample_SRS_plots, c("x_terr", "y_terr")],
  col = "red",
  pch = 16
)

# Proposal sample 2
# find average plots per cluster
ave_plots_per_cluster <- zberg_forest |>
  group_by(cluster) |>
  summarise(ave_clust = n()) |>
  ungroup() |>
  summarise(mean(ave_clust)) |>
  as.numeric()

sample_SRS_cluster <- sample(unique(zberg_forest$cluster), 50/ave_plots_per_cluster)
zberg_forest$in_sample_SRS_cluster <- zberg_forest$cluster %in% sample_SRS_cluster

plot(
  zberg_forest |> select(x_terr, y_terr),
  col = "lightgray",
  pch = 16
)

points(
  zberg_forest[zberg_forest$in_sample_SRS_cluster, c("x_terr", "y_terr")],
  col = "red",
  pch = 16
)


# Proposal sample 3
# We can add a random number (float) between 0 and 1 
# and take the minimum by selected clusters.  This is cluster by SRS and 1 plot by SRS.
sampled_clusters <- sample(unique(zberg_forest$cluster), 50)
sample_2stage <- zberg_forest |>
  filter(cluster %in% sampled_clusters) |>
  mutate(rand.number = runif(n())) |>      # Add runif to all plots in selected clusters
  group_by(cluster) |>
  filter(rand.number == min(rand.number)) |>
  ungroup() |>
  select(plot_id) |>
  unlist()

zberg_forest$in_sample_2stage <- zberg_forest$plot_id %in% sample_2stage

plot(
  zberg_forest |> select(x_terr, y_terr),
  col = "lightgray",
  pch = 16
)

points(
  zberg_forest[zberg_forest$in_sample_2stage, c("x_terr", "y_terr")],
  col = "red",
  pch = 16
)


# Proposal sample 4
#  ... Try to adapt proposal 3 to include stratification by categories
#  you create like clusters with high and low mean stem density 
#  and estimate p_h from the full inventory from last year.  Be creative.


########################
# Estimation Exercise  #
########################

# Function that takes coordinates of sample as input and give Y(x) as output
Y <- function(x1, x2, true_beta = c(30, 13, -6, -4, 3, 2)) {
  Z_x <- cbind(
    1,
    x1,
    x2,
    x1^2,
    x1 * x2,
    x2^2
  )
  
  noise <- 6 * cos(pi * x1) * sin(2 * pi * x2)
  
  as.numeric(Z_x %*% true_beta + noise)
}


# Function that takes coordinates of sample as input and give stratum id as output
get_stratum <- function(x1, x2) {
  ifelse(x2 < 2,
         ifelse(x1 < 1.5, "S1", "S2"), # lower half of x2
         ifelse(x1 < 1.5, "S3", "S4") # upper half of x2
  )
}



library(cubature)
result <- adaptIntegrate(
  f = function(x) {Y(x[1], x[2])},
  lowerLimit = c(0, 0),
  upperLimit = c(3, 4)
)
result$integral

# Total = 542
# True Mean = 542/(3*4) = 542/12 = 45.16667
true_mean <- 45.16667




# Inputs
# Note: strategy2() function only takes values n = 12*(x^2) (x = 1,2,3,etc)
n <- 768  # try 12, 48, 108, 192, 300, 432, 768, etc...
n_sim <- 1000


# Simple uniform iid sampling with one phase estimator
strategy1 <- function(){
  df <- data.frame(x1 = runif(n, 0, 3), x2 = runif(n, 0, 4))
  df$Y <- Y(df$x1, df$x2)
  n <- nrow(df)
  
  est <- mean(df$Y)
  var <- var(df$Y)/n
  lower <- est - 1.959964*sqrt(var)
  upper <- est + 1.959964*sqrt(var)
  in_CI <- ifelse(true_mean<upper & true_mean > lower, 1, 0)
  
  out <- list(type = "Uniform iid", est = est, var = var, in_CI = in_CI)
}

#initialize data.frame
uniform_iid_df <- data.frame(type = "Uniform iid",
                             est = rep(NA, n_sim), 
                             var = rep(NA, n_sim), in_CI = rep(NA, n_sim))

for(i in 1:n_sim){
  input <- strategy1()
  uniform_iid_df[i,] <- input
}



# Systematic sampling
strategy2 <- function(){
  
  # factor to find number of square grid cells in each direction based on n
  k <- sqrt(n / 12)
  
  # Randomly select reference point in first grid cell
  # to construct the systematic grid
  generate_grid <- function(){
    x1 <- runif(1, 0, 1/k) + (1/k) * (0:(3*k - 1))
    x2 <- runif(1, 0, 1/k) + (1/k) * (0:(4*k - 1))
    
    expand.grid(x1 = x1, x2 = x2)
  }
  
  df <- generate_grid()
  df$Y <- Y(df$x1, df$x2)
  n <- nrow(df)
  
  est <- mean(df$Y)
  var <- var(df$Y) / n
  
  lower <- est - 1.959964 * sqrt(var)
  upper <- est + 1.959964 * sqrt(var)
  in_CI <- ifelse(true_mean < upper & true_mean > lower, 1, 0)
  
  return(list(
    type = "Systematic",
    est = est,
    var = var,
    in_CI = in_CI
  ))
}
  
 

#initialize data.frame
grid_df <- data.frame(type = "Systematic",
                             est = rep(NA, n_sim), 
                             var = rep(NA, n_sim), in_CI = rep(NA, n_sim))

for(i in 1:n_sim){
  input <- strategy2()
  grid_df[i,] <- input
}


# Uniform iid sampling with stratification (4 equal sized quadrants...)
strategy3 <- function(){
  df <- rbind(
    data.frame(x1 = runif(n/4, 0, 1.5), x2 = runif(n/4, 0, 2)), # first quadrant
    data.frame(x1 = runif(n/4, 1.5, 3), x2 = runif(n/4, 0, 2)),
    data.frame(x1 = runif(n/4, 0, 1.5), x2 = runif(n/4, 2, 4)),
    data.frame(x1 = runif(n/4, 1.5, 3), x2 = runif(n/4, 2, 4))
  )
  df$stratum <- get_stratum(df$x1, df$x2)
  df$Y <- Y(df$x1, df$x2)
  df$strat_weight <- 1/4
  n <- nrow(df)
  df_strat_level <- df |> 
    group_by(stratum) |> 
    summarise(strat_weight = unique(strat_weight),
              strat_mean = strat_weight*mean(Y),      # p_h*mean(Y)
              strat_var = strat_weight^2*var(Y)/n())  # p_h^2*(1/n_h)*var(Y) 

  est <- sum(df_strat_level$strat_mean)
  var <- sum(df_strat_level$strat_var)
  lower <- est - 1.959964*sqrt(var)
  upper <- est + 1.959964*sqrt(var)
  in_CI <- ifelse(true_mean<upper & true_mean > lower, 1, 0)
  
  out <- list(type = "Stratified", est = est, var = var, in_CI = in_CI)
}

#initialize data.frame
stratified_df <- data.frame(type = "Stratified",
                             est = rep(NA, n_sim), 
                             var = rep(NA, n_sim), in_CI = rep(NA, n_sim))

for(i in 1:n_sim){
  input <- strategy3()
  stratified_df[i,] <- input
}


##  Post-stratification with external variance formula (i.e. "plug-in")
strategy4 <- function(){
  df <- data.frame(x1 = runif(n, 0, 3), x2 = runif(n, 0, 4))
  df$Z_x <- get_stratum(df$x1, df$x2)
  df$Y_x <- Y(df$x1, df$x2)
  n <- nrow(df)
  
  # Post-stratification IS anova (lm with a categorical variable)
  mod <- lm(Y_x ~ -1 + Z_x, data = df)
  df$R_x <- residuals(mod)

  est <- as.numeric(c(.25,.25,.25,.25)%*%coef(mod))
  var <- var(df$R_x)/n
  lower <- est - 1.959964*sqrt(var)
  upper <- est + 1.959964*sqrt(var)
  in_CI <- ifelse(true_mean<upper & true_mean > lower, 1, 0)
  
  out <- list(type = "Uniform iid", est = est, var = var, in_CI = in_CI)
}


#initialize data.frame
post_df <- data.frame(type = "Post-strat",
                      est = rep(NA, n_sim), 
                      var = rep(NA, n_sim), in_CI = rep(NA, n_sim))

for(i in 1:n_sim){
  input <- strategy4()
  post_df[i,] <- input
}


results <- rbind(
  data.frame(
    type = "Uniform iid",
    n = n,
    true_value = true_mean,
    E_estimate = mean(uniform_iid_df$est),
    Var_estimator = var(uniform_iid_df$est),
    E_variance_estimate = mean(uniform_iid_df$var),
    coverage = mean(uniform_iid_df$in_CI)
  ),
  
  data.frame(
    type = "Systematic",
    n = n,
    true_value = true_mean,
    E_estimate = mean(grid_df$est),
    Var_estimator = var(grid_df$est),
    E_variance_estimate = mean(grid_df$var),
    coverage = mean(grid_df$in_CI)
  ),
  
  data.frame(
    type = "Stratified",
    n = n,
    true_value = true_mean,
    E_estimate = mean(stratified_df$est),
    Var_estimator = var(stratified_df$est),
    E_variance_estimate = mean(stratified_df$var),
    coverage = mean(stratified_df$in_CI)
  ),
  
  data.frame(
    type = "Post-strat",
    n = n,
    true_value = true_mean,
    E_estimate = mean(post_df$est),
    Var_estimator = var(post_df$est),
    E_variance_estimate = mean(post_df$var),
    coverage = mean(post_df$in_CI)
  )
)

results
