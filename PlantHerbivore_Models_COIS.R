# Load required packages
library(deSolve)
library(tidyverse)
library(patchwork)
library(gtools)
library(viridis)
library(hillR)
pal <- viridis(256, option = "D")


# Trait matching: attack rate ####
# Species counts
n_plants <- 10
n_herbs <- 10

# Trait values
toughness <- seq(30, 50, length.out=n_plants)
mandible <- seq(2.5, 4.5, length.out=n_herbs)

# Normalize traits
mean_toughness <- mean(toughness)
mean_mandible <- mean(mandible)


# Trait-based attack rates using Gaussian trait matching
beta <- 5
a0 <- 0.2
T <- .05 * n_plants  # ensures total attack rate matches a generalist's base rate

A_raw <- matrix(0, nrow = n_plants, ncol = n_herbs)
for (i in 1:n_plants) {
  for (j in 1:n_herbs) {
    ratio <- (mandible[j] / mean_mandible) / (toughness[i] / mean_toughness)
    similarity <- exp(-beta * (ratio - 1)^2)
    A_raw[i, j] <- similarity
  }
}

# Normalize so each herbivore's total attack strength across plants = T
A_match <- sweep(A_raw, 2, colSums(A_raw), FUN = "/") * T
A <- A_match
colSums(A)
rowSums(A)

colSums(A_raw)
rowSums(A_raw)
sum(A)

# Fixed parameters
R <- rep(2, n_plants)
C <- matrix(0.0, nrow = n_plants, ncol = n_plants); diag(C) <- 0.04
HNDL <- matrix(runif(n_plants * n_herbs, 0, 0), nrow = n_plants)
E <- matrix(runif(n_plants * n_herbs, 0.1, 0.1), nrow = n_plants)
M <- rep(1, n_herbs)
g <- rep(0.00, n_herbs)  # herbivore self-limitation rate

# Initial state
state <- c(P = rep(15, n_plants), H = rep(3, n_herbs))
# state <- c(P = rep(15, n_plants), H = c(3,0,0,0,0,0,0,0,0,3))

times <- seq(0, 100, by = .1)

# Model function
competition_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    P <- state[1:n_plants]
    H <- state[(n_plants + 1):(n_plants + n_herbs)]
    dP <- numeric(n_plants)
    dH <- numeric(n_herbs)
    
    for (i in 1:n_plants) {
      comp_sum <- sum(C[i, ] * P)
      herb_sum <- 0
      for (j in 1:n_herbs) {
        denom <- 1 + sum(HNDL[i, j] * A[i, j] * P)
        herb_sum <- herb_sum + (A[i, j] * H[j]) / denom
      }
      dP[i] <- P[i] * ((R[i] - comp_sum) - herb_sum)
    }
    
    for (j in 1:n_herbs) {
      numer <- sum(E[, j] * A[, j] * P)
      denom <- 1 + sum(A[, j] * HNDL[, j] * P)
      dH[j] <- H[j] * (numer / denom - M[j] - g[j] * H[j])
    }
    
    list(c(dP, dH))
  })
}

# Solve ODE
out <- ode(y = state, times = times, func = competition_model, parms = list())

# Format for ggplot
out_df <- as.data.frame(out)
colnames(out_df) <- c("time", paste0("P", 1:n_plants), paste0("H", 1:n_herbs))
out_long <- out_df %>%
  pivot_longer(-time, names_to = "species", values_to = "value") %>%
  mutate(type = ifelse(str_detect(species, "^P"), "Plant", "Herbivore"),
         species=as_factor(species),
         spp=factor(species, levels=mixedsort(levels(species))))

## make plot ####
pop1 <- ggplot(data=out_long %>% filter(type=="Plant"), aes(x = time, y = value, color = spp)) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis(option="viridis", discrete=T, direction=-1, begin=.2, end=.8) +
  labs(#title = "Plant density (g/m2)",
       x = "Time",
       y = "Plant biomass (g/m2)",
       color="Plant spp") +
  ylim(0,30)+
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

pop2 <- ggplot(data=out_long %>% filter(type=="Herbivore"), 
                aes(x = time, y = value, color = spp)) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis(option="inferno", discrete=T, direction=-1, begin=.2, end=.8) +
  labs(title = "A) Trait matching",
       x = "Time",
       y = "Herbivore density (#/m2)",
       color="Herbivore spp") +
  ylim(0,13)+
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

popplot <- pop2/pop1

# Assume this is your output from the model at time = 100
# Filter for t = 100
abund_at_100 <- out_long %>%
  filter(time == 100) %>%
  mutate(
    trait = case_when(
      type == "Plant" ~ toughness[as.numeric(str_extract(species, "\\d+"))],
      type == "Herbivore" ~ mandible[as.numeric(str_extract(species, "\\d+"))]
    )
  )

# Plot
pop1.1 <- ggplot(abund_at_100 %>% filter(type=="Plant"), 
       aes(x = trait, y = value, color=trait)) +
  geom_point(size = 3) +
  #geom_line(aes(group = type), linetype = "dashed", alpha = 0.6) +
  scale_color_viridis(option="viridis", discrete=F, direction=-1, begin=.2, end=.8) +
  #facet_wrap(~type, scales = "free_x") +
  labs(
    #title = "Trait Value vs. Abundance at Time = 100",
    x = "LDMC (%)",
    y = "Plant biomass (g/m2)",
    color = "LDMC (%)"
  ) +
  theme_bw(base_size = 14)

pop1.2 <- ggplot(abund_at_100 %>% filter(type=="Herbivore"), 
       aes(x = trait, y = value, color=trait)) +
  geom_point(size = 3) +
  #geom_line(aes(group = type), linetype = "dashed", alpha = 0.6) +
  scale_color_viridis(option="inferno", discrete=F, direction=-1, begin=.2, end=.8) +
  #facet_wrap(~type, scales = "free_x") +
  labs(
    title = "Trait Value vs. Abund at T=100",
    x = "Mandible strength",
    y = "Herbivore abundance (#/m2)",
    color = "Mandible"
  ) +
  theme_bw(base_size = 14)

popplot12 <- pop1.2/pop1.1

Arate_plot <- popplot|popplot12
# ggsave(Arate_plot, file="Arate_plot.tiff", width=8, height=6, dpi=300, compression = "lzw")




# Spaak & De Laender (2020) N and F for MacArthur CRM #######
# Box 1 [3] and Appendix B [4] multispecies extensions
#
# Key growth rates needed per species i (Eq. 7 and 8 [3]):
#   fi(0, 0)              = intrinsic growth rate (no competitors)
#   fi(0, N-i*)           = invasion growth rate (species i rare, S-1 at equil.)
#   fi(sum_j cij*N-i*j, 0) = no-niche growth rate (species i alone at converted
#                             equilibrium density of competitor community)
#
# Conversion factor cij (Eq. 14 [3], Box 1):
#   cij = sqrt(sum(u_jl^2) / sum(u_il^2))
#       = sqrt(sum(A[,j]^2) / sum(A[,i]^2))
#   where A[,i] is the attack/consumption rate vector of herbivore i [1]

# ── Step 1: Compute conversion factors cij from A matrix ─────────────────────
# Following Eq. 14 Box 1 [3]: cij = ||uj|| / ||ui||
# cij converts density of species j into equivalent density of species i
compute_conversion_factors <- function(A) {
  n_herbs <- ncol(A)
  # L2 norm of each herbivore's consumption vector across plants
  norms <- sqrt(colSums(A^2))
  
  # cij matrix: converts j into i units
  # cij = ||uj|| / ||ui||  [3]
  C_conv <- outer(norms, norms, FUN = function(ui, uj) uj / ui)
  diag(C_conv) <- 1
  return(C_conv)
}

# ── Step 2: Run S-1 community to equilibrium (excluding focal species i) ──────
run_subcommunity <- function(focal_herb, parms_fixed, times_long) {
  state_sub <- c(P = rep(15, n_plants), H = rep(0, n_herbs))
  
  # All herbivores except focal species present
  for (h in 1:n_herbs) {
    if (h != focal_herb) state_sub[n_plants + h] <- 3
  }
  
  out <- tryCatch(
    ode(y = state_sub, times = times_long,
        func = competition_model, parms = parms_fixed,
        method = "lsoda"),
    error = function(e) NULL
  )
  return(out)
}

# ── Step 3: Compute the three key growth rates for each species ───────────────
# Following Eqs. 7, 8 [3] and Appendix B [4]

compute_growth_rates <- function(focal_herb, C_conv, parms_fixed,
                                 times_long, times_short, invasion_density = 0.01) {
  
  # --- (a) Intrinsic growth rate: fi(0, 0) ---
  # Growth of focal species alone at very low density, no competitors
  state_intrinsic <- c(P = rep(15, n_plants), H = rep(0, n_herbs))
  state_intrinsic[n_plants + focal_herb] <- invasion_density
  
  out_int <- tryCatch(
    ode(y = state_intrinsic, times = times_short,
        func = competition_model, parms = parms_fixed,
        method = "lsoda"),
    error = function(e) NULL
  )
  if (is.null(out_int)) return(NULL)
  
  df_int    <- as.data.frame(out_int)
  h_col     <- paste0("H", focal_herb)
  N_t0      <- df_int[1, h_col]
  N_tend    <- df_int[nrow(df_int), h_col]
  t_span    <- df_int$time[nrow(df_int)] - df_int$time[1]
  fi_00     <- log(N_tend / N_t0) / t_span
  
  # --- (b) Invasion growth rate: fi(0, N-i*) ---
  # Run S-1 community to equilibrium, then introduce focal species as rare
  sub_out <- run_subcommunity(focal_herb, parms_fixed, times_long)
  if (is.null(sub_out)) return(NULL)
  
  # Extract S-1 equilibrium state
  equil_sub <- as.numeric(sub_out[nrow(sub_out), 2:(n_plants + n_herbs + 1)])
  names(equil_sub) <- c(paste0("P", 1:n_plants), paste0("H", 1:n_herbs))
  
  # Store N-i* equilibrium densities (needed for no-niche rate)
  N_minus_i_star <- equil_sub[paste0("H", 1:n_herbs)]
  N_minus_i_star[focal_herb] <- 0   # focal species absent in subcommunity
  
  # Introduce focal species as rare into S-1 equilibrium
  state_inv <- equil_sub
  state_inv[n_plants + focal_herb] <- invasion_density
  
  out_inv <- tryCatch(
    ode(y = state_inv, times = times_short,
        func = competition_model, parms = parms_fixed,
        method = "lsoda"),
    error = function(e) NULL
  )
  if (is.null(out_inv)) return(NULL)
  
  df_inv    <- as.data.frame(out_inv)
  N_t0_inv  <- df_inv[1, h_col]
  N_tend_inv <- df_inv[nrow(df_inv), h_col]
  fi_0Nj    <- log(N_tend_inv / N_t0_inv) / t_span
  
  # --- (c) No-niche growth rate: fi(sum_j cij * N-i*_j, 0) ---
  # Growth of focal species ALONE at converted equilibrium density
  # sum_j cij * N-i*_j  [3] Eq. 7, Appendix B [4]
  converted_density <- sum(C_conv[focal_herb, ] * N_minus_i_star)
  converted_density <- max(converted_density, invasion_density)
  
  state_noniche <- c(P = equil_sub[paste0("P", 1:n_plants)],
                     H = rep(0, n_herbs))
  state_noniche[n_plants + focal_herb] <- converted_density
  
  out_nn <- tryCatch(
    ode(y = state_noniche, times = times_short,
        func = competition_model, parms = parms_fixed,
        method = "lsoda"),
    error = function(e) NULL
  )
  if (is.null(out_nn)) return(NULL)
  
  df_nn     <- as.data.frame(out_nn)
  N_t0_nn   <- df_nn[1, h_col]
  N_tend_nn <- df_nn[nrow(df_nn), h_col]
  fi_cNj    <- log(N_tend_nn / N_t0_nn) / t_span
  
  return(list(
    fi_00  = fi_00,
    fi_0Nj = fi_0Nj,
    fi_cNj = fi_cNj,
    N_minus_i_star = N_minus_i_star,
    converted_density = converted_density
  ))
}

# ── Step 4: Compute N and F for each species (Eqs. 7 and 8 [3]) ──────────────
compute_N_and_F <- function(gr) {
  # Ni = (fi(0, N-i*) - fi(cij*N-i*_j, 0)) /
  #      (fi(0, 0)    - fi(cij*N-i*_j, 0))    [3] Eq. 7
  Ni <- (gr$fi_0Nj - gr$fi_cNj) / (gr$fi_00 - gr$fi_cNj)
  
  # Fi = fi(cij*N-i*_j, 0) / fi(0, 0)          [3] Eq. 8
  Fi <- gr$fi_cNj / gr$fi_00
  
  return(list(N = Ni, F = Fi))
}

# ── Step 5: Coexistence criterion (Eq. 6 [3]) ─────────────────────────────────
# Species i persists if: -Fi < Ni / (1 - Ni)
check_coexistence_spaak <- function(Ni, Fi) {
  lhs <- -Fi
  rhs <- Ni / (1 - Ni)
  persists <- lhs < rhs
  return(list(persists = persists, lhs = lhs, rhs = rhs))
}

# ── Step 6: Run for all herbivore species ─────────────────────────────────────
times_long  <- seq(0, 1000, by = 0.5)   # long run to reach S-1 equilibrium
times_short <- seq(0, 10,   by = 0.1)   # short window to measure growth rates

# Use your fitted parameters; defaults used here for illustration
parms_fixed <- list(
  n_plants = n_plants, n_herbs = n_herbs,
  A    = A,
  R    = rep(2,   n_plants),
  C    = {cc <- matrix(0, n_plants, n_plants); diag(cc) <- 0.04; cc},
  HNDL = matrix(0,   n_plants, n_herbs),
  E    = matrix(0.1, n_plants, n_herbs),
  M    = rep(1,   n_herbs),
  g    = rep(0.001, n_herbs)
)

# Compute conversion factors from A matrix [1]
C_conv <- compute_conversion_factors(A)

# Run for all species
results_NF <- map_dfr(1:n_herbs, function(i) {
  gr <- compute_growth_rates(i, C_conv, parms_fixed, times_long, times_short)
  if (is.null(gr)) return(tibble(herb = i, N = NA, F_val = NA,
                                 persists = NA, lhs = NA, rhs = NA))
  nf  <- compute_N_and_F(gr)
  cx  <- check_coexistence_spaak(nf$N, nf$F)
  tibble(
    herb      = i,
    fi_00     = gr$fi_00,
    fi_0Nj    = gr$fi_0Nj,
    fi_cNj    = gr$fi_cNj,
    N         = nf$N,
    F_val     = nf$F,
    lhs       = cx$lhs,
    rhs       = cx$rhs,
    persists  = cx$persists
  )
})

print(results_NF)

## Step 7: Visualize N and F in coexistence space ####

# Transform F to 2021 Spaak et al. convention [5]
# Fi_new = -Fi_old / (1 - eta_i / mu_i)
# where eta_i = fi(sum cij * N*j, 0) and mu_i = fi(0,0)

results_NF <- results_NF %>%
  mutate(
    F_2021 = -F_val / (1 - (fi_cNj / fi_00))
  )

results_NF<- results_NF %>%
  mutate(
    nudge_x_val = ifelse(herb %in% c(5, 6, 7, 8), -0.03, 0.03)
  )

NF_plot_traitmatch <- ggplot(results_NF, aes(x = N, y = F_2021, 
                                                            shape = persists, 
                                                            label = herb,
                                                            color=as.factor(herb))) +
  geom_point(size = 4, stroke=2) +
  geom_text(aes(nudge_x = nudge_x_val), size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  # Coexistence boundary: -F = N/(1-N)  =>  F = -N/(1-N)
  stat_function(fun = function(x) x / (1 - x),
                color = "black", linetype = "dashed", linewidth=2, xlim = c(-0.5, 0.95)) +
  scale_color_viridis(option="inferno", discrete=T, direction=-1, begin=.2, end=.8) +
  #scale_fill_viridis(option="inferno", discrete=T, direction=-1, begin=.2, end=.8) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21),
                     labels = c("Excluded", "Persists")) +
  coord_cartesian(ylim = c(-.01, .062), xlim=c(-.55, .1)) +
  labs(x = "Niche difference (N)",
       y = "Fitness difference (F)",
       shape = "Outcome",
       title = "Trait matching",
       subtitle = "N-F coexistence space") +
  theme_bw(base_size = 14) +
  theme(legend.position = c(.4,.4))+
  guides(color="none")




# Strong mandible model no Handling Time ####
# Species counts
n_plants <- 10
n_herbs <- 10

# Trait values
toughness <- seq(30, 50, length.out=n_plants)
mandible <- seq(2.5, 4.5, length.out=n_herbs)

# Normalize traits
mean_toughness <- mean(toughness)
mean_mandible <- mean(mandible)

# Base attack rate
a0 <- 0.04
T <- .1 * n_plants  # ensures total attack rate matches a generalist's base rate

# Strong mandibles are better
A_strong <- outer(toughness, mandible, function(t, m) (m / mean(m)) / (t / mean(t))) * a0 * T
A <- A_strong

colSums(A)
rowSums(A)
sum(A)
A_fixed <- A

# Fixed parameters
R <- rep(2, n_plants)
C <- matrix(0.0, nrow = n_plants, ncol = n_plants); diag(C) <- 0.04
HNDL <- matrix(runif(n_plants * n_herbs, 0, 0), nrow = n_plants)
E <- matrix(runif(n_plants * n_herbs, 0.1, 0.1), nrow = n_plants)
M <- rep(1, n_herbs)
g <- rep(0.01, n_herbs)  # include a small herbivore self-limitation rate

# Initial state
state <- c(P = rep(15, n_plants), H = rep(3, n_herbs))

times <- seq(0, 100, by = .1)

# Model function
competition_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    P <- state[1:n_plants]
    H <- state[(n_plants + 1):(n_plants + n_herbs)]
    dP <- numeric(n_plants)
    dH <- numeric(n_herbs)
    
    for (i in 1:n_plants) {
      comp_sum <- sum(C[i, ] * P)
      herb_sum <- 0
      for (j in 1:n_herbs) {
        denom <- 1 + sum(HNDL[i, j] * A[i, j] * P)
        herb_sum <- herb_sum + (A[i, j] * H[j]) / denom
      }
      dP[i] <- P[i] * ((R[i] - comp_sum) - herb_sum)
    }
    
    for (j in 1:n_herbs) {
      numer <- sum(E[, j] * A[, j] * P)
      denom <- 1 + sum(A[, j] * HNDL[, j] * P)
      dH[j] <- H[j] * (numer / denom - M[j] - g[j] * H[j])
    }
    
    list(c(dP, dH))
  })
}

# Solve ODE
out <- ode(y = state, times = times, func = competition_model, parms = list())

# Format for ggplot
out_df <- as.data.frame(out)
colnames(out_df) <- c("time", paste0("P", 1:n_plants), paste0("H", 1:n_herbs))
out_long <- out_df %>%
  pivot_longer(-time, names_to = "species", values_to = "value") %>%
  mutate(type = ifelse(str_detect(species, "^P"), "Plant", "Herbivore"),
         species=as_factor(species),
         spp=factor(species, levels=mixedsort(levels(species))))

## make plot ####
pop1a <- ggplot(data=out_long %>% filter(type=="Plant"), aes(x = time, y = value, color = spp)) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis(option="viridis", discrete=T, direction=-1, begin=.2, end=.8) +
  labs(#title = "Plant density (g/m2)",
    x = "Time",
    y = "Plant biomass (g/m2)",
    color="Plant spp") +
  ylim(0,30)+
  theme_bw(base_size = 16)

pop2a <- ggplot(data=out_long %>% filter(type=="Herbivore"), 
               aes(x = time, y = value, color = spp)) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis(option="inferno", discrete=T, direction=-1, begin=.2, end=.8) +
  labs(title = "B) Strong mandibles",
       x = "Time",
       y = "Herbivore density (#/m2)",
       color="Herbivore spp") +
  ylim(0,13)+
  theme_bw(base_size = 16)

popplota <- pop2a/pop1a

# Assume this is your output from the model at time = 100
# Filter for t = 100
abund_at_100 <- out_long %>%
  filter(time == 100) %>%
  mutate(
    trait = case_when(
      type == "Plant" ~ toughness[as.numeric(str_extract(species, "\\d+"))],
      type == "Herbivore" ~ mandible[as.numeric(str_extract(species, "\\d+"))]
    )
  )

# Plot
pop1.1a <- ggplot(abund_at_100 %>% filter(type=="Plant"), 
                 aes(x = trait, y = value, color=spp)) +
  geom_point(size = 3) +
  #geom_line(aes(group = type), linetype = "dashed", alpha = 0.6) +
  scale_color_viridis(option="viridis", discrete=T, direction=-1, begin=.2, end=.8) +
  #facet_wrap(~type, scales = "free_x") +
  labs(
    #title = "Trait Value vs. Abundance at Time = 100",
    x = "LDMC (%)",
    y = "Plant biomass (g/m2)",
    color = "LDMC (%)"
  ) +
  theme_bw(base_size = 16)

pop1.2a <- ggplot(abund_at_100 %>% filter(type=="Herbivore"), 
                 aes(x = trait, y = value, color=spp)) +
  geom_point(size = 3) +
  #geom_line(aes(group = type), linetype = "dashed", alpha = 0.6) +
  scale_color_viridis(option="inferno", discrete=T, direction=-1, begin=.2, end=.8) +
  #facet_wrap(~type, scales = "free_x") +
  labs(
    title = "Trait Value vs. Abund at T=100",
    x = "Mandible strength",
    y = "Herbivore density (#/m2)",
    color = "Mandible"
  ) +
  theme_bw(base_size = 16)

popplot12a <- pop1.2a/pop1.1a

Strong_plot <- popplota|popplot12a
#ggsave(Strong_plot, file="Strong_plot.tiff", width=8, height=6, dpi=300, compression = "lzw")


# Spaak & De Laender (2020) N and F for MacArthur CRM #######
# Box 1 [3] and Appendix B [4] multispecies extensions
#
# Key growth rates needed per species i (Eq. 7 and 8 [3]):
#   fi(0, 0)              = intrinsic growth rate (no competitors)
#   fi(0, N-i*)           = invasion growth rate (species i rare, S-1 at equil.)
#   fi(sum_j cij*N-i*j, 0) = no-niche growth rate (species i alone at converted
#                             equilibrium density of competitor community)
#
# Conversion factor cij (Eq. 14 [3], Box 1):
#   cij = sqrt(sum(u_jl^2) / sum(u_il^2))
#       = sqrt(sum(A[,j]^2) / sum(A[,i]^2))
#   where A[,i] is the attack/consumption rate vector of herbivore i [1]

# ── Step 1: Compute conversion factors cij from A matrix ─────────────────────
# Following Eq. 14 Box 1 [3]: cij = ||uj|| / ||ui||
# cij converts density of species j into equivalent density of species i
compute_conversion_factors <- function(A) {
  n_herbs <- ncol(A)
  # L2 norm of each herbivore's consumption vector across plants
  norms <- sqrt(colSums(A^2))
  
  # cij matrix: converts j into i units
  # cij = ||uj|| / ||ui||  [3]
  C_conv <- outer(norms, norms, FUN = function(ui, uj) uj / ui)
  diag(C_conv) <- 1
  return(C_conv)
}

# ── Step 2: Run S-1 community to equilibrium (excluding focal species i) ──────
run_subcommunity <- function(focal_herb, parms_fixed, times_long) {
  state_sub <- c(P = rep(15, n_plants), H = rep(0, n_herbs))
  
  # All herbivores except focal species present
  for (h in 1:n_herbs) {
    if (h != focal_herb) state_sub[n_plants + h] <- 3
  }
  
  out <- tryCatch(
    ode(y = state_sub, times = times_long,
        func = competition_model, parms = parms_fixed,
        method = "lsoda"),
    error = function(e) NULL
  )
  return(out)
}

# ── Step 3: Compute the three key growth rates for each species ───────────────
# Following Eqs. 7, 8 [3] and Appendix B [4]

compute_growth_rates <- function(focal_herb, C_conv, parms_fixed,
                                 times_long, times_short, invasion_density = 0.01) {
  
  # --- (a) Intrinsic growth rate: fi(0, 0) ---
  # Growth of focal species alone at very low density, no competitors
  state_intrinsic <- c(P = rep(15, n_plants), H = rep(0, n_herbs))
  state_intrinsic[n_plants + focal_herb] <- invasion_density
  
  out_int <- tryCatch(
    ode(y = state_intrinsic, times = times_short,
        func = competition_model, parms = parms_fixed,
        method = "lsoda"),
    error = function(e) NULL
  )
  if (is.null(out_int)) return(NULL)
  
  df_int    <- as.data.frame(out_int)
  h_col     <- paste0("H", focal_herb)
  N_t0      <- df_int[1, h_col]
  N_tend    <- df_int[nrow(df_int), h_col]
  t_span    <- df_int$time[nrow(df_int)] - df_int$time[1]
  fi_00     <- log(N_tend / N_t0) / t_span
  
  # --- (b) Invasion growth rate: fi(0, N-i*) ---
  # Run S-1 community to equilibrium, then introduce focal species as rare
  sub_out <- run_subcommunity(focal_herb, parms_fixed, times_long)
  if (is.null(sub_out)) return(NULL)
  
  # Extract S-1 equilibrium state
  equil_sub <- as.numeric(sub_out[nrow(sub_out), 2:(n_plants + n_herbs + 1)])
  names(equil_sub) <- c(paste0("P", 1:n_plants), paste0("H", 1:n_herbs))
  
  # Store N-i* equilibrium densities (needed for no-niche rate)
  N_minus_i_star <- equil_sub[paste0("H", 1:n_herbs)]
  N_minus_i_star[focal_herb] <- 0   # focal species absent in subcommunity
  
  # Introduce focal species as rare into S-1 equilibrium
  state_inv <- equil_sub
  state_inv[n_plants + focal_herb] <- invasion_density
  
  out_inv <- tryCatch(
    ode(y = state_inv, times = times_short,
        func = competition_model, parms = parms_fixed,
        method = "lsoda"),
    error = function(e) NULL
  )
  if (is.null(out_inv)) return(NULL)
  
  df_inv    <- as.data.frame(out_inv)
  N_t0_inv  <- df_inv[1, h_col]
  N_tend_inv <- df_inv[nrow(df_inv), h_col]
  fi_0Nj    <- log(N_tend_inv / N_t0_inv) / t_span
  
  # --- (c) No-niche growth rate: fi(sum_j cij * N-i*_j, 0) ---
  # Growth of focal species ALONE at converted equilibrium density
  # sum_j cij * N-i*_j  [3] Eq. 7, Appendix B [4]
  converted_density <- sum(C_conv[focal_herb, ] * N_minus_i_star)
  converted_density <- max(converted_density, invasion_density)
  
  state_noniche <- c(P = equil_sub[paste0("P", 1:n_plants)],
                     H = rep(0, n_herbs))
  state_noniche[n_plants + focal_herb] <- converted_density
  
  out_nn <- tryCatch(
    ode(y = state_noniche, times = times_short,
        func = competition_model, parms = parms_fixed,
        method = "lsoda"),
    error = function(e) NULL
  )
  if (is.null(out_nn)) return(NULL)
  
  df_nn     <- as.data.frame(out_nn)
  N_t0_nn   <- df_nn[1, h_col]
  N_tend_nn <- df_nn[nrow(df_nn), h_col]
  fi_cNj    <- log(N_tend_nn / N_t0_nn) / t_span
  
  return(list(
    fi_00  = fi_00,
    fi_0Nj = fi_0Nj,
    fi_cNj = fi_cNj,
    N_minus_i_star = N_minus_i_star,
    converted_density = converted_density
  ))
}

# ── Step 4: Compute N and F for each species (Eqs. 7 and 8 [3]) ──────────────
compute_N_and_F <- function(gr) {
  # Ni = (fi(0, N-i*) - fi(cij*N-i*_j, 0)) /
  #      (fi(0, 0)    - fi(cij*N-i*_j, 0))    [3] Eq. 7
  Ni <- (gr$fi_0Nj - gr$fi_cNj) / (gr$fi_00 - gr$fi_cNj)
  
  # Fi = fi(cij*N-i*_j, 0) / fi(0, 0)          [3] Eq. 8
  Fi <- gr$fi_cNj / gr$fi_00
  
  return(list(N = Ni, F = Fi))
}

# ── Step 5: Coexistence criterion (Eq. 6 [3]) ─────────────────────────────────
# Species i persists if: -Fi < Ni / (1 - Ni)
check_coexistence_spaak <- function(Ni, Fi) {
  lhs <- -Fi
  rhs <- Ni / (1 - Ni)
  persists <- lhs < rhs
  return(list(persists = persists, lhs = lhs, rhs = rhs))
}

# ── Step 6: Run for all herbivore species ─────────────────────────────────────
times_long  <- seq(0, 1000, by = 0.5)   # long run to reach S-1 equilibrium
times_short <- seq(0, 10,   by = 0.1)   # short window to measure growth rates

# Use your fitted parameters; defaults used here for illustration
parms_fixed <- list(
  n_plants = n_plants, n_herbs = n_herbs,
  A    = A,
  R    = rep(2,   n_plants),
  C    = {cc <- matrix(0, n_plants, n_plants); diag(cc) <- 0.04; cc},
  HNDL = matrix(0,   n_plants, n_herbs),
  E    = matrix(0.1, n_plants, n_herbs),
  M    = rep(1,   n_herbs),
  g    = rep(0.001, n_herbs)
)

# Compute conversion factors from A matrix [1]
C_conv <- compute_conversion_factors(A)

# Run for all species
results_NF <- map_dfr(1:n_herbs, function(i) {
  gr <- compute_growth_rates(i, C_conv, parms_fixed, times_long, times_short)
  if (is.null(gr)) return(tibble(herb = i, N = NA, F_val = NA,
                                 persists = NA, lhs = NA, rhs = NA))
  nf  <- compute_N_and_F(gr)
  cx  <- check_coexistence_spaak(nf$N, nf$F)
  tibble(
    herb      = i,
    fi_00     = gr$fi_00,
    fi_0Nj    = gr$fi_0Nj,
    fi_cNj    = gr$fi_cNj,
    N         = nf$N,
    F_val     = nf$F,
    lhs       = cx$lhs,
    rhs       = cx$rhs,
    persists  = cx$persists
  )
})

print(results_NF)

## Step 7: Visualize N and F in coexistence space ####
## rescale to Spaak et al. 2021 N-F
results_NF <- results_NF %>%
  mutate(
    F_2021 = -F_val / (1 - (fi_cNj / fi_00))
  )

NF_plot_strongmand <- ggplot(results_NF, aes(x = N, y = F_2021, 
                                             shape = persists, 
                                             label = herb,
                                             color=as.factor(herb))) +
  geom_point(size = 4, stroke=2) +
  geom_text(nudge_y = 0.005, size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  # Coexistence boundary: -F = N/(1-N)  =>  F = -N/(1-N)
  stat_function(fun = function(x) x / (1 - x),
                color = "black", linetype = "dashed", linewidth=2, xlim = c(-0.5, 0.95)) +
  scale_color_viridis(option="inferno", discrete=T, direction=-1, begin=.2, end=.8) +
  #scale_fill_viridis(option="inferno", discrete=T, direction=-1, begin=.2, end=.8) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21),
                     labels = c("Excluded", "Persists")) +
  coord_cartesian(ylim = c(-.01, .062), xlim=c(-.55, .1)) +
  labs(x = "Niche difference (N)",
       y = "Fitness difference (F)",
       shape = "Outcome",
       title = "Strong mandibles",
       subtitle = "N-F coexistence space") +
  theme_bw(base_size = 14) +
  theme(legend.position = c(.2,.2))+
  guides(color="none")

NF_plot <- NF_plot_traitmatch+NF_plot_strongmand #+ plot_layout(guides="collect")
#ggsave(NF_plot, file="NF_plot.png", width=12, height=5, dpi=300)



# Matching vs strong mandibles plot ####
MatchStrongPlot <- popplot|popplota
ggsave(MatchStrongPlot, file="Fig_1_MatchStrongPlot.tiff", width=10, height=6, dpi=300, compression = "lzw")

MCR_plot <- (popplot|popplota)/NF_plot +
  plot_layout(heights = c(2,1))
#ggsave(MCR_plot, file="MCR_plot.tiff", width=10, height=10, dpi=300, compression = "lzw")

## Heatmaps of attack rate matrices ####
# ── Trait matching attack matrix heatmap ─────────────────────────────────────
A_match_long <- as_tibble(A_match) %>%
  setNames(paste0("H", 1:n_herbs)) %>%
  mutate(plant     = paste0("P", 1:n_plants),
         toughness = toughness) %>%
  pivot_longer(starts_with("H"),
               names_to  = "herbivore",
               values_to = "attack_rate") %>%
  mutate(herbivore = factor(herbivore, levels = paste0("H", 1:n_herbs)),
         plant     = factor(plant,     levels = paste0("P", 1:n_plants)))

p_heatmap_match <- ggplot(A_match_long,
                          aes(x = herbivore, y = plant, fill = attack_rate)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "plasma", limits=c(0,.09)) +
  labs(x     = "Herbivore (mandible strength →)",
       y     = "Plant (LDMC →)",
       title = "A) Trait matching attack matrix") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")

# ── Strong mandibles attack matrix heatmap ───────────────────────────────────
A_strong_long <- as_tibble(A_strong) %>%
  setNames(paste0("H", 1:n_herbs)) %>%
  mutate(plant     = paste0("P", 1:n_plants),
         toughness = toughness) %>%
  pivot_longer(starts_with("H"),
               names_to  = "herbivore",
               values_to = "attack_rate") %>%
  mutate(herbivore = factor(herbivore, levels = paste0("H", 1:n_herbs)),
         plant     = factor(plant,     levels = paste0("P", 1:n_plants)))

p_heatmap_strong <- ggplot(A_strong_long,
                           aes(x = herbivore, y = plant, fill = attack_rate)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "plasma", name = "Attack\nrate", limits=c(0,.09)) +
  labs(x     = "Herbivore (mandible strength →)",
       y     = "Plant (LDMC →)",
       title = "B) Strong mandibles attack matrix") +
  theme_bw(base_size = 16) 

# ── Combined ──────────────────────────────────────────────────────────────────
p_heatmaps <- p_heatmap_match | p_heatmap_strong
p_heatmaps

ggsave("Fig_A1_attack_matrices.png", p_heatmaps,
       width = 12, height = 5, dpi = 300)



# Trait matching: w/ handling time ####
# Species counts
n_plants <- 10
n_herbs <- 10

# Trait values
toughness <- seq(30, 50, length.out=n_plants)
mandible <- seq(2.5, 4.5, length.out=n_herbs)
#toughness <- rnorm(n_plants, 40, 5)
#mandible <- rnorm(n_herbs, 3.5, 1)

# Normalize traits
mean_toughness <- mean(toughness)
mean_mandible <- mean(mandible)


# Trait-based attack rates using Gaussian trait matching
beta <- 10
a0 <- 0.2
T <- .05 * n_plants  # ensures total attack rate matches a generalist's base rate

A_raw <- matrix(0, nrow = n_plants, ncol = n_herbs)
for (i in 1:n_plants) {
  for (j in 1:n_herbs) {
    ratio <- (mandible[j] / mean_mandible) / (toughness[i] / mean_toughness)
    similarity <- exp(-beta * (ratio - 1)^2)
    A_raw[i, j] <- similarity
  }
}

# Normalize so each herbivore's total attack strength across plants = T
A <- sweep(A_raw, 2, colSums(A_raw), FUN = "/") * T

colSums(A)
rowSums(A)

colSums(A_raw)
rowSums(A_raw)

# Fixed parameters
R <- rep(2, n_plants)
C <- matrix(0.0, nrow = n_plants, ncol = n_plants); diag(C) <- 0.04
HNDL <- matrix(runif(n_plants * n_herbs, .02, 0.02), nrow = n_plants)
E <- matrix(runif(n_plants * n_herbs, 0.1, 0.1), nrow = n_plants)
M <- rep(1, n_herbs)
g <- rep(0.00, n_herbs)  # herbivore self-limitation rate

# Initial state
state <- c(P = rep(5, n_plants), H = rep(2, n_herbs))

times <- seq(0, 100, by = .1)

# Model function
competition_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    P <- state[1:n_plants]
    H <- state[(n_plants + 1):(n_plants + n_herbs)]
    dP <- numeric(n_plants)
    dH <- numeric(n_herbs)
    
    for (i in 1:n_plants) {
      comp_sum <- sum(C[i, ] * P)
      herb_sum <- 0
      for (j in 1:n_herbs) {
        denom <- 1 + sum(HNDL[i, j] * A[i, j] * P)
        herb_sum <- herb_sum + (A[i, j] * H[j]) / denom
      }
      dP[i] <- P[i] * ((R[i] - comp_sum) - herb_sum)
    }
    
    for (j in 1:n_herbs) {
      numer <- sum(E[, j] * A[, j] * P)
      denom <- 1 + sum(A[, j] * HNDL[, j] * P)
      dH[j] <- H[j] * (numer / denom - M[j] - g[j] * H[j])
    }
    
    list(c(dP, dH))
  })
}

# Solve ODE
out <- ode(y = state, times = times, func = competition_model, parms = list())

# Format for ggplot
out_df <- as.data.frame(out)
colnames(out_df) <- c("time", paste0("P", 1:n_plants), paste0("H", 1:n_herbs))
out_long <- out_df %>%
  pivot_longer(-time, names_to = "species", values_to = "value") %>%
  mutate(type = ifelse(str_detect(species, "^P"), "Plant", "Herbivore"),
         species=as_factor(species),
         spp=factor(species, levels=mixedsort(levels(species))))

## make plot ####
pop1 <- ggplot(data=out_long %>% filter(type=="Plant"), aes(x = time, y = value, color = spp)) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis(option="viridis", discrete=T, direction=-1, begin=.2, end=.8) +
  labs(#title = "Plant density (g/m2)",
    x = "Time",
    y = "Plant biomass (g/m2)",
    color="Plant spp") +
  ylim(0,50)+
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

pop2 <- ggplot(data=out_long %>% filter(type=="Herbivore"), 
               aes(x = time, y = value, color = spp)) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis(option="inferno", discrete=T, direction=-1, begin=.2, end=.8) +
  labs(title = "A) Trait matching",
       x = "Time",
       y = "Herbivore density (#/m2)",
       color="Herbivore spp") +
  ylim(0,10)+
  theme_bw(base_size = 16)

popplot <- pop2/pop1

# Assume this is your output from the model at time = 100
# Filter for t = 100
abund_at_100 <- out_long %>%
  filter(time == 100) %>%
  mutate(
    trait = case_when(
      type == "Plant" ~ toughness[as.numeric(str_extract(species, "\\d+"))],
      type == "Herbivore" ~ mandible[as.numeric(str_extract(species, "\\d+"))]
    )
  )

# Plot
pop1.1 <- ggplot(abund_at_100 %>% filter(type=="Plant"), 
                 aes(x = trait, y = value, color=trait)) +
  geom_point(size = 3) +
  #geom_line(aes(group = type), linetype = "dashed", alpha = 0.6) +
  scale_color_viridis(option="viridis", discrete=F, direction=-1, begin=.2, end=.8) +
  #facet_wrap(~type, scales = "free_x") +
  labs(
    #title = "Trait Value vs. Abundance at Time = 100",
    x = "LDMC (%)",
    y = "Plant biomass (g/m2)",
    color = "LDMC (%)"
  ) +
  theme_bw(base_size = 14)

pop1.2 <- ggplot(abund_at_100 %>% filter(type=="Herbivore"), 
                 aes(x = trait, y = value, color=trait)) +
  geom_point(size = 3) +
  #geom_line(aes(group = type), linetype = "dashed", alpha = 0.6) +
  scale_color_viridis(option="inferno", discrete=F, direction=-1, begin=.2, end=.8) +
  #facet_wrap(~type, scales = "free_x") +
  labs(
    title = "Trait Value vs. Abund at T=100",
    x = "Mandible strength",
    y = "Herbivore abundance (#/m2)",
    color = "Mandible"
  ) +
  theme_bw(base_size = 14)

popplot12 <- pop1.2/pop1.1

Arate_plot <- popplot|popplot12
#ggsave(Arate_plot, file="Arate_plot.tiff", width=8, height=6, dpi=300, compression = "lzw")

# Strong mandible model with Handling Time ####
# Species counts
n_plants <- 20
n_herbs <- 10

# Trait values
toughness <- seq(30, 50, length.out=n_plants)
mandible <- seq(2.5, 4.5, length.out=n_herbs)

# Normalize traits
mean_toughness <- mean(toughness)
mean_mandible <- mean(mandible)

# Base attack rate
a0 <- 0.2
T <- .005 * n_plants  # ensures total attack rate matches a generalist's base rate

# Trait-based attack rates using Gaussian trait matching
A <- outer(toughness, mandible, function(t, m) (m / mean(m)) / (t / mean(t))) * a0 * T

colSums(A)
rowSums(A)


# Fixed parameters
R <- rep(2, n_plants)
C <- matrix(0.0, nrow = n_plants, ncol = n_plants); diag(C) <- 0.04
HNDL <- matrix(runif(n_plants * n_herbs, 0.02, 0.02), nrow = n_plants)
E <- matrix(runif(n_plants * n_herbs, 0.1, 0.1), nrow = n_plants)
M <- rep(1, n_herbs)
g <- rep(0.00, n_herbs)  # herbivore self-limitation rate

# Initial state
state <- c(P = rep(15, n_plants), H = rep(5, n_herbs))

times <- seq(0, 100, by = .1)

# Model function
competition_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    P <- state[1:n_plants]
    H <- state[(n_plants + 1):(n_plants + n_herbs)]
    dP <- numeric(n_plants)
    dH <- numeric(n_herbs)
    
    for (i in 1:n_plants) {
      comp_sum <- sum(C[i, ] * P)
      herb_sum <- 0
      for (j in 1:n_herbs) {
        denom <- 1 + sum(HNDL[i, j] * A[i, j] * P)
        herb_sum <- herb_sum + (A[i, j] * H[j]) / denom
      }
      dP[i] <- P[i] * ((R[i] - comp_sum) - herb_sum)
    }
    
    for (j in 1:n_herbs) {
      numer <- sum(E[, j] * A[, j] * P)
      denom <- 1 + sum(A[, j] * HNDL[, j] * P)
      dH[j] <- H[j] * (numer / denom - M[j] - g[j] * H[j])
    }
    
    list(c(dP, dH))
  })
}

# Solve ODE
out <- ode(y = state, times = times, func = competition_model, parms = list())

# Format for ggplot
out_df <- as.data.frame(out)
colnames(out_df) <- c("time", paste0("P", 1:n_plants), paste0("H", 1:n_herbs))
out_long <- out_df %>%
  pivot_longer(-time, names_to = "species", values_to = "value") %>%
  mutate(type = ifelse(str_detect(species, "^P"), "Plant", "Herbivore"),
         species=as_factor(species),
         spp=factor(species, levels=mixedsort(levels(species))))

## make plot ####
pop1b <- ggplot(data=out_long %>% filter(type=="Plant"), aes(x = time, y = value, color = spp)) +
  geom_line(linewidth = 1) +
  scale_color_viridis(option="viridis", discrete=T, direction=-1, begin=.2, end=.8) +
  labs(#title = "Plant density (g/m2)",
    x = "Time",
    y = "Plant biomass (g/m2)",
    color="Plant spp") +
  ylim(0,35)+
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

pop2b <- ggplot(data=out_long %>% filter(type=="Herbivore"), 
                aes(x = time, y = value, color = spp)) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis(option="inferno", discrete=T, direction=-1, begin=.2, end=.8) +
  labs(title = "Strong mandibles",
       x = "Time",
       y = "Herbivore abundance (#/m2)",
       color="Herbivore spp") +
  ylim(0,40)+
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

popplotb <- pop2b/pop1b

# Assume this is your output from the model at time = 100
# Filter for t = 100
abund_at_100 <- out_long %>%
  filter(time == 100) %>%
  mutate(
    trait = case_when(
      type == "Plant" ~ toughness[as.numeric(str_extract(species, "\\d+"))],
      type == "Herbivore" ~ mandible[as.numeric(str_extract(species, "\\d+"))]
    )
  )

# Plot
pop1.1b <- ggplot(abund_at_100 %>% filter(type=="Plant"), 
                  aes(x = trait, y = value, color=spp)) +
  geom_point(size = 3) +
  #geom_line(aes(group = type), linetype = "dashed", alpha = 0.6) +
  scale_color_viridis(option="viridis", discrete=T, direction=-1, begin=.2, end=.8) +
  #facet_wrap(~type, scales = "free_x") +
  labs(
    #title = "Trait Value vs. Abundance at Time = 100",
    x = "LDMC (%)",
    y = "Plant biomass (g/m2)",
    color = "LDMC (%)"
  ) +
  theme_bw(base_size = 14)

pop1.2b <- ggplot(abund_at_100 %>% filter(type=="Herbivore"), 
                  aes(x = trait, y = value, color=spp)) +
  geom_point(size = 3) +
  #geom_line(aes(group = type), linetype = "dashed", alpha = 0.6) +
  scale_color_viridis(option="inferno", discrete=T, direction=-1, begin=.2, end=.8) +
  #facet_wrap(~type, scales = "free_x") +
  labs(
    title = "Trait Value vs. Abund at T=100",
    x = "Mandible strength",
    y = "Herbivore abundance (#/m2)",
    color = "Mandible"
  ) +
  theme_bw(base_size = 14)

popplot12b <- pop1.2b/pop1.1b

StrongH_plot <- popplota|popplot12a
ggsave(StrongH_plot, file="StrongH_plot.tiff", width=8, height=6, dpi=300, compression = "lzw")




# PLANT ONLY - sweeps ####
# Matching - sweep across R1 values ####
# Base setup
n_plants <- 10
n_herbs <- 10
R_base <- rep(2, n_plants)

# Trait values
toughness <- seq(30, 50, length.out=n_plants)
mandible <- seq(2.5, 4.5, length.out=n_herbs)
mean_toughness <- mean(toughness)
mean_mandible <- mean(mandible)

# Attack rate matrix from trait matching
beta <- 10
a0 <- 0.2  # base attack rate
T <- .05 * n_plants  # ensures total attack rate matches a generalist's base rate

A_raw <- matrix(0, nrow = n_plants, ncol = n_herbs)
for (i in 1:n_plants) {
  for (j in 1:n_herbs) {
    ratio <- (mandible[j] / mean_mandible) / (toughness[i] / mean_toughness)
    A_raw[i, j] <- exp(-beta * (ratio - 1)^2)
  }
}
A <- sweep(A_raw, 2, colSums(A_raw), FUN = "/") * T

colSums(A)
rowSums(A)

# Fixed parameters
C <- matrix(0.0, nrow = n_plants, ncol = n_plants); diag(C) <- 0.04
HNDL <- matrix(0, nrow = n_plants, ncol = n_herbs)
E <- matrix(0.1, nrow = n_plants, ncol = n_herbs)
M <- rep(1, n_herbs)
g <- rep(0.00, n_herbs)

# Lower initial state
state <- c(P = rep(15, n_plants), H = rep(0, n_herbs))
times <- seq(0, 100, by = 1)

# ODE model function
competition_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    P <- state[1:n_plants]
    H <- state[(n_plants + 1):(n_plants + n_herbs)]
    dP <- numeric(n_plants)
    dH <- numeric(n_herbs)
    
    for (i in 1:n_plants) {
      comp_sum <- sum(C[i, ] * P)
      herb_sum <- 0
      for (j in 1:n_herbs) {
        denom <- 1 + sum(HNDL[i, j] * A[i, j] * P)
        herb_sum <- herb_sum + (A[i, j] * H[j]) / denom
      }
      dP[i] <- P[i] * ((R[i] - comp_sum) - herb_sum)
    }
    
    for (j in 1:n_herbs) {
      numer <- sum(E[, j] * A[, j] * P)
      denom <- 1 + sum(A[, j] * HNDL[, j] * P)
      dH[j] <- H[j] * (numer / denom - M[j] - g[j] * H[j])
    }
    
    list(c(dP, dH))
  })
}

# Sweep R1 from 1 to 5
R1_values <- seq(1, 3, by = .1)
tough_std <- scale(toughness)[,1]
gradient_values <- seq(-1, 1, length.out = 20)

results <- map_dfr(gradient_values, function(g) {
  #R <- R_base
  #R[7:10] <- R1_val
  R <- 2 + g * tough_std
  
  # Run ODE simulation
  out <- ode(y = state, times = times, func = competition_model, parms = list(R = R))
  
  # Convert to data frame and assign species names explicitly
  out <- as.data.frame(out)
  colnames(out) <- c("time", paste0("P", 1:n_plants), paste0("H", 1:n_herbs))
  
  # Extract the last row (time = 500), keep it as a data frame with proper column names
  final <- out[nrow(out), -1]  # Drop time column
  
  # Add R1 value and return as tibble row
  final_df <- tibble(R1 = g, !!!final)  # unpack final into tibble columns
  return(final_df)
})

# Long format and plot
results_long <- results %>%
  pivot_longer(cols = -R1, names_to = "species", values_to = "abundance") %>%
  mutate(type = ifelse(str_detect(species, "^P"), "Plant", "Herbivore"),
         species=as_factor(species),
         species=factor(species, levels=mixedsort(levels(species))))

## make plot ####
comm1p <- ggplot(data=results_long %>% filter(type=="Plant"), 
                aes(x = R1+2, y = abundance, color = species)) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis(option="viridis", discrete=T, direction=-1, begin=.2, end=.8) +
  labs(#title = "Plant abundances across R1",
       x = "Intrinsic Growth Rate of P1 (R1)",
       y = "Abundance at t = 100",
       color="Plant spp") +
  theme_bw(base_size = 16)

comm2p <- ggplot(data=results_long %>% filter(type=="Herbivore"), 
                aes(x = R1, y = abundance, color = species)) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis(option="inferno", discrete=T, direction=-1, begin=.2, end=.8) +
  labs(title = "Herbivore abundances across R1",
       x = "Intrinsic Growth Rate of P1 (R1)",
       y = "Abundance at t = 100",
       color="Herbivore spp") +
  theme_bw(base_size = 16)

comm1p/comm2p


## calculate CWM herbivore traits ####
resultsfd <- map_dfr(gradient_values, function(g) {
  #R <- R_base
  #R[7:10] <- R1_val
  R <- 2 + g * tough_std
  
  out <- ode(y = state, times = times, func = competition_model, parms = list(R = R))
  
  out <- as.data.frame(out)
  colnames(out) <- c("time", paste0("P", 1:n_plants), paste0("H", 1:n_herbs))
  
  final <- out[nrow(out), -1]  # drop time column
  herbivores <- final[paste0("H", 1:n_herbs)]
  
  total_H <- sum(herbivores)
  if (total_H > 0) {
    cwm_mandible <- sum(herbivores * mandible) / total_H
    fdis_mandible <- sum(herbivores * abs(mandible - cwm_mandible)) / total_H
  } else {
    cwm_mandible <- NA
    fdis <- NA
  }
  
  plants <- final[paste0("P", 1:n_plants)]
  
  total_P <- sum(plants)
  if (total_P > 0) {
    cwm_toughness <- sum(plants * toughness) / total_P
    fdis_toughness <- sum(plants * abs(toughness - cwm_toughness)) / total_P
  } else {
    cwm_toughness <- NA
    fdis_toughness <- NA
  }
  
  tibble(
    R1 = g,
    #cwm_mandible = cwm_mandible,
    #fdis_mandible = fdis_mandible,
    cwm_ldmc = cwm_toughness,
    fdis_toughness = fdis_toughness,
    !!!final
  )
})




### Matching - SES ####


#### Rao's randomization -- toughness trait ####
mandibleT <- (c('P1','P2','P3','P4','P5','P6','P7','P8','P9','P10'))
t8 <- cbind(as_tibble(toughness),mandibleT) %>% column_to_rownames("mandibleT")

comm3 <- resultsfd %>% dplyr::select(P1:P10) %>% as.matrix()

frao_BV <- hill_func(comm=comm3, traits=t8, q=1)[1,]
frao_BVr <- data.frame(data.frame(matrix(ncol = 10, nrow = 20)))

for (i in 1:999) {#For each row in the matrix (for each site)
  #select randomly, as many species as the species richness of the site:
  comm3rand <- comm3
  #colnames(comm3rand) <- sample(colnames(comm3), replace=F) 
  #comm3rand <- comm3rand %>% select(order(colnames(.)))
  
  t8rand <- t8 %>% as.data.frame()
  t8rownames <- rownames(t8)
  t8rand <- t8rand[sample(nrow(t8rand)),]%>% as.data.frame()
  rownames(t8rand) <- t8rownames
  colnames(t8rand) <- "value"
  
  frao_BVr[,i] <- hill_func(comm3rand, t8rand, q=1)[1,]
}

BV_SEStab <- frao_BVr %>% 
  cbind(.,frao_BV) %>%
  rowwise() %>% 
  mutate(rand_mean = mean(c_across(2:1000)),rand_sd = sd(c_across(2:1000))) %>% 
  select(frao_BV,rand_mean,rand_sd) %>% 
  mutate(SES_BV=(frao_BV-rand_mean)/rand_sd, pval=pnorm(q=SES_BV, lower.tail = T)) 

summary(BV_SEStab)  
hist(BV_SEStab$SES_BV)
resultsfd$SES_toughness <- BV_SEStab$SES_BV

### FD SES Plots ####
pcwm1 <- ggplot(resultsfd, aes(x = R1+2, y = cwm_ldmc)) +
  geom_line(size = 1.5, color = "darkgreen") +
  labs(
    #title = "Herbivore community metrics across R1",
    x = "Instrinsic Growth Rate of P1 (R1)",
    y = "CWM Leaf Dry Matter Cont.") +
  #ylim(2.75,4.75)+
  theme_bw(base_size = 16)

pfd1 <- ggplot(resultsfd, aes(x = R1, y = fdis_toughness)) +
  geom_line(size = 1.2, color = "#cf4446") +
  labs(
    #title = "Functional Dispersion of Mandible Strength Across R1",
    x = "Growth Rate of P7:P10\n(Plant Community LDMC)",
    y = "FD Mandible Strength"
  ) +
  theme_bw(base_size = 16)

pses1 <- ggplot(resultsfd, aes(x = R1, y = SES_toughness)) +
  geom_hline(aes(yintercept=0), color="grey", linetype="dashed", lwd=1.5)+
  geom_line(size = 1.2, color = "#cf4446") +
  labs(
    #title = "Functional Dispersion of Mandible Strength Across R1",
    x = "Growth Rate of P7:P10\n(Plant Community LDMC)",
    y = "SES Mandible Strength"
  ) +
  ylim(-3,3)+
  theme_bw(base_size = 16)

pcommplot2 <- (comm1p|pcwm1) + plot_annotation(tag_levels = "A")

R1_cwm <- resultsfd[1:2]

ggsave("Fig_B1_LDMC_gradient.png", pcommplot2, width=12, height = 5, dpi=300)


#SWEEPS START HERE ####
# Matching - sweep across R1 values ####
# Base setup
n_plants <- 10
n_herbs <- 10
R_base <- rep(2, n_plants)

# Trait values
toughness <- seq(30, 50, length.out=n_plants)
mandible <- seq(2.5, 4.5, length.out=n_herbs)
mean_toughness <- mean(toughness)
mean_mandible <- mean(mandible)

# Attack rate matrix from trait matching
beta <- 10
a0 <- 0.2  # base attack rate
T <- .05 * n_plants  # ensures total attack rate matches a generalist's base rate

A_raw <- matrix(0, nrow = n_plants, ncol = n_herbs)
for (i in 1:n_plants) {
  for (j in 1:n_herbs) {
    ratio <- (mandible[j] / mean_mandible) / (toughness[i] / mean_toughness)
    A_raw[i, j] <- exp(-beta * (ratio - 1)^2)
  }
}
A <- sweep(A_raw, 2, colSums(A_raw), FUN = "/") * T

colSums(A)
rowSums(A)

# Lower initial state
state <- c(P = rep(15, n_plants), H = rep(5, n_herbs))
times <- seq(0, 100, by = 1)

# Sweep growth rates of plants 
results <- map_dfr(gradient_values, function(g) {
  #R <- R_base
  #R[7:10] <- R1_val
  R <- 2 + g * tough_std
  # Run ODE simulation
  out <- ode(y = state, times = times, func = competition_model, parms = list(R = R))
  
  # Convert to data frame and assign species names explicitly
  out <- as.data.frame(out)
  colnames(out) <- c("time", paste0("P", 1:n_plants), paste0("H", 1:n_herbs))
  out <- ode(y = state, times = times, func = competition_model, parms = list(R = R))
  
  # Convert to data frame and assign species names explicitly
  out <- as.data.frame(out)
  colnames(out) <- c("time", paste0("P", 1:n_plants), paste0("H", 1:n_herbs))
  
  # Extract the last row (time = 500), keep it as a data frame with proper column names
  final <- out[nrow(out), -1]  # Drop time column
  
  # Add R1 value and return as tibble row
  final_df <- tibble(R1 = g, !!!final)  # unpack final into tibble columns
  return(final_df)
})

# Long format and plot
results_long <- results %>%
  pivot_longer(cols = -R1, names_to = "species", values_to = "abundance") %>%
  mutate(type = ifelse(str_detect(species, "^P"), "Plant", "Herbivore"),
         species=as_factor(species),
         species=factor(species, levels=mixedsort(levels(species)))) %>% 
  full_join(R1_cwm)

## make plot ####
comm1 <- ggplot(data=results_long %>% filter(type=="Plant"), aes(x = R1+2, y = abundance, color = species)) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis(option="viridis", discrete=T, direction=-1, begin=.2, end=.8) +
  labs(#title = "A",
       x = "Intrinsic Growth Rate of P1 (R1)",
       y = "Abundance at t = 100",
       color="Plant spp") +
  theme_bw(base_size = 16)+
  theme(legend.position = "top")+
  annotate("text", x = 1.1, y = 26, label = "A")+
  guides(color=guide_legend(title.position = "top", title.hjust=0.5))

comm2 <- ggplot(data=results_long %>% filter(type=="Herbivore"), 
                aes(x =cwm_ldmc, y = abundance, color = species)) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis(option="inferno", discrete=T, direction=-1, begin=.2, end=.8) +
  labs(#subtitle = "B",
       x = "CWM LDMC w/ no herbivores",
       y = "Abundance at t = 100",
       color="Herbivore spp") +
  theme_bw(base_size = 16)+
  theme(legend.position = "top") +
  annotate("text", x = 37.5, y = 22, label = "B")+
  guides(color=guide_legend(title.position = "top", title.hjust=0.5))

comm12 <- comm1/comm2


## calculate CWM herbivore traits ####
resultsfd <- map_dfr(gradient_values, function(g) {
  #R <- R_base
  #R[8:10] <- R1_val
  R <- 2 + g * tough_std
  out <- ode(y = state, times = times, func = competition_model, parms = list(R = R))
  
  out <- as.data.frame(out)
  colnames(out) <- c("time", paste0("P", 1:n_plants), paste0("H", 1:n_herbs))
  
  final <- out[nrow(out), -1]  # drop time column
  herbivores <- final[paste0("H", 1:n_herbs)]
  
  total_H <- sum(herbivores)
  if (total_H > 0) {
    cwm_mandible <- sum(herbivores * mandible) / total_H
    fdis_mandible <- sum(herbivores * abs(mandible - cwm_mandible)) / total_H
  } else {
    cwm_mandible <- NA
    fdis <- NA
  }
  
  plants <- final[paste0("P", 1:n_plants)]
  
  total_P <- sum(plants)
  if (total_P > 0) {
    cwm_toughness <- sum(plants * toughness) / total_P
    fdis_toughness <- sum(plants * abs(toughness - cwm_toughness)) / total_P
  } else {
    cwm_toughness <- NA
    fdis_toughness <- NA
  }
  
  tibble(
    R1 = g,
    cwm_mandible = cwm_mandible,
    fdis_mandible = fdis_mandible,
    cwm_toughness = cwm_toughness,
    fdis_toughness = fdis_toughness,
    !!!final
  )
})




### Matching - SES ####
#### Rao's randomization -- mandible trait ####
mandibleT <- (c('H1','H2','H3','H4','H5','H6','H7','H8','H9','H10'))
t8 <- cbind(as_tibble(mandible),mandibleT) %>% column_to_rownames("mandibleT")

comm3 <- resultsfd %>% dplyr::select(H1:H10) %>% as.matrix()

frao_BV <- hill_func(comm=comm3, traits=t8, q=1)[1,]
frao_BVr <- data.frame(data.frame(matrix(ncol = 10, nrow = 20)))

for (i in 1:999) {#For each row in the matrix (for each site)
  #select randomly, as many species as the species richness of the site:
  comm3rand <- comm3
  #colnames(comm3rand) <- sample(colnames(comm3), replace=F) 
  #comm3rand <- comm3rand %>% select(order(colnames(.)))
  
  t8rand <- t8 %>% as.data.frame()
  t8rownames <- rownames(t8)
  t8rand <- t8rand[sample(nrow(t8rand)),]%>% as.data.frame()
  rownames(t8rand) <- t8rownames
  colnames(t8rand) <- "value"
  
  frao_BVr[,i] <- hill_func(comm3rand, t8rand, q=1)[1,]
}

BV_SEStab <- frao_BVr %>% 
  cbind(.,frao_BV) %>%
  rowwise() %>% 
  mutate(rand_mean = mean(c_across(2:1000)),rand_sd = sd(c_across(2:1000))) %>% 
  select(frao_BV,rand_mean,rand_sd) %>% 
  mutate(SES_BV=(frao_BV-rand_mean)/rand_sd, pval=pnorm(q=SES_BV, lower.tail = T)) 

summary(BV_SEStab)  
hist(BV_SEStab$SES_BV)
resultsfd$SES_Mandible <- BV_SEStab$SES_BV

#### Rao's randomization -- toughness trait ####
mandibleT <- (c('P1','P2','P3','P4','P5','P6','P7','P8','P9','P10'))
t8 <- cbind(as_tibble(toughness),mandibleT) %>% column_to_rownames("mandibleT")

comm3 <- resultsfd %>% dplyr::select(P1:P10) %>% as.matrix()

frao_BV <- hill_func(comm=comm3, traits=t8, q=1)[1,]
frao_BVr <- data.frame(data.frame(matrix(ncol = 10, nrow = 20)))

for (i in 1:999) {#For each row in the matrix (for each site)
  #select randomly, as many species as the species richness of the site:
  comm3rand <- comm3
  #colnames(comm3rand) <- sample(colnames(comm3), replace=F) 
  #comm3rand <- comm3rand %>% select(order(colnames(.)))
  
  t8rand <- t8 %>% as.data.frame()
  t8rownames <- rownames(t8)
  t8rand <- t8rand[sample(nrow(t8rand)),]%>% as.data.frame()
  rownames(t8rand) <- t8rownames
  colnames(t8rand) <- "value"
  
  frao_BVr[,i] <- hill_func(comm3rand, t8rand, q=1)[1,]
}

BV_SEStab <- frao_BVr %>% 
  cbind(.,frao_BV) %>%
  rowwise() %>% 
  mutate(rand_mean = mean(c_across(2:1000)),rand_sd = sd(c_across(2:1000))) %>% 
  select(frao_BV,rand_mean,rand_sd) %>% 
  mutate(SES_BV=(frao_BV-rand_mean)/rand_sd, pval=pnorm(q=SES_BV, lower.tail = T)) 

summary(BV_SEStab)  
hist(BV_SEStab$SES_BV)
resultsfd$SES_toughness <- BV_SEStab$SES_BV

### FD SES Plots ####
resultsfd <- full_join(resultsfd,R1_cwm)

cwm1 <- ggplot(resultsfd, aes(x = cwm_ldmc, y = cwm_mandible)) +
  geom_line(size = 1.2, color = "#cf4446") +
  labs(
    #title = "Herbivore community metrics across R1",
    x = "CWM LDCM w/ no herbivores",
    y = "CWM Mandible Strength") +
  ylim(2.4,4.75)+
  annotate("text", x = 37, y = Inf, label = "C", vjust = 1.5)+
  theme_bw(base_size = 16)

fd1 <- ggplot(resultsfd, aes(x = cwm_ldmc, y = fdis_mandible)) +
  geom_line(size = 1.2, color = "#cf4446") +
  labs(
    #title = "Functional Dispersion of Mandible Strength Across R1",
    x = "CWM LDCM w/ no herbivores",
    y = "FD Mandible Strength"
  ) +
  annotate("text", x = 37, y = Inf, label = "D", vjust = 1.5)+
  theme_bw(base_size = 16)

ses1 <- ggplot(resultsfd, aes(x = cwm_ldmc, y = SES_Mandible)) +
  geom_hline(aes(yintercept=0), color="grey", linetype="dashed", lwd=1.5)+
  geom_line(size = 1.2, color = "#cf4446") +
  labs(
    #title = "Functional Dispersion of Mandible Strength Across R1",
    x = "CWM LDCM w/ no herbivores",
    y = "SES Mandible Strength"
  ) +
  ylim(-3,3)+
  annotate("text", x = 37, y = Inf, label = "E", vjust = 1.5)+
  theme_bw(base_size = 16)

commplot <- free(comm1/comm2)|(cwm1/fd1/ses1) + 
  plot_layout(guides="collect", axes="collect")

ggsave(commplot, file="Fig_2_TraitMatchingPlot.tiff", width=10, height=8, dpi=300, compression = "lzw")



# STRONG - sweep across R1 values ####
# Species counts
n_plants <- 10
n_herbs <- 10

# Trait values
toughness <- seq(30, 50, length.out=n_plants)
mandible <- seq(2.5, 4.5, length.out=n_herbs)

# Normalize traits
mean_toughness <- mean(toughness)
mean_mandible <- mean(mandible)

# Base attack rate
a0 <- 0.1
T <- .1 * n_plants  # ensures total attack rate matches a generalist's base rate

# Trait-based attack rates using Gaussian trait matching
A <- outer(toughness, mandible, function(t, m) (m / mean(m)) / (t / mean(t))) * a0 * T

colSums(A)
rowSums(A)


# Fixed parameters
R <- rep(2, n_plants)
C <- matrix(0.0, nrow = n_plants, ncol = n_plants); diag(C) <- 0.04
HNDL <- matrix(runif(n_plants * n_herbs, 0, 0), nrow = n_plants)
E <- matrix(runif(n_plants * n_herbs, 0.1, 0.1), nrow = n_plants)
M <- rep(1, n_herbs)
g <- rep(0.01, n_herbs)  # herbivore self-limitation rate

# Initial state
state <- c(P = rep(15, n_plants), H = rep(3, n_herbs))

times <- seq(0, 100, by = .1)

# ODE model function
competition_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    P <- state[1:n_plants]
    H <- state[(n_plants + 1):(n_plants + n_herbs)]
    dP <- numeric(n_plants)
    dH <- numeric(n_herbs)
    
    for (i in 1:n_plants) {
      comp_sum <- sum(C[i, ] * P)
      herb_sum <- 0
      for (j in 1:n_herbs) {
        denom <- 1 + sum(HNDL[i, j] * A[i, j] * P)
        herb_sum <- herb_sum + (A[i, j] * H[j]) / denom
      }
      dP[i] <- P[i] * ((R[i] - comp_sum) - herb_sum)
    }
    
    for (j in 1:n_herbs) {
      numer <- sum(E[, j] * A[, j] * P)
      denom <- 1 + sum(A[, j] * HNDL[, j] * P)
      dH[j] <- H[j] * (numer / denom - M[j] - g[j] * H[j])
    }
    
    list(c(dP, dH))
  })
}

# Sweep R1 from 1 to 5
R1_values <- seq(1, 3, by = .1)

results <- map_dfr(R1_values, function(R1_val) {
  R <- R_base
  R[7:10] <- R1_val
  
  # Run ODE simulation
  out <- ode(y = state, times = times, func = competition_model, parms = list(R = R))
  
  # Convert to data frame and assign species names explicitly
  out <- as.data.frame(out)
  colnames(out) <- c("time", paste0("P", 1:n_plants), paste0("H", 1:n_herbs))
  
  # Extract the last row (time = 500), keep it as a data frame with proper column names
  final <- out[nrow(out), -1]  # Drop time column
  
  # Add R1 value and return as tibble row
  final_df <- tibble(R1 = R1_val, !!!final)  # unpack final into tibble columns
  return(final_df)
})

# Long format and plot
results_long <- results %>%
  pivot_longer(cols = -R1, names_to = "species", values_to = "abundance") %>%
  mutate(type = ifelse(str_detect(species, "^P"), "Plant", "Herbivore"),
         species=as_factor(species),
         species=factor(species, levels=mixedsort(levels(species))))

## make plot ####
comm1a <- ggplot(data=results_long %>% filter(type=="Plant"), aes(x = R1, y = abundance, color = species)) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis(option="viridis", discrete=T, direction=-1, begin=.2, end=.8) +
  labs(title = "Plant abundances across R1",
       x = "Intrinsic Growth Rate of P1 (R1)",
       y = "Abundance at t = 100",
       color="Plant spp") +
  theme_bw(base_size = 14)

comm2a <- ggplot(data=results_long %>% filter(type=="Herbivore"), aes(x = R1, y = abundance, color = species)) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis(option="inferno", discrete=T, direction=-1, begin=.2, end=.8) +
  #scale_y_continuous(trans='log10', limits = c(0.00000000000000000001,20))+
  labs(title = "B) Strong mandibles",
       x = "Growth Rate of P7:P10\n(Plant Community LDMC)",
       y = "Abundance at t = 100",
       color="Herbivore spp") +
  ylim(0,15)+
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

comm1a/comm2a


## calculate CWM herbivore traits ####
resultsfd <- map_dfr(R1_values, function(R1_val) {
  R <- R_base
  R[7:10] <- R1_val
  
  out <- ode(y = state, times = times, func = competition_model, parms = list(R = R))
  
  out <- as.data.frame(out)
  colnames(out) <- c("time", paste0("P", 1:n_plants), paste0("H", 1:n_herbs))
  
  final <- out[nrow(out), -1]  # drop time column
  herbivores <- final[paste0("H", 1:n_herbs)]
  
  total_H <- sum(herbivores)
  if (total_H > 0) {
    cwm_mandible <- sum(herbivores * mandible) / total_H
    fdis_mandible <- sum(herbivores * abs(mandible - cwm_mandible)) / total_H
  } else {
    cwm_mandible <- NA
    fdis <- NA
  }
  
  plants <- final[paste0("P", 1:n_plants)]
  
  total_P <- sum(plants)
  if (total_P > 0) {
    cwm_toughness <- sum(plants * toughness) / total_P
    fdis_toughness <- sum(plants * abs(toughness - cwm_toughness)) / total_P
  } else {
    cwm_toughness <- NA
    fdis_toughness <- NA
  }
  
  tibble(
    R1 = R1_val,
    cwm_mandible = cwm_mandible,
    fdis_mandible = fdis_mandible,
    cwm_toughness = cwm_toughness,
    fdis_toughness = fdis_toughness,
    !!!final
  )
})




## Matching - SES ####
### Rao's randomization -- mandible trait ####
mandibleT <- (c('H1','H2','H3','H4','H5','H6','H7','H8','H9','H10'))
t8 <- cbind(as_tibble(mandible),mandibleT) %>% column_to_rownames("mandibleT")

comm3 <- resultsfd %>% dplyr::select(H1:H10) %>% as.matrix()

frao_BV <- hill_func(comm=comm3, traits=t8, q=1)[1,]
frao_BVr <- data.frame(data.frame(matrix(ncol = 10, nrow = 21)))

for (i in 1:999) {#For each row in the matrix (for each site)
  #select randomly, as many species as the species richness of the site:
  comm3rand <- comm3
  #colnames(comm3rand) <- sample(colnames(comm3), replace=F) 
  #comm3rand <- comm3rand %>% select(order(colnames(.)))
  
  t8rand <- t8 %>% as.data.frame()
  t8rownames <- rownames(t8)
  t8rand <- t8rand[sample(nrow(t8rand)),]%>% as.data.frame()
  rownames(t8rand) <- t8rownames
  colnames(t8rand) <- "value"
  
  frao_BVr[,i] <- hill_func(comm3rand, t8rand, q=1)[1,]
}

BV_SEStab <- frao_BVr %>% 
  cbind(.,frao_BV) %>%
  rowwise() %>% 
  mutate(rand_mean = mean(c_across(2:1000)),rand_sd = sd(c_across(2:1000))) %>% 
  select(frao_BV,rand_mean,rand_sd) %>% 
  mutate(SES_BV=(frao_BV-rand_mean)/rand_sd, pval=pnorm(q=SES_BV, lower.tail = T)) 

summary(BV_SEStab)  
hist(BV_SEStab$SES_BV)
resultsfd$SES_Mandible <- BV_SEStab$SES_BV

### Rao's randomization -- toughness trait ####
mandibleT <- (c('P1','P2','P3','P4','P5','P6','P7','P8','P9','P10'))
t8 <- cbind(as_tibble(toughness),mandibleT) %>% column_to_rownames("mandibleT")

comm3 <- resultsfd %>% dplyr::select(P1:P10) %>% as.matrix()

frao_BV <- hill_func(comm=comm3, traits=t8, q=1)[1,]
frao_BVr <- data.frame(data.frame(matrix(ncol = 10, nrow = 21)))

for (i in 1:999) {#For each row in the matrix (for each site)
  #select randomly, as many species as the species richness of the site:
  comm3rand <- comm3
  #colnames(comm3rand) <- sample(colnames(comm3), replace=F) 
  #comm3rand <- comm3rand %>% select(order(colnames(.)))
  
  t8rand <- t8 %>% as.data.frame()
  t8rownames <- rownames(t8)
  t8rand <- t8rand[sample(nrow(t8rand)),]%>% as.data.frame()
  rownames(t8rand) <- t8rownames
  colnames(t8rand) <- "value"
  
  frao_BVr[,i] <- hill_func(comm3rand, t8rand, q=1)[1,]
}

BV_SEStab <- frao_BVr %>% 
  cbind(.,frao_BV) %>%
  rowwise() %>% 
  mutate(rand_mean = mean(c_across(2:1000)),rand_sd = sd(c_across(2:1000))) %>% 
  select(frao_BV,rand_mean,rand_sd) %>% 
  mutate(SES_BV=(frao_BV-rand_mean)/rand_sd, pval=pnorm(q=SES_BV, lower.tail = T)) 

summary(BV_SEStab)  
hist(BV_SEStab$SES_BV)
resultsfd$SES_toughness <- BV_SEStab$SES_BV

### FD SES Plots ####
cwm1a <- ggplot(resultsfd, aes(x = R1, y = cwm_mandible)) +
  geom_line(size = 1.2, color = "#cf4446") +
  labs(
    #title = "Herbivore community metrics across R1",
    x = "Growth Rate of P7:P10\n(Plant Community LDMC)",
    y = "CWM Mandible Strength") +
  ylim(2.75,4.75)+
  theme_bw(base_size = 16)

fd1a <- ggplot(resultsfd, aes(x = R1, y = fdis_mandible)) +
  geom_line(size = 1.2, color = "#cf4446") +
  labs(
    #title = "Functional Dispersion of Mandible Strength Across R1",
    x = "Growth Rate of P7:P10\n(Plant Community LDMC)",
    y = "FD Mandible Strength"
  ) +
  theme_bw(base_size = 16)

ses1a <- ggplot(resultsfd, aes(x = R1, y = SES_Mandible)) +
  geom_hline(aes(yintercept = 0), color="grey", linewidth=1.5, linetype="dashed")+
  geom_line(size = 1.2, color = "#cf4446") +
  labs(
    #title = "Functional Dispersion of Mandible Strength Across R1",
    x = "Growth Rate of P7:P10\n(Plant Community LDMC)",
    y = "SES Mandible Strength"
  ) +
  ylim(-3,3)+
  theme_bw(base_size = 16)

commplot <- (comm2/cwm1/ses1)+
  plot_layout(axis_titles = "collect")
commplota <- (comm2a/cwm1a/ses1a)+
  plot_layout(axis_titles = "collect")

gradplot <- commplot|commplota  ## needs updating
#ggsave(gradplot, file="gradplot.tiff", width=8, height=10, dpi=300, compression = "lzw")


# COEXISTENCE ACROSS GRADIENTS ####
# ── Three plant community scenarios with LDMC gradient ────────────────────────
# Scenario 1: R1=1.0, R10=3.0 (low LDMC plants dominant - soft plant community)
# Scenario 2: R1=2.0, R10=2.0 (equal growth rates - null community)
# Scenario 3: R1=3.0, R10=1.0 (high LDMC plants dominant - tough plant community)
#
# Linear gradient between R1 and R10, community mean R = 2.0 in all scenarios [1]

# ── Setup ─────────────────────────────────────────────────────────────────────
n_plants  <- 10
n_herbs   <- 10

toughness <- seq(30, 50, length.out = n_plants)   # LDMC values
mandible  <- seq(2.5, 4.5, length.out = n_herbs)

# Trait-based attack rates using Gaussian trait matching
beta <- 10
a0 <- 0.2
T <- .05 * n_plants  # ensures total attack rate matches a generalist's base rate

A_raw <- matrix(0, nrow = n_plants, ncol = n_herbs)
for (i in 1:n_plants) {
  for (j in 1:n_herbs) {
    ratio <- (mandible[j] / mean_mandible) / (toughness[i] / mean_toughness)
    similarity <- exp(-beta * (ratio - 1)^2)
    A_raw[i, j] <- similarity
  }
}

# Normalize so each herbivore's total attack strength across plants = T
A_fixed <- sweep(A_raw, 2, colSums(A_raw), FUN = "/") * T

# Pre-compute fixed attack matrix (trait matching) [1]
#A_fixed <- make_attack_matrix(toughness, mandible, beta = 10, T_scale = 0.5)
C_conv  <- compute_conversion_factors(A_fixed)

# ── Define three R vectors ────────────────────────────────────────────────────
make_R_gradient <- function(R_low, R_high, n = 10) {
  # Linear gradient from R_low (soft plants, low LDMC) to R_high (tough plants)
  # Mean is always 2.0 by construction when R_low + R_high = 4 [1]
  seq(R_low, R_high, length.out = n)
}

scenarios <- list(
  scenario1 = list(
    name  = "Soft plant community\n(R1=3, R10=1)",
    R_vec = make_R_gradient(3.0, 1.0),
    color = "#35b779"
  ),
  scenario2 = list(
    name  = "Null community\n(R1=2, R10=2)",
    R_vec = make_R_gradient(2.0, 2.0),
    color = "#808080"
  ),
  scenario3 = list(
    name  = "Tough plant community\n(R1=1, R10=3)",
    R_vec = make_R_gradient(1.0, 3.0),
    color = "#d44842"
  )
)

# Verify mean R = 2.0 for all scenarios
sapply(scenarios, function(s) mean(s$R_vec))

# ── Fixed parameters (shared across scenarios) ────────────────────────────────
parms_base <- list(
  n_plants = n_plants,
  n_herbs  = n_herbs,
  A        = A_fixed,
  C        = {cc <- matrix(0, n_plants, n_plants); diag(cc) <- 0.04; cc},
  HNDL     = matrix(0,   n_plants, n_herbs),
  E        = matrix(0.1, n_plants, n_herbs),
  M        = rep(1,   n_herbs),
  g        = rep(0.001, n_herbs)
)

times_long  <- seq(0, 1000, by = 0.5)
times_short <- seq(0, 10,   by = 0.1)

# ── Run Spaak N and F for each scenario ──────────────────────────────────────
run_scenario_NF <- function(scenario) {
  
  # Update R in parameters
  parms <- parms_base
  parms$R <- scenario$R_vec
  
  # Run for all herbivore species
  results_NF <- map_dfr(1:n_herbs, function(i) {
    gr <- compute_growth_rates(
      focal_herb  = i,
      C_conv      = C_conv,
      parms_fixed = parms,
      times_long  = times_long,
      times_short = times_short
    )
    if (is.null(gr)) return(tibble(herb = i, N = NA, F_val = NA,
                                   persists = NA, lhs = NA, rhs = NA,
                                   fi_00 = NA, fi_0Nj = NA, fi_cNj = NA))
    nf <- compute_N_and_F(gr)
    cx <- check_coexistence_spaak(nf$N, nf$F)
    tibble(
      herb     = i,
      fi_00    = gr$fi_00,
      fi_0Nj   = gr$fi_0Nj,
      fi_cNj   = gr$fi_cNj,
      N        = nf$N,
      F_val    = nf$F,
      lhs      = cx$lhs,
      rhs      = cx$rhs,
      persists = cx$persists
    )
  })
  
  results_NF$scenario <- scenario$name
  results_NF$cwm_ldmc <- sum(scenario$R_vec * toughness) / sum(scenario$R_vec)
  return(results_NF)
}

# Run all three scenarios
all_results <- map_dfr(scenarios, run_scenario_NF)

# Add mandible strength for each herbivore
all_results <- all_results %>%
  mutate(mandible_str = mandible[herb],
         scenario     = factor(scenario, levels = sapply(scenarios, `[[`, "name")))

print(all_results %>% select(scenario, herb, N, F_val, persists))

## Plot 1: N-F coexistence space faceted by scenario ####
scenario_colors <- setNames(
  sapply(scenarios, `[[`, "color"),
  sapply(scenarios, `[[`, "name")
)

all_results <- all_results %>%
  mutate(
    F_2021 = -F_val / (1 - (fi_cNj / fi_00))
  )

p_NF <- ggplot(all_results, aes(x = N, y = F_2021, 
                   shape = persists, 
                   label = herb,
                   color=as.factor(herb))) +
  geom_point(size = 4, stroke=2) +
  #geom_text(aes(y = F_2021 + ifelse(herb %in% c(5,6,7,8), -0.01, 0.01)),
  #          size = 3.5, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  stat_function(fun = function(x) x / (1 - x),
                color = "black", linetype = "dashed", linewidth=2,
                xlim = c(-0.5, 0.95)) +
  scale_color_viridis(option="inferno", discrete=T, direction=-1, begin=.2, end=.8) +
  #scale_fill_viridis(option="inferno", discrete=T, direction=-1, begin=.2, end=.8) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 4),
                     labels = c("Excluded", "Persists")) +
  coord_cartesian(ylim = c(-.01, .12), xlim=c(-.01, .1)) +
  facet_wrap(~scenario, ncol = 3) +
  labs(x = "Niche difference (N)",
       y = "Fitness difference (F)",
       color = "Herbivore spp"
       #title = "Herbivore coexistence across plant community LDMC gradient"
       ) +
  theme_bw(base_size = 18) +
  theme(strip.text = element_text(size = 14))

p_NF

ggsave("Fig_3_p_NF.png", p_NF,
       width = 12, height = 4, dpi = 300)

# ── Plot 2: Number of persisting species per scenario ─────────────────────────
persist_summary <- all_results %>%
  group_by(scenario, cwm_ldmc) %>%
  summarise(
    n_persist  = sum(persists, na.rm = TRUE),
    n_excluded = sum(!persists, na.rm = TRUE),
    .groups = "drop"
  )

p_persist <- ggplot(persist_summary,
                    aes(x = cwm_ldmc, y = n_persist, fill = scenario)) +
  geom_col(width = 1.5, color = "white") +
  scale_fill_manual(values = scenario_colors) +
  labs(x = "Community weighted mean LDMC",
       y = "Number of persisting herbivore species",
       fill = "Scenario") +
  scale_y_continuous(limits=c(0, n_herbs), breaks=seq(0,10,2)) +
  theme_bw(base_size = 16)

print(p_persist)

# ── Plot 3: N and F vs mandible strength across scenarios ─────────────────────
p_N_mandible <- ggplot(all_results,
                       aes(x = mandible_str, y = N,
                           color = scenario, shape = persists)) +
  geom_point(size = 3.5, stroke=2) +
  geom_line(aes(group = scenario), linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = scenario_colors) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 4),
                     labels = c("Excluded", "Persists")) +
  labs(x = "Mandible strength",
       y = "Niche difference (N)",
       color = "Scenario", shape = "Outcome") +
  theme_bw(base_size = 13)

p_F_mandible <- ggplot(all_results,
                       aes(x = mandible_str, y = F_val,
                           color = scenario, shape = persists)) +
  geom_point(size = 3.5, stroke=2) +
  geom_line(aes(group = scenario), linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = scenario_colors) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 4),
                     labels = c("Excluded", "Persists")) +
  labs(x = "Mandible strength",
       y = "Fitness difference (F)",
       color = "Scenario", shape = "Outcome") +
  theme_bw(base_size = 13)

p_mandible_plot <- p_N_mandible / p_F_mandible +
  plot_layout(guides = "collect")

print(p_mandible_plot)

# ── Combined figure ───────────────────────────────────────────────────────────
combined_plot <- p_NF / (p_N_mandible | p_F_mandible) +
  plot_layout(heights = c(1.2, 1))



## Calculate fdis_toughness for three plant community scenarios ####
# Scenarios defined earlier:
# Scenario 1 (Soft):  R1=3.0, R10=1.0
# Scenario 2 (Null):  R1=2.0, R10=2.0
# Scenario 3 (Tough): R1=1.0, R10=3.0

scenario_fd <- data.frame()

for (s in names(scenarios)) {
  scenario <- scenarios[[s]]
  
  parms <- parms_base
  parms$R <- scenario$R_vec
  
  state_init <- c(P = rep(15, n_plants), H = rep(0, n_herbs))
  times_eq   <- seq(0, 500, by = 1)
  
  out <- ode(y = state_init, times = times_eq,
             func = competition_model, parms = parms,
             method = "lsoda")
  
  out_df <- as.data.frame(out)
  colnames(out_df) <- c("time", paste0("P", 1:n_plants), paste0("H", 1:n_herbs))
  
  final   <- out_df[nrow(out_df), ]
  plants  <- as.numeric(final[paste0("P", 1:n_plants)])
  total_P <- sum(plants)
  
  if (total_P > 0) {
    cwm_toughness  <- sum(plants * toughness) / total_P
    fdis_toughness <- sum(plants * abs(toughness - cwm_toughness)) / total_P
  } else {
    cwm_toughness  <- NA
    fdis_toughness <- NA
  }
  
  scenario_fd <- rbind(scenario_fd, data.frame(
    scenario       = scenario$name,
    R_vec_summary  = paste0("R1=", round(scenario$R_vec[1], 1),
                            ", R10=", round(scenario$R_vec[10], 1)),
    cwm_toughness  = cwm_toughness,
    fdis_toughness = fdis_toughness
  ))
}

print(scenario_fd)

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot(scenario_fd, aes(x = scenario, y = fdis_toughness, fill = scenario)) +
  geom_col(color = "white", width = 0.6) +
  geom_text(aes(label = round(fdis_toughness, 2)),
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("#35b779", "grey50", "#d44842")) +
  labs(x = "Scenario",
       y = "Functional Dispersion of LDMC",
       title = "Plant functional dispersion across LDMC gradient scenarios") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

