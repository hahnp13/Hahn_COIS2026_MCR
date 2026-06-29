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
A <- sweep(A_raw, 2, colSums(A_raw), FUN = "/") * T

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
  ylim(0,10)+
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


## calculating ND and RFD ####
# ── ND and RFD from MacArthur's CRM (Godwin et al. 2020, Eq. 12 & 13) ────────
#
# model parameters map to Godwin's notation as:
#   c_il  = A[l, i]   (attack/capture rate of herbivore i on plant l)
#   w_il  = E[l, i]   (conversion efficiency, yield per unit resource consumed)
#   K_l   = R[l] / C[l,l]  (plant carrying capacity = r / alpha_ll)
#   r_l   = R[l]      (plant intrinsic growth rate)
#   m_i   = M[i]      (herbivore mortality)

# ── Function to compute pairwise ND and RFD for all herbivore pairs ───────────
# Extract equilibrium plant abundances from ODE at t = 100
# Instead of theoretical K = R / diag(C)
equil_plants <- as.numeric(out[nrow(out), paste0("P", 1:n_plants)])

# Then substitute into ND/RFD calculations
compute_ND_RFD_empirical <- function(A, E, R, C, M, P_eq) {
  # P_eq: vector of actual plant abundances at equilibrium [n_plants]
  # replaces K in Godwin Eq. 12 and 13 [2]
  
  n_herbs  <- ncol(A)
  n_plants <- nrow(A)
  
  # Use actual plant abundances instead of K
  # K_eff still needs r_l scaling per Eq. 12: K_l / r_l = P_eq / R
  K_over_r <- P_eq / R   # replaces K_l / r_l term in Eq. 12
  K_eff    <- P_eq       # replaces K_l term in Eq. 13

  results <- tibble(
    herb_i = integer(), herb_j = integer(),
    ND = numeric(), RFD = numeric()
  )

  for (i in 1:(n_herbs - 1)) {
    for (j in (i + 1):n_herbs) {

      # Eq. 12 with P_eq substituted for K [2]
      niche_overlap <- sum(A[, i] * A[, j] * E[, i] * K_over_r)
      ND <- 1 - niche_overlap

      # Eq. 13 with P_eq substituted for K [2]
      numerator   <- sum(A[, j] * E[, i] * K_eff) - M[j]
      denominator <- sum(A[, i] * E[, j] * K_eff) - M[i]
      RFD <- numerator / denominator

      results <- results %>%
        add_row(herb_i = i, herb_j = j, ND = ND, RFD = RFD)
    }
  }
  return(results)
}

# Run with equilibrium plant abundances
nd_rfd_empirical <- compute_ND_RFD_empirical(
  A = A, E = E_mat, R = R_vec, C = C_mat, M = M_vec,
  P_eq = equil_plants
)
nd_rfd_empirical_coex <- check_coexistence(nd_rfd_empirical)

print(nd_rfd_coex, n = 45)

# ── Visualize ─────────────────────────────────────────────────────────────────
ggplot(nd_rfd_coex, aes(x = ND, y = RFD, color = outcome)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_abline(slope = -1, intercept = 1, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c(
    "Coexistence"        = "#35b779",
    "Species i excluded" = "#d44842",
    "Species j excluded" = "#3b528b"
  )) +
  labs(
    x     = "Niche Difference (ND)",
    y     = "Relative Fitness Difference (RFD)",
    color = "Outcome"
  ) +
  theme_bw(base_size = 14)

## Mutual invasibility test for all herbivore pairs ##########################
# Following Godwin et al. [2] sensitivity method logic:
# 1. Run resident to equilibrium (monoculture)
# 2. Introduce invader at low density (~1% of K)
# 3. Measure invader per capita growth rate when rare

invasion_threshold <- 0.01   # invader starting density
times_mono  <- seq(0, 500, by = 0.1)   # long enough to reach equilibrium
times_inv   <- seq(0, 50,  by = 0.1)   # short window to measure invasion growth

run_monoculture <- function(focal_herb, params_est) {
  # Run single herbivore species with all plants to equilibrium
  state_mono <- c(P = rep(15, n_plants),
                  H = rep(0,  n_herbs))
  state_mono[n_plants + focal_herb] <- 3   # only focal herbivore present
  
  out <- tryCatch(
    ode(y = state_mono, times = times_mono,
        func = competition_model, parms = params_est,
        method = "lsoda"),
    error = function(e) NULL
  )
  return(out)
}

run_invasion <- function(resident_herb, invader_herb, resident_equil, params_est) {
  # Start from resident equilibrium, add rare invader
  state_inv <- resident_equil
  state_inv[n_plants + invader_herb] <- invasion_threshold
  
  out <- tryCatch(
    ode(y = state_inv, times = times_inv,
        func = competition_model, parms = params_est,
        method = "lsoda"),
    error = function(e) NULL
  )
  return(out)
}

compute_invasion_growth_rate <- function(out, invader_herb) {
  # Per capita growth rate of invader over early invasion window [2]
  df     <- as.data.frame(out)
  h_col  <- paste0("H", invader_herb)
  N_t0   <- df[1,   h_col]
  N_tend <- df[nrow(df), h_col]
  t_span <- df$time[nrow(df)] - df$time[1]
  
  # Exponential growth rate: log(N_t / N_0) / t
  igr <- log(N_tend / N_t0) / t_span
  return(igr)
}

# ── Run all pairwise invasions ────────────────────────────────────────────────
# Use your fitted parameters here; using fixed defaults for illustration
#parms_fixed <- list(
#  n_plants = n_plants, n_herbs = n_herbs,
#  A    = A,
#  R    = rep(2,   n_plants),
#  C    = {cc <- matrix(0, n_plants, n_plants); diag(cc) <- 0.04; cc},
#  HNDL = matrix(0,   n_plants, n_herbs),
#  E    = matrix(0.1, n_plants, n_herbs),
#  M    = rep(1,   n_herbs),
#  g    = rep(0.001, n_herbs)
#)

#invasibility_results <- tibble(
#  resident = integer(), invader = integer(),
#  igr      = numeric(), can_invade = logical()
#)

for (res in 1:n_herbs) {
  
  # Run resident to equilibrium
  mono_out <- run_monoculture(res, parms_fixed)
  if (is.null(mono_out)) next
  
  # Extract equilibrium state (last time point)
  equil_state <- as.numeric(mono_out[nrow(mono_out),
                                     2:(n_plants + n_herbs + 1)])
  names(equil_state) <- c(paste0("P", 1:n_plants),
                          paste0("H", 1:n_herbs))
  
  for (inv in 1:n_herbs) {
    if (inv == res) next
    
    inv_out <- run_invasion(res, inv, equil_state, parms_fixed)
    if (is.null(inv_out)) next
    
    igr <- compute_invasion_growth_rate(inv_out, inv)
    
    invasibility_results <- invasibility_results %>%
      add_row(resident   = res,
              invader    = inv,
              igr        = igr,
              can_invade = igr > 0)
  }
}

# ── Mutual invasibility: both species must be able to invade each other ───────
mutual_invasibility <- invasibility_results %>%
  left_join(invasibility_results %>%
              rename(resident2 = resident, invader2 = invader,
                     igr2 = igr, can_invade2 = can_invade),
            by = c("resident" = "invader2", "invader" = "resident2")) %>%
  filter(resident < invader) %>%
  mutate(mutual_invasibility = can_invade & can_invade2) %>%
  select(resident, invader, igr_res_inv = igr, igr_inv_res = igr2,
         mutual_invasibility)

print(mutual_invasibility)

# ── Visualize invasion growth rates ──────────────────────────────────────────
ggplot(invasibility_results,
       aes(x = factor(invader), y = factor(resident), fill = igr)) +
  geom_tile() +
  geom_text(aes(label = round(igr, 3)), size = 3) +
  scale_fill_gradient2(low = pal[25], mid = "white", high = pal[225],
                       midpoint = 0) +
  labs(x = "Invader", y = "Resident",
       fill = "Invasion\ngrowth rate") +
  theme_bw(base_size = 14)

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



## Resource Utilization Functions (Sakarchi & Germain 2025, Eq. 2) [6] ####
# Ui = sum_l( A[l,i] * sqrt(E[l,i] * Kl / R[l]) )
# where Kl = R[l] / diag(C)[l]  (plant carrying capacity)
# E[l,i] = conversion efficiency (nutritional weight wk in Sakarchi notation)
# A[l,i] = attack rate (aik in Sakarchi notation)
# R[l]   = plant intrinsic growth rate (rk in Sakarchi notation)

compute_utilization <- function(A, E, R_vec, C_mat) {
  n_plants <- nrow(A)
  n_herbs  <- ncol(A)
  
  # Plant carrying capacities
  K <- R_vec / diag(C_mat)   # K_l = r_l / alpha_ll [1]
  
  # Resource exploitability: sqrt(wk * Kk / rk) = sqrt(E[l,i] * K[l] / R[l])
  # Note: E is constant across species in current model [1],
  # but code handles species-specific E for generality
  exploitability <- sqrt(sweep(E, 1, K / R_vec, FUN = "*"))
  # exploitability[l, i] = sqrt(E[l,i] * K[l] / R[l])
  
  # Utilization: sum over plants for each herbivore
  U <- colSums(A * exploitability)   # length n_herbs
  names(U) <- paste0("H", 1:n_herbs)
  
  # Also compute utilization matrix for plotting
  U_mat <- A * exploitability        # [n_plants x n_herbs]
  
  return(list(U = U, U_mat = U_mat,
              exploitability = K / R_vec))  # K/r per plant [6]
}

# ── Run for trait matching model ──────────────────────────────────────────────
R_vec  <- rep(2, n_plants)
C_mat  <- diag(0.04, n_plants)
E_mat  <- matrix(0.1, n_plants, n_herbs)

util <- compute_utilization(A_fixed, E_mat, R_vec, C_mat)

# Total utilization per herbivore
print(util$U)

### Plot 1: Utilization profiles across plant resource spectrum ####
# analogous to Figure 1C/D in Sakarchi & Germain [6]
util_long <- as_tibble(util$U_mat) %>%
  setNames(paste0("H", 1:n_herbs)) %>%
  mutate(plant    = paste0("P", 1:n_plants),
         toughness = toughness,
         K_over_r  = util$exploitability) %>%
  pivot_longer(starts_with("H"),
               names_to  = "herbivore",
               values_to = "utilization") %>%
  mutate(herbivore = factor(herbivore,
                            levels = paste0("H", 1:n_herbs)))

p_util <- ggplot(util_long,
                 aes(x = toughness, y = utilization,
                     color = herbivore, group = herbivore)) +
  geom_line(linewidth = 1) +
  scale_color_viridis_d(option = "inferno", begin = 0.2, end = 0.9) +
  labs(x     = "Plant LDMC (toughness)",
       y     = "Resource utilization (Ui_k)",
       color = "Herbivore",
       title = "Utilization functions across plant resource spectrum") +
  theme_bw(base_size = 14)

print(p_util)

### Plot 2: Total utilization vs mandible strength ####
util_df <- tibble(
  herbivore    = paste0("H", 1:n_herbs),
  mandible_str = mandible,
  U_total      = util$U
)

p_U_total <- ggplot(util_df, aes(x = mandible_str, y = U_total)) +
  geom_line(linewidth = 1.2, color = "#d44842") +
  geom_point(size = 3, color = "#d44842") +
  labs(x     = "Mandible strength",
       y     = "Total utilization (Ui)",
       title = "Total resource utilization by herbivore species") +
  theme_bw(base_size = 14)

print(p_U_total)

### Plot 3: Utilization overlap between herbivore pairs ####
# cosine similarity of utilization vectors = niche overlap sensu Sakarchi [6]
U_mat_norm <- sweep(util$U_mat, 2, sqrt(colSums(util$U_mat^2)), "/")
overlap_mat <- t(U_mat_norm) %*% U_mat_norm   # [n_herbs x n_herbs]

overlap_long <- as_tibble(overlap_mat) %>%
  setNames(paste0("H", 1:n_herbs)) %>%
  mutate(herb_i = paste0("H", 1:n_herbs)) %>%
  pivot_longer(-herb_i, names_to = "herb_j", values_to = "overlap")

p_overlap <- ggplot(overlap_long,
                    aes(x = herb_i, y = herb_j, fill = overlap)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", limits = c(0, 1)) +
  labs(x    = "Herbivore i",
       y    = "Herbivore j",
       fill = "Utilization\noverlap",
       title = "Pairwise utilization overlap") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_overlap)


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
a0 <- 0.2
T <- .03 * n_plants  # ensures total attack rate matches a generalist's base rate

# Strong mandibles are better
A <- outer(toughness, mandible, function(t, m) (m / mean(m)) / (t / mean(t))) * a0 * T

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
g <- rep(0.01, n_herbs)  # herbivore self-limitation rate

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
ggsave(Strong_plot, file="Strong_plot.tiff", width=8, height=6, dpi=300, compression = "lzw")


# ND and RFD calculations ####
# Replace these with your fitted parameter values after running mle2()
A_mat <- A                                         # [n_plants x n_herbs] [1]
E_mat <- E                                               # replace with fitted e_val
R_vec <- R                                               # replace with fitted R values
C_mat <- C                            # replace with fitted c_diag
M_vec <- M                                 # replace with fitted M values

# ── Run ───────────────────────────────────────────────────────────────────────
nd_rfd <- compute_ND_RFD(A = A_mat, E = E_mat, R = R_vec, C = C_mat, M = M_vec)
nd_rfd_coex <- check_coexistence(nd_rfd)

print(nd_rfd_coex, n = 45)

# ── Visualize ─────────────────────────────────────────────────────────────────
ggplot(nd_rfd_coex, aes(x = ND, y = RFD, color = outcome)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_abline(slope = -1, intercept = 1, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c(
    "Coexistence"        = "#35b779",
    "Species i excluded" = "#d44842",
    "Species j excluded" = "#3b528b"
  )) +
  labs(
    x     = "Niche Difference (ND)",
    y     = "Relative Fitness Difference (RFD)",
    color = "Outcome"
  ) +
  theme_bw(base_size = 14)

## invasion criteria ####
invasion_threshold <- 0.01   # invader starting density
times_mono  <- seq(0, 500, by = 0.1)   # long enough to reach equilibrium
times_inv   <- seq(0, 50,  by = 0.1)   # short window to measure invasion growth

run_monoculture <- function(focal_herb, params_est) {
  # Run single herbivore species with all plants to equilibrium
  state_mono <- c(P = rep(15, n_plants),
                  H = rep(0,  n_herbs))
  state_mono[n_plants + focal_herb] <- 3   # only focal herbivore present
  
  out <- tryCatch(
    ode(y = state_mono, times = times_mono,
        func = competition_model, parms = params_est,
        method = "lsoda"),
    error = function(e) NULL
  )
  return(out)
}

run_invasion <- function(resident_herb, invader_herb, resident_equil, params_est) {
  # Start from resident equilibrium, add rare invader
  state_inv <- resident_equil
  state_inv[n_plants + invader_herb] <- invasion_threshold
  
  out <- tryCatch(
    ode(y = state_inv, times = times_inv,
        func = competition_model, parms = params_est,
        method = "lsoda"),
    error = function(e) NULL
  )
  return(out)
}

compute_invasion_growth_rate <- function(out, invader_herb) {
  # Per capita growth rate of invader over early invasion window [2]
  df     <- as.data.frame(out)
  h_col  <- paste0("H", invader_herb)
  N_t0   <- df[1,   h_col]
  N_tend <- df[nrow(df), h_col]
  t_span <- df$time[nrow(df)] - df$time[1]
  
  # Exponential growth rate: log(N_t / N_0) / t
  igr <- log(N_tend / N_t0) / t_span
  return(igr)
}

# ── Run all pairwise invasions ────────────────────────────────────────────────
# Use your fitted parameters here; using fixed defaults for illustration

parms_fixed <- list(
  n_plants = n_plants, n_herbs = n_herbs,
  A    = A,
  R    = rep(2,   n_plants),
  C    = C,
  HNDL = matrix(0,   n_plants, n_herbs),
  E    = matrix(0.1, n_plants, n_herbs),
  M    = rep(1,   n_herbs),
  g    = rep(0.001, n_herbs)
)
invasibility_results <- tibble(
  resident = integer(), invader = integer(),
  igr      = numeric(), can_invade = logical()
)

for (res in 1:n_herbs) {
  
  # Run resident to equilibrium
  mono_out <- run_monoculture(res, parms_fixed)
  if (is.null(mono_out)) next
  
  # Extract equilibrium state (last time point)
  equil_state <- as.numeric(mono_out[nrow(mono_out),
                                     2:(n_plants + n_herbs + 1)])
  names(equil_state) <- c(paste0("P", 1:n_plants),
                          paste0("H", 1:n_herbs))
  
  for (inv in 1:n_herbs) {
    if (inv == res) next
    
    inv_out <- run_invasion(res, inv, equil_state, parms_fixed)
    if (is.null(inv_out)) next
    
    igr <- compute_invasion_growth_rate(inv_out, inv)
    
    invasibility_results <- invasibility_results %>%
      add_row(resident   = res,
              invader    = inv,
              igr        = igr,
              can_invade = igr > 0)
  }
}

# ── Mutual invasibility: both species must be able to invade each other ───────
mutual_invasibility <- invasibility_results %>%
  left_join(invasibility_results %>%
              rename(resident2 = resident, invader2 = invader,
                     igr2 = igr, can_invade2 = can_invade),
            by = c("resident" = "invader2", "invader" = "resident2")) %>%
  filter(resident < invader) %>%
  mutate(mutual_invasibility = can_invade & can_invade2) %>%
  select(resident, invader, igr_res_inv = igr, igr_inv_res = igr2,
         mutual_invasibility)

print(mutual_invasibility, n = 45)

# ── Visualize invasion growth rates ──────────────────────────────────────────
ggplot(invasibility_results,
       aes(x = factor(invader), y = factor(resident), fill = igr)) +
  geom_tile() +
  geom_text(aes(label = round(igr, 3)), size = 3) +
  scale_fill_gradient2(low = "#d44842", mid = "white", high = "#35b779",
                       midpoint = 0) +
  labs(x = "Invader", y = "Resident",
       fill = "Invasion\ngrowth rate") +
  theme_bw(base_size = 14)

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
ggsave(NF_plot, file="NF_plot.png", width=12, height=5, dpi=300)


## Resource Utilization Functions (Sakarchi & Germain 2025, Eq. 2) [6] ####
# Ui = sum_l( A[l,i] * sqrt(E[l,i] * Kl / R[l]) )
# where Kl = R[l] / diag(C)[l]  (plant carrying capacity)
# E[l,i] = conversion efficiency (nutritional weight wk in Sakarchi notation)
# A[l,i] = attack rate (aik in Sakarchi notation)
# R[l]   = plant intrinsic growth rate (rk in Sakarchi notation)

compute_utilization <- function(A, E, R_vec, C_mat) {
  n_plants <- nrow(A)
  n_herbs  <- ncol(A)
  
  # Plant carrying capacities
  K <- R_vec / diag(C_mat)   # K_l = r_l / alpha_ll [1]
  
  # Resource exploitability: sqrt(wk * Kk / rk) = sqrt(E[l,i] * K[l] / R[l])
  # Note: E is constant across species in current model [1],
  # but code handles species-specific E for generality
  exploitability <- sqrt(sweep(E, 1, K / R_vec, FUN = "*"))
  # exploitability[l, i] = sqrt(E[l,i] * K[l] / R[l])
  
  # Utilization: sum over plants for each herbivore
  U <- colSums(A * exploitability)   # length n_herbs
  names(U) <- paste0("H", 1:n_herbs)
  
  # Also compute utilization matrix for plotting
  U_mat <- A * exploitability        # [n_plants x n_herbs]
  
  return(list(U = U, U_mat = U_mat,
              exploitability = K / R_vec))  # K/r per plant [6]
}

# ── Run for trait matching model ──────────────────────────────────────────────
R_vec  <- rep(2, n_plants)
C_mat  <- diag(0.04, n_plants)
E_mat  <- matrix(0.1, n_plants, n_herbs)

util <- compute_utilization(A_fixed, E_mat, R_vec, C_mat)

# Total utilization per herbivore
print(util$U)

### Plot 1: Utilization profiles across plant resource spectrum ####
# analogous to Figure 1C/D in Sakarchi & Germain [6]
util_long <- as_tibble(util$U_mat) %>%
  setNames(paste0("H", 1:n_herbs)) %>%
  mutate(plant    = paste0("P", 1:n_plants),
         toughness = toughness,
         K_over_r  = util$exploitability) %>%
  pivot_longer(starts_with("H"),
               names_to  = "herbivore",
               values_to = "utilization") %>%
  mutate(herbivore = factor(herbivore,
                            levels = paste0("H", 1:n_herbs)))

p_util <- ggplot(util_long,
                 aes(x = toughness, y = utilization,
                     color = herbivore, group = herbivore)) +
  geom_line(linewidth = 1) +
  scale_color_viridis_d(option = "inferno", begin = 0.2, end = 0.9, direction=-1) +
  labs(x     = "Plant LDMC (toughness)",
       y     = "Resource utilization (Ui_k)",
       color = "Herbivore",
       title = "Utilization functions across plant resource spectrum") +
  theme_bw(base_size = 14)

print(p_util)

### Plot 2: Total utilization vs mandible strength ####
util_df <- tibble(
  herbivore    = paste0("H", 1:n_herbs),
  mandible_str = mandible,
  U_total      = util$U
)

p_U_total <- ggplot(util_df, aes(x = mandible_str, y = U_total)) +
  geom_line(linewidth = 1.2, color = "#d44842") +
  geom_point(size = 3, color = "#d44842") +
  labs(x     = "Mandible strength",
       y     = "Total utilization (Ui)",
       title = "Total resource utilization by herbivore species") +
  theme_bw(base_size = 14)

print(p_U_total)

### Plot 3: Utilization overlap between herbivore pairs ####
# cosine similarity of utilization vectors = niche overlap sensu Sakarchi [6]
U_mat_norm <- sweep(util$U_mat, 2, sqrt(colSums(util$U_mat^2)), "/")
overlap_mat <- t(U_mat_norm) %*% U_mat_norm   # [n_herbs x n_herbs]

overlap_long <- as_tibble(overlap_mat) %>%
  setNames(paste0("H", 1:n_herbs)) %>%
  mutate(herb_i = paste0("H", 1:n_herbs)) %>%
  pivot_longer(-herb_i, names_to = "herb_j", values_to = "overlap")

p_overlap <- ggplot(overlap_long,
                    aes(x = herb_i, y = herb_j, fill = overlap)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", limits = c(0, 1)) +
  labs(x    = "Herbivore i",
       y    = "Herbivore j",
       fill = "Utilization\noverlap",
       title = "Pairwise utilization overlap") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_overlap)


# Matching vs strong mandibles plot ####
MatchStrongPlot <- popplot|popplota
ggsave(MatchStrongPlot, file="MatchStrongPlot.tiff", width=10, height=6, dpi=300, compression = "lzw")

MCR_plot <- (popplot|popplota)/NF_plot +
  plot_layout(heights = c(2,1))
ggsave(MCR_plot, file="MCR_plot.tiff", width=10, height=10, dpi=300, compression = "lzw")


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
g <- rep(0.001, n_herbs)  # herbivore self-limitation rate

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
ggsave(Arate_plot, file="Arate_plot.tiff", width=8, height=6, dpi=300, compression = "lzw")

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
HNDL <- matrix(runif(n_plants * n_herbs, 0.01, 0.01), nrow = n_plants)
E <- matrix(runif(n_plants * n_herbs, 0.1, 0.1), nrow = n_plants)
M <- rep(1, n_herbs)
g <- rep(0.001, n_herbs)  # herbivore self-limitation rate

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



# Handling time model ####
# Species counts
n_plants <- 30
n_herbs <- 10

# Trait values
toughness <- seq(30, 50, length.out=n_plants)
mandible <- seq(2.5, 4.5, length.out=n_herbs)

# Normalize traits
norm_toughness <- toughness / mean(toughness)
norm_mandible  <- mandible / mean(mandible)

# Parameters
beta <- 5        # Trait sensitivity
h_max <- .1      # Maximum handling time (for poorest match)

HNDL_raw <- matrix(0, nrow = n_plants, ncol = n_herbs)
for (i in 1:n_plants) {
  for (j in 1:n_herbs) {
    ratio <- (mandible[j] / mean_mandible) / (toughness[i] / mean_toughness)
    similarity <- exp(-beta * (ratio - 1)^2)
    HNDL_raw[i, j] <- similarity
  }
}

# Normalize so each herbivore's total attack strength across plants = T
#HNDL_max <- sweep(HNDL_raw, 2, colSums(HNDL_raw), FUN = "/") 
#HNDL <- h_max * (1-HNDL_max)
HNDL <- sweep(HNDL_raw, 2, colSums(HNDL_raw), FUN = "/") 

# Label rows and columns
rownames(HNDL) <- paste0("P", 1:3)
colnames(HNDL) <- paste0("H", 1:3)

# View the matrix
HNDL

rowSums(HNDL)
colSums(HNDL)

# Fixed parameters
R <- rep(2, n_plants)
C <- matrix(0.0, nrow = n_plants, ncol = n_plants); diag(C) <- 0.04
A <- matrix(0.5/n_plants, nrow = n_plants, ncol = n_plants)

#HNDL <- matrix(runif(n_plants * n_herbs, 0, 0), nrow = n_plants)
E <- matrix(runif(n_plants * n_herbs, 0.1, 0.1), nrow = n_plants)
M <- rep(1, n_herbs)
g <- rep(0.001, n_herbs)  # herbivore self-limitation rate

# Initial state
state <- c(P = rep(25, n_plants), H = rep(5, n_herbs))

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
pop1b <- ggplot(data=out_long %>% filter(type=="Plant"), 
                aes(x = time, y = value, color = spp)) +
  geom_line(linewidth = .8) +
  scale_color_viridis(option="viridis", discrete=T, direction=-1, begin=.2, end=.8) +
  labs(#title = "Plant density (g/m2)",
    x = "Time",
    y = "Plant biomass (g/m2)",
    color="Plant spp") +
  ylim(0,45)+
  theme_bw(base_size = 14)

pop2b <- ggplot(data=out_long %>% filter(type=="Herbivore"), 
                aes(x = time, y = value, color = spp)) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis(option="inferno", discrete=T, direction=-1, begin=.2, end=.8) +
  labs(title = "Trait-matching: handling time",
       x = "Time",
       y = "Herbivore abundance (#/m2)",
       color="Herbivore spp") +
  ylim(0,10)+
  theme_bw(base_size = 14)

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
ggplot(abund_at_100 %>% filter(type=="Plant"), 
       aes(x = trait, y = value, color=spp)) +
  geom_point(size = 3) +
  #geom_line(aes(group = type), linetype = "dashed", alpha = 0.6) +
  scale_color_viridis(option="viridis", discrete=T, direction=-1, begin=.2, end=.8) +
  facet_wrap(~type, scales = "free_x") +
  labs(
    title = "Trait Value vs. Abundance at Time = 100",
    x = "LDMC (%)",
    y = "Abundance",
    color = "Species"
  ) +
  theme_bw(base_size = 14)

ggplot(abund_at_100 %>% filter(type=="Herbivore"), 
       aes(x = trait, y = value, color=spp)) +
  geom_point(size = 3) +
  #geom_line(aes(group = type), linetype = "dashed", alpha = 0.6) +
  scale_color_viridis(option="inferno", discrete=T, direction=-1, begin=.2, end=.8) +
  facet_wrap(~type, scales = "free_x") +
  labs(
    title = "Trait Value vs. Abundance at Time = 100",
    x = "Mandible strength",
    y = "Abundance",
    color = "Species"
  ) +
  theme_bw(base_size = 14)

## heatmap of attack rates ####

# Total attack pressure on each plant
row_sums <- rowSums(A_raw)
barplot(row_sums, names.arg = paste0("P", 1:4), ylab = "Total attack pressure", main = "Plant vulnerability (∑ A_ij)")

# A heatmap
heatmap(A, Rowv=NA, Colv=NA, col = hcl.colors(100, "YlOrRd", rev=TRUE),
        labRow = paste0("P", 1:4), labCol = paste0("H", 1:4), main = "Trait-Based Attack Matrix A")


## make interaction web based on A ####
library(tidyverse)
library(tidygraph)
library(ggraph)

## effect of beta ####
x <- seq(0.5, 1.5, length.out = 200)
plot(x, exp(-0.5 * (x - 1)^2), type="l", col="blue", ylim=c(0,1),
     ylab="Similarity", xlab="Trait Ratio (Mandible / Toughness)", main="Effect of Beta")
lines(x, exp(-5 * (x - 1)^2), col="darkgreen")
lines(x, exp(-10 * (x - 1)^2), col="red")
lines(x, exp(-20 * (x - 1)^2), col="purple")
legend("topright", legend=c("β = 0.5", "β = 5", "β = 10"), col=c("blue", "darkgreen", "red"), lty=1)


# Create edge list from A matrix
interaction_df <- as_tibble(A, rownames = NULL) %>%
  mutate(plant = paste0("P", row_number())) %>%
  pivot_longer(cols = starts_with("V"), names_to = "herb_idx", values_to = "weight") %>%
  mutate(herbivore = paste0("H", parse_number(herb_idx))) %>%
  select(from = plant, to = herbivore, weight)

# Get list of all unique node names
all_species <- unique(c(interaction_df$from, interaction_df$to))

# Define node types: TRUE = Plant (bottom row), FALSE = Herbivore (top row)
nodes <- tibble(
  name = all_species,
  type = str_starts(all_species, "P")  # TRUE for plants, FALSE for herbivores
)

# Build graph
graph <- tbl_graph(nodes = nodes, edges = interaction_df, directed = TRUE)

# Plot bipartite interaction web
intplot <- ggraph(graph, layout = "bipartite") +
  geom_edge_link(aes(width = weight), alpha = 0.6, color = "gray40") +
  geom_node_point(aes(color = type), size = 12, pch=21, stroke=1.5, fill="white") +
  geom_node_text(aes(label = name, vjust = ifelse(type, .5, .5)), size = 5) +
  scale_edge_width(range = c(0.2, 3)) +
  scale_color_manual(
    values = c("TRUE" = "#35b779", "FALSE" = "#d44842"),
    labels = c("Herbivore", "Plant"),
    name = "Type") +
  theme_void(base_size = 16) 

ggsave(intplot, file="intplot.tiff", width=8, height=6, dpi=300, compression = "lzw")


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
g <- rep(0.001, n_herbs)

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

ggsave("LDMC_gradient.png", pcommplot2, width=12, height = 5, dpi=300)


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

ggsave(commplot, file="TraitMatchingPlot.tiff", width=10, height=8, dpi=300, compression = "lzw")



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

gradplot <- commplot|commplota 
ggsave(gradplot, file="gradplot.tiff", width=8, height=10, dpi=300, compression = "lzw")


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

ggsave("p_NF.png", p_NF,
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

ggsave("LDMC_coexistence_gradient.png", combined_plot,
       width = 12, height = 10, dpi = 300)


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


# COEXIST PLANTS ####
# ── N and F for plant species across three scenarios ─────────────────────────
# The invasion analysis for plants follows the same logic [3]:
#   1. Run S-1 plant community + all herbivores to equilibrium
#   2. Introduce focal plant as rare, measure invasion growth rate
#   3. Run focal plant alone at converted equilibrium density, measure no-niche rate
#   4. Run focal plant alone from low density, measure intrinsic growth rate

# ── Spaak et al. 2021 N and F definitions (Box 2) [5] ────────────────────────
# Three key growth rates per species i:
#   mu_i  = fi(0, 0)                  intrinsic growth rate
#   r_i   = fi(0, N-i*)               invasion growth rate
#   eta_i = fi(sum cij * N-i*_j, 0)   no-niche growth rate
#
# Eq. 1 [5]: N_i = (r_i - eta_i) / (mu_i - eta_i)
# Eq. 2 [5]: F_i = (eta_i/mu_i - 1) / (1 - eta_i/mu_i)
#
# Note: F_i definition differs slightly from Spaak & De Laender 2020 [3]
# to ensure continuous dependence on mu_i and positive-slope boundary [5]

compute_N_and_F_2021 <- function(gr) {
  mu_i  <- gr$fi_00    # intrinsic growth rate
  r_i   <- gr$fi_0Nj   # invasion growth rate
  eta_i <- gr$fi_cNj   # no-niche growth rate
  
  # Eq. 1 [5]: niche difference
  Ni <- (r_i - eta_i) / (mu_i - eta_i)
  
  # Eq. 2 [5]: fitness difference
  Fi <- (eta_i / mu_i - 1) / (1 - eta_i / mu_i)
  
  return(list(N = Ni, F = Fi))
}

# ── Coexistence criterion: species persists if N > F [5] ─────────────────────
# Persistence line is the diagonal N = F with positive slope
check_coexistence_2021 <- function(Ni, Fi) {
  persists <- Ni > Fi
  return(list(persists = persists))
}

compute_growth_rates_plants <- function(focal_plant, C_conv_plants,
                                        parms_fixed, times_long, times_short,
                                        invasion_density = 0.01) {
  
  n_plants <- parms_fixed$n_plants
  n_herbs  <- parms_fixed$n_herbs
  
  # --- (a) Intrinsic growth rate: fi(0, 0) ---
  # Focal plant alone at low density, no herbivores, no other plants
  state_int <- c(P = rep(0, n_plants), H = rep(0, n_herbs))
  state_int[focal_plant] <- invasion_density
  
  out_int <- tryCatch(
    ode(y = state_int, times = times_short,
        func = competition_model, parms = parms_fixed, method = "lsoda"),
    error = function(e) NULL
  )
  if (is.null(out_int)) return(NULL)
  
  df_int  <- as.data.frame(out_int)
  p_col   <- paste0("P", focal_plant)
  N_t0    <- df_int[1, p_col]
  N_tend  <- df_int[nrow(df_int), p_col]
  t_span  <- df_int$time[nrow(df_int)] - df_int$time[1]
  fi_00   <- log(N_tend / N_t0) / t_span
  
  # --- (b) Run S-1 plant community + all herbivores to equilibrium ---
  state_sub <- c(P = rep(15, n_plants), H = rep(3, n_herbs))
  state_sub[focal_plant] <- 0   # focal plant absent
  
  out_sub <- tryCatch(
    ode(y = state_sub, times = times_long,
        func = competition_model, parms = parms_fixed, method = "lsoda"),
    error = function(e) NULL
  )
  if (is.null(out_sub)) return(NULL)
  
  equil_sub <- as.numeric(out_sub[nrow(out_sub), 2:(n_plants + n_herbs + 1)])
  names(equil_sub) <- c(paste0("P", 1:n_plants), paste0("H", 1:n_herbs))
  
  # Store S-1 equilibrium plant densities (excluding focal)
  N_minus_i_star <- equil_sub[paste0("P", 1:n_plants)]
  N_minus_i_star[focal_plant] <- 0
  
  # --- (c) Invasion growth rate: fi(0, N-i*) ---
  state_inv <- equil_sub
  state_inv[focal_plant] <- invasion_density
  
  out_inv <- tryCatch(
    ode(y = state_inv, times = times_short,
        func = competition_model, parms = parms_fixed, method = "lsoda"),
    error = function(e) NULL
  )
  if (is.null(out_inv)) return(NULL)
  
  df_inv   <- as.data.frame(out_inv)
  N_t0_inv <- df_inv[1, p_col]
  N_tend_inv <- df_inv[nrow(df_inv), p_col]
  fi_0Nj   <- log(N_tend_inv / N_t0_inv) / t_span
  
  # --- (d) No-niche growth rate: fi(sum cij * N-i*_j, 0) ---
  # Conversion factors for plants from plant competition matrix C [3]
  # cij = sqrt(C[j,j] / C[i,i]) for diagonal competition matrix
  converted_density <- sum(C_conv_plants[focal_plant, ] * N_minus_i_star)
  converted_density <- max(converted_density, invasion_density)
  
  # Focal plant alone at converted density, no other plants, no herbivores
  state_nn <- c(P = rep(0, n_plants), H = rep(0, n_herbs))
  state_nn[focal_plant] <- converted_density
  
  out_nn <- tryCatch(
    ode(y = state_nn, times = times_short,
        func = competition_model, parms = parms_fixed, method = "lsoda"),
    error = function(e) NULL
  )
  if (is.null(out_nn)) return(NULL)
  
  df_nn     <- as.data.frame(out_nn)
  N_t0_nn   <- df_nn[1, p_col]
  N_tend_nn <- df_nn[nrow(df_nn), p_col]
  fi_cNj    <- log(N_tend_nn / N_t0_nn) / t_span
  
  return(list(
    fi_00  = fi_00,
    fi_0Nj = fi_0Nj,
    fi_cNj = fi_cNj,
    N_minus_i_star = N_minus_i_star,
    converted_density = converted_density
  ))
}

# ── Analytical no-niche and intrinsic growth rates for plants ─────────────────
# For a plant alone, dP/dt = P * (R[l] - C[l,l] * P)
# Per capita rate at density N: R[l] - C[l,l] * N
# So:
#   fi(0, 0)          = R[l]                   (intrinsic rate at zero density)
#   fi(converted, 0)  = R[l] - C[l,l] * converted_density  (no-niche rate)

compute_growth_rates_plants_analytical <- function(focal_plant, C_conv_plants,
                                                   parms_fixed, times_long,
                                                   invasion_density = 0.01) {
  n_plants <- parms_fixed$n_plants
  n_herbs  <- parms_fixed$n_herbs
  R_vec    <- parms_fixed$R
  C_mat    <- parms_fixed$C
  
  # (a) Intrinsic growth rate: analytical from logistic equation [3]
  fi_00 <- R_vec[focal_plant]
  
  # (b) Run S-1 plant community + all herbivores to equilibrium
  sub_out <- run_subcommunity_plants(focal_plant, parms_fixed, times_long)
  if (is.null(sub_out)) return(NULL)
  
  equil_sub <- as.numeric(sub_out[nrow(sub_out), 2:(n_plants + n_herbs + 1)])
  names(equil_sub) <- c(paste0("P", 1:n_plants), paste0("H", 1:n_herbs))
  
  N_minus_i_star <- equil_sub[paste0("P", 1:n_plants)]
  N_minus_i_star[focal_plant] <- 0
  
  # (c) Invasion growth rate: ODE-based (plant invading S-1 equilibrium)
  state_inv <- equil_sub
  state_inv[focal_plant] <- invasion_density
  
  times_short <- seq(0, 10, by = 0.1)
  out_inv <- tryCatch(
    ode(y = state_inv, times = times_short,
        func = competition_model, parms = parms_fixed, method = "lsoda"),
    error = function(e) NULL
  )
  if (is.null(out_inv)) return(NULL)
  
  df_inv     <- as.data.frame(out_inv)
  p_col      <- paste0("P", focal_plant)
  N_t0_inv   <- df_inv[1, p_col]
  N_tend_inv <- df_inv[nrow(df_inv), p_col]
  t_span     <- df_inv$time[nrow(df_inv)] - df_inv$time[1]
  fi_0Nj     <- log(N_tend_inv / N_t0_inv) / t_span
  
  # (d) No-niche growth rate: analytical from logistic equation [3]
  # Per capita growth at converted density: R[l] - C[l,l] * converted_density
  converted_density <- sum(C_conv_plants[focal_plant, ] * N_minus_i_star)
  fi_cNj <- R_vec[focal_plant] - diag(C_mat)[focal_plant] * converted_density
  
  cat(sprintf("P%d: fi_00 = %.4f, fi_0Nj = %.4f, fi_cNj = %.4f, converted = %.2f, K = %.2f\n",
              focal_plant, fi_00, fi_0Nj, fi_cNj, converted_density,
              R_vec[focal_plant] / diag(C_mat)[focal_plant]))
  
  return(list(
    fi_00  = fi_00,
    fi_0Nj = fi_0Nj,
    fi_cNj = fi_cNj,
    N_minus_i_star    = N_minus_i_star,
    converted_density = converted_density
  ))
}

# ── Numerical conversion factor for plants ────────────────────────────────────
# Solves |1 - N_i(c)| = |1 - N_j(c)| numerically [3]
# For plants, limiting factors include both plant competition AND herbivory
# so c cannot be derived from C matrix alone [1]

compute_plant_conversion_factor <- function(focal_plant, other_plant,
                                            equil_sub, parms_fixed,
                                            times_short = seq(0, 5, by = 0.1),
                                            invasion_density = 0.01) {
  n_plants <- parms_fixed$n_plants
  n_herbs  <- parms_fixed$n_herbs
  R_vec    <- parms_fixed$R
  C_mat    <- parms_fixed$C
  
  # fi(0, N*_j): invasion growth rate of focal into other's subcommunity
  # Already computed upstream; pass in directly
  
  # N_i(c): growth rate of focal plant alone at converted density c * N*_j
  # We sweep c and find where |1 - N_i(c)| = |1 - N_j(c)|
  N_star_j <- equil_sub[paste0("P", other_plant)]
  
  # Per capita growth of focal plant alone at density c * N_star_j
  # using analytical logistic: fi(c * N*_j, 0) = R[i] - C[i,i] * c * N*_j
  # BUT this ignores herbivore effects — use ODE instead
  growth_at_c <- function(c_val) {
    density_c <- c_val * N_star_j
    state_c <- c(P = rep(0, n_plants), H = rep(0, n_herbs))
    state_c[focal_plant] <- max(density_c, invasion_density)
    
    out <- tryCatch(
      ode(y = state_c, times = times_short,
          func = competition_model, parms = parms_fixed, method = "lsoda"),
      error = function(e) NULL
    )
    if (is.null(out)) return(NA)
    df   <- as.data.frame(out)
    p_col <- paste0("P", focal_plant)
    log(df[nrow(df), p_col] / df[1, p_col]) / (tail(times_short, 1))
  }
  
  # Similarly for species j at density (1/c) * N*_i
  N_star_i <- equil_sub[paste0("P", focal_plant)]
  
  growth_j_at_c <- function(c_val) {
    density_c <- (1 / c_val) * N_star_i
    state_c <- c(P = rep(0, n_plants), H = rep(0, n_herbs))
    state_c[other_plant] <- max(density_c, invasion_density)
    
    out <- tryCatch(
      ode(y = state_c, times = times_short,
          func = competition_model, parms = parms_fixed, method = "lsoda"),
      error = function(e) NULL
    )
    if (is.null(out)) return(NA)
    df    <- as.data.frame(out)
    p_col <- paste0("P", other_plant)
    log(df[nrow(df), p_col] / df[1, p_col]) / (tail(times_short, 1))
  }
  
  fi_00_i <- R_vec[focal_plant]
  fi_00_j <- R_vec[other_plant]
  
  # Objective: |1 - N_i(c)| = |1 - N_j(c)|
  # N_i(c) = (fi_inv_i - fi_c_i) / (fi_00_i - fi_c_i)
  objective <- function(log_c) {
    c_val   <- exp(log_c)
    fi_c_i  <- growth_at_c(c_val)
    fi_c_j  <- growth_j_at_c(c_val)
    if (is.na(fi_c_i) | is.na(fi_c_j)) return(NA)
    # Use monoculture invasion growth rates (fi(0, N*_j) etc.)
    # approximated here as fi_00 for simplicity in the solver
    Ni <- (fi_00_i - fi_c_i) / (fi_00_i - fi_c_i + 1e-10)
    Nj <- (fi_00_j - fi_c_j) / (fi_00_j - fi_c_j + 1e-10)
    abs(1 - Ni) - abs(1 - Nj)
  }
  
  tryCatch(
    exp(uniroot(objective, interval = c(-5, 5))$root),
    error = function(e) 1.0   # fallback to 1 if solver fails
  )
}

# ── Updated compute_growth_rates_plants_analytical ────────────────────────────
# Use the multispecies formula from Appendix B [4]:
# converted_density = sum_j( c_ij * N_-i*_j )
# where c_ij is solved numerically per pair
compute_growth_rates_plants_v2 <- function(focal_plant, parms_fixed,
                                           times_long,
                                           times_short = seq(0, 5, by = 0.1),
                                           invasion_density = 0.01) {
  n_plants <- parms_fixed$n_plants
  n_herbs  <- parms_fixed$n_herbs
  R_vec    <- parms_fixed$R
  C_mat    <- parms_fixed$C
  
  # (a) Intrinsic growth rate
  fi_00 <- R_vec[focal_plant]
  
  # (b) S-1 subcommunity to equilibrium
  sub_out <- run_subcommunity_plants(focal_plant, parms_fixed, times_long)
  if (is.null(sub_out)) return(NULL)
  
  equil_sub <- as.numeric(sub_out[nrow(sub_out), 2:(n_plants + n_herbs + 1)])
  names(equil_sub) <- c(paste0("P", 1:n_plants), paste0("H", 1:n_herbs))
  
  N_minus_i_star <- equil_sub[paste0("P", 1:n_plants)]
  N_minus_i_star[focal_plant] <- 0
  
  # (c) Invasion growth rate
  state_inv <- equil_sub
  state_inv[focal_plant] <- invasion_density
  
  out_inv <- tryCatch(
    ode(y = state_inv, times = times_short,
        func = competition_model, parms = parms_fixed, method = "lsoda"),
    error = function(e) NULL
  )
  if (is.null(out_inv)) return(NULL)
  
  df_inv     <- as.data.frame(out_inv)
  p_col      <- paste0("P", focal_plant)
  fi_0Nj     <- log(df_inv[nrow(df_inv), p_col] / df_inv[1, p_col]) /
    (df_inv$time[nrow(df_inv)] - df_inv$time[1])
  
  # (d) Pairwise conversion factors solved numerically [3][4]
  other_plants <- setdiff(1:n_plants, focal_plant)
  c_ij <- rep(1, n_plants)
  
  for (j in other_plants) {
    c_ij[j] <- compute_plant_conversion_factor(
      focal_plant, j, equil_sub, parms_fixed, times_short, invasion_density
    )
  }
  
  # (e) Converted density: sum_j( c_ij * N_-i*_j ) [4]
  converted_density <- sum(c_ij * N_minus_i_star)
  
  # (f) No-niche growth rate: focal plant alone at converted density
  state_nn <- c(P = rep(0, n_plants), H = rep(0, n_herbs))
  state_nn[focal_plant] <- max(converted_density, invasion_density)
  
  out_nn <- tryCatch(
    ode(y = state_nn, times = times_short,
        func = competition_model, parms = parms_fixed, method = "lsoda"),
    error = function(e) NULL
  )
  if (is.null(out_nn)) return(NULL)
  
  df_nn  <- as.data.frame(out_nn)
  fi_cNj <- log(df_nn[nrow(df_nn), p_col] / df_nn[1, p_col]) /
    (df_nn$time[nrow(df_nn)] - df_nn$time[1])
  
  cat(sprintf("P%d: fi_00=%.4f fi_0Nj=%.4f fi_cNj=%.4f converted=%.2f K=%.2f\n",
              focal_plant, fi_00, fi_0Nj, fi_cNj, converted_density,
              R_vec[focal_plant] / diag(C_mat)[focal_plant]))
  
  return(list(
    fi_00 = fi_00, fi_0Nj = fi_0Nj, fi_cNj = fi_cNj,
    N_minus_i_star = N_minus_i_star, converted_density = converted_density,
    c_ij = c_ij
  ))
}

# ── Update run_scenario_NF_plants to use v2 ───────────────────────────────────
run_scenario_NF_plants <- function(scenario) {
  parms   <- parms_base
  parms$R <- scenario$R_vec
  
  results_NF_plants <- map_dfr(1:n_plants, function(i) {
    gr <- compute_growth_rates_plants_v2(
      focal_plant = i, parms_fixed = parms, times_long = times_long
    )
    if (is.null(gr)) return(tibble(species = paste0("P", i),
                                   N = NA, F_val = NA, persists = NA,
                                   fi_00 = NA, fi_0Nj = NA, fi_cNj = NA))
    nf <- compute_N_and_F_2021(gr)
    cx <- check_coexistence_2021(nf$N, nf$F)
    tibble(
      species   = paste0("P", i),
      toughness = toughness[i],
      fi_00     = gr$fi_00,
      fi_0Nj    = gr$fi_0Nj,
      fi_cNj    = gr$fi_cNj,
      N         = nf$N,
      F_val     = nf$F,
      persists  = cx$persists
    )
  })
  
  results_NF_plants$scenario <- scenario$name
  results_NF_plants$type     <- "Plant"
  return(results_NF_plants)
}

# Run for all scenarios
plant_NF_all <- map_dfr(scenarios, run_scenario_NF_plants) %>%
  mutate(scenario = factor(scenario, levels = sapply(scenarios, `[[`, "name")))

# ── Combine plant and herbivore N-F results ───────────────────────────────────
herb_NF_labelled <- all_results %>%
  mutate(species = paste0("H", herb), type = "Herbivore") %>%
  select(species, N, F_val, persists, scenario, type)

plant_NF_labelled <- plant_NF_all %>%
  select(species, N, F_val, persists, scenario, type)

combined_NF <- bind_rows(herb_NF_labelled, plant_NF_labelled) %>%
  mutate(scenario = factor(scenario, levels = sapply(scenarios, `[[`, "name")))

# ── Plot: combined N-F map with plants and herbivores ─────────────────────────
ggplot(combined_NF,
       aes(x = N, y = F_val, color = type, shape = persists, label = species)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(aes(y = F_val + 0.005), size = 3, show.legend = FALSE) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Plant" = "#35b779", "Herbivore" = "#d44842")) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 4),
                     labels = c("Excluded", "Persists")) +
  facet_wrap(~scenario, ncol = 3) +
  coord_cartesian(ylim = c(-1.15, 0.15)) +
  labs(x     = "Niche difference (N)",
       y     = "Fitness difference (F)",
       color = "Trophic level",
       shape = "Outcome",
       title = "Plant and herbivore coexistence across LDMC gradient [3]") +
  theme_bw(base_size = 13)

# Utilization functions across three plant community scenarios ####
# Following Sakarchi & Germain (2025) Eq. 2 [6]:
# Ui = sum_k( a_ik * sqrt(wk * Kk / rk) )
# where:
#   a_ik  = A[k, i]   attack rate of herbivore i on plant k  [1]
#   wk    = E[k, i]   nutritive weight / conversion efficiency
#   Kk    = R[k] / diag(C)[k]  plant carrying capacity
#   rk    = R[k]      plant intrinsic growth rate
#   Kk/rk = resource exploitability [6]

compute_utilization_scenarios <- function(A, E, C_mat, scenarios) {
  
  n_plants <- nrow(A)
  n_herbs  <- ncol(A)
  c_diag   <- diag(C_mat)
  
  results <- map_dfr(names(scenarios), function(s) {
    scenario <- scenarios[[s]]
    R_vec    <- scenario$R_vec
    
    # Plant carrying capacities and exploitability
    K              <- R_vec / c_diag          # Kk = rk / alpha_kk [1]
    exploitability <- K / R_vec               # Kk/rk = 1/alpha_kk [6]
    # Note: since K = R/c_diag, Kk/rk simplifies to 1/c_diag here,
    # but will differ across scenarios when R_vec varies
    
    # Exploitability-weighted attack matrix
    # U_ik = a_ik * sqrt(E[k,i] * Kk/rk)  [6]
    U_mat <- matrix(0, nrow = n_plants, ncol = n_herbs)
    for (i in 1:n_plants) {
      for (j in 1:n_herbs) {
        U_mat[i, j] <- A[i, j] * sqrt(E[i, j] * exploitability[i])
      }
    }
    
    # Total utilization per herbivore: Ui = sum_k(U_ik)  [6]
    U_total <- colSums(U_mat)
    
    # Pairwise niche overlap from normalized utilization vectors
    # cosine similarity = dot product of unit-normalized columns [6]
    U_norm    <- sweep(U_mat, 2, sqrt(colSums(U_mat^2)), "/")
    overlap   <- t(U_norm) %*% U_norm
    
    # Long format for U_mat
    U_long <- as_tibble(U_mat) %>%
      setNames(paste0("H", 1:n_herbs)) %>%
      mutate(plant        = paste0("P", 1:n_plants),
             toughness    = toughness,
             exploitability = exploitability,
             R_val        = R_vec) %>%
      pivot_longer(starts_with("H"),
                   names_to  = "herbivore",
                   values_to = "utilization") %>%
      mutate(scenario     = scenario$name,
             herbivore_id = as.integer(str_extract(herbivore, "\\d+")),
             mandible_str = mandible[herbivore_id])
    
    # Summary per herbivore
    U_summary <- tibble(
      scenario     = scenario$name,
      herbivore    = paste0("H", 1:n_herbs),
      mandible_str = mandible,
      U_total      = U_total,
      cwm_ldmc     = sum(R_vec * toughness) / sum(R_vec)
    )
    
    list(U_long = U_long, U_summary = U_summary,
         U_mat = U_mat, overlap = overlap,
         scenario = scenario$name)
  })
}

# ── Run for all three scenarios ───────────────────────────────────────────────
C_mat  <- diag(0.04, n_plants)
E_mat  <- matrix(0.1, n_plants, n_herbs)

util_scenarios <- map(names(scenarios), function(s) {
  scenario <- scenarios[[s]]
  R_vec    <- scenario$R_vec
  c_diag_v <- diag(C_mat)
  K        <- R_vec / c_diag_v
  exploit  <- K / R_vec
  
  U_mat <- A_fixed * sweep(sqrt(E_mat), 1,
                           sqrt(exploit), FUN = "*")  # [n_plants x n_herbs]
  
  U_total <- colSums(U_mat)
  
  U_norm  <- sweep(U_mat, 2, sqrt(colSums(U_mat^2)), "/")
  overlap <- t(U_norm) %*% U_norm
  
  list(name      = scenario$name,
       R_vec     = R_vec,
       exploit   = exploit,
       U_mat     = U_mat,
       U_total   = U_total,
       overlap   = overlap)
})
names(util_scenarios) <- names(scenarios)

# ── Combined long format for plotting ─────────────────────────────────────────
U_long_all <- map_dfr(names(util_scenarios), function(s) {
  u <- util_scenarios[[s]]
  as_tibble(u$U_mat) %>%
    setNames(paste0("H", 1:n_herbs)) %>%
    mutate(plant        = paste0("P", 1:n_plants),
           toughness    = toughness,
           exploitability = u$exploit) %>%
    pivot_longer(starts_with("H"),
                 names_to  = "herbivore",
                 values_to = "utilization") %>%
    mutate(scenario     = u$name,
           herbivore_id = as.integer(str_extract(herbivore, "\\d+")),
           mandible_str = mandible[herbivore_id],
           scenario     = factor(scenario,
                                 levels = sapply(scenarios, `[[`, "name")))
})

U_summary_all <- map_dfr(names(util_scenarios), function(s) {
  u <- util_scenarios[[s]]
  tibble(
    scenario     = factor(u$name, levels = sapply(scenarios, `[[`, "name")),
    herbivore    = paste0("H", 1:n_herbs),
    mandible_str = mandible,
    U_total      = u$U_total
  )
})

# ── Plot 1: Utilization profiles across plant LDMC spectrum ──────────────────
# Analogous to Figure 1C/D in Sakarchi & Germain [6]
p_util_profile <- ggplot(U_long_all,
                         aes(x = toughness, y = utilization,
                             color = herbivore, group = herbivore)) +
  geom_line(linewidth = 0.8) +
  scale_color_viridis_d(option = "inferno", begin = 0.1, end = 0.9) +
  facet_wrap(~scenario, ncol = 3) +
  labs(x     = "Plant LDMC (toughness)",
       y     = "Utilization (Uik)",
       color = "Herbivore",
       title = "Resource utilization profiles across LDMC gradient") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

print(p_util_profile)

# ── Plot 2: Total utilization vs mandible strength across scenarios ───────────
scenario_colors <- setNames(
  sapply(scenarios, `[[`, "color"),
  sapply(scenarios, `[[`, "name")
)

p_util_total <- ggplot(U_summary_all,
                       aes(x = mandible_str, y = U_total,
                           color = scenario, group = scenario)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = scenario_colors) +
  labs(x     = "Mandible strength",
       y     = "Total utilization (Ui)",
       color = "Scenario",
       title = "Total herbivore utilization across plant community scenarios") +
  theme_bw(base_size = 13)

print(p_util_total)

# ── Plot 3: Resource exploitability across scenarios ──────────────────────────
# Shows how Kk/rk shifts across scenarios — the key driver of utilization [6]
exploit_df <- map_dfr(names(util_scenarios), function(s) {
  u <- util_scenarios[[s]]
  tibble(
    scenario       = factor(u$name, levels = sapply(scenarios, `[[`, "name")),
    toughness      = toughness,
    exploitability = u$exploit
  )
})

p_exploit <- ggplot(exploit_df,
                    aes(x = toughness, y = exploitability,
                        color = scenario, group = scenario)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_color_manual(values = scenario_colors) +
  labs(x     = "Plant LDMC (toughness)",
       y     = "Resource exploitability (K/r)",
       color = "Scenario",
       title = "Plant resource exploitability across LDMC scenarios") +
  theme_bw(base_size = 13)

print(p_exploit)

# ── Plot 4: Pairwise niche overlap heatmaps faceted by scenario ───────────────
overlap_long_all <- map_dfr(names(util_scenarios), function(s) {
  u <- util_scenarios[[s]]
  as_tibble(u$overlap) %>%
    setNames(paste0("H", 1:n_herbs)) %>%
    mutate(herb_i   = paste0("H", 1:n_herbs),
           scenario = factor(u$name,
                             levels = sapply(scenarios, `[[`, "name"))) %>%
    pivot_longer(-c(herb_i, scenario),
                 names_to  = "herb_j",
                 values_to = "overlap")
})

p_overlap <- ggplot(overlap_long_all,
                    aes(x = herb_i, y = herb_j, fill = overlap)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", limits = c(0, 1)) +
  facet_wrap(~scenario, ncol = 3) +
  labs(x    = "Herbivore i",
       y    = "Herbivore j",
       fill = "Niche\noverlap",
       title = "Pairwise utilization overlap across scenarios") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_overlap)

# ── Combined figure ───────────────────────────────────────────────────────────
combined_util <- (p_exploit / p_util_total) | p_overlap +
  plot_layout(widths = c(1, 1.5))

ggsave("utilization_scenarios.tiff", combined_util,
       width = 14, height = 10, dpi = 300, compression = "lzw")


# Community Utilization ######
# ── Community Utilization Q from Sakarchi & Germain (2025) Eq. 3 [6] ─────────
# Q = U + B
# U = unutilized productivity = sum_k[ wk * Kk/rk * (rk - sum_j(ajk * Xj))^2 ]
#   = sum_k[ Kk/rk * (rk - sum_j(A[k,j] * E[k,j] * Xj))^2 ]   (at equilibrium)
# B = metabolic maintenance = sum_j( T_j * X_j ) = sum_j( M_j * H_j )
#
# Note: T_i = M_i in your model [1] (per capita resource requirement at low density)
# g only affects equilibrium density, not T [6]
#
# At consumer equilibrium, R_k = K_k - sum_j(A[k,j] * X_j * K_k/r_k)
# so unutilized productivity is based on how well consumers track resource supply

compute_community_utilization <- function(A, E, R_vec, C_mat, M_vec,
                                          H_eq, P_eq) {
  # A:     [n_plants x n_herbs]   attack rate matrix [1]
  # E:     [n_plants x n_herbs]   conversion efficiency
  # R_vec: vector [n_plants]      plant intrinsic growth rates
  # C_mat: matrix [n_plants x n_plants] plant competition
  # M_vec: vector [n_herbs]       herbivore mortality = T_i [6]
  # H_eq:  vector [n_herbs]       herbivore equilibrium densities
  # P_eq:  vector [n_plants]      plant equilibrium densities
  
  n_plants <- length(R_vec)
  n_herbs  <- length(M_vec)
  K        <- R_vec / diag(C_mat)   # plant carrying capacities
  
  # ── Resource exploitability: K_k / r_k [6] ────────────────────────────────
  exploitability <- K / R_vec       # = 1 / diag(C_mat)
  
  # ── Per-plant utilization by all herbivores: sum_j( A[k,j] * H_j ) ─────────
  # This is the total consumption rate per plant species at equilibrium [6]
  total_consumption_per_plant <- A %*% H_eq   # [n_plants x 1]
  
  # ── U: unutilized productivity ────────────────────────────────────────────
  # U = sum_k[ Kk/rk * (rk - sum_j(ajk * Xj))^2 / rk ]
  # Simplified: at equilibrium plant production = consumption
  # so U captures how much plant productivity escapes herbivory [6]
  residual_production <- R_vec - as.vector(total_consumption_per_plant)
  U <- sum(exploitability * residual_production^2 / R_vec)
  
  # ── B: metabolic maintenance cost ─────────────────────────────────────────
  # B = sum_j( T_j * X_j ) = sum_j( M_j * H_j )  [6]
  B <- sum(M_vec * H_eq)
  
  # ── Q: total community inefficiency ───────────────────────────────────────
  Q <- U + B
  
  # ── Per-plant utilization functions Uik ───────────────────────────────────
  # Uik = A[k,i] * sqrt(E[k,i] * Kk/rk)  [6]
  U_mat <- A * sweep(sqrt(E), 1, sqrt(exploitability), FUN = "*")
  
  # ── Total utilization per herbivore Ui = sum_k(Uik) ───────────────────────
  U_total <- colSums(U_mat)
  
  # ── Summed utilization across all herbivores per plant ────────────────────
  # sum_j(Uj_k) = how "packed" each plant resource is [6]
  U_sum_per_plant <- rowSums(U_mat)
  
  return(list(
    Q               = Q,
    U               = U,
    B               = B,
    U_total         = U_total,
    U_mat           = U_mat,
    U_sum_per_plant = U_sum_per_plant,
    exploitability  = exploitability,
    residual_production = residual_production
  ))
}

# ── Run ODE to equilibrium and extract herbivore/plant densities ──────────────
get_equilibrium <- function(scenario, parms_base, n_plants, n_herbs,
                            times_eq = seq(0, 1000, by = 1)) {
  parms      <- parms_base
  parms$R    <- scenario$R_vec
  
  state_init <- c(P = rep(15, n_plants), H = rep(3, n_herbs))
  
  out <- ode(y = state_init, times = times_eq,
             func = competition_model, parms = parms,
             method = "lsoda")
  
  out_df     <- as.data.frame(out)
  final      <- out_df[nrow(out_df), ]
  H_eq       <- as.numeric(final[paste0("H", 1:n_herbs)])
  P_eq       <- as.numeric(final[paste0("P", 1:n_plants)])
  
  return(list(H_eq = H_eq, P_eq = P_eq))
}

# ── Fixed parameters ──────────────────────────────────────────────────────────
C_mat  <- diag(0.04, n_plants)
E_mat  <- matrix(0.1,  n_plants, n_herbs)
M_vec  <- rep(1,       n_herbs)

# ── Compute Q for all three scenarios ─────────────────────────────────────────
Q_results <- data.frame()

for (s in names(scenarios)) {
  scenario  <- scenarios[[s]]
  eq        <- get_equilibrium(scenario, parms_base, n_plants, n_herbs)
  util      <- compute_community_utilization(
    A     = A_fixed,
    E     = E_mat,
    R_vec = scenario$R_vec,
    C_mat = C_mat,
    M_vec = M_vec,
    H_eq  = eq$H_eq,
    P_eq  = eq$P_eq
  )
  
  Q_results <- rbind(Q_results, data.frame(
    scenario       = scenario$name,
    Q_total        = util$Q,
    U_unutilized   = util$U,
    B_maintenance  = util$B,
    cwm_ldmc       = sum(scenario$R_vec * toughness) / sum(scenario$R_vec)
  ))
}

print(Q_results)

# ── Plot: Q components across scenarios ───────────────────────────────────────
Q_long <- Q_results %>%
  pivot_longer(cols = c(U_unutilized, B_maintenance),
               names_to  = "component",
               values_to = "value") %>%
  mutate(
    component = recode(component,
                       U_unutilized  = "Unutilized productivity (U)",
                       B_maintenance = "Metabolic maintenance (B)"),
    scenario  = factor(scenario, levels = sapply(scenarios, `[[`, "name"))
  )

scenario_colors <- setNames(
  sapply(scenarios, `[[`, "color"),
  sapply(scenarios, `[[`, "name")
)

ggplot(Q_long, aes(x = scenario, y = value, fill = scenario)) +
  geom_col(color = "white", width = 0.6) +
  facet_wrap(~component, scales = "free_y") +
  scale_fill_manual(values = scenario_colors) +
  labs(x = "Scenario",
       y = "Community utilization component",
       title = "Community utilization Q across plant community scenarios [6]") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))
