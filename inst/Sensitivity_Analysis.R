#Sensitivity analysis code
params<-readRDS("/Users/jessicawestworth/Desktop/Status Quo Starting Models/final7.rds")

library(multisensi)
library(sensitivity)

w(params)
# Example: create a matrix of min/max values for mulsisensi
param_ranges <- data.frame(
    Eiw_Blue_whiting = c(0.9*species_params(params)[1,]$Eiw, 1.10*species_params(params)[1,]$Eiw),
    mu_mat_Blue_whiting = c(0.9*species_params(params)[1,]$mu_mat, 1.10*species_params(params)[1,]$mu_mat),
    d_over_g_Blue_whiting = c(0.9*species_params(params)[1,]$d_over_g, 1.10*species_params(params)[1,]$d_over_g),
    l_mat_Blue_whiting = c(0.9*species_params(params)[1,]$l_mat, 1.10*species_params(params)[1,]$l_mat),
    m_Blue_whiting = c(max(0.7,0.9*species_params(params)[1,]$m), min(1.5,1.10*species_params(params)[1,]$m)),
    w_max_Blue_whiting = c(0.9*species_params(params)[1,]$w_max, min(max(w(params)),1.10*species_params(params)[1,]$w_max)),

    Eiw_Boarfish = c(0.9*species_params(params)[2,]$Eiw, 1.10*species_params(params)[2,]$Eiw),
    mu_mat_Boarfish = c(0.9*species_params(params)[2,]$mu_mat, 1.10*species_params(params)[2,]$mu_mat),
    d_over_g_Boarfish = c(0.9*species_params(params)[2,]$d_over_g, 1.10*species_params(params)[2,]$d_over_g),
    l_mat_Boarfish = c(0.9*species_params(params)[2,]$l_mat, 1.10*species_params(params)[2,]$l_mat),
    m_Boarfish = c(max(0.7,0.9*species_params(params)[2,]$m), min(1.5,1.10*species_params(params)[2,]$m)),
    w_max_Boarfish = c(0.9*species_params(params)[2,]$w_max, min(max(w(params)),1.10*species_params(params)[2,]$w_max)),

    Eiw_Cod = c(0.9*species_params(params)[3,]$Eiw, 1.10*species_params(params)[3,]$Eiw),
    mu_mat_Cod = c(0.9*species_params(params)[3,]$mu_mat, 1.10*species_params(params)[3,]$mu_mat),
    d_over_g_Cod = c(0.9*species_params(params)[3,]$d_over_g, 1.10*species_params(params)[3,]$d_over_g),
    l_mat_Cod = c(0.9*species_params(params)[3,]$l_mat, 1.10*species_params(params)[3,]$l_mat),
    m_Cod = c(max(0.7,0.9*species_params(params)[3,]$m), min(1.5,1.10*species_params(params)[3,]$m)),
    w_max_Cod = c(0.9*species_params(params)[3,]$w_max, min(max(w(params)),1.10*species_params(params)[3,]$w_max)),

    Eiw_Haddock = c(0.9*species_params(params)[4,]$Eiw, 1.10*species_params(params)[4,]$Eiw),
    mu_mat_Haddock = c(0.9*species_params(params)[4,]$mu_mat, 1.10*species_params(params)[4,]$mu_mat),
    d_over_g_Haddock = c(0.9*species_params(params)[4,]$d_over_g, 1.10*species_params(params)[4,]$d_over_g),
    l_mat_Haddock = c(0.9*species_params(params)[4,]$l_mat, 1.10*species_params(params)[4,]$l_mat),
    m_Haddock = c(max(0.7,0.9*species_params(params)[4,]$m), min(1.5,1.10*species_params(params)[4,]$m)),
    w_max_Haddock = c(0.9*species_params(params)[4,]$w_max, min(max(w(params)),1.10*species_params(params)[4,]$w_max)),

    Eiw_Hake = c(0.9*species_params(params)[5,]$Eiw, 1.10*species_params(params)[5,]$Eiw),
    mu_mat_Hake = c(0.9*species_params(params)[5,]$mu_mat, 1.10*species_params(params)[5,]$mu_mat),
    d_over_g_Hake = c(0.9*species_params(params)[5,]$d_over_g, 1.10*species_params(params)[5,]$d_over_g),
    l_mat_Hake = c(0.9*species_params(params)[5,]$l_mat, 1.10*species_params(params)[5,]$l_mat),
    m_Hake = c(max(0.7,0.9*species_params(params)[5,]$m), min(1.5,1.10*species_params(params)[5,]$m)),
    w_max_Hake = c(0.9*species_params(params)[5,]$w_max, min(max(w(params)),1.10*species_params(params)[5,]$w_max)),

    Eiw_Herring = c(0.9*species_params(params)[6,]$Eiw, 1.10*species_params(params)[6,]$Eiw),
    mu_mat_Herring = c(0.9*species_params(params)[6,]$mu_mat, 1.10*species_params(params)[6,]$mu_mat),
    d_over_g_Herring = c(0.9*species_params(params)[6,]$d_over_g, 1.10*species_params(params)[6,]$d_over_g),
    l_mat_Herring = c(0.9*species_params(params)[6,]$l_mat, 1.10*species_params(params)[6,]$l_mat),
    m_Herring = c(max(0.7,0.9*species_params(params)[6,]$m), min(1.5,1.10*species_params(params)[6,]$m)),
    w_max_Herring = c(0.9*species_params(params)[6,]$w_max, min(max(w(params)),1.10*species_params(params)[6,]$w_max)),

    Eiw_Horse_mackerel = c(0.9*species_params(params)[7,]$Eiw, 1.10*species_params(params)[7,]$Eiw),
    mu_mat_Horse_mackerel = c(0.9*species_params(params)[7,]$mu_mat, 1.10*species_params(params)[7,]$mu_mat),
    d_over_g_Horse_mackerel = c(0.9*species_params(params)[7,]$d_over_g, 1.10*species_params(params)[7,]$d_over_g),
    l_mat_Horse_mackerel = c(0.9*species_params(params)[7,]$l_mat, 1.10*species_params(params)[7,]$l_mat),
    m_Horse_mackerel = c(max(0.7,0.9*species_params(params)[7,]$m), min(1.5,1.10*species_params(params)[7,]$m)),
    w_max_Horse_mackerel = c(0.9*species_params(params)[7,]$w_max, min(max(w(params)),1.10*species_params(params)[7,]$w_max)),

    Eiw_Mackerel = c(0.9*species_params(params)[8,]$Eiw, 1.10*species_params(params)[8,]$Eiw),
    mu_mat_Mackerel = c(0.9*species_params(params)[8,]$mu_mat, 1.10*species_params(params)[8,]$mu_mat),
    d_over_g_Mackerel = c(0.9*species_params(params)[8,]$d_over_g, 1.10*species_params(params)[8,]$d_over_g),
    l_mat_Mackerel = c(0.9*species_params(params)[8,]$l_mat, 1.10*species_params(params)[8,]$l_mat),
    m_Mackerel = c(max(0.7,0.9*species_params(params)[8,]$m), min(1.5,1.10*species_params(params)[8,]$m)),
    w_max_Mackerel = c(0.9*species_params(params)[8,]$w_max, min(max(w(params)),1.10*species_params(params)[8,]$w_max)),

    Eiw_Megrim = c(0.9*species_params(params)[9,]$Eiw, 1.10*species_params(params)[9,]$Eiw),
    mu_mat_Megrim = c(0.9*species_params(params)[9,]$mu_mat, 1.10*species_params(params)[9,]$mu_mat),
    d_over_g_Megrim = c(0.9*species_params(params)[9,]$d_over_g, 1.10*species_params(params)[9,]$d_over_g),
    l_mat_Megrim = c(0.9*species_params(params)[9,]$l_mat, 1.10*species_params(params)[9,]$l_mat),
    m_Megrim = c(max(0.7,0.9*species_params(params)[9,]$m), min(1.5,1.10*species_params(params)[9,]$m)),
    w_max_Megrim = c(0.9*species_params(params)[9,]$w_max, min(max(w(params)),1.10*species_params(params)[9,]$w_max)),

    Eiw_Monkfish = c(0.9*species_params(params)[10,]$Eiw, 1.10*species_params(params)[10,]$Eiw),
    mu_mat_Monkfish = c(0.9*species_params(params)[10,]$mu_mat, 1.10*species_params(params)[10,]$mu_mat),
    d_over_g_Monkfish = c(0.9*species_params(params)[10,]$d_over_g, 1.10*species_params(params)[10,]$d_over_g),
    l_mat_Monkfish = c(0.9*species_params(params)[10,]$l_mat, 1.10*species_params(params)[10,]$l_mat),
    m_Monkfish = c(max(0.7,0.9*species_params(params)[10,]$m), min(1.5,1.10*species_params(params)[10,]$m)),
    w_max_Monkfish = c(0.9*species_params(params)[10,]$w_max, min(max(w(params)),1.10*species_params(params)[10,]$w_max)),

    Eiw_Plaice = c(0.9*species_params(params)[11,]$Eiw, 1.10*species_params(params)[11,]$Eiw),
    mu_mat_Plaice = c(0.9*species_params(params)[11,]$mu_mat, 1.10*species_params(params)[11,]$mu_mat),
    d_over_g_Plaice = c(0.9*species_params(params)[11,]$d_over_g, 1.10*species_params(params)[11,]$d_over_g),
    l_mat_Plaice = c(0.9*species_params(params)[11,]$l_mat, 1.10*species_params(params)[11,]$l_mat),
    m_Plaice = c(max(0.7,0.9*species_params(params)[11,]$m), min(1.5,1.10*species_params(params)[11,]$m)),
    w_max_Plaice = c(0.9*species_params(params)[11,]$w_max, min(max(w(params)),1.10*species_params(params)[11,]$w_max)),

    Eiw_Red_gurnard = c(0.9*species_params(params)[12,]$Eiw, 1.10*species_params(params)[12,]$Eiw),
    mu_mat_Red_gurnard = c(0.9*species_params(params)[12,]$mu_mat, 1.10*species_params(params)[12,]$mu_mat),
    d_over_g_Red_gurnard = c(0.9*species_params(params)[12,]$d_over_g, 1.10*species_params(params)[12,]$d_over_g),
    l_mat_Red_gurnard = c(0.9*species_params(params)[12,]$l_mat, 1.10*species_params(params)[12,]$l_mat),
    m_Red_gurnard = c(max(0.7,0.9*species_params(params)[12,]$m), min(1.5,1.10*species_params(params)[12,]$m)),
    w_max_Red_gurnard = c(0.9*species_params(params)[12,]$w_max, min(max(w(params)),1.10*species_params(params)[12,]$w_max)),

    Eiw_Sole = c(0.9*species_params(params)[13,]$Eiw, 1.10*species_params(params)[13,]$Eiw),
    mu_mat_Sole = c(0.9*species_params(params)[13,]$mu_mat, 1.10*species_params(params)[13,]$mu_mat),
    d_over_g_Sole = c(0.9*species_params(params)[13,]$d_over_g, 1.10*species_params(params)[13,]$d_over_g),
    l_mat_Sole = c(0.9*species_params(params)[13,]$l_mat, 1.10*species_params(params)[13,]$l_mat),
    m_Sole = c(max(0.7,0.9*species_params(params)[13,]$m), min(1.5,1.10*species_params(params)[13,]$m)),
    w_max_Sole = c(0.9*species_params(params)[13,]$w_max, min(max(w(params)),1.10*species_params(params)[13,]$w_max)),

    Eiw_Whiting = c(0.9*species_params(params)[14,]$Eiw, 1.10*species_params(params)[14,]$Eiw),
    mu_mat_Whiting = c(0.9*species_params(params)[14,]$mu_mat, 1.10*species_params(params)[14,]$mu_mat),
    d_over_g_Whiting = c(0.9*species_params(params)[14,]$d_over_g, 1.10*species_params(params)[14,]$d_over_g),
    l_mat_Whiting = c(0.9*species_params(params)[14,]$l_mat, 1.10*species_params(params)[14,]$l_mat),
    m_Whiting = c(max(0.7,0.9*species_params(params)[14,]$m), min(1.5,1.10*species_params(params)[14,]$m)),
    w_max_Whiting = c(0.9*species_params(params)[14,]$w_max, min(max(w(params)),1.10*species_params(params)[14,]$w_max))
)

# w_repro_max commented will mean changing w_max doesn't do anything...
# figure out a fix for this.
# SetBevertonHolt everytime after steady single species.


# Parameter explored evenly across ranges
m <- morris(model = NULL, factors = length(param_ranges),
            r = 10, design = list(type = "oat", levels = 10, grid.jump = 2))

# Rename columns of X to match your parameters
X <-as.data.frame(m$X)
names(X) <- names(param_ranges)
X_norm<-X

X_real <- X_norm
for (j in 1:ncol(X_norm)) {
    min_val <- param_ranges[1, j]
    max_val <- param_ranges[2, j]
    X_real[, j] <- min_val + (X_norm[, j] * (max_val - min_val))
}

input_design<-X_real

sensitivity_sim_list<-list()
final_outputs <- data.frame(total_yield=NA,total_biomass=NA,total_SSB=NA, total_N=NA)

for (i in 1:nrow(input_design)) {
    params_i <- params  # copy original model

    # Update species-specific parameters for this run
    species_params(params_i )[1,]$Eiw<-input_design[i, "Eiw_Blue_whiting"]
    species_params(params_i)[1,]$mu_mat<-input_design[i, "mu_mat_Blue_whiting"]
    species_params(params_i)[1,]$d_over_g<-input_design[i, "d_over_g_Blue_whiting"]
    species_params(params_i)[1,]$l_mat<-input_design[i, "l_mat_Blue_whiting"]
    species_params(params_i)[1,]$m<-input_design[i, "m_Blue_whiting"]
    species_params(params_i)[1,]$w_max<-input_design[i, "w_max_Blue_whiting"]

    species_params(params_i)[2,]$Eiw<-input_design[i, "Eiw_Boarfish"]
    species_params(params_i)[2,]$mu_mat<-input_design[i, "mu_mat_Boarfish"]
    species_params(params_i)[2,]$d_over_g<-input_design[i, "d_over_g_Boarfish"]
    species_params(params_i)[2,]$l_mat<-input_design[i, "l_mat_Boarfish"]
    species_params(params_i)[2,]$m<-input_design[i, "m_Boarfish"]
    species_params(params_i)[2,]$w_max<-input_design[i, "w_max_Boarfish"]

    species_params(params_i)[3,]$Eiw<-input_design[i, "Eiw_Cod"]
    species_params(params_i)[3,]$mu_mat<-input_design[i, "mu_mat_Cod"]
    species_params(params_i)[3,]$d_over_g<-input_design[i, "d_over_g_Cod"]
    species_params(params_i)[3,]$l_mat<-input_design[i, "l_mat_Cod"]
    species_params(params_i)[3,]$m<-input_design[i, "m_Cod"]
    species_params(params_i)[3,]$w_max<-input_design[i, "w_max_Cod"]

    species_params(params_i)[4,]$Eiw<-input_design[i, "Eiw_Haddock"]
    species_params(params_i)[4,]$mu_mat<-input_design[i, "mu_mat_Haddock"]
    species_params(params_i)[4,]$d_over_g<-input_design[i, "d_over_g_Haddock"]
    species_params(params_i)[4,]$l_mat<-input_design[i, "l_mat_Haddock"]
    species_params(params_i)[4,]$m<-input_design[i, "m_Haddock"]
    species_params(params_i)[4,]$w_max<-input_design[i, "w_max_Haddock"]

    species_params(params_i)[5,]$Eiw<-input_design[i, "Eiw_Hake"]
    species_params(params_i)[5,]$mu_mat<-input_design[i, "mu_mat_Hake"]
    species_params(params_i)[5,]$d_over_g<-input_design[i, "d_over_g_Hake"]
    species_params(params_i)[5,]$l_mat<-input_design[i, "l_mat_Hake"]
    species_params(params_i)[5,]$m<-input_design[i, "m_Hake"]
    species_params(params_i)[5,]$w_max<-input_design[i, "w_max_Hake"]

    species_params(params_i)[6,]$Eiw<-input_design[i, "Eiw_Herring"]
    species_params(params_i)[6,]$mu_mat<-input_design[i, "mu_mat_Herring"]
    species_params(params_i)[6,]$d_over_g<-input_design[i, "d_over_g_Herring"]
    species_params(params_i)[6,]$l_mat<-input_design[i, "l_mat_Herring"]
    species_params(params_i)[6,]$m<-input_design[i, "m_Herring"]
    species_params(params_i)[6,]$w_max<-input_design[i, "w_max_Herring"]

    species_params(params_i)[7,]$Eiw<-input_design[i, "Eiw_Horse_mackerel"]
    species_params(params_i)[7,]$mu_mat<-input_design[i, "mu_mat_Horse_mackerel"]
    species_params(params_i)[7,]$d_over_g<-input_design[i, "d_over_g_Horse_mackerel"]
    species_params(params_i)[7,]$l_mat<-input_design[i, "l_mat_Horse_mackerel"]
    species_params(params_i)[7,]$m<-input_design[i, "m_Horse_mackerel"]
    species_params(params_i)[7,]$w_max<-input_design[i, "w_max_Horse_mackerel"]

    species_params(params_i)[8,]$Eiw<-input_design[i, "Eiw_Mackerel"]
    species_params(params_i)[8,]$mu_mat<-input_design[i, "mu_mat_Mackerel"]
    species_params(params_i)[8,]$d_over_g<-input_design[i, "d_over_g_Mackerel"]
    species_params(params_i)[8,]$l_mat<-input_design[i, "l_mat_Mackerel"]
    species_params(params_i)[8,]$m<-input_design[i, "m_Mackerel"]
    species_params(params_i)[8,]$w_max<-input_design[i, "w_max_Mackerel"]

    species_params(params_i)[9,]$Eiw<-input_design[i, "Eiw_Megrim"]
    species_params(params_i)[9,]$mu_mat<-input_design[i, "mu_mat_Megrim"]
    species_params(params_i)[9,]$d_over_g<-input_design[i, "d_over_g_Megrim"]
    species_params(params_i)[9,]$l_mat<-input_design[i, "l_mat_Megrim"]
    species_params(params_i)[9,]$m<-input_design[i, "m_Megrim"]
    species_params(params_i)[9,]$w_max<-input_design[i, "w_max_Megrim"]

    species_params(params_i)[10,]$Eiw<-input_design[i, "Eiw_Monkfish"]
    species_params(params_i)[10,]$mu_mat<-input_design[i, "mu_mat_Monkfish"]
    species_params(params_i)[10,]$d_over_g<-input_design[i, "d_over_g_Monkfish"]
    species_params(params_i)[10,]$l_mat<-input_design[i, "l_mat_Monkfish"]
    species_params(params_i)[10,]$m<-input_design[i, "m_Monkfish"]
    species_params(params_i)[10,]$w_max<-input_design[i, "w_max_Monkfish"]

    species_params(params_i)[11,]$Eiw<-input_design[i, "Eiw_Plaice"]
    species_params(params_i)[11,]$mu_mat<-input_design[i, "mu_mat_Plaice"]
    species_params(params_i)[11,]$d_over_g<-input_design[i, "d_over_g_Plaice"]
    species_params(params_i)[11,]$l_mat<-input_design[i, "l_mat_Plaice"]
    species_params(params_i)[11,]$m<-input_design[i, "m_Plaice"]
    species_params(params_i)[11,]$w_max<-input_design[i, "w_max_Plaice"]

    species_params(params_i)[12,]$Eiw<-input_design[i, "Eiw_Red_gurnard"]
    species_params(params_i)[12,]$mu_mat<-input_design[i, "mu_mat_Red_gurnard"]
    species_params(params_i)[12,]$d_over_g<-input_design[i, "d_over_g_Red_gurnard"]
    species_params(params_i)[12,]$l_mat<-input_design[i, "l_mat_Red_gurnard"]
    species_params(params_i)[12,]$m<-input_design[i, "m_Red_gurnard"]
    species_params(params_i)[12,]$w_max<-input_design[i, "w_max_Red_gurnard"]

    species_params(params_i)[13,]$Eiw<-input_design[i, "Eiw_Sole"]
    species_params(params_i)[13,]$mu_mat<-input_design[i, "mu_mat_Sole"]
    species_params(params_i)[13,]$d_over_g<-input_design[i, "d_over_g_Sole"]
    species_params(params_i)[13,]$l_mat<-input_design[i, "l_mat_Sole"]
    species_params(params_i)[13,]$m<-input_design[i, "m_Sole"]
    species_params(params_i)[13,]$w_max<-input_design[i, "w_max_Sole"]

    species_params(params_i)[14,]$Eiw<-input_design[i, "Eiw_Whiting"]
    species_params(params_i)[14,]$mu_mat<-input_design[i, "mu_mat_Whiting"]
    species_params(params_i)[14,]$d_over_g<-input_design[i, "d_over_g_Whiting"]
    species_params(params_i)[14,]$l_mat<-input_design[i, "l_mat_Whiting"]
    species_params(params_i)[14,]$m<-input_design[i, "m_Whiting"]
    species_params(params_i)[14,]$w_max<-input_design[i, "w_max_Whiting"]

    params_i<-setParams(params_i)

    #insert code for recalculating Eriw
    ext_enc <- getExtEncounter(params)
    sps <- species_params(params_i)
    w_bins <- w(params_i)
    ext_enc[1, ] <- species_params(params_i)[1,]$Eiw * w_bins^sps$n[1]
    ext_enc[2, ] <- species_params(params_i)[2,]$Eiw * w_bins^sps$n[2]
    ext_enc[3, ] <- species_params(params_i)[3,]$Eiw * w_bins^sps$n[3]
    ext_enc[4, ] <- species_params(params_i)[4,]$Eiw * w_bins^sps$n[4]
    ext_enc[5, ] <- species_params(params_i)[5,]$Eiw * w_bins^sps$n[5]
    ext_enc[6, ] <- species_params(params_i)[6,]$Eiw * w_bins^sps$n[6]
    ext_enc[7, ] <- species_params(params_i)[7,]$Eiw * w_bins^sps$n[7]
    ext_enc[8, ] <- species_params(params_i)[8,]$Eiw * w_bins^sps$n[8]
    ext_enc[9, ] <- species_params(params_i)[9,]$Eiw * w_bins^sps$n[9]
    ext_enc[10, ] <- species_params(params_i)[10,]$Eiw * w_bins^sps$n[10]
    ext_enc[11, ] <- species_params(params_i)[11,]$Eiw * w_bins^sps$n[11]
    ext_enc[12, ] <- species_params(params_i)[12,]$Eiw * w_bins^sps$n[12]
    ext_enc[13, ] <- species_params(params_i)[13,]$Eiw * w_bins^sps$n[13]
    ext_enc[14, ] <- species_params(params_i)[14,]$Eiw * w_bins^sps$n[14]
    params_i <- setExtEncounter(params_i, ext_encounter = ext_enc)
    params_i<-matchBiomasses(params_i)

    #Model set up for projection
    params_i<-mizer::steadySingleSpecies(params_i)
    params_i<-setBevertonHolt(params_i)
    params_i<-setFeedingLevels(params=params_i, f=0.6, f_c=0.2)
    params_i<-mizer::steadySingleSpecies(params_i)

    dm <- dm

    # Set diffusion from d_over_g
    # d(w) = d_over_g * g(w) * w
    for (species in species_params(params_i)$species) {
        d_over_g <- species_params(params_i)[species, "d_over_g"]
        w <- params_i@w
        growth <- getEGrowth(params_i)[species, ]
        n <- params_i@species_params[species, "n"]
        g_0 <- growth[1] / w[1]^n
        d_0 <- d_over_g * g_0
        diffusion(params_i)[species, ] <- d_0 * w^(n + 1)
    }

    params_i<-setParams(params_i)
    params_i<-mizer::steadySingleSpecies(params_i)
    params_i <- validParams(params_i)

    resource_params(params_i)$w_pp_cutoff <- 1

    params_i <- params_i |>
        mizer::steadySingleSpecies() |>
        alignResource() |>
        setResourceInteraction(resource_dynamics = "resource_semichemostat") |>
        matchDiet(dm) |>
        setBevertonHolt(reproduction_level = 0.8)

    #setParams(params)?
    sim_i <- project(params_i, t_max = 500)

    sensitivity_sim_list[[paste0(i)]] <- sim_i

    # Extract final 5 outputs
    final_outputs[i, ] <- c(
        total_yield = sum(getYield(sim_i)[500,]),
        total_biomass = sum(getBiomass(sim_i)[500,]),
        total_SSB = sum(getSSB(sim_i)[500,]),
        total_N = sum(getN(sim_i)[500,])
    )
}
#save models
#saveRDS(sensitivity_sim_list, "/Users/jessicawestworth/Desktop/Sensitivity/sensitivity_sim_list.rds")
#saveRDS(final_outputs, "/Users/jessicawestworth/Desktop/Sensitivity/sensitivity_final_outputs.rds")
final_outputs<-readRDS("/Users/jessicawestworth/Desktop/Sensitivity/sensitivity_final_outputs.rds")
sensitivity_sim_list<-readRDS("/Users/jessicawestworth/Desktop/Sensitivity/sensitivity_sim_list.rds")

#compute results
# Create a copy of your morris object for each output measure
m_Yield <- m
m_Biomass <- m
m_SSB <- m
m_N <- m

# "Tell" the morris object what the results were
# This calculates the Mu and Sigma (elementary effects)
tell(m_Yield, final_outputs$total_yield)
tell(m_Biomass, final_outputs$total_biomass)
tell(m_SSB, final_outputs$total_SSB)
tell(m_N, final_outputs$total_N)

# plot

results_summary_yield <- data.frame(
    parameter = names(param_ranges),
    mu_star = apply(m_Yield$ee, 2, function(x) mean(abs(x))),
    sigma = apply(m_Yield$ee, 2, sd)
    )

ggplot(results_summary_yield, aes(x = mu_star, y = sigma, label = parameter)) +
    geom_point() +
    geom_text(vjust = -0.5) +
    labs(title = "Morris Sensitivity Yield (Benoit et al. Style)",
                x = "Mean Absolute Effect (mu*)",
                y = "Standard Deviation (sigma)") +
    theme_minimal()

results_summary_biomass <- data.frame(
    parameter = names(param_ranges),
    mu_star = apply(m_Biomass$ee, 2, function(x) mean(abs(x))),
    sigma = apply(m_Biomass$ee, 2, sd)
)

ggplot(results_summary_biomass, aes(x = mu_star, y = sigma, label = parameter)) +
    geom_point() +
    geom_text(vjust = -0.5) +
    labs(title = "Morris Sensitivity Biomass (Benoit et al. Style)",
         x = "Mean Absolute Effect (mu*)",
         y = "Standard Deviation (sigma)") +
    theme_minimal()

results_summary_SSB <- data.frame(
    parameter = names(param_ranges),
    mu_star = apply(m_SSB$ee, 2, function(x) mean(abs(x))),
    sigma = apply(m_SSB$ee, 2, sd)
)

ggplot(results_summary_SSB, aes(x = mu_star, y = sigma, label = parameter)) +
    geom_point() +
    geom_text(vjust = -0.5) +
    labs(title = "Morris Sensitivity SSB (Benoit et al. Style)",
         x = "Mean Absolute Effect (mu*)",
         y = "Standard Deviation (sigma)") +
    theme_minimal()

results_summary_N <- data.frame(
    parameter = names(param_ranges),
    mu_star = apply(m_N$ee, 2, function(x) mean(abs(x))),
    sigma = apply(m_N$ee, 2, sd)
)

ggplot(results_summary_N, aes(x = mu_star, y = sigma, label = parameter)) +
    geom_point() +
    geom_text(vjust = -0.5) +
    labs(title = "Morris Sensitivity N (Benoit et al. Style)",
         x = "Mean Absolute Effect (mu*)",
         y = "Standard Deviation (sigma)") +
    theme_minimal()

#check that all simulations ran to steady (including extinctions)
for (i in 1:nrow(input_design)) {
    sim<-sensitivity_sim_list[[i]]
    Biomass<-getBiomass(sim)
    for (s in 1:nrow(species_params(params))) {
        near<-Biomass[498,s]
        end<-Biomass[500,s]
        dif<-end-near
        if(dif>1e-4){
            warning(paste(s,i))
        }
    }
}

#Elementary effects test: tells us which parameters have large effects
#for Y top 5 are m_Horse_mackerel, d/g blue whiting, w_max whiting, mu_mat Blue whiting, w_max Cod, m_Cod
#For B top 5 are m_Horse_mackerel, w_max Whiting, d/g Blue whiting, w_max cod, m_Boarfish
#For SSB top 5 are mu_mat blue whiting, m_horse mackerel, w_max whiting, d/g blue whiting, d/g sole
#For N top 5 are eiw_red_gurnard, l_mat mackerel, d/g sole, mu_mat blue whiting, l_mat hake, m_Cod

#Regional sensitivity analysis

#What we would like to check:
#Is it that BH yield simulations are always higher than status quo simulation values
BH_LFY>status_LFY

#Is it always that biomass of larger fish is higher than status quo simulation values
BH_LFB>status_LFB

#Is it always that yield of larger fish is lower than status quo simulation values
BH_Y>status_Y

#At what c value does BH have lower Biomass than status quo simulations
#Is it always the case that at some c values BH has a higher and at other c values
#BH has a lower B.
BH_B<status_B

#At what c value does BH have lower SSB than status quo simulations
#Is it always the case that at some c values BH has a higher and at other c values
#BH has a lower SSB.
BH_SSB<status_SSB

#At what c value does BH have lower N than status quo simulations
#Is it always the case that at some c values BH has a higher and at other c values
#BH has a lower N.
BH_B<status_B

#latin hypercube sampling
library(lhs)

#number of runs
n_runs <- 100

#model
params<-readRDS("/Users/jessicawestworth/Desktop/Status Quo Starting Models/final7.rds")

#RSA outputs
RSA_sim_list<-list()
RSA_outputs <- expand.grid(
    c = c(0.1),
    alpha = c(1),
    run = seq(n_runs),
    status_LFY_t = NA,
    BH_LFY_t = NA,
    status_LFY = NA,
    BH_LFY = NA,
    status_Y = NA,
    BH_Y = NA,
    status_LFB_t = NA,
    BH_LFB_t = NA,
    status_LFB = NA,
    BH_LFB = NA,
    status_B = NA,
    BH_B = NA,
    staus_SSB = NA,
    BH_SSB = NA,
    status_N = NA,
    BH_N = NA)


#LHS design for the general top 9 of the top 5 Elementary effects
factors <- c("l_mat_Mackerel","Eiw_Red_gurnard","m_Horse_mackerel", "d_over_g_Sole","d_over_g_Blue_whiting", "w_max_Whiting", "mu_mat_Blue_whiting", "m_Cod", "w_max_Cod")
set.seed(123)
lhs_design <- randomLHS(n_runs, length(factors))

#scale LHS (0 to 1) from the 10% range
sim_params_prop <- as.data.frame(lhs_design)
colnames(sim_params_prop) <- factors

# Replace proportions with your actual baseline values
RSA_param_ranges<-param_ranges[,factors]
sim_params<-sim_params_prop
for (j in 1:ncol(sim_params_prop)) {
    min_val <- RSA_param_ranges[1, j]
    max_val <- RSA_param_ranges[2, j]
    sim_params[, j] <- min_val + (sim_params_prop[, j] * (max_val - min_val))
}

for(i in 1:n_runs) {
    params_i <- params

    # update species-specific params with the specified values from sim_params
    #(monte-carlo values) for this run
    species_params(params_i)[1,]$mu_mat<-sim_params[i, "mu_mat_Blue_whiting"]
    species_params(params_i)[1,]$d_over_g<-sim_params[i, "d_over_g_Blue_whiting"]
    species_params(params_i)[3,]$m<-sim_params[i, "m_Cod"]
    species_params(params_i)[3,]$w_max<-sim_params[i, "w_max_Cod"]
    species_params(params_i)[7,]$m<-sim_params[i, "m_Horse_mackerel"]
    species_params(params_i)[8,]$l_mat<-sim_params[i, "l_mat_Mackerel"]
    species_params(params_i)[12,]$Eiw<-sim_params[i, "Eiw_Red_gurnard"]
    species_params(params_i)[13,]$d_over_g<-sim_params[i, "d_over_g_Sole"]
    species_params(params_i)[14,]$w_max<-sim_params[i, "w_max_Whiting"]

    params_i<-setParams(params_i)

    #insert code for recalculating Eriw
    ext_enc <- getExtEncounter(params)
    sps <- species_params(params_i)
    w_bins <- w(params_i)
    ext_enc[1, ] <- species_params(params_i)[1,]$Eiw * w_bins^sps$n[1]
    ext_enc[2, ] <- species_params(params_i)[2,]$Eiw * w_bins^sps$n[2]
    ext_enc[3, ] <- species_params(params_i)[3,]$Eiw * w_bins^sps$n[3]
    ext_enc[4, ] <- species_params(params_i)[4,]$Eiw * w_bins^sps$n[4]
    ext_enc[5, ] <- species_params(params_i)[5,]$Eiw * w_bins^sps$n[5]
    ext_enc[6, ] <- species_params(params_i)[6,]$Eiw * w_bins^sps$n[6]
    ext_enc[7, ] <- species_params(params_i)[7,]$Eiw * w_bins^sps$n[7]
    ext_enc[8, ] <- species_params(params_i)[8,]$Eiw * w_bins^sps$n[8]
    ext_enc[9, ] <- species_params(params_i)[9,]$Eiw * w_bins^sps$n[9]
    ext_enc[10, ] <- species_params(params_i)[10,]$Eiw * w_bins^sps$n[10]
    ext_enc[11, ] <- species_params(params_i)[11,]$Eiw * w_bins^sps$n[11]
    ext_enc[12, ] <- species_params(params_i)[12,]$Eiw * w_bins^sps$n[12]
    ext_enc[13, ] <- species_params(params_i)[13,]$Eiw * w_bins^sps$n[13]
    ext_enc[14, ] <- species_params(params_i)[14,]$Eiw * w_bins^sps$n[14]
    params_i <- setExtEncounter(params_i, ext_encounter = ext_enc)
    params_i<-matchBiomasses(params_i)

    #Model set up for projection
    params_i<-mizer::steadySingleSpecies(params_i)
    params_i<-setBevertonHolt(params_i)
    params_i<-setFeedingLevels(params=params_i, f=0.6, f_c=0.2)
    params_i<-mizer::steadySingleSpecies(params_i)

    dm <- dm

    # Set diffusion from d_over_g
    # d(w) = d_over_g * g(w) * w
    for (species in species_params(params_i)$species) {
        d_over_g <- species_params(params_i)[species, "d_over_g"]
        w <- params_i@w
        growth <- getEGrowth(params_i)[species, ]
        n <- params_i@species_params[species, "n"]
        g_0 <- growth[1] / w[1]^n
        d_0 <- d_over_g * g_0
        diffusion(params_i)[species, ] <- d_0 * w^(n + 1)
    }

    params_i<-setParams(params_i)
    params_i<-mizer::steadySingleSpecies(params_i)
    params_i <- validParams(params_i)

    resource_params(params_i)$w_pp_cutoff <- 1

    params_i <- params_i |>
        mizer::steadySingleSpecies() |>
        alignResource() |>
        setResourceInteraction(resource_dynamics = "resource_semichemostat") |>
        matchDiet(dm) |>
        setBevertonHolt(reproduction_level = 0.8)

    params_i <- setParams(params_i)

    #Simulation
    #run status quo to steady for 500 years, blend for 100 years, run for
    #another 500 years with BH values
    t_duration <- 1100
    t_steadied <- 500
    t_blended <- 600

    for(c_val in unique(RSA_outputs$c)){
        for(alpha_val in unique(RSA_outputs$alpha)){
            my_blender <- make_blended_ssBH_FMort(
            t_max_blend = t_blended,
            target_c = c_val,
            t_steady = t_steadied,
            alpha_max = alpha_val
        )

        sim_BH_start <- setRateFunction(params_i, "FMort", "my_blender")
        sim_BH <- project(sim_BH_start, t_max = t_duration, effort = 1)
        #save in case crash
        RSA_sim_list[[paste0("run_",i,"_c_",c_val,"_a_",alpha_val)]] <- sim_BH

        #Size-species specific Fishing Mortality, Biomass to get eventual
        #Biomass, Yield, LFY, and LFB calculations
        f<-getFMort(sim_BH, drop = FALSE)
        biomass <- biomass <- sweep(sim_BH@n, 3, sim_BH@params@w * sim_BH@params@dw, "*")
        yield<-f * biomass

        status_LFY_s<-0
        status_Y_s<-0
        status_LFY<-0
        status_Y<-0

        BH_LFY_s<-0
        BH_Y_s<-0
        BH_LFY<-0
        BH_Y<-0

        status_LFB_s<-0
        status_B_s<-0
        status_LFB<-0
        status_B<-0

        BH_LFB_s<-0
        BH_B_s<-0
        BH_LFB<-0
        BH_B<-0

        for(s in 1:nrow(sim_BH@params@species_params)){
            l<-50
            a<-sim_BH@params@species_params$a[s]
            b<-sim_BH@params@species_params$b[s]
            w_cut<-a*(l^b)

            #yields
            L_yield<-yield[,s,w>w_cut]

            status_LFY_s<-sum(L_yield[499,])
            status_Y_s<-sum(yield[499,s,])
            status_LFY<-status_LFY+status_LFY_s
            status_Y<-status_Y+status_Y_s

            BH_LFY_s<-sum(L_yield[1099,])
            BH_Y_s<-sum(yield[1099,s,])
            BH_LFY<-BH_LFY+BH_LFY_s
            BH_Y<-BH_Y+BH_Y_s

            #Biomass
            L_biomass<-biomass[,s,w>w_cut]

            status_LFB_s<-sum(L_biomass[499,])
            status_B_s<-sum(biomass[499,s,])
            status_LFB<-status_LFB+status_LFB_s
            status_B<-status_B+status_B_s

            BH_LFB_s<-sum(L_biomass[1099,])
            BH_B_s<-sum(biomass[1099,s,])
            BH_LFB<-BH_LFB+BH_LFB_s
            BH_B<-BH_B+BH_B_s
        }

        status_LFY_t<-status_LFY
        status_LFY<-status_LFY/status_Y
        BH_LFY_t<-BH_LFY
        BH_LFY<-BH_LFY/BH_Y

        status_LFB_t<-status_LFB
        status_LFB<-status_LFB/status_B
        BH_LFB_t<-BH_LFB
        BH_LFB<-BH_LFB/BH_B

        #insert values into RSA_outputs
        rows <- RSA_outputs$c == c_val & RSA_outputs$alpha == alpha_val & RSA_outputs$run ==i

        #large fish proportion in yield
        RSA_outputs$status_LFY_t[rows]<-status_LFY_t
        RSA_outputs$BH_LFY_t[rows]<-BH_LFY_t
        RSA_outputs$status_LFY[rows]<-status_LFY
        RSA_outputs$BH_LFY[rows]<-BH_LFY

        #Yield
        RSA_outputs$status_Y[rows]<-status_Y
        RSA_outputs$BH_Y[rows]<-BH_Y

        #large fish proportion in system biomass
        RSA_outputs$status_LFB_t[rows]<-status_LFB_t
        RSA_outputs$BH_LFB_t[rows]<-BH_LFB_t
        RSA_outputs$status_LFB[rows]<-status_LFB
        RSA_outputs$BH_LFB[rows]<-BH_LFB

        #Biomass
        RSA_outputs$status_B[rows]<-status_B
        RSA_outputs$BH_B[rows]<-BH_B

        RSA_outputs$staus_SSB[rows] <- sum(getSSB(sim_BH)[499,])
        RSA_outputs$BH_SSB[rows] <- sum(getSSB(sim_BH)[1099,])
        RSA_outputs$status_N[rows] <- sum(getN(sim_BH)[499,])
        RSA_outputs$BH_N[rows] <- sum(getN(sim_BH)[1099,])
        }
    }
}

saveRDS(RSA_outputs, "/Users/jessicawestworth/Desktop/Sensitivity/RSA_outputs.rds")
saveRDS(RSA_sim_list, "/Users/jessicawestworth/Desktop/Sensitivity/RSA_sim_list.rds")


#library(ggplot2)

ggplot(results_comparison, aes(x = yield_sq, y = yield_bh)) +
    geom_point(aes(color = bh_wins), alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = "Is Balanced Harvest Yield Robustly Higher?",
         x = "Status Quo Yield", y = "Balanced Harvest Yield") +
    theme_minimal()


