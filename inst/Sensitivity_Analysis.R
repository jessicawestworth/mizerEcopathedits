#Sensitivity analysis code
params<-readRDS("/Users/jessicawestworth/Desktop/Status Quo Starting Models/final7.rds")
dm<-readRDS((here("inst","extdata","dm.rds")))

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
    w_repro_max_Blue_whiting = c(0.9*species_params(params)[1,]$w_repro_max, min(max(w(params)),1.10*species_params(params)[1,]$w_repro_max)),

    Eiw_Boarfish = c(0.9*species_params(params)[2,]$Eiw, 1.10*species_params(params)[2,]$Eiw),
    mu_mat_Boarfish = c(0.9*species_params(params)[2,]$mu_mat, 1.10*species_params(params)[2,]$mu_mat),
    d_over_g_Boarfish = c(0.9*species_params(params)[2,]$d_over_g, 1.10*species_params(params)[2,]$d_over_g),
    l_mat_Boarfish = c(0.9*species_params(params)[2,]$l_mat, 1.10*species_params(params)[2,]$l_mat),
    m_Boarfish = c(max(0.7,0.9*species_params(params)[2,]$m), min(1.5,1.10*species_params(params)[2,]$m)),
    w_repro_max_Boarfish = c(0.9*species_params(params)[2,]$w_repro_max, min(max(w(params)),1.10*species_params(params)[2,]$w_repro_max)),

    Eiw_Cod = c(0.9*species_params(params)[3,]$Eiw, 1.10*species_params(params)[3,]$Eiw),
    mu_mat_Cod = c(0.9*species_params(params)[3,]$mu_mat, 1.10*species_params(params)[3,]$mu_mat),
    d_over_g_Cod = c(0.9*species_params(params)[3,]$d_over_g, 1.10*species_params(params)[3,]$d_over_g),
    l_mat_Cod = c(0.9*species_params(params)[3,]$l_mat, 1.10*species_params(params)[3,]$l_mat),
    m_Cod = c(max(0.7,0.9*species_params(params)[3,]$m), min(1.5,1.10*species_params(params)[3,]$m)),
    w_repro_max_Cod = c(0.9*species_params(params)[3,]$w_repro_max, min(max(w(params)),1.10*species_params(params)[3,]$w_repro_max)),

    Eiw_Haddock = c(0.9*species_params(params)[4,]$Eiw, 1.10*species_params(params)[4,]$Eiw),
    mu_mat_Haddock = c(0.9*species_params(params)[4,]$mu_mat, 1.10*species_params(params)[4,]$mu_mat),
    d_over_g_Haddock = c(0.9*species_params(params)[4,]$d_over_g, 1.10*species_params(params)[4,]$d_over_g),
    l_mat_Haddock = c(0.9*species_params(params)[4,]$l_mat, 1.10*species_params(params)[4,]$l_mat),
    m_Haddock = c(max(0.7,0.9*species_params(params)[4,]$m), min(1.5,1.10*species_params(params)[4,]$m)),
    w_repro_max_Haddock = c(0.9*species_params(params)[4,]$w_repro_max, min(max(w(params)),1.10*species_params(params)[4,]$w_repro_max)),

    Eiw_Hake = c(0.9*species_params(params)[5,]$Eiw, 1.10*species_params(params)[5,]$Eiw),
    mu_mat_Hake = c(0.9*species_params(params)[5,]$mu_mat, 1.10*species_params(params)[5,]$mu_mat),
    d_over_g_Hake = c(0.9*species_params(params)[5,]$d_over_g, 1.10*species_params(params)[5,]$d_over_g),
    l_mat_Hake = c(0.9*species_params(params)[5,]$l_mat, 1.10*species_params(params)[5,]$l_mat),
    m_Hake = c(max(0.7,0.9*species_params(params)[5,]$m), min(1.5,1.10*species_params(params)[5,]$m)),
    w_repro_max_Hake = c(0.9*species_params(params)[5,]$w_repro_max, min(max(w(params)),1.10*species_params(params)[5,]$w_repro_max)),

    Eiw_Herring = c(0.9*species_params(params)[6,]$Eiw, 1.10*species_params(params)[6,]$Eiw),
    mu_mat_Herring = c(0.9*species_params(params)[6,]$mu_mat, 1.10*species_params(params)[6,]$mu_mat),
    d_over_g_Herring = c(0.9*species_params(params)[6,]$d_over_g, 1.10*species_params(params)[6,]$d_over_g),
    l_mat_Herring = c(0.9*species_params(params)[6,]$l_mat, 1.10*species_params(params)[6,]$l_mat),
    m_Herring = c(max(0.7,0.9*species_params(params)[6,]$m), min(1.5,1.10*species_params(params)[6,]$m)),
    w_repro_max_Herring = c(0.9*species_params(params)[6,]$w_repro_max, min(max(w(params)),1.10*species_params(params)[6,]$w_repro_max)),

    Eiw_Horse_mackerel = c(0.9*species_params(params)[7,]$Eiw, 1.10*species_params(params)[7,]$Eiw),
    mu_mat_Horse_mackerel = c(0.9*species_params(params)[7,]$mu_mat, 1.10*species_params(params)[7,]$mu_mat),
    d_over_g_Horse_mackerel = c(0.9*species_params(params)[7,]$d_over_g, 1.10*species_params(params)[7,]$d_over_g),
    l_mat_Horse_mackerel = c(0.9*species_params(params)[7,]$l_mat, 1.10*species_params(params)[7,]$l_mat),
    m_Horse_mackerel = c(max(0.7,0.9*species_params(params)[7,]$m), min(1.5,1.10*species_params(params)[7,]$m)),
    w_repro_max_Horse_mackerel = c(0.9*species_params(params)[7,]$w_repro_max, min(max(w(params)),1.10*species_params(params)[7,]$w_repro_max)),

    Eiw_Mackerel = c(0.9*species_params(params)[8,]$Eiw, 1.10*species_params(params)[8,]$Eiw),
    mu_mat_Mackerel = c(0.9*species_params(params)[8,]$mu_mat, 1.10*species_params(params)[8,]$mu_mat),
    d_over_g_Mackerel = c(0.9*species_params(params)[8,]$d_over_g, 1.10*species_params(params)[8,]$d_over_g),
    l_mat_Mackerel = c(0.9*species_params(params)[8,]$l_mat, 1.10*species_params(params)[8,]$l_mat),
    m_Mackerel = c(max(0.7,0.9*species_params(params)[8,]$m), min(1.5,1.10*species_params(params)[8,]$m)),
    w_repro_max_Mackerel = c(0.9*species_params(params)[8,]$w_repro_max, min(max(w(params)),1.10*species_params(params)[8,]$w_repro_max)),

    Eiw_Megrim = c(0.9*species_params(params)[9,]$Eiw, 1.10*species_params(params)[9,]$Eiw),
    mu_mat_Megrim = c(0.9*species_params(params)[9,]$mu_mat, 1.10*species_params(params)[9,]$mu_mat),
    d_over_g_Megrim = c(0.9*species_params(params)[9,]$d_over_g, 1.10*species_params(params)[9,]$d_over_g),
    l_mat_Megrim = c(0.9*species_params(params)[9,]$l_mat, 1.10*species_params(params)[9,]$l_mat),
    m_Megrim = c(max(0.7,0.9*species_params(params)[9,]$m), min(1.5,1.10*species_params(params)[9,]$m)),
    w_repro_max_Megrim = c(0.9*species_params(params)[9,]$w_repro_max, min(max(w(params)),1.10*species_params(params)[9,]$w_repro_max)),

    Eiw_Monkfish = c(0.9*species_params(params)[10,]$Eiw, 1.10*species_params(params)[10,]$Eiw),
    mu_mat_Monkfish = c(0.9*species_params(params)[10,]$mu_mat, 1.10*species_params(params)[10,]$mu_mat),
    d_over_g_Monkfish = c(0.9*species_params(params)[10,]$d_over_g, 1.10*species_params(params)[10,]$d_over_g),
    l_mat_Monkfish = c(0.9*species_params(params)[10,]$l_mat, 1.10*species_params(params)[10,]$l_mat),
    m_Monkfish = c(max(0.7,0.9*species_params(params)[10,]$m), min(1.5,1.10*species_params(params)[10,]$m)),
    w_repro_max_Monkfish = c(0.9*species_params(params)[10,]$w_repro_max, min(max(w(params)),1.10*species_params(params)[10,]$w_repro_max)),

    Eiw_Plaice = c(0.9*species_params(params)[11,]$Eiw, 1.10*species_params(params)[11,]$Eiw),
    mu_mat_Plaice = c(0.9*species_params(params)[11,]$mu_mat, 1.10*species_params(params)[11,]$mu_mat),
    d_over_g_Plaice = c(0.9*species_params(params)[11,]$d_over_g, 1.10*species_params(params)[11,]$d_over_g),
    l_mat_Plaice = c(0.9*species_params(params)[11,]$l_mat, 1.10*species_params(params)[11,]$l_mat),
    m_Plaice = c(max(0.7,0.9*species_params(params)[11,]$m), min(1.5,1.10*species_params(params)[11,]$m)),
    w_repro_max_Plaice = c(0.9*species_params(params)[11,]$w_repro_max, min(max(w(params)),1.10*species_params(params)[11,]$w_repro_max)),

    Eiw_Red_gurnard = c(0.9*species_params(params)[12,]$Eiw, 1.10*species_params(params)[12,]$Eiw),
    mu_mat_Red_gurnard = c(0.9*species_params(params)[12,]$mu_mat, 1.10*species_params(params)[12,]$mu_mat),
    d_over_g_Red_gurnard = c(0.9*species_params(params)[12,]$d_over_g, 1.10*species_params(params)[12,]$d_over_g),
    l_mat_Red_gurnard = c(0.9*species_params(params)[12,]$l_mat, 1.10*species_params(params)[12,]$l_mat),
    m_Red_gurnard = c(max(0.7,0.9*species_params(params)[12,]$m), min(1.5,1.10*species_params(params)[12,]$m)),
    w_repro_max_Red_gurnard = c(0.9*species_params(params)[12,]$w_repro_max, min(max(w(params)),1.10*species_params(params)[12,]$w_repro_max)),

    Eiw_Sole = c(0.9*species_params(params)[13,]$Eiw, 1.10*species_params(params)[13,]$Eiw),
    mu_mat_Sole = c(0.9*species_params(params)[13,]$mu_mat, 1.10*species_params(params)[13,]$mu_mat),
    d_over_g_Sole = c(0.9*species_params(params)[13,]$d_over_g, 1.10*species_params(params)[13,]$d_over_g),
    l_mat_Sole = c(0.9*species_params(params)[13,]$l_mat, 1.10*species_params(params)[13,]$l_mat),
    m_Sole = c(max(0.7,0.9*species_params(params)[13,]$m), min(1.5,1.10*species_params(params)[13,]$m)),
    w_repro_max_Sole = c(0.9*species_params(params)[13,]$w_repro_max, min(max(w(params)),1.10*species_params(params)[13,]$w_repro_max)),

    Eiw_Whiting = c(0.9*species_params(params)[14,]$Eiw, 1.10*species_params(params)[14,]$Eiw),
    mu_mat_Whiting = c(0.9*species_params(params)[14,]$mu_mat, 1.10*species_params(params)[14,]$mu_mat),
    d_over_g_Whiting = c(0.9*species_params(params)[14,]$d_over_g, 1.10*species_params(params)[14,]$d_over_g),
    l_mat_Whiting = c(0.9*species_params(params)[14,]$l_mat, 1.10*species_params(params)[14,]$l_mat),
    m_Whiting = c(max(0.7,0.9*species_params(params)[14,]$m), min(1.5,1.10*species_params(params)[14,]$m)),
    w_repro_max_Whiting = c(0.9*species_params(params)[14,]$w_repro_max, min(max(w(params)),1.10*species_params(params)[14,]$w_repro_max))
)

# w_repro_max commented will mean changing w_repro_max doesn't do anything...
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
    species_params(params_i)[1,]$w_repro_max<-input_design[i, "w_repro_max_Blue_whiting"]

    species_params(params_i)[2,]$Eiw<-input_design[i, "Eiw_Boarfish"]
    species_params(params_i)[2,]$mu_mat<-input_design[i, "mu_mat_Boarfish"]
    species_params(params_i)[2,]$d_over_g<-input_design[i, "d_over_g_Boarfish"]
    species_params(params_i)[2,]$l_mat<-input_design[i, "l_mat_Boarfish"]
    species_params(params_i)[2,]$m<-input_design[i, "m_Boarfish"]
    species_params(params_i)[2,]$w_repro_max<-input_design[i, "w_repro_max_Boarfish"]

    species_params(params_i)[3,]$Eiw<-input_design[i, "Eiw_Cod"]
    species_params(params_i)[3,]$mu_mat<-input_design[i, "mu_mat_Cod"]
    species_params(params_i)[3,]$d_over_g<-input_design[i, "d_over_g_Cod"]
    species_params(params_i)[3,]$l_mat<-input_design[i, "l_mat_Cod"]
    species_params(params_i)[3,]$m<-input_design[i, "m_Cod"]
    species_params(params_i)[3,]$w_repro_max<-input_design[i, "w_repro_max_Cod"]

    species_params(params_i)[4,]$Eiw<-input_design[i, "Eiw_Haddock"]
    species_params(params_i)[4,]$mu_mat<-input_design[i, "mu_mat_Haddock"]
    species_params(params_i)[4,]$d_over_g<-input_design[i, "d_over_g_Haddock"]
    species_params(params_i)[4,]$l_mat<-input_design[i, "l_mat_Haddock"]
    species_params(params_i)[4,]$m<-input_design[i, "m_Haddock"]
    species_params(params_i)[4,]$w_repro_max<-input_design[i, "w_repro_max_Haddock"]

    species_params(params_i)[5,]$Eiw<-input_design[i, "Eiw_Hake"]
    species_params(params_i)[5,]$mu_mat<-input_design[i, "mu_mat_Hake"]
    species_params(params_i)[5,]$d_over_g<-input_design[i, "d_over_g_Hake"]
    species_params(params_i)[5,]$l_mat<-input_design[i, "l_mat_Hake"]
    species_params(params_i)[5,]$m<-input_design[i, "m_Hake"]
    species_params(params_i)[5,]$w_repro_max<-input_design[i, "w_repro_max_Hake"]

    species_params(params_i)[6,]$Eiw<-input_design[i, "Eiw_Herring"]
    species_params(params_i)[6,]$mu_mat<-input_design[i, "mu_mat_Herring"]
    species_params(params_i)[6,]$d_over_g<-input_design[i, "d_over_g_Herring"]
    species_params(params_i)[6,]$l_mat<-input_design[i, "l_mat_Herring"]
    species_params(params_i)[6,]$m<-input_design[i, "m_Herring"]
    species_params(params_i)[6,]$w_repro_max<-input_design[i, "w_repro_max_Herring"]

    species_params(params_i)[7,]$Eiw<-input_design[i, "Eiw_Horse_mackerel"]
    species_params(params_i)[7,]$mu_mat<-input_design[i, "mu_mat_Horse_mackerel"]
    species_params(params_i)[7,]$d_over_g<-input_design[i, "d_over_g_Horse_mackerel"]
    species_params(params_i)[7,]$l_mat<-input_design[i, "l_mat_Horse_mackerel"]
    species_params(params_i)[7,]$m<-input_design[i, "m_Horse_mackerel"]
    species_params(params_i)[7,]$w_repro_max<-input_design[i, "w_repro_max_Horse_mackerel"]

    species_params(params_i)[8,]$Eiw<-input_design[i, "Eiw_Mackerel"]
    species_params(params_i)[8,]$mu_mat<-input_design[i, "mu_mat_Mackerel"]
    species_params(params_i)[8,]$d_over_g<-input_design[i, "d_over_g_Mackerel"]
    species_params(params_i)[8,]$l_mat<-input_design[i, "l_mat_Mackerel"]
    species_params(params_i)[8,]$m<-input_design[i, "m_Mackerel"]
    species_params(params_i)[8,]$w_repro_max<-input_design[i, "w_repro_max_Mackerel"]

    species_params(params_i)[9,]$Eiw<-input_design[i, "Eiw_Megrim"]
    species_params(params_i)[9,]$mu_mat<-input_design[i, "mu_mat_Megrim"]
    species_params(params_i)[9,]$d_over_g<-input_design[i, "d_over_g_Megrim"]
    species_params(params_i)[9,]$l_mat<-input_design[i, "l_mat_Megrim"]
    species_params(params_i)[9,]$m<-input_design[i, "m_Megrim"]
    species_params(params_i)[9,]$w_repro_max<-input_design[i, "w_repro_max_Megrim"]

    species_params(params_i)[10,]$Eiw<-input_design[i, "Eiw_Monkfish"]
    species_params(params_i)[10,]$mu_mat<-input_design[i, "mu_mat_Monkfish"]
    species_params(params_i)[10,]$d_over_g<-input_design[i, "d_over_g_Monkfish"]
    species_params(params_i)[10,]$l_mat<-input_design[i, "l_mat_Monkfish"]
    species_params(params_i)[10,]$m<-input_design[i, "m_Monkfish"]
    species_params(params_i)[10,]$w_repro_max<-input_design[i, "w_repro_max_Monkfish"]

    species_params(params_i)[11,]$Eiw<-input_design[i, "Eiw_Plaice"]
    species_params(params_i)[11,]$mu_mat<-input_design[i, "mu_mat_Plaice"]
    species_params(params_i)[11,]$d_over_g<-input_design[i, "d_over_g_Plaice"]
    species_params(params_i)[11,]$l_mat<-input_design[i, "l_mat_Plaice"]
    species_params(params_i)[11,]$m<-input_design[i, "m_Plaice"]
    species_params(params_i)[11,]$w_repro_max<-input_design[i, "w_repro_max_Plaice"]

    species_params(params_i)[12,]$Eiw<-input_design[i, "Eiw_Red_gurnard"]
    species_params(params_i)[12,]$mu_mat<-input_design[i, "mu_mat_Red_gurnard"]
    species_params(params_i)[12,]$d_over_g<-input_design[i, "d_over_g_Red_gurnard"]
    species_params(params_i)[12,]$l_mat<-input_design[i, "l_mat_Red_gurnard"]
    species_params(params_i)[12,]$m<-input_design[i, "m_Red_gurnard"]
    species_params(params_i)[12,]$w_repro_max<-input_design[i, "w_repro_max_Red_gurnard"]

    species_params(params_i)[13,]$Eiw<-input_design[i, "Eiw_Sole"]
    species_params(params_i)[13,]$mu_mat<-input_design[i, "mu_mat_Sole"]
    species_params(params_i)[13,]$d_over_g<-input_design[i, "d_over_g_Sole"]
    species_params(params_i)[13,]$l_mat<-input_design[i, "l_mat_Sole"]
    species_params(params_i)[13,]$m<-input_design[i, "m_Sole"]
    species_params(params_i)[13,]$w_repro_max<-input_design[i, "w_repro_max_Sole"]

    species_params(params_i)[14,]$Eiw<-input_design[i, "Eiw_Whiting"]
    species_params(params_i)[14,]$mu_mat<-input_design[i, "mu_mat_Whiting"]
    species_params(params_i)[14,]$d_over_g<-input_design[i, "d_over_g_Whiting"]
    species_params(params_i)[14,]$l_mat<-input_design[i, "l_mat_Whiting"]
    species_params(params_i)[14,]$m<-input_design[i, "m_Whiting"]
    species_params(params_i)[14,]$w_repro_max<-input_design[i, "w_repro_max_Whiting"]

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




for (i in 1:nrow(input_design)) {
    sim_i<-sensitivity_sim_list[[paste0(i)]]
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

# calculate the Mu and Sigma (elementary effects)
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
#for Y top 5 are m_Horse_mackerel, d/g blue whiting, w_repro_max whiting, mu_mat Blue whiting, w_repro_max Cod, m_Cod
#For B top 5 are m_Horse_mackerel, w_repro_max Whiting, d/g Blue whiting, w_repro_max cod, m_Boarfish
#For SSB top 5 are mu_mat blue whiting, m_horse mackerel, w_repro_max whiting, d/g blue whiting, d/g sole
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
    c = c(10),
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
factors <- c("Eiw_Mackerel","Eiw_Horse_mackerel","Eiw_Boarfish", "Eiw_Herring","Eiw_Whiting", "Eiw_Cod", "Eiw_Hake")
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
    # Update species-specific parameters for this run
    species_params(params_i)[2,]$Eiw<-sim_params[i, "Eiw_Boarfish"]

    species_params(params_i)[3,]$Eiw<-sim_params[i, "Eiw_Cod"]

    species_params(params_i)[5,]$Eiw<-sim_params[i, "Eiw_Hake"]

    species_params(params_i)[6,]$Eiw<-sim_params[i, "Eiw_Herring"]

    species_params(params_i)[7,]$Eiw<-sim_params[i, "Eiw_Horse_mackerel"]

    species_params(params_i)[8,]$Eiw<-sim_params[i, "Eiw_Mackerel"]

    species_params(params_i)[14,]$Eiw<-sim_params[i, "Eiw_Whiting"]

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

#saveRDS(RSA_outputs, "/Users/jessicawestworth/Desktop/Sensitivity/RSA_outputs_c_10.rds")
#saveRDS(RSA_sim_list, "/Users/jessicawestworth/Desktop/Sensitivity/RSA_sim_list_c_10.rds")


RSA_outputs_0.1<-readRDS("/Users/jessicawestworth/Desktop/Sensitivity/RSA_outputs_c_0.1.rds")
RSA_outputs_0.2<-readRDS("/Users/jessicawestworth/Desktop/Sensitivity/RSA_outputs_c_0.2.rds")
RSA_outputs_0.3<-readRDS("/Users/jessicawestworth/Desktop/Sensitivity/RSA_outputs_c_0.3.rds")
RSA_outputs_0.4<-readRDS("/Users/jessicawestworth/Desktop/Sensitivity/RSA_outputs_c_0.4.rds")
RSA_outputs_1<-readRDS("/Users/jessicawestworth/Desktop/Sensitivity/RSA_outputs_c_1.rds")
RSA_outputs_10<-readRDS("/Users/jessicawestworth/Desktop/Sensitivity/RSA_outputs_c_10.rds")
RSA_outputs_0.4_0.5<-readRDS("/Users/jessicawestworth/Desktop/Sensitivity/RSA_outputs_c_0.4_a_0.5.rds")

RSA_outputs_0.1$type<-"Whole_BH_0.1"
RSA_outputs_0.2$type<-"Whole_BH_0.2"
RSA_outputs_0.3$type<-"Whole_BH_0.3"
RSA_outputs_0.4$type<-"Whole_BH_0.4"
RSA_outputs_1$type<-"Whole_BH_1"
RSA_outputs_10$type<-"Whole_BH_10"
RSA_outputs_0.4_0.5$type<-"Hybrid_BH_0.4_0.5"
RSA_outputs_status<-RSA_outputs_0.4_0.5
RSA_outputs_status$BH_Y<-RSA_outputs_0.4_0.5$status_Y
RSA_outputs_status$BH_B<-RSA_outputs_0.4_0.5$status_B
RSA_outputs_status$BH_SSB<-RSA_outputs_0.4_0.5$staus_SSB
RSA_outputs_status$BH_N<-RSA_outputs_0.4_0.5$status_N
RSA_outputs_status$BH_LFY<-RSA_outputs_0.4_0.5$status_LFY
RSA_outputs_status$BH_LFY_t<-RSA_outputs_0.4_0.5$status_LFY_t
RSA_outputs_status$BH_LFB<-RSA_outputs_0.4_0.5$status_LFB
RSA_outputs_status$BH_LFB_t<-RSA_outputs_0.4_0.5$status_LFB_t
RSA_outputs_status$type<-"Status"
RSA_outputs_status$c<-0
RSA_outputs_status$alpha<-0


RSA_outputs<-rbind(RSA_outputs_0.1, RSA_outputs_0.2,RSA_outputs_0.3,RSA_outputs_0.4,RSA_outputs_0.4_0.5, RSA_outputs_1,RSA_outputs_10,RSA_outputs_status)

#library(ggplot2)
#is it always the case that yields are higher for BH scenarios

RSA_outputs$BH_Y_wins<-RSA_outputs$BH_Y>RSA_outputs$status_Y

ggplot(RSA_outputs, aes(x = status_Y, y = BH_Y)) +
    geom_point(aes(color = type, shape = BH_Y_wins), alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = "Is Balanced Harvest Yield Robustly Higher?",
         x = "Status Quo Yield", y = "Balanced Harvest Yield") +
    theme_minimal()

RSA_outputs$Y_dif<-RSA_outputs$BH_Y-RSA_outputs$status_Y

#for whole BH scenario 0.1 (2% were less, 98% were higher)
fr<-RSA_outputs%>%filter(type =="Whole_BH_0.1")
sum(fr$BH_Y_wins == FALSE, na.rm=TRUE)
min(fr$Y_dif) #-0.9462719
max(fr$Y_dif) #34.39567
mean(fr$Y_dif) # 9.876035

#for whole BH scenario 0.2, 0.3, 0.4 all were higher
fr<-RSA_outputs%>%filter(type =="Whole_BH_0.2")
sum(fr$BH_Y_wins == FALSE, na.rm=TRUE)
min(fr$Y_dif) #0.945
max(fr$Y_dif) #50.99154
mean(fr$Y_dif) #18.31468

fr<-RSA_outputs%>%filter(type =="Whole_BH_0.3")
sum(fr$BH_Y_wins == FALSE, na.rm=TRUE)
min(fr$Y_dif) #2.090008
max(fr$Y_dif) #53.20824
mean(fr$Y_dif) #20.61162

fr<-RSA_outputs%>%filter(type =="Whole_BH_0.4")
sum(fr$BH_Y_wins == FALSE, na.rm=TRUE)
min(fr$Y_dif) #2.719592
max(fr$Y_dif) #51.19552
mean(fr$Y_dif) #20.68017

fr<-RSA_outputs%>%filter(type =="Whole_BH_1")
sum(fr$BH_Y_wins == FALSE, na.rm=TRUE)
min(fr$Y_dif) #2.214762
max(fr$Y_dif) #34.81972
mean(fr$Y_dif) #15.26847

#for whole BH scenario 10 (9% were less, 91% were higher)
fr<-RSA_outputs%>%filter(type =="Whole_BH_10")
sum(fr$BH_Y_wins == FALSE, na.rm=TRUE)
min(fr$Y_dif) #-0.4210242
max(fr$Y_dif) #4.60189
mean(fr$Y_dif) #1.707069

#for hybrid previous best model, all were higher
fr<-RSA_outputs%>%filter(type =="Hybrid_BH_0.4_0.5")
sum(fr$BH_Y_wins == FALSE, na.rm=TRUE)
min(fr$Y_dif) #2.479561
max(fr$Y_dif) #50.88756
mean(fr$Y_dif) #19.30815

RSA_outputs$BH_B_wins<-RSA_outputs$BH_B>RSA_outputs$status_B

ggplot(RSA_outputs, aes(x = status_B, y = BH_B)) +
    geom_point(aes(color = type, shape = BH_B_wins), alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = "Is Balanced Harvest Biomass Robustly Higher?",
         x = "Status Quo Biomass", y = "Balanced Harvest Biomass") +
    theme_minimal()

RSA_outputs$B_dif<-RSA_outputs$BH_B-RSA_outputs$status_B

#c10 and c1 was worse in every model
fr<-RSA_outputs%>%filter(type =="Whole_BH_10")
sum(fr$BH_B_wins == FALSE, na.rm=TRUE)
min(fr$B_dif) #-57.2192
max(fr$B_dif) #-15.23675
mean(fr$B_dif) #-30.45759

fr<-RSA_outputs%>%filter(type =="Whole_BH_1")
sum(fr$BH_B_wins == FALSE, na.rm=TRUE)
min(fr$B_dif) #-38.99268
max(fr$B_dif) #-6.139575
mean(fr$B_dif) #-17.38153

#c0.4
fr<-RSA_outputs%>%filter(type =="Whole_BH_0.4")
sum(fr$BH_B_wins == FALSE, na.rm=TRUE)
min(fr$B_dif) #-23.07788
max(fr$B_dif) #1.411962
mean(fr$B_dif) #-7.011809
#95% of the time 0.4 c had worse B than status
#should check the range of how much lower it went
fr<-RSA_outputs%>%filter(type =="Whole_BH_0.3")
sum(fr$BH_B_wins == FALSE, na.rm=TRUE)
min(fr$B_dif) #-17.30343
max(fr$B_dif) #4.835626
mean(fr$B_dif) #-3.511683
#73% of the time 0.4 c had worse B than status
#should check the range of how much lower it went
fr<-RSA_outputs%>%filter(type =="Whole_BH_0.2")
sum(fr$BH_B_wins == FALSE, na.rm=TRUE)
min(fr$B_dif) #-8.835737
max(fr$B_dif) #8.594493
mean(fr$B_dif) #0.9787319
#43% of the time 0.4 c had worse B than status
#should check the range of how much lower it went
fr<-RSA_outputs%>%filter(type =="Whole_BH_0.1")
sum(fr$BH_B_wins == FALSE, na.rm=TRUE)
min(fr$B_dif) #-3.880961
max(fr$B_dif) #14.1626
mean(fr$B_dif) #5.375199
#16% of the time 0.4 c had worse B than status
#should check the range of how much lower it went
fr<-RSA_outputs%>%filter(type =="Hybrid_BH_0.4_0.5")
sum(fr$BH_B_wins == FALSE, na.rm=TRUE)
min(fr$B_dif) #-11.14606
max(fr$B_dif) #6.264971
mean(fr$B_dif) #-0.9748913
#60% of the time 0.4 c had worse B than status
#should check the range of how much lower it went

RSA_outputs$BH_Y_wins<-RSA_outputs$BH_Y>RSA_outputs$status_Y

ggplot(RSA_outputs, aes(x = status_Y, y = BH_Y)) +
    geom_point(aes(color = type, shape = BH_Y_wins), alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = "Is Balanced Harvest Yield Robustly Higher?",
         x = "Status Quo Yield", y = "Balanced Harvest Yield") +
    theme_minimal()

#SSB
RSA_outputs$BH_SSB_wins<-RSA_outputs$BH_SSB>RSA_outputs$staus_SSB

ggplot(RSA_outputs, aes(x = staus_SSB, y = BH_SSB)) +
    geom_point(aes(color = type, shape = BH_SSB_wins), alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = "Is Balanced Harvest Yield Robustly Higher?",
         x = "Status Quo Yield", y = "Balanced Harvest Yield") +
    theme_minimal()

RSA_outputs$SSB_dif<-RSA_outputs$BH_SSB-RSA_outputs$staus_SSB

#c10 and c1 was worse in every model
fr<-RSA_outputs%>%filter(type =="Whole_BH_10")
sum(fr$BH_SSB_wins == FALSE, na.rm=TRUE)
min(fr$SSB_dif) #-27.05369
max(fr$SSB_dif) #-5.947652
mean(fr$SSB_dif) #-13.59939

fr<-RSA_outputs%>%filter(type =="Whole_BH_1")
sum(fr$BH_SSB_wins == FALSE, na.rm=TRUE)
min(fr$SSB_dif) #-21.0008
max(fr$SSB_dif) #-1.766703
mean(fr$SSB_dif) #-8.467218

fr<-RSA_outputs%>%filter(type =="Whole_BH_0.4")
sum(fr$BH_SSB_wins == FALSE, na.rm=TRUE)
min(fr$SSB_dif) #-16.08825
max(fr$SSB_dif) #0.09748805
mean(fr$SSB_dif) #-4.973325

#98% of the time 0.4 c had worse SSB than status
#should check the range of how much lower it went
fr<-RSA_outputs%>%filter(type =="Whole_BH_0.3")
sum(fr$BH_SSB_wins == FALSE, na.rm=TRUE)
min(fr$SSB_dif) #-14.12081
max(fr$SSB_dif) # 0.739524
mean(fr$SSB_dif) #-3.607708

#92% of the time 0.3 c had worse SSB than status
#should check the range of how much lower it went
fr<-RSA_outputs%>%filter(type =="Whole_BH_0.2")
sum(fr$BH_SSB_wins == FALSE, na.rm=TRUE)
min(fr$SSB_dif) #-10.59178
max(fr$SSB_dif) # 2.110734
mean(fr$SSB_dif) #-1.550402

#71% of the time 0.2 c had worse SSB than status
#should check the range of how much lower it went
fr<-RSA_outputs%>%filter(type =="Whole_BH_0.1")
sum(fr$BH_SSB_wins == FALSE, na.rm=TRUE)
min(fr$SSB_dif) #-3.659336
max(fr$SSB_dif) # 4.92953
mean(fr$SSB_dif) #1.510525

#14% of the time 0.1 c had worse SSB than status
#should check the range of how much lower it went
fr<-RSA_outputs%>%filter(type =="Hybrid_BH_0.4_0.5")
sum(fr$BH_SSB_wins == FALSE, na.rm=TRUE)
min(fr$SSB_dif) #-12.79578
max(fr$SSB_dif) # -0.01919633
mean(fr$SSB_dif) #-3.449908
#100% of the time hybrid had worse B than status

#N
RSA_outputs$BH_N_wins<-RSA_outputs$BH_N>RSA_outputs$status_N


ggplot(RSA_outputs, aes(x = status_N, y = BH_N)) +
    geom_point(aes(color = type, shape = BH_N_wins), alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = "Is Balanced Harvest N Robustly Higher?",
         x = "Status Quo N", y = "Balanced Harvest N") +
    theme_minimal()

RSA_outputs$N_dif<-RSA_outputs$BH_N-RSA_outputs$status_N

fr<-RSA_outputs%>%filter(type =="Whole_BH_10")
sum(fr$BH_N_wins == FALSE, na.rm=TRUE)
min(fr$N_dif) #-62.02089
max(fr$N_dif) # 6.531034
mean(fr$N_dif) #-40.45637
#97%
fr<-RSA_outputs%>%filter(type =="Whole_BH_1")
sum(fr$BH_N_wins == FALSE, na.rm=TRUE)
min(fr$N_dif) #-35.43654
max(fr$N_dif) #81.14447
mean(fr$N_dif) #3.239781
#64%
fr<-RSA_outputs%>%filter(type =="Whole_BH_0.4")
sum(fr$BH_N_wins == FALSE, na.rm=TRUE)
min(fr$N_dif) #-18.52049
max(fr$N_dif) #94.99513
mean(fr$N_dif) #20.14019
#26% of the time 0.4 c had worse SSB than status
#should check the range of how much lower it went
fr<-RSA_outputs%>%filter(type =="Whole_BH_0.3")
sum(fr$BH_N_wins == FALSE, na.rm=TRUE)
min(fr$N_dif) #-13.77172
max(fr$N_dif) #93.29163
mean(fr$N_dif) #22.54266
#22% of the time 0.4 c had worse SSB than status
#should check the range of how much lower it went
fr<-RSA_outputs%>%filter(type =="Whole_BH_0.2")
sum(fr$BH_N_wins == FALSE, na.rm=TRUE)
min(fr$N_dif) #-8.197448
max(fr$N_dif) # 87.59062
mean(fr$N_dif) #23.23397
#14% of the time 0.4 c had worse SSB than status
#should check the range of how much lower it went
fr<-RSA_outputs%>%filter(type =="Whole_BH_0.1")
sum(fr$BH_N_wins == FALSE, na.rm=TRUE)
min(fr$N_dif) # -2.297293
max(fr$N_dif) # 67.17026
mean(fr$N_dif) # 19.559
#5% of the time 0.4 c had worse SSB than status
#should check the range of how much lower it went
fr<-RSA_outputs%>%filter(type =="Hybrid_BH_0.4_0.5")
sum(fr$BH_N_wins == FALSE, na.rm=TRUE)
min(fr$N_dif) # -8.475009
max(fr$N_dif) # 85.91678
mean(fr$N_dif) # 22.04883
#16% of the time 0.4 c had worse B than status

#LFY
RSA_outputs$BH_LFY_wins<-RSA_outputs$BH_LFY>RSA_outputs$status_LFY

ggplot(RSA_outputs, aes(x = status_LFY, y = BH_LFY)) +
    geom_point(aes(color = type, shape = BH_LFY_wins), alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = "Is Balanced Harvest N Robustly Higher?",
         x = "Status Quo PBLFY", y = "Balanced Harvest PBLFY") +
    theme_minimal()

#worse in every model
RSA_outputs$LFY_dif<-RSA_outputs$BH_LFY-RSA_outputs$status_LFY

fr<-RSA_outputs%>%filter(type =="Whole_BH_10")
sum(fr$BH_LFY_wins == FALSE, na.rm=TRUE)
min(fr$LFY_dif) #-0.2595394
max(fr$LFY_dif) #-0.07548432
mean(fr$LFY_dif) #-0.1552398

fr<-RSA_outputs%>%filter(type =="Whole_BH_1")
sum(fr$BH_LFY_wins == FALSE, na.rm=TRUE)
min(fr$LFY_dif) #-0.3104226
max(fr$LFY_dif) #-0.1189156
mean(fr$LFY_dif) #-0.2025145

fr<-RSA_outputs%>%filter(type =="Whole_BH_0.4")
sum(fr$BH_LFY_wins == FALSE, na.rm=TRUE)
min(fr$LFY_dif) #-0.3283499
max(fr$LFY_dif) #-0.1253051
mean(fr$LFY_dif) #-0.2170616

fr<-RSA_outputs%>%filter(type =="Whole_BH_0.3")
sum(fr$BH_LFY_wins == FALSE, na.rm=TRUE)
min(fr$LFY_dif) #-0.3310667
max(fr$LFY_dif) #-0.1258017
mean(fr$LFY_dif) #-0.2193648

fr<-RSA_outputs%>%filter(type =="Whole_BH_0.2")
sum(fr$BH_LFY_wins == FALSE, na.rm=TRUE)
min(fr$LFY_dif) #-0.3333066
max(fr$LFY_dif) #-0.126159
mean(fr$LFY_dif) #-0.2215168

fr<-RSA_outputs%>%filter(type =="Whole_BH_0.1")
sum(fr$BH_LFY_wins == FALSE, na.rm=TRUE)
min(fr$LFY_dif) #-0.3331243
max(fr$LFY_dif) #-0.1259294
mean(fr$LFY_dif) #-0.2231339

fr<-RSA_outputs%>%filter(type =="Hybrid_BH_0.4_0.5")
sum(fr$BH_LFY_wins == FALSE, na.rm=TRUE)
min(fr$LFY_dif) #-0.2731015
max(fr$LFY_dif) #-0.09676441
mean(fr$LFY_dif) #-0.1698722

#BLFY
RSA_outputs$BH_LFY_t_wins<-RSA_outputs$BH_LFY_t>RSA_outputs$status_LFY_t

ggplot(RSA_outputs, aes(x = status_LFY_t, y = BH_LFY_t)) +
    geom_point(aes(color = type, shape = BH_LFY_t_wins), alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = "Is Balanced Harvest N Robustly Higher?",
         x = "Status Quo BLFY", y = "Balanced Harvest BLFY") +
    theme_minimal()

#Status val
min(fr$status_LFY_t)
#0.5052038
max(fr$status_LFY_t)
#1.493969
mean(fr$status_LFY_t)
#0.9094909

#Worse in all whole Bh models
fr<-RSA_outputs%>%filter(type =="Whole_BH_10")
sum(fr$BH_LFY_t_wins == FALSE, na.rm=TRUE)
min(fr$LFY_t_dif) #-1.082762
max(fr$LFY_t_dif) #-0.2271678
mean(fr$LFY_t_dif) #-0.5118275

fr<-RSA_outputs%>%filter(type =="Whole_BH_1")
sum(fr$BH_LFY_t_wins == FALSE, na.rm=TRUE)
min(fr$LFY_t_dif) #-1.072078
max(fr$LFY_t_dif) #-0.2520238
mean(fr$LFY_t_dif) #-0.4674014

fr<-RSA_outputs%>%filter(type =="Whole_BH_0.4")
sum(fr$BH_LFY_t_wins == FALSE, na.rm=TRUE)
min(fr$LFY_t_dif) #-1.175242
max(fr$LFY_t_dif) #-0.3421992
mean(fr$LFY_t_dif) #-0.6240567

fr<-RSA_outputs%>%filter(type =="Whole_BH_0.3")
sum(fr$BH_LFY_t_wins == FALSE, na.rm=TRUE)
min(fr$LFY_t_dif) #-1.208075
max(fr$LFY_t_dif) #-0.3667287
mean(fr$LFY_t_dif) #-0.6683467

fr<-RSA_outputs%>%filter(type =="Whole_BH_0.2")
sum(fr$BH_LFY_t_wins == FALSE, na.rm=TRUE)
min(fr$LFY_t_dif) #-1.250066
max(fr$LFY_t_dif) #-0.3992449
mean(fr$LFY_t_dif) #-0.7269588

fr<-RSA_outputs%>%filter(type =="Whole_BH_0.1")
sum(fr$BH_LFY_t_wins == FALSE, na.rm=TRUE)
min(fr$LFY_t_dif) #-1.332519
max(fr$LFY_t_dif) #-0.451468
mean(fr$LFY_t_dif) #-0.8124832

fr<-RSA_outputs%>%filter(type =="Hybrid_BH_0.4_0.5")
sum(fr$BH_LFY_t_wins == FALSE, na.rm=TRUE)
min(fr$LFY_t_dif) #-0.09173648
max(fr$LFY_t_dif) #0.5436264
mean(fr$LFY_t_dif) #0.1388908
#15% of models worse

#BLFS
RSA_outputs$BH_LFB_t_wins<-RSA_outputs$BH_LFB_t>RSA_outputs$status_LFB_t

ggplot(RSA_outputs, aes(x = status_LFB_t, y = BH_LFB_t)) +
    geom_point(aes(color = type, shape = BH_LFB_t_wins), alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = "Is Balanced Harvest N Robustly Higher?",
         x = "Status Quo BLFS", y = "Balanced Harvest BLFS") +
    theme_minimal()
#Better in whole BH 0.1,0.2,0.3,0.4,1, and hybrid
fr<-RSA_outputs%>%filter(type =="Whole_BH_10")
sum(fr$BH_LFB_t_wins == FALSE, na.rm=TRUE)
#46% of models worse

#LFS
RSA_outputs$BH_LFB_wins<-RSA_outputs$BH_LFB>RSA_outputs$status_LFB

ggplot(RSA_outputs, aes(x = status_LFB, y = BH_LFB)) +
    geom_point(aes(color = type, shape = BH_LFB_wins), alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = "Is Balanced Harvest N Robustly Higher?",
         x = "Status Quo PLFS", y = "Balanced Harvest PLFS") +
    theme_minimal()
#Better in all alternative regimes

#points
RSA_outputs_trial<-RSA_outputs%>%
    #filter(type != "Status")%>%
    group_by(run)%>%
    mutate(Y=(BH_Y-min(BH_Y))/(max(BH_Y)-min(BH_Y)),
           N=(BH_N-min(BH_N))/(max(BH_N)-min(BH_N)),
           B=(BH_B-min(BH_B))/(max(BH_B)-min(BH_B)),
           SSB=(BH_SSB-min(BH_SSB))/(max(BH_SSB)-min(BH_SSB)),
           LFB=(BH_LFB-min(BH_LFB))/(max(BH_LFB)-min(BH_LFB)),
           LFB_t=(BH_LFB_t-min(BH_LFB_t))/(max(BH_LFB_t)-min(BH_LFB_t)),
           LFY=(BH_LFY-min(BH_LFY))/(max(BH_LFY)-min(BH_LFY)),
           LFY_t=(BH_LFY_t-min(BH_LFY_t))/(max(BH_LFY_t)-min(BH_LFY_t)),
           conservation=N+B+SSB+LFB+LFB_t,
           economic=Y+LFY+LFY_t,
           MCCS=conservation+economic,
           beat_status=case_when(MCCS > MCCS[type == "Status"]~ TRUE,
                                 TRUE ~ FALSE),
           hybrid_wins=case_when(type == "Hybrid_BH_0.4_0.5" & MCCS > max(MCCS[type != "Hybrid_BH_0.4_0.5"], na.rm = TRUE) ~ TRUE,
                                 TRUE ~ FALSE),
           status_wins=case_when(type == "Status" & MCCS > max(MCCS[type != "Status"], na.rm = TRUE) ~ TRUE,
                                 TRUE ~ FALSE),
           c_0.1_wins=case_when(type == "Whole_BH_0.1" & MCCS > max(MCCS[type != "Whole_BH_0.1"], na.rm = TRUE) ~ TRUE,
                                 TRUE ~ FALSE),
           c_0.2_wins=case_when(type == "Whole_BH_0.2" & MCCS > max(MCCS[type != "Whole_BH_0.2"], na.rm = TRUE) ~ TRUE,
                                TRUE ~ FALSE),
           c_0.3_wins=case_when(type == "Whole_BH_0.3" & MCCS > max(MCCS[type != "Whole_BH_0.3"], na.rm = TRUE) ~ TRUE,
                                TRUE ~ FALSE),
           c_0.4_wins=case_when(type == "Whole_BH_0.4" & MCCS > max(MCCS[type != "Whole_BH_0.4"], na.rm = TRUE) ~ TRUE,
                                TRUE ~ FALSE),
           c_1_wins=case_when(type == "Whole_BH_1" & MCCS > max(MCCS[type != "Whole_BH_1"], na.rm = TRUE) ~ TRUE,
                                TRUE ~ FALSE),
           c_10_wins=case_when(type == "Whole_BH_10" & MCCS > max(MCCS[type != "Whole_BH_10"], na.rm = TRUE) ~ TRUE,
                                TRUE ~ FALSE))%>%
        ungroup()

ggplot(RSA_outputs_trial, aes(x = type, y = MCCS)) +
    geom_point(aes(color=run), alpha = 0.5) +
    labs(title = "Is Balanced Harvest Yield Robustly Higher?",
         x = "type", y = "MCCS") +
    theme_minimal()


ggplot(RSA_outputs_trial, aes(x = type, y = MCCS)) +
    geom_point(aes(color=run), alpha = 0.5) +
    labs(title = "Is Balanced Harvest Yield Robustly Higher?",
         x = "type", y = "MCCS") +
    theme_minimal()

RSA_outputs_trial_1<-RSA_outputs_trial%>%filter(type== "Hybrid_BH_0.4_0.5")
sum(RSA_outputs_trial_1$hybrid_wins == FALSE, na.rm = TRUE)
mean(RSA_outputs_trial_1$MCCS) #4.928283
#doesn't win 75% of the time

#hybrid not necessarily the best ranking scenario, potential that other hybrid
#scenarios rank better under each regime, we did 6 BH scenarios versus 1 hybrid
#what this is saying is that a regime of 0.4 and 0.5 wansn't always better than
#the rest in fact in most cases it wasn't 75% of the time in fact
#this could mean that different hybrid scenarios in each run would have been
#better but the computing power to create 220 models for 100 runs with different
#parameter values, was 22,000 models each taking 10 min with a total run time of
#220,000 minutes, or 3333 hours, or 138.888 days.

RSA_outputs_trial_1<-RSA_outputs_trial%>%filter(type== "Status")
sum(RSA_outputs_trial_1$status_wins == FALSE, na.rm = TRUE)
mean(RSA_outputs_trial_1$MCCS) #4.282
#doesn't win in 94% of cases

RSA_outputs_trial_1<-RSA_outputs_trial%>%filter(type== "Whole_BH_0.1")
sum(RSA_outputs_trial_1$c_0.1_wins == FALSE, na.rm = TRUE)
mean(RSA_outputs_trial_1$MCCS) #4.688
#doesn't win in 95% of cases

RSA_outputs_trial_1<-RSA_outputs_trial%>%filter(type== "Whole_BH_0.2")
sum(RSA_outputs_trial_1$c_0.2_wins == FALSE, na.rm = TRUE)
mean(RSA_outputs_trial_1$MCCS) #4.966
#doesn't win in 66% of cases

RSA_outputs_trial_1<-RSA_outputs_trial%>%filter(type== "Whole_BH_0.3")
sum(RSA_outputs_trial_1$c_0.3_wins == FALSE, na.rm = TRUE)
mean(RSA_outputs_trial_1$MCCS) #4.93677
#doesn't win in 77% of cases

RSA_outputs_trial_1<-RSA_outputs_trial%>%filter(type== "Whole_BH_0.4")
sum(RSA_outputs_trial_1$c_0.4_wins == FALSE, na.rm = TRUE)
mean(RSA_outputs_trial_1$MCCS)# 4.825315
#doesn't win in 93% of cases

RSA_outputs_trial_1<-RSA_outputs_trial%>%filter(type== "Whole_BH_1")
sum(RSA_outputs_trial_1$c_1_wins == FALSE, na.rm = TRUE)
mean(RSA_outputs_trial_1$MCCS) #4.16
#doesn't win in 100% of cases

RSA_outputs_trial_1<-RSA_outputs_trial%>%filter(type== "Whole_BH_10")
sum(RSA_outputs_trial_1$c_10_wins == FALSE, na.rm = TRUE)
mean(RSA_outputs_trial_1$MCCS) #1.734
#doesn't win in 100% of cases


#BEAT status MCCS
RSA_outputs_trial_1<-RSA_outputs_trial%>%filter(type== "Hybrid_BH_0.4_0.5")
sum(RSA_outputs_trial_1$beat_status == FALSE, na.rm = TRUE)
#6% of the time is worse

RSA_outputs_trial_1<-RSA_outputs_trial%>%filter(type== "Whole_BH_0.1")
sum(RSA_outputs_trial_1$beat_status == FALSE, na.rm = TRUE)
#21% of the time is worse

RSA_outputs_trial_1<-RSA_outputs_trial%>%filter(type== "Whole_BH_0.2")
sum(RSA_outputs_trial_1$beat_status == FALSE, na.rm = TRUE)
#9% of the time is worse

RSA_outputs_trial_1<-RSA_outputs_trial%>%filter(type== "Whole_BH_0.3")
sum(RSA_outputs_trial_1$beat_status == FALSE, na.rm = TRUE)
#9% of the time is worse

RSA_outputs_trial_1<-RSA_outputs_trial%>%filter(type== "Whole_BH_0.4")
sum(RSA_outputs_trial_1$beat_status == FALSE, na.rm = TRUE)
#16% of the time is worse

RSA_outputs_trial_1<-RSA_outputs_trial%>%filter(type== "Whole_BH_1")
sum(RSA_outputs_trial_1$beat_status == FALSE, na.rm = TRUE)
#57% of the time is worse

RSA_outputs_trial_1<-RSA_outputs_trial%>%filter(type== "Whole_BH_10")
sum(RSA_outputs_trial_1$beat_status == FALSE, na.rm = TRUE)
#100% of the time is worse
