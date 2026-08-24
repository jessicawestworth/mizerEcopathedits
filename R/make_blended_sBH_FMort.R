#' Create a Blended Fishing Mortality Function for Balanced Harvesting Transitions
#' @description
#' This function returns a customized mortality function for use within
#' a \code{mizer} simulation. It creates a linear transition from
#' traditional gear-based fishing mortality to species-level Balanced
#' Harvesting (sBH) based on species production.
#'
#' @param params A \code{MizerParams} object.
#' @param size Optional size cut-offs passed to precompute_cutoff.
#' @param length Optional length cut-offs passed to precompute_cutoff.
#' @param t_max_blend Numeric. Time point at which maximum transition is reached.
#' @param target_c Numeric. Proportionality constant for Balanced Harvest (\eqn{F = c \cdot P}).
#' @param alpha_max Numeric (0 to 1). Maximum extent of transition (1.0 = total shift to sBH).
#' @param t_steady Numeric. Burn-in time before blending begins.
#' @param Type_Fmort Integer (1 or 2). Mode of calculating f_sBH:
#' \itemize{
#'   \item 1: Rescale f_standard by current yield/biomass ratio.
#'   \item 2: Rescale normalised gear selectivity x catchability profile directly by target production.
#' }
#'
#' @return A function compatible with \code{mizer}'s \code{FMort} rate slot.
#' @export
make_blended_sBH_FMort <- function(params, size = NULL, length = NULL, t_max_blend, target_c,
                                   alpha_max = 1, t_steady, Type_Fmort = 1) {
    if (t_max_blend <= t_steady) {
        stop("t_max_blend must be greater than t_steady")
    }

    denom <- t_max_blend - t_steady

    # Pre-calculate cutoff matrices once in the closure to save simulation runtime
    precalc <- precompute_cutoff(params, size = size, length = length)

    # For Type_Fmort == 2, compute the normalized selectivity profile ONCE in factory scope
    if (Type_Fmort == 2) {
        selectivity_array <- mizer::getSelectivity(params) # [Gear x Species x Size]
        catchability      <- mizer::getCatchability(params)  # [Gear] or [Gear x Species]

        num_sp   <- nrow(params@species_params)
        num_w    <- length(params@w)
        combined_selectivity <- matrix(0, nrow = num_sp, ncol = num_w)

        gears <- dimnames(selectivity_array)$gear

        for (gear in gears) {
            # Safely extract gear catchability vector across species
            if (is.matrix(catchability)) {
                c_gear <- catchability[gear, ]
            } else {
                c_gear <- catchability[gear]
            }

            gear_F_profile <- selectivity_array[gear, , ] * c_gear
            combined_selectivity <- combined_selectivity + gear_F_profile
        }

        # Normalize relative selectivity profile across size classes per species
        spec_sel_sum <- rowSums(combined_selectivity)
        spec_sel_sum[spec_sel_sum == 0] <- 1
        norm_selectivity <- combined_selectivity / spec_sel_sum
    }

    # Returned closure for mizer FMort rate slot
    function(params, n, n_pp, n_other, t, effort, ...) {
        # Calculate standard mortality matrix [Species x Size]
        f_standard <- mizerFMort(params, n = n, n_pp = n_pp, n_other = n_other,
                                 t = t, effort = effort)

        if (t <= t_steady) {
            return(f_standard)
        }

        # Compute growth & ecosystem production flux (Ps = N * G * dw)
        G <- getEGrowth(params, n = n, n_pp = n_pp, n_other = n_other, t = t)
        Ps <- n * G * precalc$dw_mat
        Ps[!precalc$w_mask] <- 0
        Prod <- rowSums(Ps) # Vector [num_species]

        target_production <- target_c * Prod # Target F flux vector

        if (Type_Fmort == 1) {
            # Biomass: B = N * dw * w
            B <- n * precalc$dw_mat * precalc$w_mat
            B[!precalc$w_mask] <- 0
            Biomass <- rowSums(B)

            # Yield: Y = N * dw * w * f_standard
            Y <- n * precalc$dw_mat * precalc$w_mat * f_standard
            Y[!precalc$w_mask] <- 0
            Yield <- rowSums(Y)

            f_standard_total <- Yield / Biomass
            f_standard_total[is.na(f_standard_total) | is.infinite(f_standard_total)] <- 0

            multiplier <- ifelse(f_standard_total > 0, target_production / f_standard_total, 0)

            # Scale f_standard by target/current ratio
            f_sBH <- f_standard * multiplier

        } else if (Type_Fmort == 2) {
            # Scale pre-computed normalized selectivity profile by target production vector
            f_sBH <- norm_selectivity * target_production
        }

        # Compute alpha transition rate
        if (t < t_max_blend) {
            alpha <- alpha_max * ((t - t_steady) / denom)
        } else {
            alpha <- alpha_max
        }

        # Return linearly blended mortality matrix
        return((1 - alpha) * f_standard + alpha * f_sBH)
    }
}




#P_a calculate production based fishing mortality as
#' @export
make_blended_P_y<- function(params, size = NULL, length = NULL,
                                   t_max_blend, target_c,
                                   alpha_max = 1, t_steady, Type_Fmort = 1,
                             y=1.8) {
    if (t_max_blend <= t_steady) {
        stop("t_max_blend must be greater than t_steady")
    }

    denom <- t_max_blend - t_steady

    # Pre-calculate cutoff matrices once in the closure to save simulation runtime
    precalc <- precompute_cutoff(params, size = size, length = length)

    # For Type_Fmort == 2, compute the normalized selectivity profile ONCE in factory scope
    if (Type_Fmort == 2) {
        selectivity_array <- mizer::getSelectivity(params) # [Gear x Species x Size]
        catchability      <- mizer::getCatchability(params)  # [Gear] or [Gear x Species]

        num_sp   <- nrow(params@species_params)
        num_w    <- length(params@w)
        combined_selectivity <- matrix(0, nrow = num_sp, ncol = num_w)

        gears <- dimnames(selectivity_array)$gear

        for (gear in gears) {
            # Safely extract gear catchability vector across species
            if (is.matrix(catchability)) {
                c_gear <- catchability[gear, ]
            } else {
                c_gear <- catchability[gear]
            }

            gear_F_profile <- selectivity_array[gear, , ] * c_gear
            combined_selectivity <- combined_selectivity + gear_F_profile
        }

        # Normalize relative selectivity profile across size classes per species
        spec_sel_sum <- rowSums(combined_selectivity)
        spec_sel_sum[spec_sel_sum == 0] <- 1
        norm_selectivity <- combined_selectivity / spec_sel_sum
    }

    # Returned closure for mizer FMort rate slot
    function(params, n, n_pp, n_other, t, effort, ...) {
        # Calculate standard mortality matrix [Species x Size]
        f_standard <- mizerFMort(params, n = n, n_pp = n_pp, n_other = n_other,
                                 t = t, effort = effort)

        if (t <= t_steady) {
            return(f_standard)
        }

        # Compute growth & ecosystem production flux (Ps = N * G * dw)
        G <- getEGrowth(params, n = n, n_pp = n_pp, n_other = n_other, t = t)
        Ps <- n * G * precalc$dw_mat
        Ps[!precalc$w_mask] <- 0
        Prod <- rowSums(Ps) # Vector [num_species]

        B <- n * precalc$dw_mat * precalc$w_mat
        B[!precalc$w_mask] <- 0
        Biomass <- rowSums(B)

        target_production <- target_c * (Prod^(y)/Biomass) # Target F flux vector

        if (Type_Fmort == 1) {

            # Yield: Y = N * dw * w * f_standard
            Y <- n * precalc$dw_mat * precalc$w_mat * f_standard
            Y[!precalc$w_mask] <- 0
            Yield <- rowSums(Y)

            f_standard_total <- Yield / Biomass
            f_standard_total[is.na(f_standard_total) | is.infinite(f_standard_total)] <- 0

            multiplier <- ifelse(f_standard_total > 0, target_production / f_standard_total, 0)

            # Scale f_standard by target/current ratio
            f_sBH <- f_standard * multiplier

        } else if (Type_Fmort == 2) {
            # Scale pre-computed normalized selectivity profile by target production vector
            f_sBH <- norm_selectivity * target_production
        }

        # Compute alpha transition rate
        if (t < t_max_blend) {
            alpha <- alpha_max * ((t - t_steady) / denom)
        } else {
            alpha <- alpha_max
        }

        # Return linearly blended mortality matrix
        return((1 - alpha) * f_standard + alpha * f_sBH)
    }
}


