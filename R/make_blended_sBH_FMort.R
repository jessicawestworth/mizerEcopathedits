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
#'   \item 1: Keep the shape of the current gear-based F and rescale it so that
#'     the biomass-weighted mean F above the cutoff equals the target.
#'   \item 2: Rescale the peak-normalised gear selectivity x catchability
#'     profile so that fully-selected individuals experience the target F.
#' }
#'
#' @section Limitation of `Type_Fmort = 1`:
#' Because Type 1 rescales the existing fishing mortality, a species with zero
#' status-quo F stays unfished for the whole projection, and a species with a
#' very small status-quo F receives a very large multiplier. Balanced harvest as
#' a policy asks that every species be harvested in proportion to its
#' productivity regardless of what the current fleet does, so Type 1 cannot
#' express the rule for species the fleet does not currently target. Use
#' `Type_Fmort = 2` if that matters for the question being asked.
#'
#' @section Definition of production:
#' Production is the somatic production \eqn{\int N_i(w) g_i(w) dw} above the
#' size cutoff. This matches the Ecopath production of [getProduction()]
#' restricted to that size range, because the other component of the Ecopath
#' production — offspring biomass, see [getOffspringProduction()] — enters the
#' spectrum at `w_min`, far below any cutoff of interest.
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

        # Normalise each species' profile to a peak of 1, so that it is a
        # selectivity in [0, 1] and the target F is the mortality experienced by
        # fully-selected individuals. Normalising by rowSums() instead would make
        # the result depend on the number of size bins and would rescale each
        # species by a different, arbitrary factor.
        spec_sel_max <- apply(combined_selectivity, 1, max)
        spec_sel_max[spec_sel_max == 0] <- 1
        norm_selectivity <- combined_selectivity / spec_sel_max
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




#' Create a Blended Fishing Mortality Function for a Power-Law Harvest Rule
#'
#' @description
#' As [make_blended_sBH_FMort()], but implements the more general harvest rule
#' \deqn{F_i = c\, P_i^{\gamma} / B_i,} equivalently \eqn{Y_i = c\,P_i^{\gamma}},
#' where \eqn{P_i} and \eqn{B_i} are the production and biomass of species
#' \eqn{i} above the size cutoff. \eqn{\gamma = 1} is the fixed exploitation
#' ratio \eqn{Y_i = cP_i}.
#'
#' @section Relation to balanced harvest:
#' If the model satisfies \eqn{B_i = k P_i^{\alpha}} then setting
#' \eqn{\gamma = 1 + \alpha} gives \eqn{F_i = (c/k) P_i}, the balanced harvest
#' rule of [make_blended_sBH_FMort()]. That equivalence holds only at the state
#' where the \eqn{B}-\eqn{P} fit was made: \eqn{B_i} in the denominator here is
#' the *dynamic* biomass, so once the community responds to fishing the two
#' rules diverge. Note also that the fitted constant \eqn{k} is absorbed into
#' \eqn{c}, so `target_c` is not comparable across models or across values of
#' `gamma` unless it is recalibrated for each; see
#' [compute_LogP_LogB_slope()] and the "Species Level Balanced Harvest
#' Implementation" vignette.
#'
#' @inheritParams make_blended_sBH_FMort
#' @param gamma Numeric. Exponent on production in the harvest rule. `gamma = 1`
#'   is the fixed exploitation ratio; `gamma = 1 + alpha` reproduces balanced
#'   harvest at the reference state.
#' @param y Deprecated. Former name of `gamma`.
#'
#' @return A function compatible with \code{mizer}'s \code{FMort} rate slot.
#' @seealso [make_blended_sBH_FMort()]
#' @export
make_blended_P_y <- function(params, size = NULL, length = NULL,
                             t_max_blend, target_c,
                             alpha_max = 1, t_steady, Type_Fmort = 1,
                             gamma = 1.8, y = NULL) {
    if (!is.null(y)) {
        warning("The 'y' argument of make_blended_P_y() has been renamed to ",
                "'gamma'. Please use 'gamma' instead.")
        gamma <- y
    }
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

        # Normalise each species' profile to a peak of 1, so that it is a
        # selectivity in [0, 1] and the target F is the mortality experienced by
        # fully-selected individuals. Normalising by rowSums() instead would make
        # the result depend on the number of size bins and would rescale each
        # species by a different, arbitrary factor.
        spec_sel_max <- apply(combined_selectivity, 1, max)
        spec_sel_max[spec_sel_max == 0] <- 1
        norm_selectivity <- combined_selectivity / spec_sel_max
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

        target_production <- target_c * (Prod^gamma / Biomass) # Target F vector

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


