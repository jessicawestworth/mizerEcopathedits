#Plot Log P, Log B and determine slope
#20 cm for the length cutoff for all species, due to trawl survey reliability.
#Check the P vs B relationship and see if it changes as simulations are run
# Explore possible c and alpha values
# Find the values which keep the yields the same and look at the changes in the other metrics
# Also test the changes in metrics across a range of values and identify tradeoffs
# Metrics: Yield, large fish above 50 cm
# Review Hill Numbers
# Need to be able to plot the production values of species overtime


#For the function we should be able to compute the log P log B relationship dynamically for anytime in the simulation or for a params object
#For the function we should be able to compute the log P log Y relationship dynamically for anytime in the simulation or for a params object
#Selectivity should stay the same but the logP-log Y relationship should change dynamically


#' Parse the `log` argument of a time-series plot
#'
#' Mirrors the `log` argument of base [graphics::plot()]: a string naming the
#' axes to put on a log10 scale, e.g. `"x"`, `"y"`, `"xy"` or `""`. When
#' supplied it overrides `log_x` and `log_y`.
#'
#' This reimplements a helper that older versions of mizer provided internally
#' as `parseTimePlotLog()` but which is not present in mizer 3.3.0.
#'
#' @param log Character string, or NULL to use `log_x`/`log_y`.
#' @param log_x,log_y Logical fallbacks used when `log` is NULL.
#' @return A list with logical elements `log_x` and `log_y`.
#' @noRd
parse_time_plot_log <- function(log = NULL, log_x = FALSE, log_y = FALSE) {
    if (is.null(log)) {
        return(list(log_x = log_x, log_y = log_y))
    }
    if (!is.character(log) || base::length(log) != 1) {
        stop("'log' must be a single character string such as \"x\", \"y\", ",
             "\"xy\" or \"\".")
    }
    chars <- strsplit(log, "")[[1]]
    if (!all(chars %in% c("x", "y"))) {
        stop("'log' may only contain the characters \"x\" and \"y\", but was ",
             deparse(log), ".")
    }
    list(log_x = "x" %in% chars, log_y = "y" %in% chars)
}

#' Plot log10(Biomass) vs log10(Production) Slope Over Time
#'
#' @export
plotBiomassProductionSlope <- function(object, sim2 = NULL, species = NULL,
                                       tlim = c(NA, NA), return_data = FALSE,
                                       size = NULL, length = NULL, ...) {
    UseMethod("plotBiomassProductionSlope")
}

#' @rdname plotBiomassProductionSlope
#' @export
plotBiomassProductionSlope.MizerSim <- function(object, sim2 = NULL,
                                                species = NULL,
                                                tlim = c(NA, NA),
                                                return_data = FALSE,
                                                size = NULL, length = NULL,
                                                ...) {
    assert_that(is(object, "MizerSim"), is.flag(return_data))

    species <- mizer:::valid_species_arg(object, species, error_on_empty = TRUE)

    calc_slope_time <- function(sim_obj) {
        precalc <- precompute_cutoff(sim_obj@params, size = size, length = length)
        times <- as.numeric(dimnames(sim_obj@n)$time)
        sp_names <- dimnames(sim_obj@n)$sp
        sp_indices <- which(sp_names %in% species)

        slopes <- numeric(length(times))

        for (i in seq_along(times)) {
            t_i <- times[i]
            n_i <- sim_obj@n[i, , ]
            n_pp_i <- sim_obj@n_pp[i, ]
            n_other_i <- sim_obj@n_other[i, ]

            # Growth [Species x Size]
            G_i <- getEGrowth(sim_obj@params, n = n_i, n_pp = n_pp_i, n_other = n_other_i, t = t_i)

            # Production flux: Ps = N * G * dw
            Ps_i <- n_i * G_i * precalc$dw_mat
            Ps_i[!precalc$w_mask] <- 0
            Prod <- rowSums(Ps_i)[sp_indices]

            # Biomass: B = N * w * dw
            B_i <- n_i * precalc$w_mat * precalc$dw_mat
            B_i[!precalc$w_mask] <- 0
            Biomass <- rowSums(B_i)[sp_indices]

            # Filter valid non-zero pairs
            valid <- Biomass > 0 & Prod > 0
            if (sum(valid) >= 2) {
                # MATCHING YOUR plot_LogB_LogP FORMULA: log10(Biomass) ~ log10(Production)
                fit <- lm(log10(Biomass[valid]) ~ log10(Prod[valid]))
                slopes[i] <- coef(fit)[2]
            } else {
                slopes[i] <- NA
            }
        }

        dat <- data.frame(Year = times, Slope = slopes)

        if (!is.na(tlim[1])) dat <- subset(dat, Year >= tlim[1])
        if (!is.na(tlim[2])) dat <- subset(dat, Year <= tlim[2])

        return(dat)
    }

    if (is.null(sim2)) {
        plot_dat <- calc_slope_time(object)
        if (return_data) return(plot_dat)

        ggplot2::ggplot(plot_dat, ggplot2::aes(x = Year, y = Slope)) +
            ggplot2::geom_line(size = 1, color = "#0072B2") +
            ggplot2::labs(x = "Year", y = "log10(Biomass) vs log10(Production) Slope") +
            ggplot2::theme_bw()

    } else {
        d1 <- calc_slope_time(object)
        d2 <- calc_slope_time(sim2)

        d1$Simulation <- "Sim 1"
        d2$Simulation <- "Sim 2"
        plot_dat <- rbind(d1, d2)

        if (return_data) return(plot_dat)

        ggplot2::ggplot(plot_dat, ggplot2::aes(x = Year, y = Slope, color = Simulation)) +
            ggplot2::geom_line(size = 1) +
            ggplot2::labs(x = "Year", y = "log10(Biomass) vs log10(Production) Slope") +
            ggplot2::theme_bw()
    }
}


#' Plot Production over time
#'
#' @param object A \code{MizerSim} object.
#' @param sim2 Optional second \code{MizerSim} object to compare.
#' @param species Species to include. Defaults to all valid species.
#' @param total Logical. If TRUE, includes total ecosystem production across species.
#' @param log_x Logical. Log scale x axis? Defaults to FALSE.
#' @param log_y Logical. Log scale y axis? Defaults to TRUE.
#' @param log Character. Legacy log setting string.
#' @param ylim Numeric vector of length 2 for y-axis limits.
#' @param tlim Numeric vector of length 2 for time limits.
#' @param highlight Species name to highlight.
#' @param return_data Logical. If TRUE, returns data frame instead of ggplot.
#' @param size Optional size cut-off vector passed to precompute_cutoff.
#' @param length Optional length cut-off vector passed to precompute_cutoff.
#' @param ... Additional arguments.
#'
#' @export
plotProduction <- function(object, sim2 = NULL, species = NULL, total = FALSE,
                           log_x = FALSE, log_y = TRUE, log = NULL,
                           ylim = c(NA, NA), tlim = c(NA, NA), highlight = NULL,
                           return_data = FALSE, size = NULL, length = NULL, ...) {
    UseMethod("plotProduction")
}

#' @rdname plotProduction
#' @export
plotProduction.MizerSim <- function(object, sim2 = NULL,
                                    species = NULL,
                                    total = FALSE,
                                    log_x = FALSE, log_y = TRUE, log = NULL,
                                    ylim = c(NA, NA), tlim = c(NA, NA),
                                    highlight = NULL, return_data = FALSE,
                                    size = NULL, length = NULL,
                                    ...) {
    log_axes <- parse_time_plot_log(log, log_x = log_x, log_y = log_y)

    assert_that(is(object, "MizerSim"),
                is.flag(total),
                is.flag(log_axes$log_x),
                is.flag(log_axes$log_y),
                is.flag(return_data))

    params <- object@params
    species <- mizer:::valid_species_arg(object, species, error_on_empty = TRUE)

    if (is.null(sim2)) {
        # Precompute size/length masks and dw_mat once
        precalc <- precompute_cutoff(params, size = size, length = length)

        times <- as.numeric(dimnames(object@n)$time)
        sp_names <- dimnames(object@n)$sp
        num_times <- length(times)
        num_sp <- length(sp_names)

        # Matrix to hold total production per species per time step [Time x Species]
        prod_mat <- matrix(0, nrow = num_times, ncol = num_sp,
                           dimnames = list(Year = times, Species = sp_names))

        # Calculate production at each saved time step
        for (i in seq_along(times)) {
            t_i <- times[i]
            n_i <- object@n[i, , ]
            n_pp_i <- object@n_pp[i, ]
            n_other_i <- object@n_other[i, ]

            # Growth rate [Species x Size]
            G_i <- getEGrowth(params, n = n_i, n_pp = n_pp_i, n_other = n_other_i, t = t_i)

            # Production flux: Ps = N * G * dw
            Ps_i <- n_i * G_i * precalc$dw_mat
            Ps_i[!precalc$w_mask] <- 0

            # Sum production across size classes for each species
            prod_mat[i, ] <- rowSums(Ps_i)
        }

        # Filter time bounds (tlim)
        if (!is.na(tlim[1])) {
            prod_mat <- prod_mat[times >= tlim[1], , drop = FALSE]
            times <- as.numeric(rownames(prod_mat))
        }
        if (!is.na(tlim[2])) {
            prod_mat <- prod_mat[times <= tlim[2], , drop = FALSE]
        }

        # Calculate total across species if requested
        prod_total <- rowSums(prod_mat)

        # Subset requested species
        prod_mat <- prod_mat[, (as.character(dimnames(prod_mat)[[2]]) %in% species), drop = FALSE]

        if (total) {
            prod_mat <- cbind(prod_mat, "Total" = prod_total)
        }

        # Reshape data for ggplot using plotDataFrame
        plot_dat <- reshape2::melt(prod_mat, varnames = c("Year", "Species"),
                                   value.name = "Production")
        plot_dat <- subset(plot_dat, plot_dat$Production > 0)

        # Reorder columns to expected format: Year, Production, Species
        plot_dat <- plot_dat[, c(1, 3, 2)]

        if (nrow(plot_dat) == 0) {
            warning("There is no positive production to plot.")
        }
        if (return_data) return(plot_dat)

        mizer:::plotDataFrame(plot_dat, params,
                              ylab = "Production [g/year]",
                              xtrans = ifelse(log_axes$log_x, "log10", "identity"),
                              ytrans = ifelse(log_axes$log_y, "log10", "identity"),
                              ylim = ylim,
                              highlight = highlight)
    } else {
        # Compare two simulations side-by-side
        if (!all(dimnames(object@n)$time == dimnames(sim2@n)$time)) {
            stop("The two simulations do not have the same time steps")
        }
        pm1 <- plotProduction(object, species = species, tlim = tlim, total = total,
                              log_x = log_axes$log_x, log_y = log_axes$log_y,
                              ylim = ylim, highlight = highlight, size = size,
                              length = length, return_data = TRUE, ...)

        pm2 <- plotProduction(sim2, species = species, tlim = tlim, total = total,
                              log_x = log_axes$log_x, log_y = log_axes$log_y,
                              ylim = ylim, highlight = highlight, size = size,
                              length = length, return_data = TRUE, ...)

        pm1$Simulation <- rep(1, nrow(pm1))
        pm2$Simulation <- rep(2, nrow(pm2))
        pm <- rbind(pm1, pm2)

        if (return_data) return(pm)

        mizer:::plotDataFrame(pm, params,
                              ylab = "Production [g/year]",
                              xtrans = ifelse(log_axes$log_x, "log10", "identity"),
                              ytrans = ifelse(log_axes$log_y, "log10", "identity"),
                              ylim = ylim,
                              highlight = highlight, wrap_var = "Simulation")
    }
}


#' Plot fishing mortality
#'
#' Generic version of [mizer::plotFMort()]. The `MizerSim` method plots a
#' summary of the fishing mortality over time, which mizer's function does not
#' provide. Every other class, in particular `MizerParams`, is passed straight
#' through to [mizer::plotFMort()], so that attaching this package does not
#' change the behaviour of `plotFMort()` for existing code.
#'
#' @param object A `MizerSim` or `MizerParams` object.
#' @param sim2 Optional second `MizerSim` object to compare.
#' @param species Species to include. Defaults to all valid species.
#' @param total Logical. If TRUE, includes the total across species.
#' @param log_x,log_y Logical. Use log scales on the x and y axes?
#' @param log Character. Legacy log setting string, as in base `plot()`.
#' @param ylim Numeric vector of length 2 for y-axis limits.
#' @param tlim Numeric vector of length 2 for time limits.
#' @param highlight Species name to highlight.
#' @param return_data Logical. If TRUE, returns the data frame instead of a plot.
#' @param summary_fn Function used to reduce fishing mortality over sizes to one
#'   value per species and time. Defaults to `max`.
#' @param ... Additional arguments passed on to the method.
#'
#' @return A ggplot2 object, or a data frame if `return_data = TRUE`.
#' @export
plotFMort <- function(object, sim2 = NULL, species = NULL, total = FALSE,
                      log_x = FALSE, log_y = FALSE, log = NULL,
                      ylim = c(NA, NA), tlim = c(NA, NA), highlight = NULL,
                      return_data = FALSE, summary_fn = max, ...) {
    UseMethod("plotFMort")
}

#' @rdname plotFMort
#' @usage NULL
#' @export
plotFMort.default <- function(object, sim2 = NULL, species = NULL, total = FALSE,
                              log_x = FALSE, log_y = FALSE, log = NULL,
                              ylim = c(NA, NA), tlim = c(NA, NA),
                              highlight = NULL, return_data = FALSE,
                              summary_fn = max, ...) {
    mizer::plotFMort(object, species = species, return_data = return_data,
                     highlight = highlight, ...)
}

#' @rdname plotFMort
#' @usage NULL
#' @export
plotFMort.MizerSim <- function(object, sim2 = NULL,
                               species = NULL,
                               total = FALSE,
                               log_x = FALSE, log_y = FALSE, log = NULL,
                               ylim = c(NA, NA), tlim = c(NA, NA),
                               highlight = NULL, return_data = FALSE,
                               summary_fn = max,
                               ...) {
    log_axes <- parse_time_plot_log(log, log_x = log_x, log_y = log_y)

    assert_that(is(object, "MizerSim"),
                is.flag(total),
                is.flag(log_axes$log_x),
                is.flag(log_axes$log_y),
                is.flag(return_data))

    params <- object@params
    species <- mizer:::valid_species_arg(object, species, error_on_empty = TRUE)

    if (is.null(sim2)) {
        times <- as.numeric(dimnames(object@n)$time)
        sp_names <- dimnames(object@n)$sp
        num_times <- length(times)
        num_sp <- length(sp_names)
        num_w <- length(params@w)

        # Array to hold F across time x species x size
        f_mort_array <- array(0, dim = c(num_times, num_sp, num_w),
                              dimnames = list(time = times, sp = sp_names, w = params@w))

        # Check if a custom FMort function is registered in the params rate slots
        custom_fmort_name <- params@rates_funcs[["FMort"]]

        # If custom rate function exists, evaluate it dynamically at each saved time point
        if (!is.null(custom_fmort_name) && custom_fmort_name != "mizerFMort") {
            fmort_fn <- match.fun(custom_fmort_name)
            effort_mat <- object@effort

            for (i in seq_along(times)) {
                t_i <- times[i]
                n_i <- object@n[i, , ]
                n_pp_i <- object@n_pp[i, ]
                n_other_i <- object@n_other[i, ]

                # Extract effort for time step i if present
                eff_i <- if (!is.null(effort_mat) && nrow(effort_mat) >= i) effort_mat[i, ] else 1

                f_mort_array[i, , ] <- fmort_fn(params = params, n = n_i, n_pp = n_pp_i,
                                                n_other = n_other_i, t = t_i, effort = eff_i)
            }
        } else {
            # Standard mizer fallback
            f_mort_array <- getFMort(object, ...)
        }

        # Collapse size dimension using summary_fn (default: max F across sizes)
        y <- apply(f_mort_array, c(1, 2), summary_fn)

        times <- as.numeric(rownames(y))
        if (!is.na(tlim[1])) {
            y <- y[times >= tlim[1], , drop = FALSE]
            times <- as.numeric(rownames(y))
        }
        if (!is.na(tlim[2])) {
            y <- y[times <= tlim[2], , drop = FALSE]
        }

        y_total <- rowMeans(y)
        y <- y[, (as.character(dimnames(y)[[2]]) %in% species), drop = FALSE]

        if (total) {
            y <- cbind(y, "Total" = y_total)
        }

        plot_dat <- reshape2::melt(y, varnames = c("Year", "Species"),
                                   value.name = "FMort")
        plot_dat <- subset(plot_dat, plot_dat$FMort > 0)

        # plotDataFrame expects: Year, FMort, Species
        plot_dat <- plot_dat[, c(1, 3, 2)]

        if (nrow(plot_dat) == 0) {
            warning("There is no fishing mortality to include.")
        }
        if (return_data) return(plot_dat)

        mizer:::plotDataFrame(plot_dat, params,
                              ylab = "Fishing mortality [1/year]",
                              xtrans = ifelse(log_axes$log_x, "log10", "identity"),
                              ytrans = ifelse(log_axes$log_y, "log10", "identity"),
                              ylim = ylim,
                              highlight = highlight)
    } else {
        if (!all(dimnames(object@n)$time == dimnames(sim2@n)$time)) {
            stop("The two simulations do not have the same times")
        }
        ym <- plotFMort(object, species = species,
                        tlim = tlim, total = total,
                        log_x = log_axes$log_x, log_y = log_axes$log_y,
                        ylim = ylim, highlight = highlight,
                        return_data = TRUE, summary_fn = summary_fn, ...)

        ym2 <- plotFMort(sim2, species = species,
                         tlim = tlim, total = total,
                         log_x = log_axes$log_x, log_y = log_axes$log_y,
                         ylim = ylim, highlight = highlight,
                         return_data = TRUE, summary_fn = summary_fn, ...)

        ym$Simulation <- rep(1, nrow(ym))
        ym2$Simulation <- rep(2, nrow(ym2))
        ym <- rbind(ym, ym2)

        if (return_data) return(ym)

        mizer:::plotDataFrame(ym, params,
                              ylab = "Fishing mortality [1/year]",
                              xtrans = ifelse(log_axes$log_x, "log10", "identity"),
                              ytrans = ifelse(log_axes$log_y, "log10", "identity"),
                              ylim = ylim,
                              highlight = highlight, wrap_var = "Simulation")
    }
}

#' Helper function that precomputes the weights corresponding to the desired
#' size cutoff for the metrics such as Yield, Biomass and Production
#' Precompute the params for species-selective balanced harvest
#' @description This function precomputes static properties ONCE before running
#' simulations
#' @export
precompute_cutoff <- function(object, size = NULL, length = NULL,...) {
    if (is(object, "MizerSim")){
        params <- object@params
    } else if (is(object, "MizerParams")){
        params <- object
    } else {
        stop("Input 'object' must be of class 'MizerParams' or 'MizerSim'.")
    }

    sp_params <- species_params(params)
    gp <- gear_params(params)
    no_sp <- nrow(sp_params)

    # Initialize weight vector
    weight <- numeric(no_sp)

    # Helper: recycle a length-1 or length-no_sp vector to one value per species
    recycle_to_species <- function(x, what) {
        n <- base::length(x)
        if (n == 1) return(rep(x, no_sp))
        if (n == no_sp) return(x)
        stop("'", what, "' must have length 1 or match the number of species (",
             no_sp, "), but has length ", n, ".")
    }

    if (!is.null(size) && is.numeric(size)) {
        # Cutoff supplied directly as a weight in grams
        if (!is.null(length)) {
            stop("Supply either 'size' (weights in g) or 'length' (cm), not both.")
        }
        weight <- recycle_to_species(size, "size")

    } else if (!is.null(length) ||
               (!is.null(size) && identical(as.character(size), "length"))) {
        # Cutoff supplied as a length in cm; convert with the length-weight relation
        if (!is.null(size) && !identical(as.character(size), "length")) {
            stop("'size' must be numeric weights or the string \"length\", not ",
                 deparse(size), ".")
        }
        if (is.null(length)) {
            stop("The 'length' argument must be provided when size == \"length\".")
        }
        cm <- recycle_to_species(length, "length")
        weight <- sp_params$a * (cm ^ sp_params$b)

    } else if (!is.null(size)) {
        stop("'size' must be numeric weights in g or the string \"length\", not ",
             deparse(size), ".")

    } else {
        # Precompute weight per species using gear parameters
        for (s in seq_len(no_sp)) {
            sp_name <- sp_params$species[s]
            g_sp <- gp[gp$species == sp_name, ]

            if (nrow(g_sp) > 0 && sum(g_sp$catchability) > 0) {
                total_catchability <- sum(g_sp$catchability)
                weighting <- g_sp$catchability / total_catchability
                chosen_length <- sum(g_sp$l50 * weighting)
                weight[s] <- sp_params$a[s] * (chosen_length ^ sp_params$b[s])
            } else {
                weight[s] <- 0
            }
        }
    }

    # Create static matrix mask (Species x Size) where w >= weight
    w_vec <- w(params)
    dw_vec <- dw(params)

    w_mask <- outer(weight, w_vec, "<=")
    w_mat <- matrix(w_vec, nrow = no_sp, ncol = length(w_vec), byrow = TRUE)
    dw_mat <- matrix(dw_vec, nrow = no_sp, ncol = length(w_vec), byrow = TRUE)

    return(list(
        weight = weight,
        w_mask = w_mask,
        w_mat = w_mat,
        dw_mat = dw_mat
    ))
}

#' Species-specific production, yield and biomass beyond a given size cutoff
#' this cutoff should be specified in weight and can either be a
#' species-specific list or a single value for all species, if a size cutoff is
#' not supplied the cutoff will be taken as the l50 of the gears using a
#' weighted average by the catchabilities.
#' @export
Compute_Yield_Biomass_Production <- function(object,
                                             size = NULL,
                                             length = NULL,
                                             n = NULL,
                                             n_pp = NULL,
                                             n_other = NULL,
                                             effort = NULL,
                                             t = NULL,
                                             ...) {
    # Resolve MizerParams, time, and state inputs dynamically
    if (is(object, "MizerSim")) {
        params <- object@params

        # Use final time step by default if not specified
        if (is.null(t)) {
            idx <- idxFinalT(object)
            t_val <- as.numeric(dimnames(object@n)$time[idx])
        } else {
            t_val <- t
            idx <- which(as.numeric(dimnames(object@n)$time) == t_val)
            if (length(idx) == 0) stop("Specified time 't' not found in MizerSim object.")
        }

        # Extract states at time step 't'
        if (is.null(n)) n <- object@n[idx, , ]
        if (is.null(n_pp)) n_pp <- object@n_pp[idx, ]
        if (is.null(n_other)) {
            n_other <- if (length(object@n_other) > 0) object@n_other[idx, ] else initialNOther(params)
        }

        # Extract dynamic Fishing Mortality matrix directly from the simulation
        FMort <- getFMort(object, time_range = t_val)
        if (length(dim(FMort)) == 3) FMort <- FMort[1, , ]

    } else if (is(object, "MizerParams")) {
        params <- object
        t_val <- if (is.null(t)) 0 else t

        if (is.null(n)) n <- initialN(params)
        if (is.null(n_pp)) n_pp <- initialNResource(params)
        if (is.null(n_other)) n_other <- initialNOther(params)
        if (is.null(effort)) effort <- 1

        # Compute static params Fishing Mortality
        FMort <- getFMort(params, effort = effort, t = t_val)
    } else {
        stop("Input 'object' must be of class 'MizerParams' or 'MizerSim'.")
    }

    # Extract static threshold parameters
    precalc <- precompute_cutoff(params, size = size, length = length)

    # Compute single time step Growth [Species x Size]
    G <- getEGrowth(params, n = n, n_pp = n_pp, n_other = n_other, t = t_val)

    # Calculate Production (Ps), Yield (Y), and Biomass (B)
    Ps <- n * G * precalc$dw_mat
    Ps[!precalc$w_mask] <- 0
    Prod <- rowSums(Ps)

    Y <- n * precalc$dw_mat * precalc$w_mat * FMort
    Y[!precalc$w_mask] <- 0
    Yield <- rowSums(Y)

    B <- n * precalc$dw_mat * precalc$w_mat
    B[!precalc$w_mask] <- 0
    Biomass <- rowSums(B)

    # Return standard single time step data frame
    result <- data.frame(
        Species    = names(Biomass),
        Biomass    = as.numeric(Biomass),
        Production = as.numeric(Prod[names(Biomass)]),
        Yield      = as.numeric(Yield[names(Biomass)])
    )

    return(result)
}

#' Plot Log Y Log P relationship
#' @export
plot_LogY_LogP <- function(object, ...) {
    UseMethod("plot_LogY_LogP")
}

#' @rdname plot_LogY_LogP
#' @export
plot_LogY_LogP.MizerParams <- function(object, size = NULL, length = NULL,
                                       n = initialN(object),
                                       n_pp = initialNResource(object),
                                       n_other = initialNOther(object),
                                       ...) {
    params <- object
    params_name <- deparse(substitute(object))

    result <- Compute_Yield_Biomass_Production(params, size = size, length = length,
                                               n = n, n_pp = n_pp, n_other = n_other)

    result <- subset(result, Yield > 0 & Production > 0)

    model <- lm(log(Yield) ~ log(Production), data = result)
    intercept <- unname(coef(model)[1])
    slope <- unname(coef(model)[2])

    subtitle_text <- sprintf("Y = %.3f * P^%.3f", exp(intercept), slope)
    if (!is.null(length)) {
        subtitle_text <- paste(subtitle_text, paste0("Species Length >= ", length, " cm"), sep = "\n")
    }

    ggplot(result, aes(x = Production, y = Yield)) +
        stat_function(
            fun = function(x) exp(intercept) * x^slope,
            colour = "black",
            linewidth = 0.6,
            linetype = "dashed"
        ) +
        geom_point(aes(colour = Species)) +
        theme_cowplot(12) +
        scale_colour_manual(values = params@linecolour) +
        scale_x_log10() +
        scale_y_log10() +
        labs(
            title = paste("Yield-Production Relationship of", params_name),
            subtitle = subtitle_text,
            x = "Production",
            y = "Yield",
            colour = "Species"
        )
}


#' @rdname plot_LogY_LogP
#' @export
plot_LogY_LogP.MizerSim <- function(object, size = NULL, length = NULL, ...) {
    params <- object@params

    # Run Compute_Yield_Biomass_Production passing the MizerSim object directly
    result <- Compute_Yield_Biomass_Production(object, size = size, length = length, ...)

    result <- subset(result, Yield > 0 & Production > 0)

    model <- lm(log(Yield) ~ log(Production), data = result)
    intercept <- unname(coef(model)[1])
    slope <- unname(coef(model)[2])

    subtitle_text <- sprintf("Y = %.3f * P^%.3f", exp(intercept), slope)
    if (!is.null(length)) {
        subtitle_text <- paste(subtitle_text, paste0("Species Length >= ", length, " cm"), sep = "\n")
    }

    ggplot(result, aes(x = Production, y = Yield)) +
        stat_function(
            fun = function(x) exp(intercept) * x^slope,
            colour = "black",
            linewidth = 0.6,
            linetype = "dashed"
        ) +
        geom_point(aes(colour = Species)) +
        theme_cowplot(12) +
        scale_colour_manual(values = params@linecolour) +
        scale_x_log10() +
        scale_y_log10() +
        labs(
            title = "Yield-Production Relationship (Final Time Step)",
            subtitle = subtitle_text,
            x = "Production",
            y = "Yield",
            colour = "Species"
        )
}

#' Plot of the Biomass Production relationship
#' @export
plot_LogB_LogP <- function(object, ...) {
    UseMethod("plot_LogB_LogP")
}

#' @rdname plot_LogB_LogP
#' @export
plot_LogB_LogP.MizerParams <- function(object, size = NULL, length = NULL,
                                       n = initialN(object),
                                       n_pp = initialNResource(object),
                                       n_other = initialNOther(object), ...) {
    params <- object
    params_name <- deparse(substitute(object))

    # Compute production, biomass, and yields
    result <- Compute_Yield_Biomass_Production(params, size = size, length = length,
                                               n = n, n_pp = n_pp, n_other = n_other)

    # Filter out non-positive values
    result <- subset(result, Biomass > 0 & Production > 0)

    # Fit linear model
    model <- lm(log(Biomass) ~ log(Production), data = result)

    intercept <- unname(coef(model)[1])
    slope <- unname(coef(model)[2])

    subtitle_text <- sprintf("B = %.3f * P^%.3f", exp(intercept), slope)
    if (!is.null(length)) {
        subtitle_text <- paste(subtitle_text, paste0("Species Length >= ", length, " cm"), sep = "\n")
    }

    ggplot(result, aes(x = Production, y = Biomass)) +
        stat_function(
            fun = function(x) exp(intercept) * x^slope,
            colour = "black",
            linewidth = 0.6,
            linetype = "dashed"
        ) +
        geom_point(aes(colour = Species)) +
        theme_cowplot(12) +
        scale_colour_manual(values = params@linecolour) +
        scale_x_log10() +
        scale_y_log10() +
        labs(
            title = paste("Biomass-Production Relationship of", params_name),
            subtitle = subtitle_text,
            x = "Production",
            y = "Biomass",
            colour = "Species"
        )
}

#' @rdname plot_LogB_LogP
#' @export
plot_LogB_LogP.MizerSim <- function(object, size = NULL, length = NULL, ...) {
    # Extract final time step indices
    idx <- idxFinalT(object)

    n_final <- object@n[idx, , ]
    n_pp_final <- object@n_pp[idx, ]
    n_other_final <- if (length(object@n_other) > 0) object@n_other[idx, ] else initialNOther(object@params)

    plot_LogB_LogP.MizerParams(
        object = object@params,
        size = size,
        length = length,
        n = n_final,
        n_pp = n_pp_final,
        n_other = n_other_final,
        ...
    )
}

#'Function that will compute the Biomass Production slope relationship
#'may be used in future for detecting changes in this relationship as the
#'fishing mortality is adjusted
#'
#' The fit is \eqn{\log B = \log k + \alpha \log P} on natural logs, so the
#' returned `slope` is \eqn{\alpha} and the fitted proportionality constant is
#' \eqn{k = \exp(\mathrm{intercept})}. Note that \eqn{k} is *not* the harvest
#' rule constant \eqn{c}: see the "Species Level Balanced Harvest
#' Implementation" vignette.
#'
#' @return A list with the fitted `model`, its `intercept` and `slope`, and the
#'   per-species `data` the fit was made from.
#' @export
compute_LogP_LogB_slope<-function(params, size=NULL, length=NULL, n=initialN(params), n_pp=initialNResource(params), n_other=initialNOther(params)){
    result <- Compute_Yield_Biomass_Production(params, size = size, length = length,
                                               n = n, n_pp = n_pp, n_other = n_other)
    valid <- result$Biomass > 0 & result$Production > 0
    if (sum(valid) < 2) {
        stop("Need at least two species with positive biomass and production ",
             "above the size cutoff to fit a slope.")
    }
    model <- lm(log(Biomass) ~ log(Production), data = result[valid, ])
    list(
        model = model,
        intercept = unname(coef(model)[1]),
        slope = unname(coef(model)[2]),
        data = result
    )
}
