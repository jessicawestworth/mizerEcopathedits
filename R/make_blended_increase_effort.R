#' Create a Blended Fishing Mortality Function for Slowly adjusting fishing effort
#' @description
#' This "Function Factory" returns a customized mortality function for use within
#' a \code{mizer} simulation. It creates a linear transition from
#' traditional gear-based fishing mortality to changed effort traditional
#' gear-based fishing
#'
#' @details
#' ## How it Works
#' The function uses a blending parameter \eqn{\alpha} that is tied to the
#' simulation time (\eqn{t}).
#' \itemize{
#'   \item At \eqn{t = 0}, \eqn{\alpha = 0} and fishing effort is 0
#'   \item As \eqn{t} approaches \eqn{t_{max\_sim}}, \eqn{\alpha} moves toward
#'   \code{alpha_max}, phasing in the desired effort rate.
#'   \item The transition is calculated as: \eqn{F = (1-\alpha)f_{standard} + \alpha f_{ssBH}}.
#' }
#'
#' ## Implementation Note
#' This function must be registered using \code{setRateFunction()} before
#' running \code{project()}. Because it is a closure, it "captures" the
#' \code{t_max_blend}, and \code{alpha_max} values.
#'
#' @param t_max_blend Numeric. The total duration of the simulation (in years).
#' This defines the "speed" of the transition.
#' @param t_steady Numeric. The time it takes for the simulation to reach steady
#'  state with normal fishing employed
#' @param alpha_max Numeric desired effort value.
#'
#' @return A function compatible with \code{mizer}'s \code{FMort} rate slot.
#' @export


make_blended_increase_effort <- function(t_max_blend, alpha_max, t_steady){
    function(params, n, n_pp, n_other, t, effort, ...) {
    # Calculate the 'Standard' Mizer Mortality (Gear-based)
    # We assume a base effort of 1 for the underlying gears for this calculation
    f_standard <- mizerFMort(params, n = n, n_pp = n_pp, n_other = n_other,
                             t = t, effort = 1, ...)

    #Check that t_steady is less than t_max_blend
    if(t_max_blend <= t_steady) stop("t_max_blend must be greater
                                         than t_steady")

    # 3. Blend them based on the 'effort' parameter passed to project()
    # 'effort' here acts as our transition alpha (0 to 1)
    if(t<t_max_blend){
        alpha <- alpha_max*((t-t_steady)/(t_max_blend-t_steady))
    }

    if(t>=t_max_blend){
        alpha<-alpha_max
    }

    if(t>t_steady){
        f_combined <- f_standard * alpha
    }
    else {
        f_combined <- f_standard
    }

    return(f_combined)
}
}
