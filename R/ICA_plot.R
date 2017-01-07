#  roxygen
# Plot ICA object
#'
#' Plot method for an object of class "ICA".
#' @param x an object of class "ICA".
#' @param iter iteration. if \code{NULL}, is set to the last iteration.
#' @param sensitivity logical; if \code{TRUE} the best design in iteration \code{iter} is checked by equivalence theorem.
#' @param print_equivalence logical; if \code{TRUE} the result of checking by equivalence theorem will be printed.
#' Only applicable when \code{sensitivity = TRUE}.
#' @param ... argument with no further use.
#' @seealso \code{\link{mica}}, \code{\link{multica_4pl}} and \code{\link{ave}}
#' @export

plot.ICA <- function(x, iter = NULL,  sensitivity = TRUE, print_equivalence = FALSE, ...){

  if (any(class(x) != c("list", "ICA")))
    stop("'x' must be of class 'ICA'")


  ## to not get confused with design points
  object <- x
  if(is.null(iter))
    totaliter <- length(object$evol) else
      totaliter <- iter

  if (totaliter > length(x$evol))
    stop("'iter' is larger than the maximum number of iterations")

  ################################################################################
  ## in this block we extract all the values and variables for plotting to use the same names as 'iterate.ICA'
  ## ectracting all min_cost and mean_cost up to 'totaliter'
  mean_cost <- sapply(1:totaliter, FUN = function(i)object$evol[[i]]$mean_cost)
  min_cost <- sapply(1:totaliter, FUN = function(i)object$evol[[i]]$min_cost)

  if( grepl("on_average",  object$arg$type))
    type = "optim_on_average" else
      type <- object$arg$type

  ################################################################################

  #################################################################################################
  ###################################################################################### plot cost

  ##############################################################
  # plot setting
  legend_place <- "topright"
  legend_text <- c( "Best Imperialist", "Mean of Imperialists")
  line_col <- c("firebrick3", "blue4")
  if (type == "minimax")
    title1 <- "cost value"
  if (type == "standardized")
    title1 <- "minimum efficiency"
  if (type == "locally" || type == "optim_on_average")
    title1 <- "log determinant of inverse of FIM"
  if (type == "multiple_locally")
    title1 <- "criterion value"

  ###############################################################

  ylim = switch(type,
                "minimax" = c(min(min_cost) - .07, max(mean_cost[1:(totaliter)]) + .2),
                "standardized" = c(min(mean_cost[1:(totaliter)]) - .07, max( min_cost) + .2),
                "locally" = c(min(min_cost) - .07, max(mean_cost[1:(totaliter)]) + .2),
                "multiple_locally" = c(min(min_cost) - .07, max(mean_cost[1:(totaliter)]) + .2))


  PlotEffCost(from = 1,
              to = (totaliter),
              AllCost = min_cost, ##all criterion up to now (all cost function)
              UserCost = NULL,
              DesignType = type,
              col = line_col[1],
              xlab = "Iteration",
              ylim = ylim,
              lty = 1,
              title1 = title1,
              plot_main = TRUE)
  ICAMean_line <-  mean_cost[1:(totaliter)]
  lines(x = 1:(totaliter),
        y = ICAMean_line,
        col = line_col[2], type = "s", lty = 5)
  legend(x = legend_place,  legend = legend_text,lty = c(1,5, 3), col = line_col, xpd = TRUE, bty = "n")

  #################################################################################################
  ###################################################################################### end of plot cost

  #################################################################################################
  ###################################################################################### checking equivalence theory
  if (sensitivity){
    #if (type == "locally" || type == "minimax" || type == "standardized")
    if (type == "locally" || type == "minimax" || type == "standardized")
      out <- equivalence(fimfunc = object$arg$FIM,
                         locally = object$arg$locally,
                         x = object$evol[[totaliter]]$x,
                         w = object$evol[[totaliter]]$w,
                         lp = object$arg$lp,
                         up = object$arg$up,
                         lx = object$arg$lx,
                         ux = object$arg$ux,
                         type = type,
                         n.seg = object$arg$control$n.seg,
                         maxeval_equivalence = object$arg$control$maxeval_equivalence)
    if (type == "optim_on_average")
      out <- equivalence_ave(fimfunc = object$arg$FIM,
                                    x = object$evol[[totaliter]]$x,
                                    w = object$evol[[totaliter]]$w,
                                    param = object$arg$param,
                                    prior = object$arg$prior,
                                    lx = object$arg$lx,
                                    ux = object$arg$ux,
                                    maxeval_equivalence = object$arg$control$maxeval_equivalence)

    if (type == "multiple_locally")
      out <- equivalence_multiple(x = object$evol[[totaliter]]$x,
                                  w = object$evol[[totaliter]]$w,
                                  param = object$arg$param,
                                  lambda = object$arg$lambda,
                                  delta = object$arg$delta,
                                  lx = object$arg$lx,
                                  ux = object$arg$ux,
                                  maxeval_equivalence = object$arg$control$maxeval_equivalence)
    #Psi_Mu = object$arg$Psi_Mu)

    out$iter <- totaliter
  }else out <- NULL
  #################################################################################################
  ###################################################################################### end of equivalence theory


  if(sensitivity && print_equivalence)
    return(out) else
      return(invisible(out))
}
