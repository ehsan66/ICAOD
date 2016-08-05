PlotEffCost <- function(from,
                        to ,
                        AllCost, ##all criterion up to now (all costfunction)
                        UserCost,
                        title1,
                        DesignType,
                        plot_main = TRUE,
                        ...){

  ## plot the evolution per each iteration
  # from: start iteration
  # to: last iteration
  # AllCost: all the cost values that should be plotted
  # title1: what title do u want? for example ICA, FWA or multiple
  # DesignType: 'standardized' or 'minimax'or 'locally'
  # UserCost: the cost value that user has. if not NULL, then the relative efficiency is plotted.
  cex.main <- .9
  prec <- 8
  if (is.null(UserCost)){
    ##cost function plot
    if (plot_main)
      main1 <- paste(title1, ": ", round( AllCost[length(AllCost)], prec), paste = "") else
        main1 = NULL

      plot(x = from:to, y = AllCost ,
           xlim = c(from, to),  ylab = "Imperialist Cost", type = "s",
           main = main1,
           cex.main = cex.main,
           xaxt = "n",...)

  }else{
    ##Efficiency
    Efficiency <- switch(DesignType, "minimax" = UserCost/AllCost, "standardized" = AllCost/UserCost, "locally" =  UserCost/AllCost)
    if (plot_main)
      main1 <- paste(title1, ": ", round(Efficiency[length(Efficiency)], prec), paste = "") else
        main1 = NULL

      plot(x = from:to, y = Efficiency ,
           xlim = c(from, to),  ylab = "Efficiency", type = "s",
           main = main1,
           cex.main = cex.main,
           xaxt = "n",...)

  }
  ## here we plot the axis to control everything.
  if (to < 5)
    axis(side = 1, at = c(from:to),
         labels = c(from:to)) else{
           pretty_plot <- pretty(from:to)
           axis(side = 1, at = pretty_plot,
                labels = pretty_plot)
         }
}




