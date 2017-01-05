
PlotELB <- function(Iter ,
                         ELB,
                         plot_main = TRUE,
                         ...){
  # Iter: the iterations of ELB
  # ELB: a vector of Efficiency lower bounds

  cex.main <- .9
  prec <- 8

  if(plot_main)
    main1 <- paste("D-Efficiency Lower Bound (ELB): ", round(ELB[length(ELB)], prec), paste = "") else
      main1 = NULL

  plot(x = Iter, y = ELB ,
       xlim = c(Iter[1], Iter[length(Iter)]), xlab = "Iteration", ylab = "D-Efficiency Lower Bound (ELB)", type = "s",
       main = main1,
       cex.main = cex.main,
       xaxt = "n",...)


  if(Iter[length(Iter)] < 5)
    axis(side = 1, at = c(Iter),
         labels = c(Iter)) else{
           pretty_plot <- pretty(Iter)
           axis(side = 1, at = pretty_plot,
                labels = pretty_plot)
         }
}


