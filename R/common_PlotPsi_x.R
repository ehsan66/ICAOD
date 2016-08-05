#' @importFrom grDevices  rainbow
#' @importFrom graphics abline axis legend lines mtext plot points


PlotPsi_x <- function(x, w, lower, upper, Psi_x, FIM, answering, mu, plot_3d  = "lattice"){
  # plot the equivalanece theorem, for both one and two dimensional.
  # x: vector of the point of the design
  # w: vector of the weights of the design
  # lower: lower bound
  # psi_x psi_x as a function of x. in multiobjective optimal design ('multi_ICA'), 'PsiMulti_x' is 'given as Psi_x'
  # FIM: fisher information matrix function
  # AsnweringSet: answering set matrix, each row is one element of the A(\xi) or N(\xi)
  # vector of the measure. For locally optimal design mu = 1.
  # plot_3d: which package should be used to plot the 3d plots, 'rgl' or 'lattice'


  if(length(lower) == 1){
    xPlot <- seq(lower,upper, length.out = 1000)

    PsiPlot<- sapply(1:length(xPlot), FUN = function(j) Psi_x(x1 = xPlot[j], mu = mu, FIM = FIM, x = x,
                                                              w = w, answering = answering))
    plot(xPlot, PsiPlot, type = "l",
         col = "blue",   xlab = "Design Interval",
         ylab = "Sensitivity Function", xaxt = "n", main = "Sensitivity Plot")
    abline(h = 0, v = c(x) ,col = "grey", lty = 3)

    Point_y<- sapply(1:length(x), FUN = function(j) Psi_x(x1 = x[j], mu = mu, FIM = FIM, x = x,
                                                              w = w, answering = answering))

    points(x = x,  y = Point_y, col = "red" ,pch = 16, cex = 1)

    axis(side = 1, at = c(lower, x, upper, 0),
         labels = round(c(lower, x, upper, 0), 4))
  }

  if(length(lower) == 2){
    xlab_3d <- "S"
    ylab_3d <- "I"
    len <- 40
    xPlot <- seq(lower[1],upper[1], length.out = len)
    yPlot <- seq(lower[2],upper[2], length.out = len)
    ncol = 100
    color1 <- rev(rainbow(ncol, start = 0/6, end = 4/6))
    ## Psi_x here is Psi_xy
    Psi_xy1 <- Vectorize(FUN = Psi_x, vectorize.args=c("x1", "y1"))
    ###does not work, we can njot vectorize
    z <- -outer(X = xPlot,
                Y =yPlot,
                FUN = Psi_xy1,
                answering = answering,
                mu=mu,
                FIM = FIM,
                w = w,
                x = x)

    if (plot_3d  == "lattice"){
      if (requireNamespace("lattice", quietly = TRUE)){
        wireframe_dat <-expand.grid(x = xPlot, y = yPlot)
        wireframe_dat$z <- as.vector(z)
        p1 <- lattice::wireframe( z ~ x * y, data = wireframe_dat,
                                  col.regions=  color1,
                                  drap = TRUE,
                                  scales = list(arrows = FALSE),
                                  shade = FALSE,
                                  xlab = xlab_3d, ylab = ylab_3d,
                                  pretty = TRUE,
                                  zlab = "c",
                                  #zlab = expression(c(x, xi, mu)),
                                  #main = "Optimality Verification Plot",
                                  colorkey = list(col = color1, tick.number = 16))
        #lattice::print.trellis(p1)
        print(p1)
      }else {
        warning("Package 'lattice' is not installed in your system. The 3D derivation plot can not be plotted unless this packages is installed.")
      }
    }
    if(plot_3d  == "rgl"){

      if (requireNamespace("rgl", quietly = TRUE)){
        rgl::.check3d()



        zcol  = cut(z, ncol)

        rgl::persp3d(x = xPlot, y = yPlot, z = z, col = color1[zcol],
                     smooth = FALSE,
                     alpha = .8,
                     shininess    = 128,
                     xlab = xlab_3d, ylab = ylab_3d, zlab = "c(x, \\xi, \\mu)")


        Point_mat <- matrix(round(x, 3), ncol = length(lower), nrow = length(x)/length(lower))
        ##adding poin
        Point_y <- sapply(1:dim(Point_mat)[1], FUN = function(j) Psi_x(x1 = Point_mat[j, 1], y1 = Point_mat[j, 2], mu=mu, FIM=FIM, x=x,
                                                                       w=w, answering=answering))
        rgl::points3d(x = Point_mat[, 1],
                      y = Point_mat[, 2],
                      z = Point_y,
                      col = "darkred",
                      size = 11,
                      point_antialias = TRUE)

        text_point <- paste("(", xlab_3d, "=", Point_mat[, 1], ", ", ylab_3d, "=", Point_mat[, 2], "; c=",  Point_y, ")", sep ="")


        rgl::text3d(x = Point_mat[, 1],
                    y = Point_mat[, 2],
                    z = Point_y + .1,
                    texts = text_point,
                    col = "darkred",
                    font = 4)

        rgl::grid3d("z")
      }else{
        warning("Package 'rgl' is not installed in your system. The 3D derivation plot can not be plotted unless this packages is installed.")
      }
    }

    ###contour plot
    xyz <- expand.grid(x = xPlot, y = yPlot)
    xyz$z <- as.vector(z)

    if (requireNamespace("lattice", quietly = TRUE)){
      p2 <- lattice::contourplot( z ~ x * y, data = xyz,
                                 xlab = xlab_3d,
                                 ylab = ylab_3d,
                                 #main = "Optimality Verification Contour Plot",
                                 region = TRUE,
                                 col.regions = color1,
                                 colorkey = list(col = color1, tick.number = 16),
                                 cuts = 13)
      #lattice::print.trellis(p2)
      print(p2)
    }else{
      warning("Package 'lattice' is not installed in your system. The countor plot can not be plotted.")
    }
  }


}



