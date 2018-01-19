############################
# Exponential growth model
############################
# Verifying optimailty of a locally D-optimal design
senslocally(formula = ~a + exp(-b*x),
            predvars = "x", parvars = c("a", "b"),
            x = c(.1, 1), w = c(.5, .5),
            lx = 0, ux = 1, inipars = c(1, 10))


##############################
# A model with two predictors
##############################
x0 <- c(30, 3.861406, 30, 4.600633, 0, 0, 5.111376, 4.168798)
w0 <- rep(.25, 4)
senslocally(formula =  ~ V*S/(Km * (1 + I/Kic)+ S * (1 + I/Kiu)),
            predvars = c("S", "I"),
            parvars = c("V", "Km", "Kic", "Kiu"),
            x = x0, w = w0,
            lx = c(0, 0), ux = c(30, 60),
            inipars = c(1.5, 5.2, 3.4, 5.6))

# using package rgl for 3d plot:
res<- senslocally(formula =  ~ V*S/(Km * (1 + I/Kic)+ S * (1 + I/Kiu)),
                  predvars = c("S", "I"),
                  parvars = c("V", "Km", "Kic", "Kiu"),
                  x = x0, w = w0,
                  lx = c(0, 0), ux = c(30, 60),
                  inipars = c(1.5, 5.2, 3.4, 5.6),
                  plot_3d = "rgl")

