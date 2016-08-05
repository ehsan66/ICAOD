Imperialist Competitive Algorithm (ICA) to find optimal designs for nonlinear models.

Example: locally D-optimal design for the exponential model
mica(fimfunc = "FIM_exp_2par", lx = 0, ux = 1, lp = c(2, 3), up = c(2, 3),
    iter = 40, k = 2, type = "locally", control = list(seed = 215))
    
    
How to isntall:

install.packages("ICAOD")
require(ICAOD)

The most important function is mica that finds locally, minimax and standardized maximin D-optimal design for nonlinear models.
on_average_ica also finds optim on the average optimal designs.
