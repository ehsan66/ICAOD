# Finding a robust design for the two-parameter logistic model
# See how we set a stopping rule.
# The ELB is computed every checkfreq = 30 iterations
# The optimization stops when the ELB is larger than stoptol = .95
res <- robust(formula = ~1/(1 + exp(-b *(x - a))),
              predvars = c("x"), parvars = c("a", "b"),
              family = binomial(),
              lx = -5, ux = 5, prob = rep(1/4, 4),
              parset = matrix(c(0.5, 1.5, 0.5, 1.5, 4.0, 4.0, 5.0, 5.0), 4, 2),
              iter = 1, k =3,
              ICA.control = list(stop_rule = "equivalence",
                                 stoptol = .95, checkfreq = 30))

\dontrun{
res <- iterate(res, 100)
# stops at iteration 51
}
