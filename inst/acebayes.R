set.seed(1)
n <- 6
p <- 5
k <- 4
C <- 10
start.d <- list()
for(i in 1:C){
  start.d[[i]] <- randomLHS(n = n, k = k) * 2 - 1
  colnames(start.d[[i]]) <- c("x1", "x2", "x3", "x4") }

a1 <- c(-3, 4, 5, -6, -2.5)
a2 <- c(3, 10, 11, 0, 3.5)
prior <- list(support = rbind(a1, a2))
ex411 <- paceglm(formula = ~ x1 + x2 + x3 + x4, family = binomial,
                 method = "quadrature",
                 B = c(3, 8),
                 start.d = start.d, prior = prior, criterion = "D")

?paceglm
ex411$utility(ex411$d)
-10.15454 #default
-10.22591 8 8
