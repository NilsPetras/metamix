# test the accuracy of the mixed model estimation for trivial cases

# assume ncp = 0 for all studies and n1 = n2 out of c(30, 300)
# without any selective publishing (can only overfit, but how much?)

nrep_global <- 100
design <- data.frame(
  ID = 1:4,
  n_studies = c(30, 30, 300, 300),
  n_per_group = c(30, 300, 30, 300)
)

design_long <- design[rep(seq_len(nrow(design)), nrep_global), ]
results <- design_long
results$p_sp <- NA
results$SE_p <- NA
results$ci_lower <- NA
results$ci_upper <- NA

for (i in seq_along(design_long$ID)) {
  temp <- metamix::metamix(
    seed = i,
    t = rt(n = design_long[i,"n_studies"],
           df = design_long[i,"n_per_group"] * 2 - 2,
           ncp = 0),
    n1 = design_long[i,"n_per_group"],
    n2 = design_long[i,"n_per_group"],
    nrep = 10
  )
  results[i,"p_sp"] <- temp$estimates$p_est
  results[i,"SE_p"] <- temp$estimates$SE_p
  results[i,"ci_lower"] <- temp$estimates$CI_p[1]
  results[i,"ci_upper"] <- temp$estimates$CI_p[2]
  
  if (i%%10 == 0) print(i)
}

tapply(results$p_sp, list(results$n_studies, results$n_per_group), mean)
tapply(results$ci_lower,
       list(results$n_studies, results$n_per_group),
       function(x) {
         mean(x > 0.01, na.rm = TRUE)
       })



library(ggplot2)

p_sp <- ggplot(results, aes(x = p_sp)) +
  geom_histogram() +
  facet_grid(list("n_per_group", "n_studies"), labeller = label_both)

save(p_sp, results, file = "bias_in_p_sp_estimates.RData")
