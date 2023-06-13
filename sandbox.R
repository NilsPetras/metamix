rm(list = ls())
devtools::load_all()

# Sigurdson et al. (2023) data
sigurdson <- read.csv(system.file("extdata", "Sigurdson2023.csv", package = "metamix"))

sigurdson <- sigurdson[1:33,c("Index", "Pln", "Trn", "t")]
n1_sigurdson <- sigurdson$Pln
n2_sigurdson <- sigurdson$Trn

# multiply t with correction factor for Hedge's g in small samples
sigurdson$tcorr <- sigurdson$t / (1 - 3 / (4 * (sigurdson$Pln + sigurdson$Trn - 2) - 1))
tcorr_sigurdson <- sigurdson$tcorr

m_sigurdson <- metamix(
  tcorr_sigurdson, n1_sigurdson, n2_sigurdson,
  nrep = 1e6)
# m_sigurdson$estimates
testplot <- mixplot(m_sigurdson) +
  ggplot2::xlim(-4, 4)  +
  # ggplot2::ylim(0, 8e4) +
  ggplot2::annotate(
    "text",
    x = 3,
    y = 50000,
    label = paste(
      "p[sp] == ",
      round(m_sigurdson$estimates$p_est, 3),
      sep = ""),
    parse = TRUE) +
  ggplot2::annotate(
    "text",
    x = 3,
    y = 40000,
    label = paste(
      "p[pub] == ",
      round(mean(m_sigurdson$theoretical_distribution$published), 3),
      sep = ""),
    parse = TRUE)

ggplot2::ggsave("testplot.pdf", testplot, width = 8, height = 5)

m_sigurdson_twosided <- metamix(
  tcorr_sigurdson,
  n1_sigurdson, n2_sigurdson,
  nrep = 1e5,
  TwoSided = TRUE)
testplot2 <- mixplot(m_sigurdson_twosided) +
  ggplot2::xlim(-4, 4)  +
  # ggplot2::ylim(0, 8e4) +
  ggplot2::annotate(
    "text",
    x = 3,
    y = 50000,
    label = paste(
      "p[sp] == ",
      round(m_sigurdson_twosided$estimates$p_est, 3),
      "(", round(m_sigurdson_twosided$estimates$CI_p[1], 3), ", ",
      round(m_sigurdson_twosided$estimates$CI_p[2], 3), ")",
      sep = ""),
    parse = TRUE) +
  ggplot2::annotate(
    "text",
    x = 3,
    y = 40000,
    label = paste(
      "p[pub] == ",
      round(mean(m_sigurdson_twosided$theoretical_distribution$published), 3),
      sep = ""),
    parse = TRUE)
ggplot2::ggsave("testplot2.pdf", testplot2, width = 8, height = 5)








# cross-validation of software implementation and paper

# Table 1 example in Ulrich et al. (2018) - all significant (replicates exactly)
if(TRUE) {
  n1_t1 <- c(34, 55, 28, 32, 88, 70, 33, 35, 81, 62, 40, 80, 25, 42, 18)
  n2_t1 <- c(41, 50, 19, 25, 90, 40, 45, 22, 93, 48, 50, 80, 38, 44, 22)
  t_t1 <- c(2.04, 2.21, 1.85, 2.01, 2.58, 2.61, 1.73, 1.79, 1.68, 1.93, 1.68, 1.94, 2.16, 1.82, 2.22)
  
  m_t1 <- metamix(t_t1, n1_t1, n2_t1, SigOnly = TRUE)
}

# Table 2 example in Ulrich et al. (2018) - mix (replicates except one minor rounding error)
if(FALSE) {n1_t2 <- c(20, 30, 35, 25, 60, 40, 45, 35, 70, 65, 30, 40, 30, 20, 90, 70, 50, 65, 50, 45)
n2_t2 <- c(30, 35, 35, 20, 50, 40, 50, 30, 80, 60, 25, 40, 25, 20, 80, 75, 50, 70, 50, 55)
t_t2 <- c(.59, 1.84, 1.72, -.4, .1, .48, 1.17, 2.18, .15, 1.76, 1.9, .5, 1.27, 1.71, 1.17, 1.91, .08, -.48, .98, 2.28)
m_t2 <- metamix(t = t_t2, n1 = n1_t2, n2 = n2_t2, nrep = 1e5)}
