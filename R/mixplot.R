#' Metamix plot
#'
#' Create an overview plot of the model-implied and empirical distribution of
#' t-values.
#'
#' @param x object of the class "metamix" as created by the `metamix` function
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' # values from Table 2 in Ulrich et al. (2018)
#' n1_t2 <- c(20, 30, 35, 25, 60, 40, 45, 35, 70, 65,
#'            30, 40, 30, 20, 90, 70, 50, 65, 50, 45)
#' n2_t2 <- c(30, 35, 35, 20, 50, 40, 50, 30, 80, 60,
#'            25, 40, 25, 20, 80, 75, 50, 70, 50, 55)
#' t_t2 <- c(.59, 1.84, 1.72, -.4, .1, .48, 1.17, 2.18, .15, 1.76,
#'           1.9, .5, 1.27, 1.71, 1.17, 1.91, .08, -.48, .98, 2.28)
#' m_t2 <- metamix(t_t2, n1_t2, n2_t2)
#' mixplot(m_t2)
mixplot <- function(x){
  
  if(!any(class(x) == "metamix")) stop("x is not a metamix object")
  
  # abbreviate
  dat <- data.frame(
    t = stats::na.omit(c(x$theoretical_distribution$t_values)),
    published = stats::na.omit(c(x$theoretical_distribution$published)))
  
  # adjust the second y-axis to overlap the theoretical and empirical distributions
  # the log adjusts for the difference in binwidth (cf. usage below)
  scale_factor <- sum(dat$published) *
    log(length(x$data$t)) /
    length(x$data$t) /
    log(length(dat$t))
    
  
  # plot creation
  gg <- ggplot2::ggplot(
    data = dat,
    mapping = ggplot2::aes(
      x = t)) +
    ggplot2::geom_vline( # natural y-axis at x = 0
      xintercept = 0,
      linetype = "dashed",
      color = "gray") +
    ggplot2::geom_area( # background: all bootstrapped t-values
      ggplot2::aes(y = ggplot2::after_stat(count)), stat = "bin",
      bins = floor(10 * log(length(dat$t))),
      color = "#6c0016",
      fill = "#6c0016",
      alpha = .2) +
    ggplot2::geom_area( # overlay 1: bootstrapped published t-values
      ggplot2::aes(y = ggplot2::after_stat(count)), stat = "bin",
      data = dat[dat$published, ],
      bins = floor(10 * log(length(dat$t))),
      color = "#004c6c",
      fill = "#004c6c") +
    ggplot2::geom_histogram( # overlay 2: empirical distribution histogram outlines
      data = data.frame(
        t = rep(
          x$data$t,
          each = scale_factor)),
      bins = floor(10 * log(length(x$data$t))),
      fill = NA,
      color = "black") +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous( # double count axis
      name = "Simulation Count",
      sec.axis = ggplot2::sec_axis(
        ~. / scale_factor,
        name = "Data Count"
      ))
  
  return(gg)
}