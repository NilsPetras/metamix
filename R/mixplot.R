mixplot <- function(x){
  
  if(!any(class(x) == "metamix")) stop("x is not a metamix object")
  
  dat <- data.frame(
    t = x$theoretical_distribution$t_values,
    published = x$theoretical_distribution$published)
  # p <- x$estimates$p_est
  coeff <- length(x$theoretical_distribution$t_values) /
    length(x$data$t) *
    mean(x$theoretical_distribution$published) ^ 2
    
  
  
  gg <- ggplot2::ggplot(
    data = dat,
    mapping = ggplot2::aes(
      x = t)) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "gray") +
    ggplot2::geom_area(
      ggplot2::aes(y = ..count..), stat = "bin",
      bins = floor(10 * log(length(dat$t))),
      color = "#6c0016",
      fill = "#6c0016",
      alpha = .2) +
    ggplot2::geom_area(
      ggplot2::aes(y = ..count..), stat = "bin",
      data = dat[dat$published, ],
      bins = floor(10 * log(length(dat$t))),
      color = "#004c6c",
      fill = "#004c6c") +
    ggplot2::theme_classic() +
    ggplot2::geom_histogram(
      data = data.frame(
        t = rep(
          x$data$t,
          each = coeff)),
      bins = floor(10 * log(length(x$data$t))),
      fill = NA,
      color = "black") +
    ggplot2::scale_y_continuous(
      name = "simulation count",
      sec.axis = ggplot2::sec_axis(
        ~. / coeff,
        name = "data count"
      ))
  
  return(gg)
}