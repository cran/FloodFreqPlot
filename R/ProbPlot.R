#' Flood Probability Plotting
#'
#' \code{ProbPlot} checks that a probability distribution fits a set of flood data. 
#'
#' This is a function for frequency analysis by a graphical method. The flood data are plotted on 
#'    an appropriate probability paper that linearizes the cumulative distribution function. Then the
#'    plotted flood data are fitted with a straight line for interpolation and extrapolation purposes.
#'    If \code{probs = NULL}, then a Weibull plotting position formula is used to calculate probability
#'    values for quantiles. If \code{PP = NULL}, then a Weibull plotting position formula is used to
#'    calculate the probabilities corresponding to the quantiles. If \code{dist = NULL}, then Normal
#'    distribution is used as the default frequency distribution. It should be noted that the distribution
#'    parameters are estimated by Method Of Moments (MOM). If \code{beta_CL = NULL}, then the
#'    confidence level is considered equal to 0.95 (that means the significance level is equal to 1-0.95=0.05).
#'
#' @param data_obs A vector, data frame or matrix containing observed data or flood quantiles.
#' @param probs Optional. The vector of plotting position probability values corresponding to
#'    the quantiles. If \code{probs = NULL}, then a Weibull plotting position formula is used to calculate
#'    probability values for quantiles.
#' @param PP Optional. A character string that represents the plotting position formula used to 
#'    calculate the empirical probability. The formula can be chosen from the list: \code{"Blom"},
#'    \code{"California_1"}, \code{"California_2"}, \code{"Chegodayev"}, \code{"Gringorten"}, 
#'    \code{"Hazen"}, \code{"Tukey"}, and \code{"Weibull"}. If \code{PP = NULL}, then
#'    \code{PP = 'Weibull'}.
#' @param dist Optional. A string that represents CDF and it can be 'Norm' for Normal distribution,
#'    'LNorm' for Log-Normal distribution, 'Gumb' for Gumbel distribution, 'Pea3' for Pearson type III
#'    distribution, and 'LPea3' for Log-Pearson type III distribution. If \code{dist = NULL}, then
#'    \code{dist = 'Norm'}.
#' @param T_rp Optional. A numeric vector including the return periods of interest for the flood
#'    quantile estimation.
#' @param beta_CL Optional. A numeric scalar that represents the confidence level for calculating 
#'    and plotting the confidence limits (bounds). If \code{beta_CL = NULL},
#'    then \code{beta_CL = 0.95}. It means that the significance level is equal to \code{0.05}.
#' @param T_lim Optional. A two-member numeric vector including the lower and upper return period
#'    limits determining the horizontal (x) axis range.
#' @param Q_lim Optional. A two-member numeric vector including the lower and upper limits
#'    determining the vertical (y) axis range to show quantile values.
#' @param main_title Optional. A character string representing the main title of the plot. The default title
#'    denotes the name of the theoretical probability distribution chosen to fit the data.
#' @param x_lab Optional. A character string representing the label of horizontal axis. The default label of
#'    the axis is \code{F(x) = P(X <= x)}.
#' @param y_lab Optional. A character string representing the label of vertical axis. The default label of
#'    the axis is \code{"Quantile"}.
#' @param Pcol Optional. A specification for the observed flood quantile points color. Defaults to
#'    \code{"black"}.
#' @param Ppch Optional. Either an integer specifying a symbol or a single character to be used as
#'    the default in plotting observed flood quantile points. See \code{points} for possible values and
#'    their interpretation. Defaults to 1. 
#' @param Pcex Optional. A numerical value giving the amount by which plotting point symbols should 
#'    be magnified relative to the default. Defaults to 1.
#' @param Lcol Optional. A specification for the theoretical probability line color. Defaults to \code{"blue"}.
#' @param Lty Optional. The theoretical probability line type. Line types can either be specified as an
#'    integer (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash) 
#'    or as one of the character strings \code{"blank"}, \code{"solid"}, \code{"dashed"}, \code{"dotted"},
#'    \code{"dotdash"}, \code{"longdash"}, or \code{"twodash"}, where \code{"blank"} uses
#'    'invisible lines' (i.e., does not draw them). Defaults to 1. 
#' @param Lwd Optional. The theoretical probability line width, a positive number, defaulting to 1.5.
#' @param CPlot Logical. If \code{CPlot = TRUE}, the confidence limits (bounds) are plotted. Defaults to \code{TRUE}.
#' @param CLcol Optional. A specification for the confidence limits (bounds) color. Defaults to \code{"red"}.
#' @param CLty Optional. The confidence limits (bounds) line type. Line types can either be specified as an
#'    integer (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash) 
#'    or as one of the character strings \code{"blank"}, \code{"solid"}, \code{"dashed"}, \code{"dotted"},
#'    \code{"dotdash"}, \code{"longdash"}, or \code{"twodash"}, where \code{"blank"} uses
#'    'invisible lines' (i.e., does not draw them). Defaults to 1. 
#' @param CLwd Optional. The confidence limits (bounds) line width, a positive number, defaulting to 1.5.
#' @param QTcol Optional. A specification for the T-year flood quantile estimate point color. Defaults to
#'    \code{"green"}.
#' @param QTpch Optional. Either an integer specifying a symbol or a single character to be used as
#'    the default in plotting the T-year flood quantile estimate points. See \code{points} for possible
#'    values and their interpretation. Defaults to 15. 
#' @param QTcex Optional. A numerical value giving the amount by which the T-year flood quantile
#'    estimate point symbols should be magnified. Defaults to 1.5.
#' @param GumbRV Logical. If \code{dist = 'Gumb'} and \code{GumbRV = 'TRUE'}, an extra horizontal
#'     axis is plotted to show Reduced Variable values.
#' @param P3SkewCheck Logical. If \code{P3SkewCheck = 'TRUE'} (default), the skewness of data is checked and
#'     if the coefficient of skewness is greater than 2.5, the confidence limits are not plotted for some data in the
#'     left tail of the dataset.
#'
#' @return The function returns a graph including the plotted flood data and the fitted distribution
#'    and the confidence limits (bounds). Also, it returns and shows the flood quantile estimates 
#'    corresponding to the return period(s) \code{T_rp}.
#'
#' @examples
#' # First Example
#' data('Harricana')
#' ProbPlot(data_obs = Harricana, PP = 'Cunnane', dist = 'LPea3', T_rp = c(100, 1000))
#'
#' # Second Example
#' data('AH_Tab12_1_1')
#' ProbPlot(data_obs = AH_Tab12_1_1, PP = 'Weibull', dist = 'Gumb', T_rp = 250, T_lim = c(2, 1000))
#'
#' @importFrom graphics abline axis grid legend lines mtext par points segments title
#' @importFrom stats ks.test pgamma pnorm qgamma qnorm qt sd
#' @importFrom grDevices dev.new
#'
#' @seealso \code{\link{PlotPos}} for the plotting position probability.
#' @export
ProbPlot = function(data_obs, probs = NULL, PP = NULL, dist = NULL, T_rp = NULL, beta_CL = NULL,
                     T_lim = NULL, Q_lim = NULL, main_title = NULL, x_lab = NULL, y_lab = NULL,
                     Pcol = 'black', Ppch = 1, Pcex = 1, Lcol = 'blue', Lty = 1, Lwd = 1.5,
                     CPlot = TRUE, CLcol = 'red', CLty = 2, CLwd = 1.5,
                     QTcol = 'green', QTpch = 15, QTcex = 1.5, GumbRV = FALSE, P3SkewCheck = TRUE){

##############################################################################
 if(!is.vector(data_obs)){
  n = dim(data_obs)[1]
  x = data_obs[1:n, 1]}else{
  n = length(data_obs)
  x = data_obs[1:n]}

 x_bar = mean(x)
 s_x = sd(x)
 y = log(x)
 y_bar = mean(y)
 s_y = sd(y)
 y10 = log10(x)
 y10_bar = mean(y10)
 s_y10 = sd(y10)

 if(is.null(T_rp)){T_r = 5}else{T_r = T_rp}
 n_T = length(T_r)

 p_vec = sort(c(0.0001, 0.001, 0.01, 0.05, seq(0.1, 0.9, 0.1), 0.95, 0.99, 0.999, 0.9999))

 squants = sort(x)
 lsquants = log(squants)

 if(is.null(x_lab)){x_lab = expression(paste(F(x) == P(X <= x),' (%)'))}
 if(is.null(y_lab)){y_lab = 'Quantile'} 

 if(is.null(PP)){PP = 'Weibull'} 

 if(is.null(probs)){
  probs = PlotPos(data_obs = x, PP = PP)
 }

 if(is.null(beta_CL)){
  beta_CL = 0.95}

 if(is.null(dist)){
  dist = 'Norm'}

##############################################################################
 if(dist == 'Norm'){
  dev.new(noRStudioGD = TRUE)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mar = c(5, 5, 8, 1))

  FXx = pnorm(squants, x_bar, s_x)
  q = qnorm(p_vec, x_bar, s_x)
  xq = qnorm(probs, x_bar, s_x)
  pT = 1 - (1/T_r)
  qT = qnorm(pT, x_bar, s_x)

  if(is.null(T_rp)){
   QuantsB = xq
   z_sn = qnorm(probs)}else{
   QuantsB = sort(c(xq, qT))
   z_sn = qnorm(sort(c(probs, pT)))
  }

  S_e = s_x * sqrt((2 + z_sn^2) / n)

  Bounds = ConfInt(quant = QuantsB, n_samp = n, SE = S_e, ConfLev = beta_CL)

  mn = min(squants, xq, qT)
  MX = max(xq, qT)

  if(is.null(T_lim)){
   x_lim = c(mn, MX)}else{
   qT_lim = qnorm((1 - (1/T_lim)), x_bar, s_x)
   if((min(T_lim) >= min(round(1/(1 - p_vec), digits = 2))) & (max(T_lim) <= max(round(1/(1 - p_vec), digits = 2)))){
    x_lim = c(min(qT_lim), max(qT_lim))}else{
    x_lim = c(mn, MX)
    warning('Inappropriate return period bounds')}
  }

 if(is.null(Q_lim)){
   y_lim = c(min(x, xq, qT), max(x, xq, qT))}else{
   y_lim = c(min(Q_lim), max(Q_lim))
  }

  plot(x = NA, y = NA, xaxt = 'n',
       xlab = x_lab, ylab = y_lab,
       xlim = x_lim,
       ylim = y_lim)

  if(is.null(main_title)){main_title = 'Normal Probability Plot'}

  title(main = main_title, line = 6)
  axis(1, at = q, labels = p_vec * 100)
  axis(3, at = q, labels = round(1/(1 - p_vec), digits = 2))
  mtext('Return Period (years)', side = 3, line = 3)

  grid(nx = NA, ny = NULL)
  for(i in 1:length(p_vec)){
   abline(v = q[i], untf = FALSE, lty = 3)}

  points(xq, squants, pch = Ppch, col = Pcol, cex = Pcex)
  if(!is.null(T_rp)){
   points(qT, qT, pch = QTpch, col = QTcol, cex = QTcex)
   segments(x0 = -10^6, y0 = qT, x1 = qT, y1 = qT, col = QTcol)
  }

  LX = sort(c(xq, qT))
  LY = sort(c(xq, qT))

  lines(LX, LY, col = Lcol, lty = Lty, lwd = Lwd)
 
  if(CPlot == TRUE){
   lines(QuantsB, Bounds[ , 1], col = CLcol, lty = CLty, lwd = CLwd)
   lines(QuantsB, Bounds[ , 2], col = CLcol, lty = CLty, lwd = CLwd)
  }

  legend('bottomright', legend = c('Observed Data', 'Theoretical Distribution', paste('Confidence Limits (', beta_CL * 100, '%)', sep = ''),
   ifelse(!is.null(T_rp), 'Quantile Estimate', NA)),
   col = c(Pcol, Lcol, CLcol, ifelse(!is.null(T_rp), QTcol, NA)),
   lty = c(0, Lty, CLty, ifelse(!is.null(T_rp), NA, NA)),
   pch =c(Ppch, NA, NA, ifelse(!is.null(T_rp), QTpch, NA)),
   bg = 'white')

  KStest = ks.test(squants, qnorm(probs, x_bar, s_x))

  if(!is.null(T_rp)){Q_T = qT; names(Q_T) = T_rp}

  return(list(GOF = KStest, if(!is.null(T_rp)){Q_Tr = Q_T}))}

##############################################################################
 if(dist == 'LNorm'){
  dev.new(noRStudioGD = TRUE)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mar = c(5, 5, 8, 1))

  FYy = pnorm(lsquants, y_bar, s_y)
  q = qnorm(p_vec, y_bar, s_y)
  yq = qnorm(probs, y_bar, s_y)
  pT = 1 - (1/T_r)
  qTlog = qnorm(pT, y_bar, s_y)

  if(is.null(T_rp)){
   QuantsB = yq
   z_sn = qnorm(probs)}else{
   QuantsB = sort(c(yq, qTlog))
   z_sn = qnorm(sort(c(probs, pT)))
  }

  S_e = s_y * sqrt((2 + z_sn^2) / n)
 
  Bounds = ConfInt(quant = QuantsB, n_samp = n, SE = S_e, ConfLev = beta_CL)

  mn = min(lsquants, yq, qTlog)
  MX = max(yq, qTlog)

  if(is.null(T_lim)){
   x_lim = c(mn, MX)}else{
   qTlog_lim = qnorm((1 - (1/T_lim)), y_bar, s_y)
   if((min(T_lim) >= min(round(1/(1 - p_vec), digits = 2))) & (max(T_lim) <= max(round(1/(1 - p_vec), digits = 2)))){
    x_lim = c(min(qTlog_lim), max(qTlog_lim))}else{
    x_lim = c(mn, MX)
    warning('Inappropriate return period bounds')}
  }

 if(is.null(Q_lim)){
   y_lim = c(min(x, exp(qTlog)), max(x, exp(yq), exp(qTlog)))}else{
   y_lim = c(min(Q_lim), max(Q_lim))
  }

  plot(x = NA, y = NA, xaxt = 'n',
       xlab = x_lab, ylab = y_lab,
       xlim = x_lim,
       ylim = y_lim, log = 'y')

  if(is.null(main_title)){main_title = 'Log-Normal Probability Plot'}

  title(main = main_title, line = 6)
  axis(1, at = q, labels = p_vec * 100)
  axis(3, at = q, labels = round(1/(1-p_vec), digits = 2))
  mtext('Return Period (years)', side = 3, line = 3)

  grid(nx = NA, ny = NULL) 
  for(i in 1:length(p_vec)){
   abline(v = q[i], untf = FALSE, lty = 3)}

  points(qnorm(probs, y_bar, s_y), squants, pch = Ppch, col = Pcol, cex = Pcex)
  if(!is.null(T_rp)){
   points(qTlog, exp(qTlog), pch = QTpch, col = QTcol, cex = QTcex)
   segments(x0 = -10^6, y0 = exp(qTlog), x1 = qTlog, y1 = exp(qTlog), col = QTcol)
  }

  LX = sort(c(yq, qTlog))
  LY = sort(c(exp(yq), exp(qTlog)))

  lines(LX, LY, col = Lcol, lty = Lty, lwd = Lwd)

  if(CPlot == TRUE){
   lines(QuantsB, exp(Bounds[ , 1]), col = CLcol, lty = CLty, lwd = CLwd)
   lines(QuantsB, exp(Bounds[ , 2]), col = CLcol, lty = CLty, lwd = CLwd)
  }

  legend('bottomright', legend = c('Observed Data', 'Theoretical Distribution', paste('Confidence Limits (', beta_CL * 100, '%)', sep = ''),
   ifelse(!is.null(T_rp), 'Quantile Estimate', NA)),
   col = c(Pcol, Lcol, CLcol, ifelse(!is.null(T_rp), QTcol, NA)),
   lty = c(0, Lty, CLty, ifelse(!is.null(T_rp), NA, NA)),
   pch =c(Ppch, NA, NA, ifelse(!is.null(T_rp), QTpch, NA)),
   bg = 'white')

  KStest = ks.test(lsquants, qnorm(probs, y_bar, s_y))

  if(!is.null(T_rp)){Q_T = exp(qTlog); names(Q_T) = T_rp}

  return(list(GOF = KStest, if(!is.null(T_rp)){Q_Tr = Q_T}))}

##############################################################################
 if(dist == 'Gumb'){
  dev.new(noRStudioGD = TRUE)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mar = c(5, 5, 8, 1))

  yrv = -log(-log(p_vec))
  yrvx = -log(-log(probs))
  yT = -log(-log(1 - (1/T_r)))

  ahat = 0.7797 * s_x
  x0 = x_bar - (0.45 * s_x)
  q = x0 + (ahat * yrv)
  xq = x0 + (ahat * yrvx)
  pT = 1 - (1/T_r)
  qTgum = x0 + (ahat * yT)

  if(is.null(T_rp)){
   QuantsB = xq
   pgum = probs
   yrvB = yrvx}else{
   QuantsB = sort(c(xq, qTgum))
   pgum = sort(c(probs, pT))
   yrvB = sort(c(yrvx, yT))
  }

  K_T = (sqrt(6)/pi) * (0.5772 + log(log(1/pgum)))
  S_e = (s_x/sqrt(n)) * sqrt(1 + 1.1396 * K_T + 1.10 * K_T^2)

  Bounds = ConfInt(quant = QuantsB, n_samp = n, SE = S_e, ConfLev = beta_CL)

  mn = min(yrvx, yT)
  MX = max(yrvx, yT)

  if(is.null(T_lim)){
   x_lim = c(mn, MX)}else{
   yT_lim = -log(-log(1 - (1/T_lim)))
   qTgum_lim = x0 + (ahat * yT_lim)
   if((min(T_lim) >= min(round(1/(1 - p_vec), digits = 2))) & (max(T_lim) <= max(round(1/(1 - p_vec), digits = 2)))){
    x_lim = c(min(yT_lim), max(yT_lim))}else{
    x_lim = c(mn, MX)
    warning('Inappropriate return period bounds')}
  }

 if(is.null(Q_lim)){
   y_lim = c(min(x, xq, qTgum), max(x, xq, qTgum))}else{
   y_lim = c(min(Q_lim), max(Q_lim))
  }

  plot(x = NA, y = NA, xaxt = 'n',
       xlab = x_lab, ylab = y_lab,
       xlim = x_lim, ylim = y_lim)

  if(is.null(main_title)){main_title = 'Gumbel Probability Plot'}

  title(main = main_title, line = 6)
  axis(1, at = yrv, labels = p_vec * 100)
  if(GumbRV == TRUE){
   axis(1, at = yrv[which(yrv <= (max(yrvx, yT) + 1))],
    labels = round(yrv[which(yrv <= (max(yrvx, yT) + 1))], digits = 2), line = 6)
   mtext('Reduced Variable', side = 1, line = 9)
  }

  axis(3, at = yrv, labels = round(1/(1 - p_vec), digits = 2))
  mtext('Return Period (years)', side = 3, line = 3)

  grid(nx = NA, ny = NULL) 
  for(i in 1:length(p_vec)){
   abline(v = yrv[i], untf = FALSE, lty = 3)}

  FXx = exp(-exp(-((squants - x0)/ahat)))
  FXxq = exp(-exp(-((xq - x0)/ahat)))
  FXqT = exp(-exp(-((qTgum - x0)/ahat)))

  points(yrvx, squants, pch = Ppch, col = Pcol, cex = Pcex)

  if(!is.null(T_rp)){
   points(yT, qTgum, pch = QTpch, col = QTcol, cex = QTcex)
   segments(x0 = -10^6, y0 = qTgum, x1 = yT, y1 = qTgum, col = QTcol)
  }

  LX = sort(-log(-log(c(FXxq, FXqT))))
  LY = sort(c(xq, qTgum))

  lines(LX, LY, col = Lcol, lty = Lty, lwd = Lwd)

  if(CPlot == TRUE){
   lines(yrvB, Bounds[ , 1], col = CLcol, lty = CLty, lwd = CLwd)
   lines(yrvB, Bounds[ , 2], col = CLcol, lty = CLty, lwd = CLwd)
  }

  legend('bottomright', legend = c('Observed Data', 'Theoretical Distribution', paste('Confidence Limits (', beta_CL * 100, '%)', sep = ''),
   ifelse(!is.null(T_rp), 'Quantile Estimate', NA)),
   col = c(Pcol, Lcol, CLcol, ifelse(!is.null(T_rp), QTcol, NA)),
   lty = c(0, Lty, CLty, ifelse(!is.null(T_rp), NA, NA)),
   pch =c(Ppch, NA, NA, ifelse(!is.null(T_rp), QTpch, NA)),
   bg = 'white')

  KStest = ks.test(squants, x0 + (ahat * yrv))

  if(!is.null(T_rp)){Q_T = qTgum; names(Q_T) = T_rp}

  return(list(GOF = KStest, if(!is.null(T_rp)){Q_Tr = Q_T}))}

##############################################################################
 if(dist == 'Pea3'){
  dev.new(noRStudioGD = TRUE)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mar = c(5, 5, 8, 1))
  
  if(n < 30){
   C_s = (n * sum((x - x_bar)^3))/((n - 1) * (n - 2) * s_x^3)}else{
   C_s = sum((x - x_bar)^3)/(n * s_x^3)}

  alpha = 4/(C_s^2)
  beta = (s_x * C_s) / 2
  tau =  sign(C_s) * (x_bar - 2 * (s_x / C_s))

  FXx = rep(0, n)
  qFXx = rep(0, n)
  for(i in 1:n){FXx[i] = cdfP3(squants[i], mu_p = x_bar, sigma_p = s_x, gamma_p = C_s)
   qFXx[i] = quaP3(P = FXx[i], mu_p = x_bar, sigma_p = s_x, gamma_p = C_s)
  }

  q = rep(0, length(p_vec))
  for(i in 1:length(p_vec)){
   q[i] = quaP3(P = p_vec[i], mu_p = x_bar, sigma_p = s_x, gamma_p = C_s)
  }

  xq = rep(0, n)
  for(i in 1:n){
   xq[i] = quaP3(P = probs[i], mu_p = x_bar, sigma_p = s_x, gamma_p = C_s)
  }

  pT = 1 - 1/T_r
  qT = rep(0, n_T)
  for(i in 1:n_T){
   qT[i] = quaP3(P = pT[i], mu_p = x_bar, sigma_p = s_x, gamma_p = C_s)
  }

  p2 = 1 - 1/2
  q2 = quaP3(P = p2, mu_p = x_bar, sigma_p = s_x, gamma_p = C_s)

  if(is.null(T_rp)){
   QuantsB = xq
   z_p = qnorm(probs)}else{
   QuantsB = sort(c(xq, qT))
   z_p = qnorm(sort(c(probs, pT)))
  }

  K_p_prime = (1/6) * (z_p^2 - 1) + (1/9) * (z_p^3 - 6 * z_p) * C_s/6 - (1/2) * (z_p^2 - 1) * (C_s/6)^2 + (2 * z_p / 3) * (C_s/6)^3 - (5/18) * (C_s/6)^4
  K_p = (2/C_s) * ((z_p - C_s/6) * (C_s/6) + 1)^3 - (2/C_s)
  delta2 = 1 + C_s * K_p + (1/2) * (1 + (3/4) * C_s^2) * K_p^2 + 6 * (1 + (1/4) * C_s^2) * K_p_prime * ((1/2) * C_s * K_p + (1 + (5/4) * C_s^2) * K_p_prime)
  delta = sqrt(delta2)
  S_e = s_x * delta / sqrt(n)

  Bounds = ConfInt(quant = QuantsB, n_samp = n, SE = S_e, ConfLev = beta_CL)

  mn = min(squants, xq, qT)
  MX = max(xq, qT)

  if(is.null(T_lim)){
   x_lim = c(mn, MX)}else{
   qT_lim =   rep(0, 2)
   qT_lim[1] = quaP3(P = 1 - 1/T_lim[1], mu_p = x_bar, sigma_p = s_x, gamma_p = C_s)
   qT_lim[2] = quaP3(P = 1 - 1/T_lim[2], mu_p = x_bar, sigma_p = s_x, gamma_p = C_s)
   if((min(T_lim) >= min(round(1/(1 - p_vec), digits = 2))) & (max(T_lim) <= max(round(1/(1 - p_vec), digits = 2)))){
    x_lim = c(min(qT_lim), max(qT_lim))}else{
    x_lim = c(mn, MX)
    warning('Inappropriate return period bounds')}
  }

 if(is.null(Q_lim)){
   y_lim = c(min(x, xq, qT), max(x, xq, qT))}else{
   y_lim = c(min(Q_lim), max(Q_lim))
  }

  plot(x = NA, y = NA, xaxt = 'n',
       xlab = x_lab, ylab = y_lab,
       xlim = x_lim,
       ylim = y_lim)

  if(is.null(main_title)){main_title = 'Pearson type III Probability Plot'}

  title(main = main_title, line = 6)
  axis(1, at = q, labels = p_vec * 100)
  axis(3, at = q, labels = round(1/(1 - p_vec), digits = 2))
  mtext('Return Period (years)', side = 3, line = 3)

  grid(nx = NA, ny = NULL) 
  for(i in 1:length(p_vec)){
   abline(v = q[i], untf = FALSE, lty = 3)}

  points(xq, squants, pch = Ppch, col = Pcol, cex = Pcex)

  if(!is.null(T_rp)){
   points(qT, qT, pch = QTpch, col = QTcol, cex = QTcex)
   segments(x0 = -10^6, y0 = qT, x1 = qT, y1 = qT, col = QTcol)
  }

  LX = sort(c(xq, qT))
  LY = sort(c(xq, qT))

  lines(LX, LY, col = Lcol, lty = Lty, lwd = Lwd)

  if(CPlot == TRUE){
   if((P3SkewCheck == TRUE) & (C_s > 2.5)){
    warning('High skewness in the left tail of dataset!')
    lines(QuantsB[which.min(Bounds[ , 2] - Bounds[ , 1]):length(QuantsB)], Bounds[which.min(Bounds[ , 2] - Bounds[ , 1]):length(QuantsB), 1], col = CLcol, lty = CLty, lwd = CLwd)
    lines(QuantsB[which.min(Bounds[ , 2] - Bounds[ , 1]):length(QuantsB)], Bounds[which.min(Bounds[ , 2] - Bounds[ , 1]):length(QuantsB), 2], col = CLcol, lty = CLty, lwd = CLwd)
   }else{
    lines(QuantsB, Bounds[ , 1], col = CLcol, lty = CLty, lwd = CLwd)
    lines(QuantsB, Bounds[ , 2], col = CLcol, lty = CLty, lwd = CLwd)
   }
  }

  legend('bottomright', legend = c('Observed Data', 'Theoretical Distribution', paste('Confidence Limits (', beta_CL * 100, '%)', sep = ''),
   ifelse(!is.null(T_rp), 'Quantile Estimate', NA)),
   col = c(Pcol, Lcol, CLcol, ifelse(!is.null(T_rp), QTcol, NA)),
   lty = c(0, Lty, CLty, ifelse(!is.null(T_rp), NA, NA)),
   pch =c(Ppch, NA, NA, ifelse(!is.null(T_rp), QTpch, NA)),
   bg = 'white')

  KStest = ks.test(squants, xq)

  if(!is.null(T_rp)){Q_T = qT; names(Q_T) = T_rp}

  return(list(GOF = KStest, if(!is.null(T_rp)){Q_Tr = Q_T}))}

##############################################################################
 if(dist == 'LPea3'){
  dev.new(noRStudioGD = TRUE)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mar = c(5, 5, 8, 1))

  if(n < 30){
   C_sy = (n * sum((y - y_bar)^3))/((n - 1) * (n - 2) * s_y^3)}else{
   C_sy = sum((y - y_bar)^3)/(n * s_y^3)}

  alpha = 4/(C_sy^2)
  beta = (s_y * C_sy) / 2
  tau =  sign(C_sy) * (y_bar - 2 * (s_y / C_sy))

  FYy = rep(0, n)
  qFYy = rep(0, n)
  for(i in 1:n){FYy[i] = cdflP3(lsquants[i], mu_p = y_bar, sigma_p = s_y, gamma_p = C_sy)
   qFYy[i] = qualP3(P = FYy[i], mu_p = y_bar, sigma_p = s_y, gamma_p = C_sy)
  }

  q = rep(0, length(p_vec))
  for(i in 1:length(p_vec)){
   q[i] = qualP3(P = p_vec[i], mu_p = y_bar, sigma_p = s_y, gamma_p = C_sy)
  }

  yq = rep(0, n)
  for(i in 1:n){
   yq[i] = qualP3(P = probs[i], mu_p = y_bar, sigma_p = s_y, gamma_p = C_sy)
  }

  pT = 1 - 1/T_r
  qTlP3 = rep(0, n_T)
  for(i in 1:n_T){
   qTlP3[i] = qualP3(P = pT[i], mu_p = y_bar, sigma_p = s_y, gamma_p = C_sy)
  }

  if(is.null(T_rp)){
   QuantsB = yq
   z_p = qnorm(probs)}else{
   QuantsB = sort(c(yq, qTlP3))
   z_p = qnorm(sort(c(probs, pT)))
  }

  K_p_prime = (1/6) * (z_p^2 - 1) + (1/9) * (z_p^3 - 6 * z_p) * C_sy/6 - (1/2) * (z_p^2 - 1) * (C_sy/6)^2 + (2 * z_p / 3) * (C_sy/6)^3 - (5/18) * (C_sy/6)^4
  K_p = (2/C_sy) * ((z_p - C_sy/6) * (C_sy/6) + 1)^3 - (2/C_sy)
  delta2 = 1 + C_sy * K_p + (1/2) * (1 + (3/4) * C_sy^2) * K_p^2 + 6 * (1 + (1/4) * C_sy^2) * K_p_prime * ((1/2) * C_sy * K_p + (1 + (5/4) * C_sy^2) * K_p_prime)
  delta = sqrt(delta2)
  S_e = s_y * delta / sqrt(n)

  Bounds = ConfInt(quant = QuantsB, n_samp = n, SE = S_e, ConfLev = beta_CL)

  mn = min(yq, qTlP3)
  MX = max(yq, qTlP3)

  if(is.null(T_lim)){
   x_lim = c(mn, MX)}else{
   qTlP3_lim =   rep(0, 2)
   qTlP3_lim[1] = qualP3(P = 1 - 1/T_lim[1], mu_p = y_bar, sigma_p = s_y, gamma_p = C_sy)
   qTlP3_lim[2] = qualP3(P = 1 - 1/T_lim[2], mu_p = y_bar, sigma_p = s_y, gamma_p = C_sy)
   if((min(T_lim) >= min(round(1/(1 - p_vec), digits = 2))) & (max(T_lim) <= max(round(1/(1 - p_vec), digits = 2)))){
    x_lim = c(min(qTlP3_lim), max(qTlP3_lim))}else{
    x_lim = c(mn, MX)
    warning('Inappropriate return period bounds')}
  }

 if(is.null(Q_lim)){
   y_lim = c(min(x, exp(yq), exp(qTlP3)), max(x, exp(yq), exp(qTlP3)))}else{
   y_lim = c(min(Q_lim), max(Q_lim))
  }

  plot(x = NA, y = NA, xaxt = 'n',
       xlab = x_lab, ylab = y_lab,
       xlim = x_lim,
       ylim = y_lim, log = 'y')

  if(is.null(main_title)){main_title = 'Log-Pearson type III Probability Plot'}

  title(main = main_title, line = 6)
  axis(1, at = q, labels = p_vec * 100)
  axis(3, at = q, labels = round(1/(1 - p_vec), digits = 2))
  mtext('Return Period (years)', side = 3, line = 3)

  grid(nx = NA, ny = NULL) 
  for(i in 1:length(p_vec)){
   abline(v = q[i], untf = FALSE, lty = 3)}

  points(yq, squants, pch = Ppch, col = Pcol, cex = Pcex)

  if(!is.null(T_rp)){
   points(qTlP3, exp(qTlP3), pch = QTpch, col = QTcol, cex = QTcex)
   segments(x0 = -10^6, y0 = exp(qTlP3), x1 = qTlP3, y1 = exp(qTlP3), col = QTcol)
  }

  LX = sort(c(yq, qTlP3))#qFYy
  LY = sort(c(exp(yq), exp(qTlP3)))

  lines(LX, LY, col = Lcol, lty = Lty, lwd = Lwd)

  if(CPlot == TRUE){
   lines(QuantsB, exp(Bounds[ , 1]), col = CLcol, lty = CLty, lwd = CLwd)
   lines(QuantsB, exp(Bounds[ , 2]), col = CLcol, lty = CLty, lwd = CLwd)
  }

  legend('bottomright', legend = c('Observed Data', 'Theoretical Distribution', paste('Confidence Limits (', beta_CL * 100, '%)', sep = ''),
   ifelse(!is.null(T_rp), 'Quantile Estimate', NA)),
   col = c(Pcol, Lcol, CLcol, ifelse(!is.null(T_rp), QTcol, NA)),
   lty = c(0, Lty, CLty, ifelse(!is.null(T_rp), NA, NA)),
   pch =c(Ppch, NA, NA, ifelse(!is.null(T_rp), QTpch, NA)),
   bg = 'white')

  KStest = ks.test(lsquants, yq)

  if(!is.null(T_rp)){Q_T = exp(qTlP3); names(Q_T) = T_rp}

  return(list(GOF = KStest, if(!is.null(T_rp)){Q_Tr = Q_T}))}
}