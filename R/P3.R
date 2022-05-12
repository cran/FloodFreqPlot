pdfP3 = function(xp3, alpha, beta, tau){
 expon = ifelse((((xp3 - tau)/beta) < 0) & (((alpha - 1) %% 1) != 0), round(alpha - 1), alpha - 1)
 density = (1/(abs(beta) * gamma(alpha))) * ((xp3 - tau)/beta)^expon * exp(-(xp3 - tau)/beta)
 return(density)
}

###########################################################
#' @importFrom stats pnorm pgamma

cdfP3 = function(xp3, mu_p = x_bar, sigma_p = s_x, gamma_p = C_s){
 if(xp3 <= 0){PP3 = 0; return(PP3)}
 if(sigma_p <= 0){stop('invalid parameters')}
 if(gamma_p <= 1e-06){PP3 = pnorm(xp3, mu_p, sigma_p); return(PP3)}else{
  alpha = 4/gamma_p^2
  z = 2 * (xp3 - mu_p)/(sigma_p * gamma_p) + alpha
  PP3 = pgamma(pmax(0, z), alpha)
  if(gamma_p < 0){PP3 = 1 - PP3}
  return(PP3)
 }
}

###########################################################

P3P = function(xp3, P){cdfP3(xp3) - P}

###########################################################

#' @importFrom stats qnorm qgamma

quaP3 = function(P, mu_p = x_bar, sigma_p = s_x, gamma_p = C_s){
 if(P < 0){stop('invalid probability')}
 if(sigma_p <= 0){stop('invalid parameters')}
 if(abs(gamma_p) <= 1e-08){QP3 = qnorm(P, mu_p, sigma_p); return(QP3)}else{
  alpha = 4/gamma_p^2
  beta = abs(0.5 * sigma_p * gamma_p)
  QP3 = ifelse(gamma_p > 0,  mu_p - alpha * beta + beta * pmax(0, qgamma(P, alpha)),
   mu_p + alpha * beta - beta * pmax(0, qgamma(1 - P, alpha)))
  return(QP3)}
}
