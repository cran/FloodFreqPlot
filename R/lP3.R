pdflP3 = function(yp3, alpha, beta, tau){
 expon = ifelse((((yp3 - tau)/beta) < 0) & (((alpha - 1) %% 1) != 0), round(alpha - 1), alpha - 1)
 density = (1/(abs(beta) * gamma(alpha))) * ((yp3 - tau)/beta)^expon * exp(-(yp3 - tau)/beta)
 return(density)
}

#############################################################
#' @importFrom stats pnorm pgamma

cdflP3 = function(yp3, mu_p = y_bar, sigma_p = s_y, gamma_p = C_sy){
 if(yp3 <= 0){PlP3 = 0; return(PlP3)}
 if(sigma_p <= 0){stop('invalid parameters')}
 if(gamma_p <= 1e-06){PlP3 = pnorm(yp3, mu_p, sigma_p); return(PlP3)}else{
  alpha = 4/gamma_p^2
  zlP3 = 2 * (yp3 - mu_p)/(sigma_p * gamma_p) + alpha
  PlP3 = pgamma(pmax(0, zlP3), alpha)
  if(gamma_p < 0){PlP3 = 1 - PlP3}
  return(PlP3)
 }
}

#############################################################

lP3P = function(yp3, P){cdflP3(yp3) - P}

#############################################################

#' @importFrom stats qnorm qgamma

qualP3 = function(P, mu_p = y_bar, sigma_p = s_y, gamma_p = C_sy){
 if(P < 0){stop('invalid probability')}
 if(sigma_p <= 0){stop('invalid parameters')}
 if(abs(gamma_p) <= 1e-08){QlP3 = qnorm(P, mu_p, sigma_p); return(QlP3)}else{
  alpha = 4/gamma_p^2
  beta = abs(0.5 * sigma_p * gamma_p)
  QlP3 = ifelse(gamma_p > 0,  mu_p - alpha * beta + beta * pmax(0, qgamma(P, alpha)),
  mu_p + alpha * beta - beta * pmax(0, qgamma(1 - P, alpha)))
  return(QlP3)}
}