#' Plotting Position Probability
#'
#' \code{PlotPos} returns the empirical probability values corresponding to 
#'    the observed data of hydrological extreme events as a vector of numerics.
#'
#' This is a function to calculate the emprical probability values assigned to the observed data
#'    of hydrological extreme events to be plotted.
#'
#' @param PP A character string that determines the empirical formula used to calculate the probability.
#'     The formula can be chosen from the list: \code{"Blom"}, \code{"Chegodayev"}, \code{"California"}
#'     \code{"Gringorten"}, \code{"Hazen"}, \code{"Tukey"}, and \code{"Weibull"}.
#' @param data_obs A vector, data frame or matrix containing observed data or flood quantiles.
#'
#' @return The function returns the probabilities assigned to the observed data as a vector of numerics.
#'
#' @examples
#' # First Example
#' data('Harricana')
#' PlotPos(data_obs = Harricana, PP = 'Weibull')
#'
#' # Second Example
#' data('B17C_Tab8_1')
#' PlotPos(data_obs = B17C_Tab8_1, PP = 'Cunnane')
#'
#' @section Reference:
#' Chow, V. T., Maidment, D. R., & Mays, L. W. (1988). \emph{Applied Hydrology}. McGraw-Hill, New York, U.S.
#'
#' @seealso \code{\link{ProbPlot}} for graphical frequency analysis.
#' @export
PlotPos = function(data_obs, PP){

 if(!is.vector(data_obs)){
  n = dim(data_obs)[1]
  x = data_obs[1:n, 1]}else{
  n = length(data_obs)
  x = data_obs[1:n]}

 data_rnk = rank(x, ties.method = 'max')

 i_dr = sort(data_rnk)

 probs = rep(0, n)

 if(PP == 'Blom'){
  for(i in 1:n){
   probs[i] = (i_dr[i] - (3/8))/(n + 1 - 2 * (3/8))}}

#  if(PP == 'California_1'){
#   for(i in 1:n){
#    probs[i] = i_dr[i]/n}}

  if(PP == 'California'){
   for(i in 1:n){
    probs[i] = (i_dr[i] - 1)/n}}

  if(PP == 'Chegodayev'){
   for(i in 1:n){
    probs[i] = (i_dr[i] - 0.3)/(n + 0.4)}}

  if(PP == 'Cunnane'){
   for(i in 1:n){
    probs[i] = (i_dr[i] - 0.4)/(n + 0.2)}}

  if(PP == 'Gringorten'){
   for(i in 1:n){
    probs[i] = (i_dr[i] - 0.44)/(n + 1 - 2 * 0.44)}}

  if(PP == 'Hazen'){
   for(i in 1:n){
    probs[i] = (i_dr[i] - 0.5)/n}}

  if(PP == 'Tukey'){
   for(i in 1:n){
    probs[i] = (i_dr[i] - (1/3))/(n + 1 - 2 * (1/3))}}

  if(PP == 'Weibull'){
   for(i in 1:n){
    probs[i] = i_dr[i]/(n + 1)}}

 return(probs)
}