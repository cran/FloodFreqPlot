#' @importFrom stats qt qnorm

ConfInt = function(quant, n_samp, SE, ConfLev){

  if(n_samp <= 30){
   L_Talpha = quant - qt(p = 1 - (1 - ConfLev)/2, df = n_samp - 1) * SE
   U_Talpha = quant + qt(p = 1 - (1 - ConfLev)/2, df = n_samp - 1) * SE}else{

   L_Talpha = quant - qnorm(p = 1 - (1 - ConfLev)/2) * SE
   U_Talpha = quant + qnorm(p = 1 - (1 - ConfLev)/2) * SE
  }

 return(ConfLim = data.frame(Lower_Limit = L_Talpha, Upper_Limit = U_Talpha))
}