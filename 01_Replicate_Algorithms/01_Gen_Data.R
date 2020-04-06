library(dplyr)
library(mvtnorm)
library(rlang)
library(copula)

#Function To Generate Multivariate Data From a Given Distribution

Gen.MV.Dat = function(size, correlation, 
                      x_transform, y_transform, 
                      x_param1, x_param2 = NULL,   # Accomodates 1 or 2-parameter functions
                      y_param1, y_param2 = NULL,
                      x_quantiles, y_quantiles){
  DATA = rmvnorm(n = size, mean = c(0,0),
                 sigma = matrix(c(1,correlation,correlation,1),
                                nrow = 2,byrow = TRUE)) %>%                      #S1: Gen Normal Data
    data.frame() %>%
    rename(x = X1, y = X2) %>%
    mutate_all(function(x) pnorm(x)) %>%                                         #S2: Transform to Uniform
    mutate(x_cat = jitter(x_transform(x, x_param1, x_param2),amount = 0.0001),
           y_cat = jitter(y_transform(y, y_param1, y_param2),amount = 0.0001))   #S3: Transform to Desired
  x_breaks = c(-Inf, as.numeric(quantile(DATA$x_cat, x_quantiles)),Inf)  
  y_breaks = c(-Inf, as.numeric(quantile(DATA$y_cat, y_quantiles)),Inf)          #S4: Extract Quantles
  DATA = DATA %>%
    mutate(x_cat = cut(x_cat, breaks = x_breaks, 
                       labels = as.character(1:(length(x_breaks)-1))),
           y_cat = cut(y_cat, breaks = y_breaks, 
                       labels = as.character(1:(length(y_breaks)-1))))           #S5:  Cut into Factor
  #S6:  Return DATA
  return(DATA)}   


Gen.Copula.Dat = function(size, parameter, CopulaType,
                          x_quantiles, y_quantiles, equal = TRUE){
  #Define Copula
  
  REAL_COP = eval(parse_expr(CopulaType))
  if(str_detect(CopulaType, "Nelsen")){
    #Implements Algorithm From Nelsen Exercise 4.16, or from paper directly
    DATA = REAL_COP$Gen_Dat(size, parameter) 
  }
  else{
    #Implement method available in "copula" package
    COPULA = REAL_COP(param  = parameter, dim = 2)
    DATA   = rCopula(size, COPULA) %>%
      data.frame()
    colnames(DATA) = c("x","y")
  }
  #Uniform Discretizing -- Evenly-Spaced Intervals
  if(equal == TRUE){
      x_breaks = c(0, as.numeric(quantile(DATA$x, x_quantiles)),1)  
      y_breaks = c(0, as.numeric(quantile(DATA$y, y_quantiles)),1) 
  }
  else{
      x_breaks = c(0, as.numeric(x_quantiles), 1)
      y_breaks  = c(0, as.numeric(y_quantiles), 1)
  }
  DATA = DATA %>%
    mutate(x_cat = cut(x, breaks = x_breaks, 
                       labels = as.character(1:(length(x_breaks)-1))),
           y_cat = cut(y, breaks = y_breaks, 
                       labels = as.character(1:(length(y_breaks)-1)))) 
  return(DATA)
}

Test = Gen.Copula.Dat(100, 0.5,"claytonCopula",x_quantiles = c(0.1,0.3,0.7,0.9), y_quantiles =  c(0.25,0.75), equal = TRUE)
head(Test)

aggregate(Test$x, by = list(Test$x_cat), FUN = count)

Test %>%
  group_by(x_cat) %>%
  tally()


COPULA = normalCopula(param = 0.5, dim = 2)
DAT = rCopula(100, COPULA) %>%
  data.frame() 
colnames(DAT) = c("x","y")

head(DAT)

