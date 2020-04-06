library(mvtnorm)
library(dplyr)
library(polycor)
library(purrr)
library(ggplot2)
library(stargazer)

set.seed(18)

#--------------------Test 2: Bootstrapped Samples------------#
#Simulating Data
###Simulated data with cor = 0.5
###Sizes Are 20,40,60,100,250
List_of_Normals = map(.x = c(20,40,60,100,250),
                      .f =  ~ Gen.MV.Dat(size = .x, correlation = 0.5,
                                         x_transform = qnorm, y_transform = qnorm,
                                         x_param1 = 0, x_param2 = 1,
                                         y_param1 = 0, y_param2 = 1,
                                         x_quantiles =  c(0.2,0.5,0.8),
                                         y_quantiles =  c(0.1,0.6,0.9)))
#Bootstrapping Estimates for Normal PCC
Normal_Res      = matrix(nrow = 100,ncol = length(List_of_Normals))
Generalized_Res = matrix(nrow = 100,ncol = length(List_of_Normals))
Empirical_Res   = matrix(nrow = 100,ncol = length(List_of_Normals))
for(i in 1:length(List_of_Normals)){
  Normal_Res[,i] = map_dbl(.x = 1:100,
                           .f = function(.x){
                             set.seed(.x)
                             Indeces = sample(1:nrow(List_of_Normals[[i]]), 
                                              nrow(List_of_Normals[[i]]),replace = TRUE)
                             DATA_BS = (List_of_Normals[[i]])[Indeces,]
                             return(optimize(Normal_PCC_Two_Step, lower = -1, upper = 1, 
                                             maximum = TRUE,
                                             DATA = DATA_BS)$maximum)
                           })
  Generalized_Res[,i] = map_dbl(.x = 1:100,
                            .f = function(.x){
                              set.seed(.x)
                              Indeces = sample(1:nrow(List_of_Normals[[i]]), 
                                               nrow(List_of_Normals[[i]]),replace = TRUE)
                              DATA_BS = (List_of_Normals[[i]])[Indeces,]
                              return(Optimize.Copula(CopulaType = normalCopula,
                                                     DATA = DATA_BS,
                                                     NORM = 2,
                                                     Neg_Support = c(-1,1), Pos_Support = 0))
                            })
  Empirical_Res[,i]  = map_dbl(.x = 1:100,
                               .f = function(.x){
                               set.seed(.x)
                               Indeces = sample(1:nrow(List_of_Normals[[i]]), 
                                                nrow(List_of_Normals[[i]]),replace = TRUE)
                               DATA_BS = (List_of_Normals[[i]])[Indeces,]
                               return(2 * sin(Empirical_PCC(DATA_BS) * pi/6))})
}

Results = data.frame(Normal_Res,Generalized_Res,Empirical_Res) 
colnames(Results) = c("Normal_20", "Normal_40","Normal_60","Normal_100","Normal_250",
                      "Generalized_20", "Generalized_40","Generalized_60","Generalized_100","Generalized_250",
                      "Empirical_20", "Empirical_40","Empirical_60","Empirical_100","Empirical_250")

head(Results)
                      
                     
List_of_Copulas = map(.x = c(20,40,60,100,250),
                      .f = ~ Gen.Copula.Dat(size = .x,
                                            parameter = )
                      
                      size, parameter, CopulaType,
                      x_quantiles, y_quantiles
  
  
  .x = c(20,40,60,100,250),
                      .f =  ~ Gen.MV.Dat(size = .x, correlation = 0.5,
                                         x_transform = qnorm, y_transform = qnorm,
                                         x_param1 = 0, x_param2 = 1,
                                         y_param1 = 0, y_param2 = 1,
                                         x_quantiles =  c(0.2,0.5,0.8),
                                         y_quantiles =  c(0.1,0.6,0.9)))







