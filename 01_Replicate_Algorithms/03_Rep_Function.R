#load libraries
library(polycor)
library(purrr)
library(tidyverse)
library(rlang)
library(mvtnorm)
library(copula)
library(scatterplot3d)
library(ggplot2)

#--------------------------------------------Replicating Normal Polychoric---------------------------------#
#-------METHOD 1: Two-Step Estimation (Olson, 1979)----------#
#Note: This method uses quantiles obtained through qnorm
#Maximization Function -- Normal Two-Step Polychoric Function
Normal_PCC_Two_Step = function(rho, DATA){
  #Obtaining Table and Qnorm Thresholds
  X_cuts = Gen_Thresholds(DATA, quote(x_cat))
  Y_cuts = Gen_Thresholds(DATA, quote(y_cat))
  Table = Cont.Table(DATA)
  #Generating Results Table
  Cell_Likelihoods = matrix(data = NA, nrow = length(X_cuts), 
                          ncol = length(Y_cuts))
  #Populating Results Table
  for(i in 1:length(X_cuts)){
    for(j in 1:length(Y_cuts)){
      lower.x = ifelse(i==1, -Inf, X_cuts[i-1])
      lower.y = ifelse(j==1, -Inf, Y_cuts[j-1])
      upper.x = X_cuts[i]
      upper.y = Y_cuts[j]
      Cell_Likelihoods[i,j] = Table[i,j] * 
         log(pmvnorm(lower =  c(lower.x,lower.y),
                     upper =  c(upper.x,upper.y),
                     mean  =  c(0,0),
                     corr  =  matrix(c(1,rho,rho,1),nrow=2,byrow = TRUE)))
    }
  }
  return(sum(Cell_Likelihoods))
}
#-------Method 2: Generalize with Lp Penalty-------#
Generalized_PCC_Norm = function(parameter, CopulaType, DATA, NORM, RMSE = FALSE){
  #First Define the Copula
  real_cop = eval(parse_expr(CopulaType))
  if(!str_detect(CopulaType,"Nelsen")){
  if(CopulaType == "plackettCopula"){
    COPULA = real_cop(parameter)
  }
  else{
  COPULA  = real_cop(param = parameter, dim = 2)
  }
  }
  #Contingency Table Option
  if(is.table(DATA) | is.matrix(DATA)){
    X_cuts = Gen.Cont.Table.Margins(DATA, TRUE)
    Y_cuts = Gen.Cont.Table.Margins(DATA, FALSE)
    Table = DATA
  }
  else{
    X_cuts = Gen_Emp_Thresholds(DATA, quote(x_cat))
    Y_cuts = Gen_Emp_Thresholds(DATA, quote(y_cat))
    Table  = Cont.Table(DATA) 
  }
  #Now Define the Cells
  Cell_Likelihoods = matrix(data = NA, nrow = length(X_cuts), 
                            ncol = length(Y_cuts))
  for(i in 1:length(X_cuts)){
    for(j in 1:length(Y_cuts)){
      lower.x = ifelse(i==1, -Inf, X_cuts[i-1])
      lower.y = ifelse(j==1, -Inf, Y_cuts[j-1])
      upper.x = X_cuts[i]
      upper.y = Y_cuts[j]
      if(str_detect(CopulaType, "Nelsen")){
        Cell_Likelihoods[i,j] = abs((real_cop$Main_Function(parameter, upper.x,upper.y) +
                                     real_cop$Main_Function(parameter, lower.x,lower.y) - 
                                     real_cop$Main_Function(parameter, upper.x,lower.y) - 
                                     real_cop$Main_Function(parameter, lower.x,upper.y)) -  (Table[i,j]/sum(Table)))^(NORM)
      }
      else{
        Cell_Likelihoods[i,j] = 
          abs((pCopula(c(upper.x,upper.y), COPULA) +    #Why was this abs()???
                 pCopula(c(lower.x, lower.y),COPULA) - 
                 pCopula(c(upper.x, lower.y),COPULA) - 
                 pCopula(c(lower.x, upper.y),COPULA)) - (Table[i,j]/sum(Table)))^(NORM)
      }
      
    }
  }
  if(RMSE == TRUE){
    Cell_Likelihoods = Cell_Likelihoods/(length(X_cuts) * length(Y_cuts))
    Final = (sum(Cell_Likelihoods))^(1/NORM)
  }
  else{
    Final = (sum(Cell_Likelihoods))^(1/NORM)
  }
  return(Final)
} 

#--------------------Method 2A: Optimizing With Non-Continuous Support------------------------#
#Note: Given How The non-continuous corresponds to negative and positive correlation, 
#the comparison of 
Optimize.Copula = function(CopulaType, DATA, NORM = 2, Neg_Support, Pos_Support = 0){
  X1 = optimize(Generalized_PCC_Norm, interval = Neg_Support, maximum = FALSE,
                CopulaType = CopulaType, DATA = DATA, NORM = 2)
  if(length(Pos_Support) == 2){
    X2 = optimize(Generalized_PCC_Norm, interval = Pos_Support, maximum = FALSE,
                  CopulaType = CopulaType, DATA = DATA, NORM = 2)
  }
  else{
    X2 = X1$objective + 1
  }
  if(!is.list(X2)){
    Result = Convert_To_Correlation(CopulaType, optimal_param = X1$minimum)
    }
  else{
    if(X1$objective < X2$objective){
      Result = Convert_To_Correlation(CopulaType, optimal_param = X1$minimum)
      }
    else{
      Result = Convert_To_Correlation(CopulaType, optimal_param = X2$minimum)}
  }
  return(Result) }

Test_Dat = Gen.Copula.Dat(100,3,"frankCopula", map_dbl(1:(4-1), ~ .x*1/4),map_dbl(1:(4-1), ~ .x*1/4))
cor(Test_Dat$x,Test_Dat$y)
Optimize.Copula("frankCopula",Test_Dat, 2,c(-1000,-0.0001),c(0.01,10))

#-------Method 4: Empirical Copula Estimation-------#
Empirical_PCC = function(DATA){
  #Generating Necessary Data Structures
  if(is.table(DATA) | is.matrix(DATA)){
    Empirical_Props = Gen_Empirical_Perc_V2(DATA)
    Cumulative_Props = Gen_Empirical_Cumulative(DATA)
  }
  else{
  Table = Cont.Table(DATA)
  Empirical_Props = Gen_Empirical_Perc_V2(Table)
  Cumulative_Props = Gen_Empirical_Cumulative(Table)
  }
  #Creating Results Matrix
  Results_Table = matrix(nrow = nrow(Table), ncol = ncol(Table))
  rownames(Results_Table) = c(1:nrow(Table))
  colnames(Results_Table) = c(1:ncol(Table))
  for (i in 1:nrow(Results_Table)){
    for(j in 1:ncol(Results_Table)){
      Results_Table[i,j] = Empirical_Props[i,j] * (1/4)*(sum(Cumulative_Props[i:(i+1),j:(j+1)]))
    } 
  }
  FINAL = 12 * sum(Results_Table) - 3
  return(FINAL)
}





#--------------------------What's Wrong With the Empirical-------------------------#
#A: We know that 2sin(E * pi/6) must be between -1,1
Spearmans = seq(-1.5,1.5, by = 0.01)
EPC = 2 * sin(Spearmans*pi/6)
ggplot() + 
  geom_line(aes(x=Spearmans, y = EPC)) + 
  geom_hline(yintercept = c(-1,1),col= "red") +
  geom_vline(xintercept = c(-1,1), col= "red")
#B: Therefore, E must be between -1 and 1
#C: This means that 12*sum(stuff) -3 is between -1,1
#D: 12 * sum(stuff) must be between -2,4
#E: sum(stuff) must be between -1/6 and 1/3
