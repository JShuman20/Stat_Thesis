library(polycor)
library(purrr)
library(dplyr)
library(rlang)
library(mvtnorm)
library(copula)
library(scatterplot3d)
library(ggplot2)
library(tidyverse)



#Pre-Cursor: Building a List of FilePaths

PATH_BUILDER = list(
  DATA = "~/Desktop/THESIS/Thesis_Computing/DATA/",
  OUTPUTS = "~/Desktop/THESIS/Thesis_Computing/Outputs/"
)


#------01. Generate Uniform Thresholds From Inverse Normal-------#
Gen_Thresholds = function(DATA,grouping){
  DATA = DATA %>%
    group_by(!!grouping) %>%
    summarise(Count = n()/nrow(DATA)) %>%
    mutate(Cumulative = cumsum(Count),
           Thresholds = qnorm(Cumulative)) %>%
    pull(Thresholds) #Extracts as a vector
  return(DATA)
}
#-------02. Generate Empirical Thresholds-------#
#Note: Returns a Vector of Marginal Cumulative Thresholds Based on Empirical Proportions
Gen_Emp_Thresholds = function(DATA, grouping){
  DATA = DATA %>%
    group_by(!!grouping) %>%
    summarise(Count = n()/nrow(DATA)) %>%
    mutate(Cumulative = cumsum(Count)) %>%
    pull(Cumulative) #Extracts as a vector
  return(DATA)
}
#-------03. Generate Contingency Table-------#
Cont.Table = function(DATA){
  DATA = DATA %>% 
    dplyr::select(ends_with("cat")) %>% 
    table()
  return(DATA)
}
#-------0.4. Margins on Contingency Table----------------
Gen.Cont.Table.Margins = function(Table, Row){
  RES = vector()
  if(Row == TRUE){
    RES[1] = sum(Table[1,])/ sum(Table)
    for(i in 2:nrow(Table)){
      RES[i] = sum(Table[i,])/ sum(Table) + RES[i-1]
    }
  }
  else{
    RES[1] = sum(Table[,1])/sum(Table)
    for(i in 2:ncol(Table)){
      RES[i] = sum(Table[,i])/sum(Table) + RES[i-1] 
    }
  }
  return(RES)
}
#-------04. Generate Empirical Cumulative Thresholds-------#
#Notes: This Function Returns With First Row and Column Filled With Zeros
#This is for use on empirical PCC, where need to average at corners -- might be way around this with if statement
Gen_Empirical_Cumulative = function(TABLE){
  #Define Results Object
  RES = matrix(nrow = nrow(TABLE)+1,ncol = ncol(TABLE)+1)
  colnames(RES) = c(0:ncol(TABLE))
  rownames(RES) = c(0:nrow(TABLE))
  RES[1,] = 0; RES[,1] = 0
  for(i in 1:nrow(TABLE)){
    for(j in 1:ncol(TABLE)){
      RES[i+1,j+1] =  sum(TABLE[1:i,1:j])/sum(TABLE)
    }
  }
  return(RES)
}
#-------05. Generate Empirical Percentages of Contingency Table-------#
#Note: These are the non-cumulative values in each cell of the contingency table
Gen_Empirical_Perc = function(TABLE){
  RES = matrix(nrow = nrow(TABLE),ncol = ncol(TABLE))
  colnames(RES) = c(1:nrow(TABLE))
  rownames(RES) = c(1:nrow(TABLE))
  for(i in 1:nrow(TABLE)){
    for(j in 1:ncol(TABLE)){
      RES[i,j] =   TABLE[i,j]/sum(TABLE)
    }
  }
  return(RES)
}

#----------05a. Generate Empirical Percentages Under Null Hypothesis-------------#
Gen_Empirical_Perc_V2 = function(TABLE){
  RES = matrix(nrow = nrow(TABLE),ncol = ncol(TABLE))
  colnames(RES) = c(1:nrow(TABLE))
  rownames(RES) = c(1:nrow(TABLE))
  for(i in 1:nrow(TABLE)){
    for(j in 1:ncol(TABLE)){
      RES[i,j] =   (sum(TABLE[i,])/sum(TABLE)) * (sum(TABLE[,j])/sum(TABLE))
    }
  }
  return(RES)
}

#-------06. Extract Correlation From Copula Estimate-------#
Convert_To_Correlation = function(CopulaType, optimal_param){
  if(str_detect(CopulaType, "Nelsen")){
    Spearman = 12 * integral2(Vectorize(eval(parse_expr(CopulaType))$Main_Function), 
                              xmin = 0, xmax = 1,
                              ymin = 0, ymax = 1, theta = optimal_param)$Q - 3
  }
  else{
    real_cop = eval(parse_expr(CopulaType))
    Spearman = rho(real_cop(param = optimal_param))
  }
  Result = 2 * sin(Spearman * pi/6)
  return(Result)
}

