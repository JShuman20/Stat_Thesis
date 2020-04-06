#This Script Will Replicate Results for the Polychoric Correlation Using 
#The Normal Definition, the generalized normal definition, and the empirical 
#definition under a variety of simulation parameters

#Loading Libraries
library(copula)
library(tidyverse)
library(rlang)
library(rlist)
library(mvtnorm)
library(readxl)
library(purrr)
library(pracma)
library(cubature)


#----------------#NOTE: COMMANDS BELOW EXECUTED ON HPCC: ~ 10hrs with 8 processors each----------#
options("future.fork.enable" = TRUE) # Allow parallel processing
future::plan(multicore(workers = 5)) # Set Number of Cores
argList = list(
  SAMPLE_SIZES = c(50,100,500),
  COPULAS      = c("normalCopula", "frankCopula", "claytonCopula", "Nelsen_2", "Nelsen_15_GG")
)
crossedArgs <- cross_df(argList)
COR_VALS     = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)
#Reading Data Frame
SUPPORTS = read_xlsx("~/Desktop/THESIS/Thesis_Computing/DATA/02_Copula_Supports.xlsx")

#Bootstrap
for(i in c(3,5,7)){
   future_map2_dfr(.x = crossedArgs$SAMPLE_SIZES,.y = crossedArgs$COPULAS, 
              .f = function(.x,.y){
                        print(.y)
                        return(Bootstrap_Test_Copulas(nsamp = .x,copulaType = .y,nquantiles = i,
                                                    PPCs = COR_VALS, iterations = 1000))
                         }) %>%
    write.table(str_c("~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/01_Test_Normal/Bootstrap/SIM_BS_"), #Clean Paths
                as.character(i))
    
}

#Independent Tables
for(i in c(3,5,7)){
  future_map2_dfr(.x = crossedArgs$SAMPLE_SIZES,.y = crossedArgs$COPULAS, 
              .f = function(.x,.y){
                print(.y)
                return(Test_Indep_Tables(nsamp = .x,copulaType = .y,nquantiles = i,
                                         PPCs = COR_VALS, iterations = 1000))
              }) %>%
    write.table(str_c("~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/01_Test_Normal/Bootstrap/SIM_ID_"),
                as.character(i))
  
}

                                  
