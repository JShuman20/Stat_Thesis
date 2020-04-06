library(purrr)
library(rlist)
library(rlang)
library(copula)
library(mvtnorm)
library(dplyr)
library(readxl)
library(stringr)
library(furrr)

#-----------------------Test 1: Normal(rho = 0.4)--------------------------------#
DAT_Norm = Gen.MV.Dat(size = 100, correlation = 0.4,
                      x_transform = qnorm, y_transform = qnorm,
                      x_param1 = 0, x_param2 = 1,
                      y_param1 = 0, y_param2 = 1,
                      x_quantiles =  c(0.25,0.7,0.9),
                      y_quantiles =  c(0.3, 0.55,0.8))
#First Test With Polycor Package
polychor(Cont.Table(DAT_Norm))
#With Normal PCC Function
optimize(Normal_PCC_Two_Step, lower = -1, upper = 1, maximum = TRUE,
         DATA = DAT_Norm)
#With Generalized PCC Function
X = optimize(Generalized_PCC_Norm, lower = -0.999, upper = 0.999,
             maximum = FALSE,
             DATA = DAT_Norm,
             CopulaType = normalCopula, NORM = 2)
#With Empirical PCC
2 * sin(Empirical_PCC(DAT_Norm)*pi/6)

#-----------------Replicating Results From Ekstrom Generalized Table 2---------------#
#List of Tables From Table 2 of Generalized Paper
Tables_List = list(
  matrix(c(0.25,0.25,0.25,0.25), nrow = 2, byrow = TRUE),
  matrix(c(0.5,0,0,0.5), nrow = 2, byrow = TRUE),
  matrix(c(0.64, 0.16, 0.16,0.04), nrow = 2, byrow = TRUE),
  matrix(c(0.16,0.64,0.04,0.16), nrow = 2, byrow = TRUE),
  matrix(c(0.6,0.2,0.2,0), nrow = 2, byrow = TRUE),
  matrix(c(0.2,0.6,0,0.2), nrow = 2, byrow = TRUE),
  matrix(c(0.08,0.12,0.12,0.68), nrow = 2, byrow = TRUE),
  matrix(c(0.68,0.12,0.12,0.08), nrow = 2, byrow = TRUE)
)

Rep_Results = data.frame(
  Table   = c(1:8),
  Normal  = map_dbl(.x = Tables_List,
                    .f = ~ round(Optimize.Copula(CopulaType = normalCopula, 
                                                 DATA = .x, NORM = 2, 
                                                 Neg_Support = c(-1,1)),2)),
  Frank   = map_dbl(.x = Tables_List, 
                    .f = ~ round(Optimize.Copula(CopulaType = frankCopula, 
                                                 DATA = .x, NORM = 2, 
                                                 Neg_Support = c(-100,-0.001), 
                                                 Pos_Support = c(0.001,100)),2)),
  Clayton = map_dbl(.x = Tables_List, 
                    .f = ~ round(Optimize.Copula(CopulaType = claytonCopula, 
                                                 DATA = .x, NORM = 2, 
                                                 Neg_Support = c(-1,-0.0001), 
                                                 Pos_Support = c(0.0001,1000)),2))
)
Rep_Results
#Normal and Frank match up to a small deviation
#Clayton mostly matches, with exception of 5 and 6

#--------------------------------------Replicating Table 3--------------------------------#
#Replicating Contingency Table
Table3 = matrix(c(0,7,0,0,3,
                  3,10,25,10,3,
                  18,84,80,47,7,
                  40,54,65,43,10,
                  43,29,29,14,10),nrow = 5,byrow = TRUE)
#Replicating Column 3 of Results
Tab3.Gaussian = round(Optimize.Copula(CopulaType = normalCopula, 
                                      DATA = Table3, NORM = 2, 
                                      Neg_Support = c(-1,1)),2)
Tab3.Frank = round(Optimize.Copula(CopulaType = frankCopula, 
                                   DATA = Table3, NORM = 2, 
                                   Neg_Support = c(-100,-0.001), 
                                   Pos_Support = c(0.001,100)),2)
Tab3.Clayton = round(Optimize.Copula(CopulaType = claytonCopula, 
                                     DATA = Table3, NORM = 2, 
                                     Neg_Support = c(-1,-0.0001), 
                                     Pos_Support = c(0.0001,1000)),2)
#All Three Match

#---------------------------------------------------Replicating Ekstrom Empirical: Table 1-------------------------------#
#----------------------------------#--Version 1: Bootstrapped Samples For Same DataSet-------#---------------------------#
#Define an Optimize Function To Pick Copula Parameter to Match Population Polychoric
Find_Param = function(param, Copula, desired_PCC){
  RES = abs(Convert_To_Correlation(Copula,param) - desired_PCC)
  return(RES)
}
#Defining a Dataset of Positive and Negative Supports
SUPPORTS = read_xlsx("~/Desktop/THESIS/Thesis_Computing/DATA/02_Copula_Supports.xlsx")
#Carries Out any of Ekstrom's Simulations With Desired Parameters
Bootstrap_Test_Copulas = function(nsamp,copulaType, nquantiles, PPCs, iterations = 100){
  if(copulaType == "normalCopula"){
    Copula_List = map(.x = PPCs,
                           .f = ~ Gen.Copula.Dat(size        = nsamp, 
                                                 parameter   = .x,
                                                 CopulaType  = copulaType,
                                                 x_quantiles = map_dbl(1:(nquantiles-1), ~ .x*1/nquantiles),
                                                 y_quantiles = map_dbl(1:(nquantiles-1), ~ .x*1/nquantiles)))
  }
  else{
    Copula_List = map(.x = PPCs,
                      .f = function(.x){
                        FILT_SUP = SUPPORTS[which(SUPPORTS$COPULA==copulaType),]
                        NEG_SUP = eval(parse_expr(FILT_SUP$NEG))
                        POS_SUP = eval(parse_expr(FILT_SUP$POS))
                        POS_OPT = optimize(Find_Param, interval = POS_SUP,
                                         copulaType,.x)
                        if(length(NEG_SUP)==1){
                          Final_Param = POS_OPT$minimum
                        }
                        else{
                          NEG_OPT = optimize(Find_Param, interval = NEG_SUP,
                                             copulaType,.x)
                          if(.x<0){
                            Final_Param = NEG_OPT$minimum
                          }
                          else{
                            Final_Param = POS_OPT$minimum
                          }
                        }
                        Gen.Copula.Dat(size        = nsamp, 
                                       parameter   = Final_Param,
                                       CopulaType  = copulaType,
                                       x_quantiles = map_dbl(1:(nquantiles-1), ~ .x*1/nquantiles),
                                       y_quantiles = map_dbl(1:(nquantiles-1), ~ .x*1/nquantiles))})
  }
  EMP_RES  = matrix(nrow = iterations, ncol = length(PPCs))
  NORM_RES = matrix(nrow = iterations, ncol = length(PPCs))
  GEN_RES = matrix(nrow = iterations, ncol = length(PPCs))
  for(i in 1:length(Copula_List)){
    NORM_RES[,i] = map_dbl(.x = 1:iterations, 
                           .f = function(.x){
                             set.seed(.x)
                             Indeces = sample(1:nrow(Copula_List[i][[1]]),
                                              nrow(Copula_List[i][[1]]), replace = TRUE)
                             DATA_BS = (Copula_List[i][[1]])[Indeces,]
                             return(optimize(Normal_PCC_Two_Step, lower = -1, upper = 1, 
                                             maximum = TRUE,
                                             DATA = DATA_BS)$maximum)})
    GEN_RES[,i] = map_dbl(.x = 1:iterations,
                          .f = function(.x){
                            set.seed(.x)
                            Indeces = sample(1:nrow(Copula_List[i][[1]]),
                                             nrow(Copula_List[i][[1]]), replace = TRUE)
                            DATA_BS = (Copula_List[i][[1]])[Indeces,]
                            return(optimize(Generalized_PCC_Norm, interval = c(-1,1),
                                            "normalCopula", DATA_BS, 2)$minimum)
                          })
    EMP_RES[,i] = map_dbl(.x = 1:iterations, 
                          .f = function(.x){
                            set.seed(.x)
                            Indeces = sample(1:nrow(Copula_List[i][[1]]),
                                             nrow(Copula_List[i][[1]]), replace = TRUE)
                            DATA_BS = (Copula_List[i][[1]])[Indeces,]
                            return(2 * sin(Empirical_PCC(DATA_BS)*pi/6))
                          })
    
  }
  Final_RES = cbind(NORM_RES, GEN_RES, EMP_RES) %>%
    data.frame()
  colnames(Final_RES) = c(map_chr(.x = PPCs,.f = ~ str_c("Normal_",.x,sep = "")),
                          map_chr(.x = PPCs,.f = ~ str_c("GEN_",.x,sep = "")),
                          map_chr(.x = PPCs,.f = ~ str_c("EMP_",.x,sep = "")))
  Final_RES = Final_RES %>%
    mutate(SIZE = nsamp,
           REAL_TYPE = copulaType,
           QUANTILES = nquantiles)
  return(Final_RES)
}


map_dfr(.x = c("Nelsen_2", "Nelsen_15_GG"), 
               .f = ~ Bootstrap_Test_Copulas(nsamp=100,copulaType = .x, nquantiles = 3, PPCs = c(-0.67,0.67),iterations = 5)) %>%
  write.table("~/Desktop/BS_Test")

SIM_BS_100_3_Nelsen = map_dfr(.x = c("Nelsen_2", "Nelsen_15_GG"),
                              .f = ~ Bootstrap_Test_Copulas(nsamp=100,copulaType = .x,nquantiles = 3,
                                                            iterations = 1000))

#This Has been Run -- took list 20 mins
SIM_BS_100_3 = map_dfr(.x = c("claytonCopula", "frankCopula","normalCopula"), 
               .f = ~ Bootstrap_Test_Copulas(nsamp = 100,copulaType = .x,nquantiles = 3,iterations = 1000))

#--------------------------------------Version 2 of Replication: Generate New Table Each Time With Random Seed-------------------------#
Test_Indep_Tables = function(nsamp,copulaType, nquantiles,PPCs, iterations = 100){
  #First Generating a bunch of tables
  TOT_RES = list()
  #Create Revised Parameter Values
  if(copulaType=="normalCopula"){
    Revised_PPCs = PPCs
  }
  else{
    Revised_PPCs = map_dbl(.x = PPCs,
                           .f = function(.x){
                             FILT_SUP = SUPPORTS[which(SUPPORTS$COPULA==copulaType),]
                             NEG_SUP = eval(parse_expr(FILT_SUP$NEG))
                             POS_SUP = eval(parse_expr(FILT_SUP$POS))
                             POS_OPT = optimize(Find_Param, interval = POS_SUP,
                                                copulaType,.x)
                             if(length(NEG_SUP)==1){
                               Final_Param = POS_OPT$minimum
                             }
                             else{
                               NEG_OPT = optimize(Find_Param, interval = NEG_SUP,
                                                  copulaType,.x)
                               if(.x<0){
                                 Final_Param = NEG_OPT$minimum
                               }
                               else{
                                 Final_Param = POS_OPT$minimum
                               }
                             }
                           })
  }
  #Initialize Result Vector
  NORM_RES = matrix(nrow = iterations, ncol = length(PPCs))
  GEN_RES  = matrix(nrow = iterations, ncol = length(PPCs))
  EMP_RES  = matrix(nrow = iterations, ncol = length(PPCs))
  for(i in 1:length(PPCs)){
    TOT_RES[[i]] = map(.x = 1:iterations,
                       .f = function(.x){
                       set.seed(.x)
                       Gen.Copula.Dat(size        = nsamp, 
                                      parameter   = Revised_PPCs[i],
                                      CopulaType  = copulaType,
                                      x_quantiles = map_dbl(1:(nquantiles-1), ~ .x*1/nquantiles),
                                      y_quantiles = map_dbl(1:(nquantiles-1), ~ .x*1/nquantiles))
                       })
    library(dplyr)
    NORM_RES[,i] = map_dbl(.x = TOT_RES[[i]],
                           .f = ~ optimize(Normal_PCC_Two_Step, interval = c(-1,1),
                                           DATA = .x, maximum = TRUE)$maximum)
    GEN_RES[,i] = map_dbl(.x = TOT_RES[[i]],
                           .f = ~ optimize(Generalized_PCC_Norm, interval = c(-1,1),
                                           "normalCopula", .x, 2)$minimum)
    EMP_RES[,i]  = map_dbl(.x = TOT_RES[[i]],
                          .f = ~ return(2 * sin(Empirical_PCC(.x)*pi/6)))
  }
  #Formatting Results Matrix
  FINAL_RES = cbind(NORM_RES,GEN_RES,EMP_RES) %>%
    data.frame()
  colnames(FINAL_RES) =  c(map_chr(.x = PPCs,.f = ~ str_c("Normal_",.x,sep = "")),
                           map_chr(.x = PPCs,.f = ~ str_c("GEN_",.x,sep = "")),
                           map_chr(.x = PPCs,.f = ~ str_c("EMP_",.x,sep = "")))
  
  FINAL_RES = FINAL_RES %>%
    mutate(SIZE      = as.character(nsamp),
           REAL_TYPE = copulaType,
           QUANTILES = as.character(nquantiles))
  return(FINAL_RES)
}

#---------------------------Summary Function---------------------------#
MySum = function(DAT){
  MEAN = DAT %>%
    group_by(SIZE, REAL_TYPE, QUANTILES) %>%
    summarise_if(is.numeric, ~mean(.)) %>% 
    data.frame() %>% 
    mutate(STAT = "MEAN")
  SD = DAT %>%
    group_by(SIZE, REAL_TYPE, QUANTILES) %>%
    summarise_if(is.numeric, ~sd(.)) %>% 
    data.frame() %>% 
    mutate(STAT = "SD")
  RES = rbind(MEAN,SD)
  return(RES)
}








#---------------------------------------------------------------------------------------------------------------#
#----------------------------------------------ARCHIVE----------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#

#------------3,5,7 Categories, N=100, 1000 Iterations---------------------#
for(i in c(3,5,7)){
  PATH = str_c("SIM_ID_", as.character(i),"_100")
  map(.x = c("normalCopula","frankCopula","claytonCopula"),
      .f = ~ Test_Indep_Tables(100,.x,i,c(-0.67,-0.33,0.33,0.67), iterations = 1000)) %>%
  map_dfr(.f = ~ MySum(.x)) %>%
    write.table(file = str_c(getwd(), "/Outputs/Simulation/01_Test_Normal/", PATH))
}

#-----------3,5,7 Categories, N=500, 1000 iterations---------------------#
for(i in c(3,5,7)){
  PATH = str_c("SIM_ID_", as.character(i),"_500")
  map(.x = c("normalCopula","frankCopula","claytonCopula"),
      .f = ~ Test_Indep_Tables(500,.x,i,c(-0.67,-0.33,0.33,0.67), iterations = 1000)) %>%
    map_dfr(.f = ~ MySum(.x)) %>%
    write.table(file = str_c(getwd(), "/Outputs/Simulation/01_Test_Normal/", PATH))
}
#Gumbel
for(i in c(3,5,7)){
  PATH = str_c("SIM_ID_Gumbel", as.character(i),"_500")
  map(.x = c("gumbelCopula"),
      .f = ~ Test_Indep_Tables(100,.x,i,c(0.25,0.5,0.75), iterations = 1000)) %>%
    map_dfr(.f = ~ MySum(.x)) %>%
    write.table(file = str_c(getwd(), "/Outputs/Simulation/01_Test_Normal/", PATH))
}
#Joe
for(i in c(3,5,7)){
  PATH = str_c("SIM_ID_Joe_", as.character(i),"_500")
  map(.x = c("joeCopula"),
      .f = ~ Test_Indep_Tables(100,.x,i,c(-0.67,-0.33,0.33,0.67), iterations = 1000)) %>%
    map_dfr(.f = ~ MySum(.x)) %>%
    write.table(file = str_c(getwd(), "/Outputs/Simulation/01_Test_Normal/", PATH))
}
Gumbels = map_dfr(.x = c("3_500","5_500","7_500"),
                  .f = ~ read.table(str_c(getwd(), "/Outputs/Simulation/01_Test_Normal/SIM_ID_Gumbel",.x)))

#Test with Nelsen- 2
X = Test_Indep_Tables(300,"Nelsen_15",3,c(-0.67,-0.33,0.33,0.67), iterations = 100) 




#Bootstrap CIs
set.seed(17)
for(i in c(3,5,7)){
  PATH = str_c("SIM_BS_", as.character(i),"_100")
  map(.x = c("normalCopula","frankCopula","claytonCopula"),
      .f = ~ Bootstrap_Test_Copulas(100,.x,i,iterations = 1000)) %>%
    map_dfr(.f = ~ MySum(.x)) %>%
    write.table(file = str_c(getwd(), "/Outputs/Simulation/01_Test_Normal/", PATH))
}
for(i in c(3,5,7)){
  PATH = str_c("SIM_BS_", as.character(i),"_500")
  map(.x = c("normalCopula","frankCopula","claytonCopula"),
      .f = ~ Bootstrap_Test_Copulas(500,.x,i,iterations = 1000)) %>%
    map_dfr(.f = ~ MySum(.x)) %>%
    write.table(file = str_c(getwd(), "/Outputs/Simulation/01_Test_Normal/", PATH))
}


