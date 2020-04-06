#The Original Function for one ordinal one continuous
#I added a section at the end to only consider points in the training set
Polyserial_LP_V2 = function(theta, DATA, COPULA, NORM = 2, CHOOSE_INDEX = FALSE){
  COP = eval(parse_expr(COPULA))
  Thresholds = Gen_Emp_Thresholds(DATA, quote(y_cat))
  DATA = DATA %>%
    mutate(y_cat = as.numeric(Thresholds[y_cat]))
  if(COPULA %in% c("frankCopula", "normalCopula", "claytonCopula", "amhCopula", "joeCopula")){
    #Define Copula Object
    COP = COP(theta, dim = 2)
    DIFF = map_dbl(.x = 1:nrow(DATA),
                   .f = function(.x){
                     COP_VAL = pCopula(c(DATA$x[.x], DATA$y_cat[.x]), COP)
                     REAL_PERC = mean(DATA$x <= DATA$x[.x] & DATA$y_cat <= DATA$y_cat[.x])
                     return((COP_VAL - REAL_PERC)^NORM)
                   })
  }
  else{
    DIFF = map_dbl(.x = 1:nrow(DATA),
                   .f = function(.x){
                     COP_VAL = COP$Main_Function(theta, DATA$x[.x], DATA$y_cat[.x])
                     REAL_PERC = mean(DATA$x <= DATA$x[.x] & DATA$y_cat <= DATA$y_cat[.x])
                     return((COP_VAL - REAL_PERC)^NORM)
                   })
  }
  if(length(CHOOSE_INDEX) > 1){ # False Has Length 1
    INDICATOR = DATA %>%
      mutate(IN_INDEX = ifelse(row_number() %in% CHOOSE_INDEX, 1,0)) %>%
      pull(IN_INDEX)
    DIFF = DIFF * INDICATOR
    RES = (sum(DIFF))^(1/NORM)
    return(RES)
  }
  else{
    RES = (sum(DIFF))^(1/NORM)
    return(RES)
  }
}

#Generate 50 Samples from Clayton

#Estimate With Given Cuts
PCC_Est_K = function(DAT,CUTS,COP,SUP, RETURN = "min"){
  k = seq(0,1,length.out = CUTS)
  NEW_DAT =  mutate(DAT, x_cat = cut(x,breaks = k,
                                     labels = as.character(1:(length(k)-1))))
  RES = optimize(Generalized_PCC_Norm, interval = SUP, maximum = FALSE,
                 CopulaType = COP, DATA = NEW_DAT, NORM = 2)
  if(RETURN == "all"){
    return(RES)
  }
  else if (RETURN == "obj"){
    return(RES$objective)
  }
  else{
    return(RES$minimum)
  }}

#Cross Validation Function
Cross_Validate_PCC = function(DATA,Cuts,COP,SUP,N_Cross, ADD_OBJ = FALSE){
  Size = nrow(DATA)/N_Cross
  RES = matrix(nrow = length(Cuts),ncol = N_Cross + 2)
  RES[,1] = Cuts
  #List of 10s
  ALL_TRAIN = map(.x = 0:(N_Cross-1), .f = ~ DATA[-((Size*.x) + 1):-(Size*(.x+1)),])
  for(i in 1:length(ALL_TRAIN)){
    Param_Ests = sapply(Cuts, PCC_Est_K, DAT = ALL_TRAIN[[i]], COP = COP, SUP = SUP, RETURN = "min")
    RES[,(i+1)] = sapply(Param_Ests, Polyserial_LP_V2, DATA, COP, 2, 
                         as.numeric(setdiff(as.character(1:nrow(DATA)),rownames(ALL_TRAIN[[i]])))) 
  }
  RES[,(N_Cross+2)] = map_dbl(.x = 1:length(Cuts), .f = ~ mean(RES[.x,2:(N_Cross+1)]))
  RES = RES[,c(1,(N_Cross+2))]
  RES = data.frame(RES) 
  colnames(RES) = c("Cuts","AVG")
  BEST_BINS = RES$Cuts[which(RES$AVG == min(RES$AVG))]
  if(ADD_OBJ == TRUE){
    RES = PCC_Est_K(DATA, BEST_BINS,COP=COP,SUP=SUP,RETURN = "all") 
    return(data.frame(
      CV_VAL  = RES$minimum,
      CV_OBJ  = RES$objective
    ))
  }
  else{
  return(data.frame(
    CV_CUTS = BEST_BINS,
    CV_VAL =  PCC_Est_K(DATA, BEST_BINS,COP=COP,SUP=SUP,RETURN = "min"))
  )
}}


MIN_L2 = function(DAT, Cuts,COP,SUP, ADD_OBJ = FALSE){
  OBJECTIVES = sapply(Cuts, PCC_Est_K, DAT = DAT,COP = COP,SUP=SUP, RETURN = "obj")
  MIN_L2_VAL = PCC_Est_K(DAT, CUTS = Cuts[which.min(OBJECTIVES)],COP,SUP, RETURN = "min")
  if(ADD_OBJ == TRUE){
    data.frame(
      MIN_L2_VAL = MIN_L2_VAL,
      MIN_L2_OBJ = min(OBJECTIVES)
    )
  }
  else{
  data.frame(
    MIN_L2_CUTS = Cuts[which.min(OBJECTIVES)],
    MIN_L2_VAL = MIN_L2_VAL
  )
  }
}



map_dbl(.x = c(5,10,15,20), ~ PCC_Est_K(DAT[[1]], .x,"claytonCopula",c(0.001,10)))


SUPPORTS = readxl::read_excel("~/Desktop/THESIS/Thesis_Computing/DATA/02A_Copula_Supports.xlsx")



#Test Concensus Method
Test_Concensus  =function(nsamp, copulaType, nquantiles,Param, n_to_test){
  FILT_SUP = SUPPORTS[which(SUPPORTS$COPULA==copulaType),]
  #For Now, we are going to only use the positive support for positive, and vice versa
  if(Param < 0 & copulaType %in% c("frankCopula","claytonCopula")){
    SUP = eval(parse_expr(FILT_SUP$NEG))
  }
  else{
    SUP = eval(parse_expr(FILT_SUP$POS))
  }
  #Allowing a larger upper range
  CUTS = unique(floor(seq(from = 5, to = nsamp/(2*nquantiles), length.out = n_to_test)))
  #Proeeding with the evenly-spaced intervals
  DAT =  Gen.Copula.Dat(          size        = nsamp, 
                                  parameter   = Param,
                                  CopulaType  = copulaType,
                                  x_quantiles = map_dbl(1:(nquantiles-1), ~ .x*1/nquantiles),
                                  y_quantiles = map_dbl(1:(nquantiles-1), ~ .x*1/nquantiles))
  CONCENSUS = mean(map_dbl(.x = CUTS, .f = ~ PCC_Est_K(DAT, .x, copulaType, SUP)))
  RES = data.frame(
    COPULA = copulaType,
    SIZE = nsamp,
    QUANTILES = nquantiles,
    TRUE_PARAM= Param,
    NUMBER_AVGED = n_to_test,
    CONCENSUS = CONCENSUS
  )
  return(RES)
}

args = list(
  BINS = c(3,5,7),
  SIZE = c(100,500),
  PARAM = c(-4,-2,2,4),
  IT = c(1:2)
)
crossedArgs = cross_df(args)
library(furrr)
options("future.fork.enable"= TRUE)
future::plan(multiprocess(workers = 5))
X = future_pmap_dfr(list(..1 = crossedArgs$SIZE, ..2 = crossedArgs$BINS, ..3 = crossedArgs$PARAM, ..4 = crossedArgs$IT),
                .f = ~ Test_Concensus(..1,"frankCopula",..2,..3,10))


map_dfr(1:10, ~Test_Concensus(100, "claytonCopula", 5,2,10))





#Generate Data Under Whatever Specs
PS_COMP_METHODS = function(nsamp, copulaType, nquantiles,Param, n_to_test,Crosses){
  FILT_SUP = SUPPORTS[which(SUPPORTS$COPULA==copulaType),]
  #For Now, we are going to only use the positive support for positive, and vice versa
  if(Param < 0 & copulaType %in% c("frankCopula","claytonCopula")){
    SUP = eval(parse_expr(FILT_SUP$NEG))
  }
  else{
    SUP = eval(parse_expr(FILT_SUP$POS))
  }
  #Determining Ad-Hoc Number of Bins
  K_ad_hoc = nsamp/(2*nquantiles)
  #Number of Cuts to try for Min_L2 and CV
  CUTS = unique(floor(seq(from = 5, to = nsamp/(2*nquantiles), length.out = n_to_test))) # Approximately n_to_test 
  #Generate Data
  DAT =  Gen.Copula.Dat(          size        = nsamp, 
                                  parameter   = Param,
                                  CopulaType  = copulaType,
                                  x_quantiles = map_dbl(1:(nquantiles-1), ~ .x*1/nquantiles),
                                  y_quantiles = map_dbl(1:(nquantiles-1), ~ .x*1/nquantiles))
  
  AH = data.frame(AH = PCC_Est_K(DAT,K_ad_hoc,copulaType,SUP, RETURN = "min"))
  L2 =  MIN_L2(DAT,CUTS,copulaType,SUP)
  CV = Cross_Validate_PCC(DAT, CUTS, copulaType, SUP, N_Cross = Crosses)
  bind_cols(AH,L2,CV) %>%
    mutate(COPULA = copulaType,
            SIZE = nsamp,
            QUANTILES = nquantiles)
}


         
#-----------Look at Preliminary Results
#Let's Read in New Files
FILES = list.files("~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/02_Polyserial/Recover_Param/New_Methods/FINAL/FINAL/FINALFINAL/")
PS_New = map_dfr(FILES, .f  = function(.x){
  library(dplyr)
  DAT = read.table(str_c("~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/02_Polyserial/Recover_Param/New_Methods/FINAL/FINAL/FINALFINAL/",.x))
  if(colnames(DAT)[1] == "V1"){
    colnames(DAT) = c("ROWNUM","AH","MIN_L2_CUTS","MIN_L2_VAL","CV_CUTS","CV_VAL","COPULA","SIZE","QUANTILES")
    DAT = select(DAT,-ROWNUM)
  }
  DAT = DAT %>%
    mutate(TRUE_PARAM = as.numeric(str_extract(.x, "[-]{0,1}[0-9]{1}\\.{0,1}[0-9]{0,1}")))
  return(DAT)
}) %>%
  pivot_longer(cols = c(AH, MIN_L2_VAL, CV_VAL),names_to = "METHOD")


#Add on New Data
PCC_Base = map_dfr(.x = c("Clayton", "Frank"),
                            .f = ~ read.table(str_c("~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/01_Test_Normal/01_Recover_Param/",.x))) %>%
  dplyr::rename(COPULA = TRUE_COP,
         value = RECOVERED_PARAM) %>%
  select(COPULA,SIZE,QUANTILES,value, TRUE_PARAM) %>%
  mutate(METHOD = "PCC")

#Add on First Method:
PS_BadIdea = map_dfr(.x = c("Clayton_Final", "Frank_Final"),
                     .f = ~ read.table(str_c("/Users/jacobshuman/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/02_Polyserial/Recover_Param/FINAL_PS_METHOD/",.x))) %>%
  dplyr::rename(value = RECOVERED_VALUE, COPULA = REAL_TYPE) %>%
  select(value, TRUE_PARAM, SIZE, COPULA, QUANTILES) %>%
  mutate(METHOD = "PS_BAD")
head(PS_New)
#Inverse Rho
PS_Irho = read.table("~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/02_Polyserial/Recover_Param/IRHO") %>%
  mutate(METHOD = "IRHO")
library(purrr)
FILES = list.files("~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/02_Polyserial/Recover_Param/Concensus/")
PS_Concensus = map_dfr(.x = FILES,
                       .f = function(.x){
                         DAT = read.table(str_c("~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/02_Polyserial/Recover_Param/Concensus/", .x))
                         colnames(DAT) = c("IT", "COPULA", "SIZE", "QUANTILES", "TRUE_PARAM", "TESTED_n","value")
                         return(DAT)
                       }) %>%
  select(-c(IT,TESTED_n)) %>%
  mutate(METHOD = "Concensus")

PS_Compare = PS_New %>%
  select(METHOD, value, SIZE, QUANTILES, COPULA, TRUE_PARAM) %>%
  rbind(PCC_Base, PS_BadIdea, PS_Irho, PS_Concensus)
head(PS_Compare)


#First, Convince People that the Two Continuous Method Doesn't Work

map(.x = c("IRHO", "PS_BAD"),
    .f = function(.x){
      MEANS = PS_Compare %>%
        filter(METHOD == .x) %>%
        dplyr::group_by(SIZE,QUANTILES,COPULA,TRUE_PARAM) %>%
        dplyr::summarize(MEAN = mean(value, na.rm = TRUE)) %>%
        pivot_wider(names_from = c(SIZE,QUANTILES), values_from = MEAN)  %>%
        as.matrix()
      SDS = PS_Compare %>%
        filter(METHOD == .x) %>%
        dplyr::group_by(SIZE,QUANTILES,COPULA,TRUE_PARAM) %>%
        dplyr::summarize(SDS = sd(value, na.rm = TRUE)) %>%
        pivot_wider(names_from = c(SIZE,QUANTILES), values_from = SDS)  %>%
        as.matrix()
      RES = matrix(nrow = nrow(MEANS), ncol = ncol(MEANS))
      RES[,1:2] = MEANS[,1:2]
      for(i in 1:nrow(RES)){
        for(j in 3:ncol(RES)){
          RES[i,j] = str_c(sprintf("%.2f", as.numeric(MEANS[i,j])), "(", sprintf("%.2f", as.numeric(SDS[i,j])), ")")
        }
      }
      colnames(RES) = c("COPULA", "THETA", "3","5","7","3","5","7")
      RES = xtable(RES, caption = str_c("Biased Parameter Estimated Using", .x, " Method"))
      print(RES, file = "~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/02_Polyserial/Recover_Param/FORMATTED/BAD_METHODS")
    })

            
RES
library(tidyverse)


PS_Compare %>%
  filter(METHOD == "IRHO" & QUANTILES == 5) %>%
  ggplot() +
  geom_boxplot(aes(x= factor(SIZE), y = value)) +
  facet_grid(~c(TRUE_PARAM,COPULA))
unique(PS_Compare$METHOD)
#Bad Methods
PLOT_VALS = data.frame(QUANTILES = rep(c(3,5,7), 4),
                       TRUE_PARAM = c(rep(-0.7,3),rep(-0.2,3),rep(2,3),rep(5,3)),
                       Z= c(rep(-0.7,3),rep(-0.2,3),rep(2,3),rep(5,3)))
PS_Compare %>%
  filter(METHOD %in% c("PS_BAD","IRHO","PCC"),
         SIZE == 100,
         COPULA == "claytonCopula") %>%
  mutate(METHOD = ifelse(METHOD == "PS_BAD", "PS",METHOD)) %>%
  ggplot() +
  geom_boxplot(aes(x = METHOD, y = value)) +
  geom_hline(data = PLOT_VALS, aes(yintercept = Z), col = "red") +
  theme_bw() +
  facet_grid(cols = vars(QUANTILES), rows = vars(TRUE_PARAM), scales = "free_y") +
  xlab("Method") + 
  ylab("Recovered Parameter Values") +
  ggtitle("Comparison of Parameter Estimates for Ordinal-Continuous")


#Now we can look at the good properties:
unique(PS_Compare$METHOD)

for(COP in c("frankCopula", "claytonCopula")){
  PS_Comp_Means = PS_Compare %>%
    filter(METHOD %in% c("AH", "MIN_L2_VAL", "CV_VAL", "PCC"),
           COPULA == COP) %>%
    group_by(SIZE, QUANTILES,COPULA, TRUE_PARAM, METHOD) %>%
    summarize(MEAN = mean(value)) %>%
    pivot_wider(names_from = c(METHOD,SIZE), values_from = MEAN) %>%
    arrange(COPULA)  %>%
    as.matrix()
  
  PS_Comp_Sd = PS_Compare %>%
    filter(METHOD %in% c("AH", "MIN_L2_VAL", "CV_VAL", "PCC")) %>%
    group_by(SIZE, QUANTILES,COPULA, TRUE_PARAM, METHOD) %>%
    summarize(MEAN = sd(value)) %>%
    pivot_wider(names_from = c(METHOD,SIZE), values_from = MEAN) %>%
    arrange(COPULA)  %>%
    as.matrix()
  
  RES = matrix(nrow = nrow(PS_Comp_Means), ncol = ncol(PS_Comp_Means))
  RES[,1:3] = PS_Comp_Means[,1:3]
  for(i in 1:nrow(RES)){
    for(j in 4:ncol(RES)){
      RES[i,j] = str_c(sprintf("%.2f", as.numeric(PS_Comp_Means[i,j])), "(", sprintf("%.2f", as.numeric(PS_Comp_Sd[i,j])), ")")
    }
  }
  RES = as.data.frame(RES) 
  colnames(RES) = colnames(PS_Comp_Means) 
  RES = select(RES, -COPULA)
  RES = xtable(RES, caption = str_c("Parameter Estimates for ",COP,"With Four Methods"), size = "small")
  print(RES, file = "~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/02_Polyserial/Recover_Param/FORMATTED/BoxMethods.tex", append = TRUE)
}


for(COP in c("frankCopula", "claytonCopula")){
  PS_Comp_Means = PS_Compare %>%
    filter(METHOD %in% c("AH", "Concensus", "PCC"),
           COPULA == COP,
           value > -29) %>%
    group_by(SIZE, QUANTILES,COPULA, TRUE_PARAM, METHOD) %>%
    summarize(MEAN = mean(value)) %>%
    pivot_wider(names_from = c(METHOD,SIZE), values_from = MEAN) %>%
    arrange(COPULA)  %>%
    as.matrix()
  
  PS_Comp_Sd = PS_Compare %>%
    filter(METHOD %in% c("AH", "Concensus", "PCC"),
           COPULA == COP, value > -29) %>%
    group_by(SIZE, QUANTILES,COPULA, TRUE_PARAM, METHOD) %>%
    summarize(MEAN = sd(value)) %>%
    pivot_wider(names_from = c(METHOD,SIZE), values_from = MEAN) %>%
    arrange(COPULA)  %>%
    as.matrix()
  
  RES = matrix(nrow = nrow(PS_Comp_Means), ncol = ncol(PS_Comp_Means))
  RES[,1:3] = PS_Comp_Means[,1:3]
  for(i in 1:nrow(RES)){
    for(j in 4:ncol(RES)){
      RES[i,j] = str_c(sprintf("%.2f", as.numeric(PS_Comp_Means[i,j])), "(", sprintf("%.2f", as.numeric(PS_Comp_Sd[i,j])), ")")
    }
  }
  RES = as.data.frame(RES) 
  colnames(RES) = colnames(PS_Comp_Means) 
  RES = select(RES, -COPULA)
  RES = xtable(RES, caption = str_c("Parameter Estimates for ",COP,"With Four Methods"), size = "small")
  print(RES, file = "~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/02_Polyserial/Recover_Param/FORMATTED/Concensus.tex", append = TRUE)
}

PS_Compare %>%
  filter(METHOD %in% c("Concensus", "PCC"), 
         COPULA == "frankCopula", 
         SIZE == 500) %>%
  ggplot() +
  geom_boxplot(aes(x = METHOD, y = value)) +
  facet_grid(rows = vars(TRUE_PARAM), cols = vars(QUANTILES), scales = "free_y")


for(i in 1:100){
  set.seed(i)
  X = (Test_Concensus(100, "frankCopula", 5,-4,10)$CONCENSUS)
  if(X < -10){
    print(i)
    print(X)
  }
}

set.seed(85)
Test_Concensus(nsamp = 100,copulaType = "frankCopula",nquantiles = 5, Param = -4, n_to_test = 10)


FILT_SUP = SUPPORTS[which(SUPPORTS$COPULA=="frankCopula"),]
#For Now, we are going to only use the positive support for positive, and vice versa
  SUP = eval(parse_expr(FILT_SUP$NEG))
CUTS = unique(floor(seq(from = 5, to = 100/(2*5), length.out = 10)))
set.seed(85)
#Proeeding with the evenly-spaced intervals
DAT =  Gen.Copula.Dat(          size        = 100, 
                                parameter   = -4,
                                CopulaType  = "frankCopula",
                                x_quantiles = map_dbl(1:(5-1), ~ .x*1/5),
                                y_quantiles = map_dbl(1:(5-1), ~ .x*1/5))

mean(map_dbl(.x = CUTS, .f = ~ PCC_Est_K(DAT,.x,"frankCopula",c(-1000,-0.001))))

head(SUPPORTS)






#More Robustness 

PATH = "~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/02_Polyserial/Recover_Param/Robustness/"
FILES = list.files(PATH)

PS_GOOD_ROB = map_dfr(.x = FILES,
                      .f = function(.x){
                        DAT = read.table(str_c(PATH, .x)) %>%
                          select(-V1)
                        if(str_detect(.x,"CONC")){
                          DAT = DAT %>%
                            select(-V6)
                          colnames(DAT) = c("COPULA","SIZE","QUANTILES","TRUE_PARAM","value")
                          DAT = mutate(DAT, method = "CONC")
                        }
                        else{
                          colnames(DAT) = c("value","COPULA","SIZE","QUANTILES","TRUE_PARAM")
                          DAT = mutate(DAT, method = "AH")
                        }
                        return(DAT)
                      })

#Create Robustness Table
PS_GOOD_ROB %>%
  filter(SIZE == 100,
         COPULA == "frankCopula",
         value > -30) %>%
  ggplot() +
  geom_boxplot(aes(x = factor(method), y = value)) +
  facet_grid(rows = vars(TRUE_PARAM), cols = vars(QUANTILES), scales = "free_y")

head(PS_GOOD_ROB)


MEANS = PS_GOOD_ROB %>%
  filter(value > -30) %>%
  group_by(COPULA,SIZE,QUANTILES,TRUE_PARAM,method) %>%
  summarize(MEAN = mean(value)) %>%
  pivot_wider(names_from = c(SIZE,method), values_from = MEAN)  %>%
  as.matrix()
SD = PS_GOOD_ROB %>%
  filter(value > -30) %>%
  group_by(COPULA,SIZE,QUANTILES,TRUE_PARAM,method) %>%
  summarize(MEAN = sd(value)) %>%
  pivot_wider(names_from = c(SIZE,method), values_from = MEAN)  %>%
  as.matrix()
RES = matrix(nrow = nrow(MEAN), ncol = ncol(MEAN))
RES[,1:3] = MEANS[,1:3]
for(i in 1:nrow(MEANS)){
  for(j in 4:7){
    RES[i,j] = str_c(sprintf("%.2f", as.numeric(MEANS[i,j])), "(", sprintf("%.2f", as.numeric(SD[i,j])), ")")
  }
}
colnames(RES) = c("COPULA","BINS","PARAM","AH","CONC","AH","CONC")
RES = xtable(RES, caption = "Sensitivity to Different Boundary Definitions",
             label = "tab:PS_rob",
             align = "llrrrrrr")
print(RES, file = "~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/02_Polyserial/Recover_Param/FORMATTED/PS_good_rob.tex")
head(RES)

  









ggplot(PS_Compare %>% 
         filter(COPULA == "frankCopula" & TRUE_PARAM == -2 & value > -400)) +
  geom_boxplot(aes(x = METHOD, y = value)) +
 
  facet_grid(cols = vars(SIZE), rows = vars(QUANTILES)) +
  geom_hline(aes(yintercept = -2))

CLAYTON = read.table("~/Desktop/CLAYTON_-0.7")
head(CLAYTON)

CLAYTON %>%
  pivot_longer(cols = c(AH,MIN_L2_VAL,CV_VAL)) %>%
  ggplot() +
  geom_boxplot(aes(x = name, y = value)) + 
  facet_grid(rows = vars(SIZE), cols = vars(QUANTILES))

FRANK %>%
  filter(CV_VAL > -100) %>%
  group_by(SIZE,QUANTILES) %>%
  summarize(SD = sd(CV_VAL))

FRANK = read.table("~/Desktop/FRANK_-4")

FRANK %>%
  pivot_longer(cols = c(AH,MIN_L2_VAL,CV_VAL)) %>%
  #filter(value >- 100) %>%
  ggplot() +
  geom_boxplot(aes(x = name, y = value)) + 
  facet_grid(rows = vars(SIZE), cols = vars(QUANTILES))
  




PS_Compare %>%
  filter(COPULA == "frankCopula") %>%
  distinct(TRUE_PARAM)











