#Load Libraries
library(tidyverse)
library(purrr)
library(mvtnorm)
library(xtable)

#--------------------------------Analyzing Bootstrap Data--------------------#

ALL_PCC_BS = map_dfr(.x = as.character(c(3,5,7)),
    .f = ~read.table(str_c("~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/01_Test_Normal/Bootstrap/SIM_BS_", .x)))
nrow(ALL_PCC_BS)

#First, Summary Statistics
ALL_PCC_BS %>%
  group_by(SIZE,REAL_TYPE,QUANTILES) %>%
  summarise_if(is.numeric, ~sd(.))
glimpse(ALL_PCC_BS)
#For Large
ALL_PCC_BS %>%
  filter(SIZE == 500, QUANTILES = 7) %>%
  ggplot() +
  geom_histogram(aes(x=))

#First, Melting into Really Long Table
ALL_PCC_ID = map_dfr(.x = c("SIM_ID_3_FINAL", "02_SIM_ID_5","02_SIM_ID_7"),
                     .f  = ~ read.table(str_c("~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/01_Test_Normal/Indep_Tables/", .x))) %>% 
  pivot_longer(-c("SIZE","REAL_TYPE","QUANTILES")) %>%
  mutate(METHOD = str_extract(name, "[A-Z]{1}"),
         TRUE_COR = str_extract(name, "_(.*)")) %>%
  mutate(TRUE_COR = str_sub(TRUE_COR,2,)) %>%
  mutate(TRUE_COR = case_when(
    TRUE_COR == ".0.75" ~ "-0.75",
    TRUE_COR == ".0.5"  ~ "-0.5",
    TRUE_COR == ".0.25" ~ "-0.25",
    TRUE_COR %in% c("0", "0.25", "0.5", "0.75") ~ TRUE_COR)) %>% # This is a mess :(
  select(-name) 


#------------------------------Full Table of Results-----------------------------#
for(i in c(3,5,7)){
  for(j in c(100,500)){
    print("\n", file = "~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/01_Test_Normal/Indep_Tables/ALL_RES.tex", 
          type = "latex", append = TRUE)
    FILT = ALL_PCC_ID %>% filter(SIZE == j, QUANTILES == i)
    ALL_MEAN = map_dfc(.x = c(-0.75,-0.25,0,0.25,0.75),
                       .f = function(.x){
                         TRUE_COR = .x
                         FILT %>%
                           mutate(REAL = ifelse(REAL_TYPE %in% c("Nelsen_15_GG","Nelsen_2"), str_c(REAL_TYPE,"Copula",sep="_"),as.character(REAL_TYPE))) %>%
                           filter(TRUE_COR == .x) %>%
                           aggregate(data = ., value ~ SIZE + QUANTILES + REAL + METHOD, mean)}) %>%
      select(c(1:5,10,15,20,25)) %>%
      arrange(REAL) %>%
      rename(`-0.75` = value, `-0.25` = value1, `0` = value2, `0.25` = value3, `0.75` = value4) %>%
      mutate_at(vars(5:9),~round(.,2))
    ALL_SD = map_dfc(.x = c(-0.75,-0.25,0,0.25,0.75),
                     .f = function(.x){
                       TRUE_COR = .x
                       FILT %>%
                         mutate(REAL = ifelse(REAL_TYPE %in% c("Nelsen_15_GG","Nelsen_2"), str_c(REAL_TYPE,"Copula",sep="_"),as.character(REAL_TYPE))) %>%
                         filter(TRUE_COR == .x) %>%
                         aggregate(data = ., value ~ SIZE + QUANTILES + REAL + METHOD, sd)}) %>%
      select(c(1:5,10,15,20,25)) %>%
      arrange(REAL) %>%
      rename(`-0.75` = value, `-0.25` = value1, `0` = value2, `0.25` = value3, `0.75` = value4) %>%
      mutate_at(vars(5:9),~round(.,2))
    FULL = matrix(nrow = 15,ncol = 7)
    FULL[,1] = str_remove_all(ALL_MEAN$REAL, paste(c("Copula","_","GG"),collapse = "|"))
    FULL[,2] = ALL_MEAN$METHOD
    for(a in 1:15){
      for(b in 3:7){
        FULL[a,b] = str_c(sprintf('%.2f',ALL_MEAN[a,b+2]),"(",sprintf('%.2f',ALL_SD[a,b+2]),")")
      }
    }
    colnames(FULL) = c("True Cop","Method","-0.75","-0.25","0","0.25","0.75")
    RES = xtable(FULL, caption = str_c("Polychoric Correlation Estimated Using Three Methods For Sample Size ", j, "\n and contingency table size ", i," x ", i)) 
    print(RES)
    print(RES, file = "~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/01_Test_Normal/Indep_Tables/ALL_RES.tex", 
          type = "latex", append = TRUE)
  }
}



levels(ALL_PCC_ID$REAL_TYPE) = c("Clayton","Frank","Nelsen 15", "Nelsen 2", "Normal")
ALL_PCC_ID$QUANTILES = as.factor(ALL_PCC_ID$QUANTILES)
levels(ALL_PCC_ID$QUANTILES) = c("3x3","5x5","7x7")

#Creating Figure With Lots of Boxplots for N=100
library(plyr)
map(.x = c(3,5,7),
    .f = function(.x){
      PLOT_VALS = data.frame(TRUE_COR = c(rep(-0.25, 1), rep(-0.75,1),rep(0,1),rep(0.25,1),rep(0.75,1)),  
                             QUANTILES = rep(c(.x),5),
                             Z = c(rep(-0.25, 1), rep(-0.75,1),rep(0,1),rep(0.25,1),rep(0.75,1)))
      ALL_PCC_ID %>%
        filter(METHOD == "N",
               SIZE == 100,
               QUANTILES == .x,
               !TRUE_COR %in% c(-0.5,0.5)) %>%
        mutate(REAL_TYPE = case_when(
          REAL_TYPE == "normalCopula" ~ "Normal", 
          REAL_TYPE ==     "frankCopula" ~ "Frank", 
          REAL_TYPE == "claytonCopula" ~ "Clayton",
          REAL_TYPE == "Nelsen_2" ~ "N_2",
          REAL_TYPE == "Nelsen_15_GG" ~ "N_15")) %>%
        mutate(TRUE_COR = as.numeric(TRUE_COR)) %>%
        ggplot() +
        geom_boxplot(aes(x = REAL_TYPE, y = value)) +
        geom_hline(data = PLOT_VALS, aes(yintercept = Z), lty = 2, col = "red") +
        facet_grid(cols = vars(TRUE_COR), rows = vars(QUANTILES)) + 
        xlab("True Copula") +
        ylab("Recovered Correlation Value") +
        ggtitle(str_c("Distribution of PCC Estimates, N = 100, Bins = ", .x)) +
        theme_bw() 
    })


#Creating Figure for 100 v. 500 at 7x7 cotingency table
PLOT_VALS = data.frame(TRUE_COR = c(rep(-0.25, 5), rep(-0.75,5),rep(0,5),rep(0.25,5),rep(0.75,5)),  
                       REAL_TYPE = rep(c("Clayton", "Frank", "N_15", "N_2", "Normal"),5),
                       Z = c(rep(-0.25, 5), rep(-0.75,5),rep(0,5),rep(0.25,5),rep(0.75,5)))
ALL_PCC_ID %>%
  filter(METHOD == "N",
         SIZE %in% c(100,500),
         !TRUE_COR %in% c(-0.5,0.5),
         QUANTILES == 7) %>%
  mutate(REAL_TYPE = case_when(
    REAL_TYPE == "normalCopula" ~ "Normal", 
    REAL_TYPE ==     "frankCopula" ~ "Frank", 
    REAL_TYPE == "claytonCopula" ~ "Clayton",
    REAL_TYPE == "Nelsen_2" ~ "N_2",
    REAL_TYPE == "Nelsen_15_GG" ~ "N_15")) %>%
  mutate(TRUE_COR = as.numeric(TRUE_COR)) %>%
  ggplot() +
  geom_boxplot(aes(x = factor(SIZE), y = value)) +
  geom_hline(data = PLOT_VALS, aes(yintercept = Z), lty = 2, col = "red") +
    facet_grid(rows = vars(TRUE_COR), cols = vars(REAL_TYPE), scales = "free_y") +
  xlab("Sample Size") +
  ylab("Recocvered Correlation Value") +
  ggtitle("Distribution of PCC Estimates By True Correlation, Copula, and Sample Size") +
  theme_bw()
  





#Plots of Normal Method -- Exp[loratory
USECOR = -0.25
map(.x = c(50,100,500),
    .f = function(.x){
      FILT = ALL_PCC_ID %>%
        filter(METHOD == "N",
               SIZE == .x,
               TRUE_COR == USECOR)
      TRUE_PARAMS = FILT %>%
        group_by(REAL_TYPE,QUANTILES) %>%
        summarize(MED = mean(value))
      FILT %>%
        ggplot() +
        geom_histogram(aes(x = value), bins = 20, fill = "black") +
        geom_vline(aes(xintercept = USECOR, color = "Population")) +
        geom_vline(data = TRUE_PARAMS, aes(xintercept = MED, color = "Empirical")) + 
        ggtitle(str_c("Distribution of Polychoric Correlation Estimates, N = ", as.character(.x))) +
        scale_color_manual(name = "LEGEND:", values = c(Population = "red", Empirical = "blue")) +
        theme_minimal()+
        theme(legend.position = "bottom") +
        xlab("Value") +
        ylab("Count") +
        facet_grid(cols = vars(REAL_TYPE), rows = vars(QUANTILES))
    } )

#Producing Actual Output
FILT = ALL_PCC_ID %>%
  filter(METHOD == "N",
         SIZE == 100,
         TRUE_COR == USECOR)
TRUE_PARAMS = FILT %>%
  group_by(REAL_TYPE,QUANTILES) %>%
  summarize(MED = mean(value))
FILT %>%
  ggplot() +
  geom_histogram(aes(x = value), bins = 20, fill = "black") +
  geom_vline(aes(xintercept = USECOR, color = "Population")) +
  geom_vline(data = TRUE_PARAMS, aes(xintercept = MED, color = "Empirical")) + 
  ggtitle(str_c("Distribution of Polychoric Correlation Estimates, N = 100")) +
  scale_color_manual(name = "LEGEND:", values = c(Population = "red", Empirical = "blue")) +
  theme_minimal()+
  theme(legend.position = "bottom") +
  xlab("Value") +
  ylab("Count") +
  facet_grid(cols = vars(REAL_TYPE), rows = vars(QUANTILES))



#Robustness
ROB3 = read.csv("~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/01_Test_Normal/Indep_Tables/Robustness/ROB/ROB_3BINS", sep="")
ROB5 = read.csv("~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/01_Test_Normal/Indep_Tables/Robustness/ROB/ROB_5BINS", sep="")

Table_Rob = function(DAT){
  INT = DAT %>%
    pivot_longer(-c("SIZE","REAL_TYPE","QUANTILES"))%>%
    mutate(METHOD = str_extract(name, "[A-Z]{1}"),
           TRUE_COR = str_extract(name, "_(.*)")) %>%
    mutate(TRUE_COR = str_sub(TRUE_COR,2,)) %>%
    mutate(TRUE_COR = case_when(
      TRUE_COR == ".0.75" ~ "-0.75",
      TRUE_COR == ".0.5"  ~ "-0.5",
      TRUE_COR == ".0.25" ~ "-0.25",
      TRUE_COR %in% c("0", "0.25", "0.5", "0.75") ~ TRUE_COR)) %>%
    mutate(REAL_TYPE = case_when(
      REAL_TYPE == "normalCopula" ~ "Normal",
      REAL_TYPE == "claytonCopula" ~"Clayton",
      REAL_TYPE == "frankCopula"  ~ "Frank",
      REAL_TYPE == "Nelsen_15_GG" ~ "N_15",
      REAL_TYPE == "Nelsen_2" ~ "N_2")) %>%
    select(-name)
MEANS = INT %>%
  filter(METHOD == "N") %>%
  aggregate(data = ., value ~ SIZE + REAL_TYPE + TRUE_COR, mean) %>%
  pivot_wider(names_from = TRUE_COR, values_from = value) %>%
  arrange(SIZE) 
MEANS  = MEANS[,c(1,2,5,4,3,6,7,8,9)]
MEANS = as.matrix(MEANS)
SDS = INT %>%
  filter(METHOD == "N") %>%
  aggregate(data = ., value ~ SIZE + REAL_TYPE + TRUE_COR, sd) %>%
  pivot_wider(names_from = TRUE_COR, values_from = value) %>%
  arrange(SIZE) 
SDS  = SDS[,c(1,2,5,4,3,6,7,8,9)]
SDS = as.matrix(SDS)
RES = matrix(nrow = nrow(MEANS), ncol = ncol(MEANS))
RES[,1:2] = MEANS[,1:2]
for(i in 1:10){
  for(j in 3:9){
    RES[i,j] =  str_c(sprintf("%.2f", as.numeric(MEANS[i,j])), "(", sprintf("%.2f", as.numeric(SDS[i,j])), ")")
  }
}
colnames(RES) = c("SIZE","COPULA", "-0.75","-0.5","-0.25,","0","0.25","0.5","0.75")
RES = xtable(RES)

print(RES, file = "~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/01_Test_Normal/Indep_Tables/Robustness/ROB.tex", append = TRUE)
}

Table_Rob(ROB3)
Table_Rob(ROB5)

lapply(X = list(ROB3,ROB5), FUN = Table_Rob(X))

    
    
    .f = function(.x){
      INT = .x %>%
       
      #Create Table
      


 

head(ROB3)



as.numeric(MEANS[10,9])

RES


PLOT_VALS = data.frame(TRUE_COR = c(rep(-0.25, 5), rep(-0.75,5),rep(0,5),rep(0.25,5),rep(0.75,5)),  
                       REAL_TYPE = rep(c("Clayton", "Frank", "N_15", "N_2", "Normal"),5),
                       Z = c(rep(-0.25, 5), rep(-0.75,5),rep(0,5),rep(0.25,5),rep(0.75,5)))
ROB3 %>%
  filter(METHOD == "N",
         SIZE %in% c(100,500),
         !TRUE_COR %in% c(-0.5,0.5)) %>%
  mutate(REAL_TYPE = case_when(
    REAL_TYPE == "normalCopula" ~ "Normal", 
    REAL_TYPE ==     "frankCopula" ~ "Frank", 
    REAL_TYPE == "claytonCopula" ~ "Clayton",
    REAL_TYPE == "Nelsen_2" ~ "N_2",
    REAL_TYPE == "Nelsen_15_GG" ~ "N_15")) %>%
  mutate(TRUE_COR = as.numeric(TRUE_COR)) %>%
  ggplot() +
  geom_boxplot(aes(x = factor(SIZE), y = value)) +
  geom_hline(data = PLOT_VALS, aes(yintercept = Z), lty = 2, col = "red") +
  facet_grid(rows = vars(TRUE_COR), cols = vars(REAL_TYPE), scales = "free_y") +
  xlab("Sample Size") +
  ylab("Recocvered Correlation Value") +
  ggtitle("Distribution of PCC Estimates By True Correlation, Copula, and Sample Size") +
  theme_bw()
  
















#Plots Comparing Normal to Generalized Method
####01. For Normal Copula
PLOT_VALS = data.frame(REAL_TYPE = rep("normalCopula", 4), REAL = c("-0.75", "-0.25", "0.25","0.75"), 
                       Z = c(-0.75,-0.25,0.25,0.75))
ALL_PCC_ID %>%
  filter(METHOD %in% c("N", "G"),
         SIZE == 100, 
         QUANTILES == 5,
         REAL_TYPE == "normalCopula",
         TRUE_COR %in% c("-0.75", "-0.25", "0.25","0.75")) %>%
  mutate(REAL = case_when(
    TRUE_COR == "-0.25" ~ -0.25, TRUE_COR == "-0.75" ~ -0.75, 
    TRUE_COR == "0.25"  ~  0.25, TRUE_COR == "0.75"  ~ 0.75)) %>%
  arrange(REAL) %>%
  ggplot() +
  geom_boxplot(aes(x = METHOD, y= value, group = METHOD, col = METHOD)) +
  scale_color_manual(values = group.colors) +
  theme_bw() +
  ggtitle("Comparison of Normal and Generalized Gaussian Definitions") +
  facet_grid(cols = vars(as.factor(REAL)), rows = vars(REAL_TYPE)) +
  geom_hline(data = PLOT_VALS, aes(yintercept = Z), lty = 2, col = "black")
  


group.colors = c("N" = "blue4", "G" = "red3", "E" = "seagreen4")
PLOT_LINES = data.frame(REAL_TYPE = c(rep("Nelsen_2",4),rep("Nelsen_15_GG",4)),
                        REAL = rep(c(-0.75, -0.25, 0.25,0.75),2),
                        VALS = rep(c(-0.75, -0.25, 0.25,0.75),2))
PLOT_LINES
ALL_PCC_ID %>%
  filter(METHOD %in% c("N", "G"),
         SIZE == 500, 
         QUANTILES == 5,
         REAL_TYPE %in% c("Nelsen_2", "Nelsen_15_GG"),
         TRUE_COR %in% c("-0.75", "-0.25", "0.25","0.75")) %>%
  mutate(REAL = case_when(
    TRUE_COR == "-0.25" ~ -0.25, TRUE_COR == "-0.75" ~ -0.75, 
    TRUE_COR == "0.25"  ~  0.25, TRUE_COR == "0.75"  ~ 0.75)) %>%
  arrange(REAL) %>%
  ggplot() +
     geom_boxplot(aes(x = METHOD, y= value, group = METHOD, col = METHOD)) +
     scale_color_manual(values = group.colors) +
     theme_bw() +
     ggtitle("Comparison of Normal and Generalized Gaussian Definitions") +
     facet_grid(cols = vars(as.factor(REAL)), rows = vars(REAL_TYPE)) +
     geom_hline(data = PLOT_LINES, aes(yintercept = VALS), lty = 2, col = "black")

#Look at the Empirical Definition
ALL_PCC_ID %>%
  filter(SIZE == 100, QUANTILES == 5,
         TRUE_COR %in% c("-0.75", "-0.25", "0.25","0.75")) %>%
  mutate(REAL = case_when(
    TRUE_COR == "-0.25" ~ -0.25, TRUE_COR == "-0.75" ~ -0.75, 
    TRUE_COR == "0.25"  ~  0.25, TRUE_COR == "0.75"  ~ 0.75)) %>%
  arrange(REAL) %>%
  ggplot() +
  geom_boxplot(aes(x = METHOD, y= value, group = METHOD, col = METHOD)) +
  scale_color_manual(values = group.colors) +
  theme_bw() +
  ggtitle("Comparison of Normal and Generalized Gaussian Definitions") +
  facet_grid(cols = vars(as.factor(REAL)), rows = vars(REAL_TYPE)) +
  geom_hline(data = PLOT_LINES, aes(yintercept = VALS), lty = 2, col = "black")
  
  

ALL_PCC_ID %>%
  group_by(METHOD,SIZE,QUANTILES,TRUE_COR,REAL_TYPE) %>%
  summarise(MEAN = mean(value)) %>%
  filter(METHOD %in% c("N","G")) %>%
  arrange(REAL_TYPE, TRUE_COR,QUANTILES,SIZE) %>%
  View()








