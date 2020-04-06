#Load Libraries
library(copula)
library(ggplot2)
library(dplyr)
library(readxl)
library(rlang)
library(stringr)
library(purrr)
library(gridExtra)
library(pracma)
library(xtable)

#Reading in Copula Supports
SUPPORTS = read_xlsx("~/Desktop/THESIS/Thesis_Computing/DATA/02_Copula_Supports.xlsx")
#Generating Plots of Archimedian Copulas
Gen_Cop_Plots = function(Cop_Type, PPC, size=10000, PARAM =FALSE, title = "DEFAULT", R.S = FALSE){
  FILT_SUP = SUPPORTS[which(SUPPORTS$COPULA==Cop_Type),]
  NEG_SUP = eval(parse_expr(FILT_SUP$NEG))
  POS_SUP = eval(parse_expr(FILT_SUP$POS))
  #Optimization to Find Final Param
  if(PARAM == TRUE){
    FINAL_PARAM = PPC
  }
  else{
  if(Cop_Type == "normalCopula"){
    FINAL_PARAM = PPC
  }
  else{
    POS_OPT = optimize(Find_Param, interval = POS_SUP,
                       Cop_Type,PPC)
    if(length(NEG_SUP)==1){
      FINAL_PARAM = POS_OPT$minimum
    }
    else{
      NEG_OPT = optimize(Find_Param, interval = NEG_SUP,
                         Cop_Type,PPC)
      if(PPC<0){
        FINAL_PARAM =  NEG_OPT$minimum
      }
      else{
        FINAL_PARAM = POS_OPT$minimum
      }
    }
  }
  }
  #Using That to Generate Desired Data
  library(dplyr)
  DAT = Gen.Copula.Dat(size = size,
                       parameter = FINAL_PARAM,
                       CopulaType = Cop_Type,
                       x_quantiles = c(0.5),y_quantiles = c(0.5)) %>%
    dplyr::select(1:2)
  if(R.S == TRUE){
    return(DAT)
  }
  else{
  #Adding a Convex Hull
  HULL = DAT %>%
    slice(chull(x,y)) 
  if(title != "DEFAULT"){
    T_FINAL =title
  }
  else{
    T_FINAL = str_c("PLOT OF ", toupper(str_remove(Cop_Type,"Copula")), " COPULA \n PPC = ",as.character(PPC))
  }
  #Creating Plot
  ggplot() +
    geom_point(aes(x=DAT$x,y=DAT$y),col="black") +
    theme_minimal() +
    #geom_polygon(aes(x=HULL$x,y=HULL$y),alpha = 0.3)+
    xlab("X") + 
    ylab("Y") +
    ggtitle(T_FINAL)
}}

set.seed(10)


Gen_Cop_Plots("normalCopula", 0.8,PARAM=TRUE)

#Generating Explanatory Plots
map(.x = c("Nelsen_2", "claytonCopula","normalCopula"),
    .f = ~ Gen_Cop_Plots(.x, 0.9,title = ""))



RES = matrix(nrow = 2, ncol = 3)
rownames(RES) = c("Lower Quadrant", "Upper Quadrant")
colnames(RES) = c("Normal", "Clayton", "Nelsen-2")

COPS = c("normalCopula","claytonCopula","Nelsen_2")
for(i in 1:3){
  X = Gen_Cop_Plots(COPS[i], 0.75,R.S =TRUE)
  LQ  = X %>% filter(x < 0.25, y < 0.25)
  UQ  = X %>% filter(x > 0.75, y > 0.75)
  RES[1,i] = cor(LQ$x,LQ$y)
  RES[2,i] = cor(UQ$x,UQ$y)
}

RES = xtable(RES)
print(RES, file = "~/Desktop/THESIS/Thesis_Computing/Outputs/Cop_Examples/Comp_Quads.tex")

write.table(RES, 
cat(RES, "\n", file = )


X = Gen_Cop_Plots("Nelsen-2", 0.75,R.S =TRUE)
LQ = X %>% filter(x < 0.25, y <0.25)
MID = X %>% filter(x >= 0.25, x <=0.75, y >= 0.25, y<= 0.75)
UQ  = X %>% filter(x > 0.75, y > 0.75)
RES[1,1] = cor(LQ$x,LQ$y)
RES[2,1] = cor(MID$x,MID$y)
RES[3,1] = cor(UQ$x,UQ$y)
RES

cor(X$x,X$y)

head(X)

