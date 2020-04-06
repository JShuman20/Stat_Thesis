#Show That Inverse Rho is Inconsistent
library(furrr)
library(ggplot2)
IRHO = function(DAT,COPULA){
  Thresholds = Gen_Emp_Thresholds(DAT, quote(y_cat))
                 DAT = DAT %>%
                   select(x,y_cat) %>%
                   mutate(y_cat = as.numeric(Thresholds[y_cat]))
                 COR = cor(DAT$x,DAT$y_cat, method = "spearman")
                 if(COPULA == "frankCopula"){
                   PARAM = iRho(frankCopula(), COR)
                 }
                 else if(COPULA == "normalCopula"){
                   PARAM = iRho(normalCopula(),COR)
                 }
                 else if(COPULA == "claytonCopula"){
                   PARAM = iRho(claytonCopula(), COR)
                 }
                 return(PARAM)
}


Test_IRHO = function(SIZE, parameter, COP, QUANTILES,iterations){
  DAT = map(.x = 1:iterations,
            .f = ~ Gen.Copula.Dat(size        = SIZE, 
                                  parameter   = parameter,
                                  CopulaType  = COP,
                                  x_quantiles = map_dbl(1:(QUANTILES-1), ~ .x*1/QUANTILES),
                                  y_quantiles = map_dbl(1:(QUANTILES-1), ~ .x*1/QUANTILES)))
  VALS = map_dbl(DAT, ~ IRHO(.,COP))
  return(data.frame(
    value = VALS,
    SIZE =SIZE,
    QUANTILES = QUANTILES, 
    COPULA = COP,
    TRUE_PARAM = parameter
  ))
}

args_clayton = list(
  SIZE = c(100,500),
  QUANTILES = c(3,5,7),
  PARAM = c(-0.7,-0.2,2,5)
)

args_frank = list(
  SIZE = c(100,500),
  QUANTILES = c(3,5,7),
  PARAM = c(-4,-2,2,4)
)

ALL_PARAMS = rbind(cross_df(args_clayton) %>% mutate(COPULA = "claytonCopula"),
                   cross_df(args_frank) %>% mutate(COPULA = "frankCopula"))

library(furrr)
options("future.fork.enable" = TRUE)
future::plan(multicore(workers = 8))
future_pmap_dfr(list(..1 = ALL_PARAMS$SIZE, ..2 = ALL_PARAMS$QUANTILES, ..3 = ALL_PARAMS$PARAM, ..4 = ALL_PARAMS$COPULA),
                ~ Test_IRHO(SIZE = ..1, parameter = ..3, COP = ..4, QUANTILES = ..2, iterations = 1000)) %>%
  write.table("~/Desktop/THESIS/Thesis_Computing/Outputs/Simulation/02_Polyserial/Recover_Param/IRHO")

ggplot(X) +
  geom_histogram(aes(x = VALS),fill = "black") +
  facet_grid(rows = vars(TRUE_P), cols = vars(SIZE)) +
  geom_vline(aes(xintercept = 5)) +
  geom_vline(aes(xintercept = 8))





