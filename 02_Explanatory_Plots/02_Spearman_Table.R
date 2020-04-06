#--------------------------------------------Replicating Figure From Pearson, 1907------------------------------#
library(xtable)
library(ggplot2)
library(dplyr)

V1 = c(-2,-1,1,2)
V2 = c(-2,-1,1,2)
G1 = c(-2,-1.99,1.99,2)
G2 = c(-2,-0.01,0,0.01)
cor(G1,G2,method = "spearman")
cor(G1,G2)
#Matrix of Variates
Var = matrix(c(V1,V2,G1,G2),nrow = 4)
rownames(Var) = c("X1","Y1","X2","Y2")
colnames(Var) = 1:4
xtable(Var,caption = "Continuous Variables With Same Grade Correlation")

V1 = c(-2,-1.99,1.99,2)
V2 = c(-1,-0.01,0.01,0.02)
cor(V1,V2)


#Making a plot to demonstrate discretization:

X = seq(-3,3,by = 0.001)
DAT = X %>%
  data.frame() %>%
  mutate(Y = dnorm(X,0,1),
         Category = case_when(
           X < -1.5 ~ 1,
           X >= -1.5 & X <0 ~ 2,
           X >= 0 & X <1.5 ~ 3,
           X >= 1.5 ~ 4)) %>%
  rename(X = ".")
head(DAT)
plot(X,Y)
ggplot() +
  geom_line(aes(x = DAT$X,y = DAT$Y)) +
  geom_segment(aes(x=c(-1.5,0,1.5),xend = c(-1.5,0,1.5),
                   y = c(0,0,0), yend= dnorm(c(-1.5,0,1.5))),col = "darkblue") +
  xlab("X Star") + ylab("Density") +
  ggtitle("Discretized Version of X Star")
  


