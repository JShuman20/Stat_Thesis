library(cubature)
library(mvtnorm)
library(pracma)
library(tidyverse)
library(rlang)

#----------------------------------All Max-Based Functions Written Independently----------------------------------#
AMHCOPULA = list(
  COPULA = function(theta,u,v){
    val = (u*v)/(1 - theta*(1-u)*(1-v))
    return(val)
  },
  D_Dx = function(theta, u,v){
    val = (theta * (v-1)*v + v)/((1 - theta*(u-1)*(v-1))^2)
    return(val)
  },
  D_Dx_Dt = function(theta,u,v){
    val = -(2*(u-1)*(v-1))/(theta*(u-1)*(v-1) -1)^3 + ((v-1)*v)/((1 - theta*(u-1)*(v-1))^2)
    return(val)
  }
)
NORMALCOPULA = list(
  COPULA = function(theta, u,v){
    val = dmvnorm(c(u,v))
    
  }
)
CLAYTONCOPULA = list(
  COPULA = function(theta,u,v){
    val = u^(-theta) + v^(-theta) - 1
    return(max(val,0)^(-1/theta))
  },
  D_Dx = function(theta,u,v){
    val = - (v^theta*(u^(-theta)+v^(-theta)-1)^(-1/theta))/(u*(u^theta*(v^theta - 1) - v^theta))
    test = isTRUE((u^(-theta) + v^(-theta) > 1))
    if(test){
      return(val)
    }
    else{
      return(0)
    }
  },
  D_Dx_Dt = function(theta,u,v){
    NUM1 = - (v^(theta)  * (u^(-theta)+v^(-theta)-1)^(-1/theta) *
                (log(u^(-theta)+v^(-theta)-1)/theta^2) - ((u^(-theta)*(-log(u)) - v^(-theta)*log(v))/(theta*(u^(-theta)+v^(-theta)-1))))
    NUM2 = v^(theta)*log(v)*(u^(-theta)+v^(-theta)-1)^(-1/theta)
    NUM3 = v^(theta)*(u^(-theta)+v^(-theta)-1)^(-1/theta) * (u^(theta)*(v^(theta)-1)*log(u) + u^(theta)*v^(theta)*log(v) - v^(theta)*log(v))
    DENOM1 = u*(u^(theta) *(v^theta -1) - v^theta)
    DENOM2 = u*(u^(theta)*(v^(theta)-1))^2
    val = (NUM1/DENOM1) - (NUM2/DENOM1) + (NUM3/DENOM2)
    test = isTRUE((u^(-theta) + v^(-theta) > 1))
    if(test){
      return(val)
    }
    else{
      return(0)
    }
  }
)

FRANKCOPULA = list(
  COPULA = function(theta, u,v){
    val = (-1/theta) * log(1 + ((exp(-theta*u)-1)*(exp(-theta*v)-1))/(exp(-theta)-1))
    return(val)
  }
)

JOECOPULA = list(
  COPULA = function(theta,u,v){
    val = 1 - ((1-u)^theta + (1-v)^theta - (((1-u)^theta) * ((1-v)^theta)))^(1/theta)
    return(val)
  }
)

NELSEN_2_COPULA = list(
  COPULA = function(theta, u,v){
    Val =  1 - ((1-u)^(theta) +  (1-v)^(theta))^(1/theta)
    return(max(Val,0))
  }
)

NELSEN_15_GG_COPULA = list(
  COPULA = function(theta, u,v){
    Val = 1 - ((1-u^theta)^(1/theta) + (1-v^theta)^(1/theta))^(theta)
    return(max(Val,0)^(1/theta))
  }
)




Nelsen_2 = list(
  Main_Function = function(theta, u,v){
  Val =  1 - ((1-u)^(theta) +  (1-v)^(theta))^(1/theta)
  return(max(Val,0))
  },
  Generator = function(t,theta){
    Val = (1-t)^theta
    return(Val)
  },
  Deriv = function(t,theta){
    Val = -theta * (1-t)^(theta-1)
    return(Val)
  },
  Inv = function(t,theta){
    Val = 1 - t^(1/theta)
    return(Val)
  },
  Inv.Deriv = function(t,theta){
    Val = 1 - (-t/theta)^(1/(theta-1))
    return(Val)
  },
  Gen_Dat = function(size,parameter){
    u = runif(size,0,1)
    t = runif(size,0,1)
    w = map_dbl(.x = 1:size,
                .f = ~ max(Nelsen_2$Inv.Deriv((Nelsen_2$Deriv(u[.x],parameter)/t[.x]),parameter),0))
    v = map_dbl(.x = 1:size,
                .f = ~ max(Nelsen_2$Inv((Nelsen_2$Generator(w[.x],parameter) - Nelsen_2$Generator(u[.x],parameter)),parameter),0))
    DATA = data.frame(x=u,
                      y=v)
    return(DATA)
  }
)

Nelsen_15_GG = list(
  Main_Function = function(theta, u,v){
    Val = 1 - ((1-u^theta)^(1/theta) + (1-v^theta)^(1/theta))^(theta)
    return(max(Val,0)^(1/theta))
  },
  Generator = function(t,theta){
    Val = (1-t^theta)^(1/theta)
    return(Val)
  },
  Gen_Dat = function(size, parameter){
    u = runif(size,0,1)
    v = runif(size,0,1)
    R1 = Nelsen_15_GG$Generator((Nelsen_15_GG$Generator((v^(1/(1-parameter))),parameter) * u), parameter)
    R2 = Nelsen_15_GG$Generator((Nelsen_15_GG$Generator((v^(1/(1-parameter))),parameter) * (1-u)), parameter)
    RES = data.frame(x=R1,
                     y=R2)
    return(RES)
  }
)


#----------------------------------Testing Bounds----------------------------------#
#Works Over Whole Range
for(i in c(1,1.01,1.1, 2, 10,100)){
  print(12 *integral2(Vectorize(Nelsen_2),xmin = 0, xmax = 1, ymin = 0, ymax = 1, theta = i)$Q - 3)
}
#Works Over Negatives
for(i in c(0.01, 0.1,0.5,0.9,0.99,0.999)){
  print(12 * integral2(Vectorize(Nelsen_7),xmin = 0, xmax = 1, ymin = 0, ymax = 1, theta = i)$Q - 3)
}
#Currently Problematic
for(i in c(1, 2, 10,100)){
  print(12 *integral2(Vectorize(Nelsen_8),xmin = 0, xmax = 1, ymin = 0, ymax = 1, theta = i)$Q - 3)
}
#The Domain is a little weird -- seems to work up to beyond the domain of theta specified in book
for(i in c(0.0001,0.1,0.25,0.4,0.49, 0.5,0.6,0.6855)){
  print(12 *integral2(Vectorize(Nelsen_11),xmin = 0, xmax = 1, ymin = 0, ymax = 1, theta = i)$Q - 3)
}
#Works Over While Range -- Practical Bounds Associated with numerical integration stuff
for(i in c(1.01, 1.1,1.3,1.5,1.6,1.8, 2, 10,50)){
  print(12 *integral2(Vectorize(Nelsen_15),xmin = 0, xmax = 1, ymin = 0, ymax = 1, theta = i)$Q - 3)
}
#This One Doesn't Work So Well
for(i in c(2,2.01,2.2,4,8.10,100)){
  print(12 *integral2(Vectorize(Nelsen_18),xmin = 0, xmax = 1, ymin = 0, ymax = 1, theta = i)$Q - 3)
}
#Nelsen 21: 
for(i in c(1,1.01,1.2,1.5,3,5,10,11,15,16)){
  print(12 *integral2(Vectorize(Nelsen_21),xmin = 0, xmax = 1, ymin = 0, ymax = 1, theta = i)$Q - 3)
}

for(i in c(0.0001,0.1, 0.15)){
  print(12 *integral2(Vectorize(Nelsen_22),xmin = 0, xmax = 1, ymin = 0, ymax = 1, theta = i)$Q - 3)
}

#For the three good ones: 7,15,21 -- verifying bounds mechanically
##01. Nelsen_2 -- c(1.001, 80)
map_dbl(.x = c(1.001,1.5,2,80),.f = ~ 12 *integral2(Vectorize(Nelsen_2),xmin = 0, xmax = 1, ymin = 0, ymax = 1, theta = .x)$Q - 3)
##02. Nelsen_15 -- c(1.0005,40)
map_dbl(.x = c(1.0005,1.5,2,40),.f = ~ 12 *integral2(Vectorize(Nelsen_15),xmin = 0, xmax = 1, ymin = 0, ymax = 1, theta = .x)$Q - 3)
##03. Nelsen_21 = c(1.0005,16.6)
map_dbl(.x = c(1.0005,1.5,2,16.6),.f = ~ 12 *integral2(Vectorize(Nelsen_21),xmin = 0, xmax = 1, ymin = 0, ymax = 1, theta = .x)$Q - 3)





