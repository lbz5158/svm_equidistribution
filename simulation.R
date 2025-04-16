rm(list=ls())
cat("\014")

library(ggplot2)
library(dplyr)
library(e1071)
library(gridExtra)

#### Exercise 1: Asymptotics ####
n <- 50 # Sample size
S <-  round(500*log(n)) # Grid size [User choice]

# Generate data using Weyl sequences
#X <-  sqrt(3)*1:S-floor(sqrt(3)*1:S)
#Y <-  sqrt(2)*1:S-floor(sqrt(2)*1:S)

# Generate data using Monte Carlo sequences
set.seed(1071)
X <- runif(S)
Y <- runif(S)


#### True regrions ####
{
xc <- 0.25 # center x_c or h
yc <- 0.75 # y_c or k
a <- 0.125 # major axis length
b <- 0.125 # minor axis length
phi <- 0 # angle of major axis with x axis phi or tau

t <- seq(0, 2*pi, 0.01) 
x1 <- xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi)
y1 <- yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi)
region1_true <- data.frame(x=x1,y=y1)

xc2 <- 0.5 # center x_c or h
yc2 <- 0.5 # y_c or k
a2 <- 0.25*sqrt(2) # major axis length
b2 <- 0.125*sqrt(2) # minor axis length
phi2 <- pi/4 # angle of major axis with x axis phi or tau
x2 <- xc2 + a2*cos(t)*cos(phi2) - b2*sin(t)*sin(phi2)
y2 <- yc2 + a2*cos(t)*sin(phi2) + b2*sin(t)*cos(phi2)
region2_true <- data.frame(x=x2,y=y2)

region3_true <- data.frame(x = c(0.25, 0.375, 0.5, 0.375),
                      y = c(0.375, 0.25, 0.375, 0.5))
region4_true <- data.frame(x = c(0.5, 0.625, 0.75, 0.625),
                      y = c(0.625, 0.5, 0.625, 0.75))
}

# True data
data_truth <- data.frame(X=X, Y=Y, label = 0)
for(i in 1:S){
  x <- data_truth[i,1]
  y <- data_truth[i,2]
  
  # region 1
  if((x-0.25)^2+(y-0.75)^2<=0.125^2){data_truth[i,3] <- 1}
  
  # region 2
  if((x+y-1)^2/0.5^2+(y-x)^2/0.25^2<=1){data_truth[i,3] <- 1}
  
  # regions 3 & 4
  if(abs(x-0.375)+abs(y-0.375)<0.125 | abs(x-0.625)+abs(y-0.625)<0.125)
    {data_truth[i,3] <- 0}
}
data_truth$label <- factor(data_truth$label)


#### Approximate Regions ####
{
  xc <- 0.25+1/n # center x_c or h
  yc <- 0.75-1/n # y_c or k
  a <- 0.125+1/n # major axis length
  b <- 0.125-1/n # minor axis length
  phi <- 10*pi/n # angle of major axis with x axis phi or tau
  
  t <- seq(0, 2*pi, 0.01) 
  x1 <- xc + a*cos(t)*cos(0.5*pi-phi) - b*sin(t)*sin(0.5*pi-phi)
  y1 <- yc + a*cos(t)*sin(0.5*pi-phi) + b*sin(t)*cos(0.5*pi-phi)
  region1 <- data.frame(x=x1,y=y1)
  
  xc2 <- 0.5+1/n # center x_c or h
  yc2 <- 0.5-1/n # y_c or k
  a2 <- 0.25*sqrt(2)+1/n # major axis length
  b2 <- 0.125*sqrt(2)-1/n # minor axis length
  phi2 <- pi/4-10*pi/n # angle of major axis with x axis phi or tau
  x2 <- xc2 + a2*cos(t)*cos(0.5*pi-phi2) - b2*sin(t)*sin(0.5*pi-phi2)
  y2 <- yc2 + a2*cos(t)*sin(0.5*pi-phi2) + b2*sin(t)*cos(0.5*pi-phi2)
  region2 <- data.frame(x=x2,y=y2)

  shrink <- 1-1/sqrt(n)
  inflate <- 1+1/sqrt(n)
  region3 <- data.frame(x = c(0.25, 0.375, 0.5, 0.375)*shrink, 
                             y = c(0.375, 0.25, 0.375, 0.5)*shrink)
  region4 <- data.frame(x = c(0.5, 0.625, 0.75, 0.625)*inflate, 
                             y = c(0.625, 0.5, 0.625, 0.75)*inflate)
}



#### Helper Functions ####
make_data <- function(df, case = 1){
  for(i in 1:S){
    x <- df[i,1]
    y <- df[i,2]
    
    # region 1
    if(((x-xc)*cos(phi)-(y-yc)*sin(phi))^2/(b)^2+
       ((x-xc)*sin(phi)+(y-yc)*cos(phi))^2/(a)^2<=1){df[i,3] <- 1}
    
    # region 2
    if(case >= 2){
      if(((x-xc2)*cos(phi2)-(y-yc2)*sin(phi2))^2/b2^2+
         ((x-xc2)*sin(phi2)+(y-yc2)*cos(phi2))^2/a2^2<=1){df[i,3] <- 1}
    }
    
    # regions 3 & 4
    if(case == 3){
      if(abs(x-0.375*shrink)+abs(y-0.375*shrink)<0.125*shrink | 
         abs(x-0.625*inflate)+abs(y-0.625*inflate)<0.125*inflate){df[i,3] <- 0}
    }
  }
  df$label <- factor(df$label)
  return(df)
}

make_plot <- function(case = 1){
  return(ggplot() +
           geom_point(data = predicted, 
                      aes(X, Y, color = factor(label)),
                      size = 0.5,
                      show.legend = FALSE) +
           scale_color_manual(values=c("blue","red")) + 
           geom_polygon(data = region1, aes(x,y), colour="black", fill = NA) +
           xlab("") + ylab("") +
           {if(case >= 2)geom_polygon(data = region2, aes(x,y), colour="black", fill = NA)} +
           {if(case == 3)geom_polygon(data = region3, aes(x = x, y = y), colour='black', fill=NA)} +
           {if(case == 3)geom_polygon(data = region4, aes(x = x, y = y), colour='black', fill=NA)} +
           geom_point(data = errors, 
                      aes(X, Y, color = factor(label)),
                      shape = 4,
                      stroke = 1.5,
                      show.legend = F)
  )
}



#### Generate Plots ####
{
data <- data.frame(X=X,Y=Y,label=0)

data <- make_data(data,case=3)

ggplot() +
  geom_point(data = data,
             aes(X, Y, color=factor(label)),
             size = 1,
             show.legend = FALSE) +
  scale_color_manual(values=c("blue","green")) +
  geom_polygon(data = region1, aes(x,y), colour="black", fill = NA) +
  geom_polygon(data = region2, aes(x,y), colour="black", fill = NA) +
  geom_polygon(data = region3, aes(x,y), colour="black", fill = NA) +
  geom_polygon(data = region4, aes(x,y), colour="black", fill = NA) +
  xlab("") + ylab("") +
  geom_polygon(data = region1_true, aes(x,y), colour="red", fill = NA) +
  geom_polygon(data = region2_true, aes(x,y), colour="red", fill = NA) +
  geom_polygon(data = region3_true, aes(x,y), colour="red", fill = NA) +
  geom_polygon(data = region4_true, aes(x,y), colour="red", fill = NA) 
}




#### Case 1: 1 region ####
# {
#   data <- data.frame(X=X,Y=Y,label=0)
#   
#   # Add class 1
#   data <- make_data(data,1)
#   
#   # Plot the grid
#   p2 <- 
#     ggplot() +
#     geom_point(data = data, 
#                aes(X, Y, color=factor(label)),
#                size = 1,
#                show.legend = FALSE) +
#     scale_color_manual(values=c("blue","green")) +
#     geom_polygon(data = region1, aes(x,y), colour="black", fill = NA) +
#     xlab("") + ylab("") +
#     geom_polygon(data = region1_true, aes(x,y), colour="red", fill = NA)
#   
#  grid.arrange(p1, p2, nrow = 1)
#
#
#   # gamma =0.45
#   svmfit = svm(label ~ ., data = data, kernel = "radial", gamma = 0.45)
#   print(svmfit)
# 
#   predicted <- data.frame(X=X,Y=Y,label=predict(svmfit,data[,1:2]))
#   errors <- predicted[which(predicted$label != data$label),]
#   p1 <- make_plot(1)
# 
#   # gamma = 45
#   svmfit = svm(label ~ ., data = data, kernel = "radial", gamma = 45)
#   print(svmfit)
# 
#   predicted <- data.frame(X=X,Y=Y,label=predict(svmfit,data[,1:2]))
#   errors <- predicted[which(predicted$label != data$label),]
#   p2 <- make_plot(1)
# 
#   grid.arrange(p1, p2, nrow = 1)
# }






#### Case 3: 4 regions ####
{
  data <- data.frame(X=X,Y=Y,label=0)
  
  # Add class 1
  data <- make_data(data,3)
  
  # gamma = 43
  svmfit = svm(label ~ ., data = data, kernel = "radial", gamma = 43)
  print(svmfit)
  
  predicted <- data.frame(X=X,Y=Y,label=predict(svmfit,data[,1:2]))
  errors <- predicted[which(predicted$label != data$label),]
  p1 <- make_plot(3)
  
  errors_truth <- predicted[which(predicted$label != data_truth$label),]
  p1_truth <- ggplot() +
    geom_point(data = predicted, 
               aes(X, Y, color = factor(label)),
               size = 0.5,
               show.legend = FALSE) +
    scale_color_manual(values=c("blue","red")) + 
    xlab("") + ylab("") + 
    geom_polygon(data = region1_true, aes(x,y), colour="green", fill = NA) +
    geom_polygon(data = region2_true, aes(x,y), colour="green", fill = NA) +
    geom_polygon(data = region3_true, aes(x,y), colour="green", fill = NA) +
    geom_polygon(data = region4_true, aes(x,y), colour="green", fill = NA) +
    geom_point(data = errors_truth, 
               aes(X, Y, color = factor(label)),
               shape = 4,
               stroke = 1.5,
               show.legend = F) 
  
  # gamma = 50000
  svmfit = svm(label ~ ., data = data, kernel = "radial", gamma = 50000)
  print(svmfit)
  
  predicted <- data.frame(X=X,Y=Y,label=predict(svmfit,data[,1:2]))
  errors <- predicted[which(predicted$label != data$label),]
  p2 <- make_plot(3)
  
  errors_truth <- predicted[which(predicted$label != data_truth$label),]
  p2_truth <- ggplot() +
    geom_point(data = predicted, 
               aes(X, Y, color = factor(label)),
               size = 0.5,
               show.legend = FALSE) +
    scale_color_manual(values=c("blue","red")) + 
    xlab("") + ylab("") + 
    geom_polygon(data = region1_true, aes(x,y), colour="green", fill = NA) +
    geom_polygon(data = region2_true, aes(x,y), colour="green", fill = NA) +
    geom_polygon(data = region3_true, aes(x,y), colour="green", fill = NA) +
    geom_polygon(data = region4_true, aes(x,y), colour="green", fill = NA) +
    geom_point(data = errors_truth, 
               aes(X, Y, color = factor(label)),
               shape = 4,
               stroke = 1.5,
               show.legend = F) 
  
  
  grid.arrange(p1, p2, nrow = 1) # library(gridExtra)
  
  # Against truth
  grid.arrange(p1_truth, p2_truth, nrow=1)
  
}










#### Exercise 2: OLS ####
rm(list=ls())
cat("\014")

library(e1071)
library(dplyr)
library(ggplot2)

# Results: n=500&L=100: 0.91, 0.88

#### Idea 1 ####

#### OLS Estimates ####
n <- 500
d <- 2
S <- round(6^(log(n)))
L <- 500

{
  set.seed(2)
  # True beta
  beta <- c(runif(d,-1,1)*5)#,runif(2,-1,1)*5)
}

captures <- numeric(L)############
svm_captures <- numeric(L)
set.seed(2)
for(i in 1:L){
  cat(i,": ") ##################
  # Generate simulation data
  x <- cbind(#runif(n)*5,
    #runif(n)*10,
    #runif(n)*10,
    runif(n,2.2,3.8))
  designX <- cbind(rep(1,n), x)
  
  sigma <- 3
  epsilon <- rnorm(n,0,sigma)
  
  y <- designX %*% beta + epsilon  
  
  # Estimation
  beta_hat <- solve(t(designX) %*% designX) %*% t(designX) %*% y
  
  V_hat <- sigma^2*solve(t(designX)%*%designX)
  
  captures[i] <- t(beta-beta_hat) %*% solve(V_hat) %*% (beta-beta_hat) <= qchisq(0.95,d)
  
  
  # Generate grid
  ols_data <- data.frame(runif(S,-10,10),
                         runif(S,-10,10),
                         #runif(S,-1,1)*5,
                         #runif(S,-1,1)*10,
                         #runif(S,-1,1)*10,
                         label = 0)
  
  # Make labels
  for(j in 1:S){
    b <- t(ols_data[j,1:2])
    if(t(b-beta_hat) %*% solve(V_hat) %*% (b-beta_hat) <= qchisq(0.95,d))
    {ols_data[j,3] <- 1}
  }
  
  ols_data$label <- factor(ols_data$label)
  
  names(ols_data)[1:2] <- c("X","Y")
  # ggplot() +
  #   geom_point(data = ols_data,
  #              aes(X, Y, color = factor(label)),
  #              size = 0.5,
  #              show.legend = FALSE) +
  #   scale_color_manual(values=c("blue","red","green")) +
  #   xlab("") + ylab("") +
  #   geom_point(aes(x=beta[1],y=beta[2],color='green'),
  #              show.legend = F)
  
  
  
  # Fit SVM
  gamma <- 1000
  while(TRUE){
    ols_svm = svm(label ~ ., data = ols_data, kernel = "radial", gamma = gamma)
    wrong_rate <- sum(predict(ols_svm,ols_data[,1:2]) == ols_data$label)/S
    #if(wrong_rate >= 1-100/n)
    {break} #############
    gamma <- gamma*2
    cat(wrong_rate,"; ",gamma,'\n')
  }
  
  svm_captures[i] <- as.numeric(as.character(predict(ols_svm, t(beta))))
  cat(svm_captures[i],'\n')
}

mean(captures)
mean(svm_captures)








#### Idea 2 ####
rm(list=ls())
cat("\014")

n <- 500
d <- 5
S <- round(6^(log(n)))

{
  set.seed(2)
  # True beta
  beta <- c(runif(d,-1,1)*5)
}

set.seed(10)
x <- cbind(runif(n)*5,
           runif(n)*10,
           runif(n)*10,
           runif(n,2.2,3.8))
designX <- cbind(rep(1,n), x)

sigma <- 3
epsilon <- rnorm(n,0,sigma)

y <- designX %*% beta + epsilon  

# Estimation
beta_hat <- solve(t(designX) %*% designX) %*% t(designX) %*% y

V_hat <- sigma^2*solve(t(designX)%*%designX)
ols_grid <- data.frame(runif(S,-4,-2.5),
                       runif(S,1.7,2.3),
                       runif(S,0.5,1),
                       runif(S,-4,-2.8),
                       runif(S,4,5),
                       label = 0)

# Make labels
for(j in 1:S){
  b <- t(ols_grid[j,1:d])
  if(t(b-beta_hat) %*% solve(V_hat) %*% (b-beta_hat) <= qchisq(0.95,d))
  {ols_grid[j,d+1] <- 1}
}
table(ols_grid$label)
ols_grid$label <- factor(ols_grid$label)

names(ols_grid)[1:d] <- c("X1","X2","X3","X4",'X5')


# ggplot() +
#   geom_point(data = ols_data,
#              aes(X, Y, color = factor(label)),
#              size = 0.5,
#              show.legend = FALSE) +
#   scale_color_manual(values=c("blue","red","green")) +
#   xlab("") + ylab("") +
#   geom_point(aes(x=beta[1],y=beta[2],color='green'),
#              show.legend = F)



# Iterations
system.time({
  means <- c()
  correct <- c()
  set.seed(10)
  for(iter in 1:2000){
    train <- sample(1:S, round(0.05*S), repl=F)
    cat(iter,'\n')############
    ols_svm <- svm(label ~ ., data = ols_grid[train,], kernel = "radial", gamma = 1.5)
    b <- as.numeric(as.character(predict(ols_svm,t(beta))))
    c <- mean(as.numeric(as.character(predict(ols_svm, ols_grid[-train,]))) == ols_grid[-train,"label"])
    means <- c(means,c)
    correct <- c(correct,b)
  }
  #plot(seq(1,10,1),means)
  cat(mean(means), "&", mean(correct))
})


# 0.80*S & 2000 iters: 0.99148 & 1; 23630 sec
# 0.42*S & 2000 iters: 0.9899284 & 0.9465; 8146 sec
# 0.05*S & 2000 iters: 0.9877293 & 0.088; 1140 sec

