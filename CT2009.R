#### Comparison with Ciliberto and Tamer (Econometrica 2009) ####
library(optimization)
library(e1071)

#########################################
load('Home/work/svm/myEnvironment.RData')
#########################################

#### Data generation ####
set.seed(100)
N <- 1000
B <- 100 # number of bootstrap subsamples
dX <- 2

eps <- rnorm(N,1)

means <- runif(dX,-1,1)*10
theta0 <- runif(dX,-5,10)
X <- matrix(0,N,dX)
for(i in 1:dX){
  X[,i] <- rnorm(N,means[i],1)
}

Y <- 3 + X%*%theta0 + eps # real DGP
designX <- cbind(rep(1,N),X) # design matrix

Q <- function(theta,x=designX,y=Y){ 
  return(norm(y-x%*%theta,type='2'))
}  
Q(c(3,theta0))
theta_hat <- optim_nm(Q,k=3)




#### CHT CS on a Large Grid ####
{
  # Generate grid (size 101*101*101) #
  theta1_grid <- seq(-10,10,length.out=101)
  theta2_grid <- seq(-10,10,length.out=101)
  theta3_grid <- seq(-10,10,length.out=101)
  Grid <- expand.grid(theta1_grid,theta2_grid,theta3_grid)
  rm(theta1_grid,theta2_grid,theta3_grid)
  
  
  # Evaluate nQ on the grid
  for(i in 1:nrow(Grid)){
    Grid[i,4] <- Q(t(as.matrix(Grid[i,1:3])))
    if(i%%10000==0){cat(i,'\n')}
  }
  
  c0 <- 900
  CS_c0 <- Grid[Grid[,4]<(c0+min(Grid[,4])),1:3] 
  
  
  # Bootstrap Subsampling #
  Q_bs <- function(theta){  # helper function: bootstrap objective
    x <- designX[bs_index,] 
    y <- Y[bs_index]
    return(norm(y-x%*%theta,type='2'))
  }  
  
  {
    Crits <- numeric(B)
    for(b in 1:B){
      # Draw bootstrap samples
      bs_index <- sample(N,N/4,repl=T)
      minQ <- optim_nm(Q_bs,k=3,start=as.vector(theta_hat$par))$function_value
      
      # Evaluate the bootstrap objective at the set of parameters given by c0
      bs_max <- -999
      for(j in 1:nrow(CS_c0)){
        v <- Q(t(as.numeric(CS_c0[j,])),designX[bs_index],Y[bs_index])-minQ
        if(v > bs_max){
          bs_max <- v
        }
      }
      
      Crits[b] <- bs_max # store the maximum
      cat(b,'\n')
    }
    c1 <- quantile(Crits,0.05) 
    CS_c1 <- Grid[Grid[,4]<(c1+min(Grid[,4])),1:3] # Q here is nQ(theta) in CT2009
  }
  
  {
    Crits <- numeric(B)
    for(b in 1:B){
      # Draw bootstrap samples
      bs_index <- sample(N,N/4,repl=T)
      minQ <- optim_nm(Q_bs,k=3,start=as.vector(theta_hat$par))$function_value
      
      # Evaluate the bootstrap objective at the set of parameters given by c0
      bs_max <- -999
      for(j in 1:nrow(CS_c1)){
        v <- Q(t(as.numeric(CS_c1[j,])),designX[bs_index],Y[bs_index])-minQ
        if(v > bs_max){
          bs_max <- v
        }
      }
      Crits[b] <- bs_max # store the maximum
      cat(b,'\n')
    }
    c2 <- quantile(Crits,0.05) 
    CS_c2 <- Grid[Grid[,4]<(c2+min(Grid[,4])),1:3] 
  }
}

true_labels <- Grid[,4]<(c2+min(Grid[,4]))




#### CHT Prediction vs Truth ####
{
  # Generate grid (size 51*51*51) #
  theta1_grid <- seq(-10,10,length.out=51)
  theta2_grid <- seq(-10,10,length.out=51)
  theta3_grid <- seq(-10,10,length.out=51)
  Grid51 <- expand.grid(theta1_grid,theta2_grid,theta3_grid)
  rm(theta1_grid,theta2_grid,theta3_grid)
  
  
  # Evaluate nQ on the grid
  for(i in 1:nrow(Grid51)){
    Grid51[i,4] <- Q(t(as.matrix(Grid51[i,1:3])))
    if(i%%10000==0){cat(i,'\n')}
  }
  
  c0_new <- 900
  CS_c0_new <- Grid51[Grid51[,4]<(c0_new+min(Grid51[,4])),1:3] 
  
  {
    Crits <- numeric(B)
    for(b in 1:B){
      # Draw bootstrap samples
      bs_index <- sample(N,N/4,repl=T)
      minQ <- optim_nm(Q_bs,k=3,start=as.vector(theta_hat$par))$function_value
      
      # Evaluate the bootstrap objective at the set of parameters given by c0
      bs_max <- -999
      for(j in 1:nrow(CS_c0_new)){
        v <- Q(t(as.numeric(CS_c0_new[j,])),designX[bs_index],Y[bs_index])-minQ
        if(v > bs_max){
          bs_max <- v
        }
      }
      
      Crits[b] <- bs_max # store the maximum
      cat(b,'\n')
    }
    c1_new <- quantile(Crits,0.05) 
    CS_c1_new <- Grid51[Grid51[,4]<(c1_new+min(Grid51[,4])),1:3] 
  }
  
  {
    Crits <- numeric(B)
    for(b in 1:B){
      # Draw bootstrap samples
      bs_index <- sample(N,N/4,repl=T)
      minQ <- optim_nm(Q_bs,k=3,start=as.vector(theta_hat$par))$function_value
      
      # Evaluate the bootstrap objective at the set of parameters given by c0
      bs_max <- -999
      for(j in 1:nrow(CS_c1_new)){
        v <- Q(t(as.numeric(CS_c1_new[j,])),designX[bs_index],Y[bs_index])-minQ
        if(v > bs_max){
          bs_max <- v
        }
      }
      Crits[b] <- bs_max # store the maximum
      cat(b,'\n')
    }
    c2_new <- quantile(Crits,0.05) 
    CS_c2_new <- Grid51[Grid51[,4]<(c2_new+min(Grid51[,4])),1:3] 
  }
}

cbind(c(3,theta0),t(apply(CS_c2_new[,1:3],2,range))) # projected confidence intervals

# % labeled CS (51^3)
CS_Grid51 <- Grid51[Grid51[,2]>min(CS_c2_new[,2]),] # accepted points by rectangle
nrow(CS_Grid51)/nrow(Grid51)

# % correct (51^3)
(nrow(CS_c2_new)+(nrow(Grid51)-nrow(CS_Grid)))/nrow(Grid51)
# (correctly accepted + corrected rejected) / total grid points

# % labeled CS (101^3)
cht_pred_Grid <- numeric(nrow(Grid))
cht_pred_Grid[Grid[,2]>min(Grid51[,2])] <- 1
table(cht_pred_Grid)/length(cht_pred_Grid)

# % correct (101^3)
sum(cht_pred_Grid == true_labels)/length(true_labels)



#### SVM Prediction vs Truth ####
# Prepare training grid
train <- Grid51
train[train[,4]<(c2_new+min(train[,4])),5] <- 1
train[train[,4]>=(c2_new+min(train[,4])),5] <- 0
names(train)[5] <- 'label'
train <- train[,-4]
train$label <- as.factor(train$label)


## 1. SVM with tuning -------
svmfit = svm(label ~ ., data = train, kernel = "radial", gamma = 500)
print(svmfit)

predicted <- predict(svmfit,train[,1:3])

# % labeled CS (51^3)
table(predicted)/length(predicted)

# % correct (51^3) (this is 100%)
sum(predicted2==train$label)/length(predicted2)

# Predict on the large grid
svm_pred_Grid <- numeric(nrow(Grid)) 
for(i in 1:nrow(Grid)){
  svm_pred_Grid[i]<- predict(svmfit,Grid[i,1:3])
  if(i%%10000==0){cat(i,'\n')}
}

# % labeled CS (101^3)
table(svm_pred_Grid)/length(svm_pred_Grid)

# % correct by tuned SVM (101^3)
sum(svm_pred_Grid == (true_labels+1))/length(true_labels) 



## 2. SVM with no tuning -------
svmfit2 = svm(label ~ ., data = train, kernel = "radial", gamma = 1)
print(svmfit2)

predicted2 <- predict(svmfit2,train[,1:3])

# % labeled CS (51^3)
table(predicted2)/length(predicted2)

# % correct (51^3)
sum(predicted2==train$label)/length(predicted2)

# Predict on the large grid
svm_pred_Grid2 <- numeric(nrow(Grid)) 
t <- Sys.time()
for(i in 1:nrow(Grid)){
  svm_pred_Grid2[i]<- predict(svmfit2,Grid[i,1:3])
  if(i%%10000==0){cat(i,Sys.time()-t,'\n')}
}

# % labeled CS (101^3)
table(svm_pred_Grid2)/length(svm_pred_Grid2)

# % correct by untuned SVM (101^3)
sum(svm_pred_Grid2 == (true_labels+1))/length(true_labels) 



#####################################################
save.image(file='/storage/home/lbz5158/work/svm/myEnvironment.RData')
#####################################################