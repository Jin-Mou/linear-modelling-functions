### the work is completed by Jin Mou, s2103110

### code overview: 
## We build our own functions for fitting linear models using QR decomposition. 
## "linmod" function does the fit and returns the corresponding objects.
## "print.linmod" function prints the summmary of fit, including formula, parameter estimates and their standard deviations.
## "plot.linmod" function plots a default residual plot.
## "predict.linmod" function prints a vector of predictions



### (1) define linmod function
linmod <- function(formula,dat) {
  ## function brief introduction
  # it takes QR decomposition to estimate the specified linear model
  # formula: a linear model
  # dat: a data frame containing the corresponding data

  ## process data
  dat <- data.frame(dat) # convert dat into regular data frame
  dat <- model.frame(formula,dat) # get the data that will be used in formula

  ## convert character variables to factor variables and store information of factor variables
  flev <- list() # create a list that contains information of factor variables
  # a loop goes through every variables
  for (d in 1:ncol(dat)) {
    if (class(dat[,d]) == "character") {
      dat[,d] <- factor(dat[,d]) # Convert character variables into factor variables
    }
    if (class(dat[,d]) == "factor") {
      factor_name <- colnames(dat)[d] # get the name of a factor variable
      flev <- c(flev, list(levels(dat[,d]))) # store levels information of factor variables
      names(flev)[length(flev)] <- factor_name # change the name
    }
  }
  
  ## QR decomposition and store response variable y
  X <- model.matrix(formula, dat) # model variable matrix
  qrx <- qr(X) # QR decomposition of matrix X
  y <- dat[[1]] # get the values of response variable
  beta <- backsolve(qr.R(qrx),qr.qty(qrx, y)[1:ncol(X)]) ## R^{-1} Q^T y
  
  ## set the format of beta
  beta_names <- c("(Intercept)", colnames(X)[2:dim(X)[2]]) # get the names of coefficient
  
  names(beta) <- beta_names # set names of beta
  
  ## calculate the vector of expected values of the response variable mu
  mu <- X %*% matrix(beta)
  
  ## set the estimated covariance matrix V, yname
  R2 <- solve(qr.R(qrx)) # get the inverse matrix of R
  yname <- all.vars(formula)[1] # store the name of response variable
  sigma_square <- sum((y - mu)**2) / (nrow(dat)-ncol(dat)+1) # calculate sigma square
  sigma <- sqrt(sigma_square) # store sigma
  V <- R2 %*% t(R2) * sigma_square # get the estimated covariance matrix of the least squares estimators
  
  ## returns an object of class "linmod"
  linmod <- list(beta=beta, V=V, mu=drop(mu), y=y, yname=yname, formula=formula, flev=flev, sigma = sigma)
  class(linmod) <- "linmod"
  linmod
}


### (2) define function print.linmod
print.linmod <- function(x, ...) {
  ## function brief introduction
  # it prints the model formula and reports the parameter estimates and their standard deviations
  # x: an object of class "linmod"

  ## show the model formula
  print(x$formula) 
  cat("\n")
  
  ## calculate and show the parameter estimates and their standard deviations
  estimate <- t(t(x$beta)) # transpose beta
  colnames(estimate) <- "Estimate" # assign column name to estimate
  se <- sqrt(matrix(diag(x$V))) # calculate standard deviations
  colnames(se) <- "s.e." # assign column name to estimate
  print(cbind(estimate, se))
} 


### (3) define function plot.linmod
plot.linmod <- function(x, ...) {
  ## function brief introduction
  # it plots the model residuals against the model fitted values
  # x: an object of class "linmod"
  
  ## calculate residuals and plot
  residuals <- x$y - x$mu # get the residuals
  plot(x$mu, residuals, xlab="fitted values", ylab="residuals", col="blue") # plot residuals vs. fitted values
  abline(h=0, lty=3) # plot lines that residuals = 0
}


### (4) define function predict.linmod
predict.linmod <- function(x, newdata, ...) {
  ## function brief introduction
  # it return a vector of predictions
  # x: an object of class "linmod"
  # newdata: a data frame containing values of the predictor variables
  
  ## pre-process newdata
  newdata <- data.frame(newdata) # convert newdata into regular data frame
  # if response variable is included in the newdata, then we don't need to modify new data
  # otherwise, we construct dummy response variables into newdata
  if (is.null(newdata[[x$yname]])) { # if response variable isn't included in the newdata
    dummy_y <- matrix(1, dim(newdata)[1],1) # build dummy response variables dummy_y
    colnames(dummy_y) <- x$yname
    newdata <- cbind(dummy_y, newdata) # add dummy_y into newdata
  }
  
  ## if original model data has factor variables
  if (! is.null(names(x$flev))) {
    # for every variable in in newdata,
    for (d in 2:dim(newdata)[2]) {
      # if it is a factor variable in original data,
      # then it is also a factor variable in newdata
      if (names(newdata)[d] %in% names(x$flev)) {
        newdata[,d] <- factor(newdata[,d], levels=x$flev[[names(newdata)[d]]])
      }
    }
  }
  
  ## do predictions
  X2 <- model.matrix(x$formula, newdata)
  prediction_y <- X2 %*% matrix(x$beta) # prediction using y = x*beta
  colnames(prediction_y) <- "predictions"
  print(prediction_y)
}





# ### test linmod using lm()
# mm <- linmod(dist ~ speed + I(speed^2), cars)
# mm2 <- lm(dist ~ speed + I(speed^2), cars)
# mm$beta
# mm2$coefficients
# mm$V
# head(mm$mu - mm2$fitted.values)
# head(mm$y)
# head(cars[2])
# mm$formula
# mm$flev
# mm2$xlevels
# mm$sigma
# sqrt(sum(mm2$residuals**2)/(dim(cars)[1]-dim(cars)[2]))
# 
# ### test print.linmod
# mm <- linmod(dist ~ speed + I(speed^2), cars)
# print(mm)
# 
# ### test plot.linmod
# mm <- linmod(dist ~ speed + I(speed^2), cars)
# mm2 <- lm(dist ~ speed + I(speed^2), cars)
# plot(mm)
# plot(mm2)
# 
# ### test plot.predict
# dd <- linmod(weight ~ group, PlantGrowth)
# dd2 <- lm(weight ~ group, PlantGrowth)
# aa <- matrix(PlantGrowth[,2][15:30])
# colnames(aa) <- "group"
# aa <- data.frame(aa)
# levels(aa[,1])
# levels(PlantGrowth[,2])
# predict(dd, aa)
# predict(dd2, aa)
