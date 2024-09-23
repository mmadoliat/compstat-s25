mu.logo <- data.frame(read.csv("mu-logo.csv"))
lm.fit <- lm(Y~X1+X2+X3+X4+X5, data= mu.logo)
plot(lm.fit$fitted.values, lm.fit$residuals)

library("neuralnet");
net1 <- neuralnet(Y~X1+X2+X3+X4+X5, data=mu.logo, hidden=0, act.fct=function(x) {x})
plot(net1)

net2 <- neuralnet(Y~X1+X2+X3+X4+X5, data=mu.logo,
                  hidden=10, act.fct=function(x) {x})
plot(net2)


net3 <- neuralnet(Y~X1+X2+X3+X4+X5, data=mu.logo, 
                  hidden=c(10,10), 
                  act.fct=function(x) {x})
plot(net3)


errors <- c(net1$result.matrix[1],
            net2$result.matrix[1],
            net3$result.matrix[1])

par(mfrow=c(2,2))
plot(lm.fit$fitted.values,lm.fit$residuals)
plot(net1$net.result[[1]], net1$data$Y-net1$net.result[[1]])
plot(net2$net.result[[1]], net1$data$Y-net2$net.result[[1]])
plot(net3$net.result[[1]], net1$data$Y-net3$net.result[[1]])

lm.fit$coefficients

net1$weights[[1]][[1]]

net2$weights[[1]][[1]]%*%net2$weights[[1]][[2]][-1,]+
  c(net2$weights[[1]][[2]][1,],rep(0,5))

net3$weights[[1]][[1]]%*%net3$weights[[1]][[2]][-1,]%*%
  net3$weights[[1]][[3]][-1,]+
  c(net3$weights[[1]][[2]][1,]%*%
      net3$weights[[1]][[3]][-1,]+
      net3$weights[[1]][[3]][1,],rep(0,5))

mu.logo <- data.frame(cbind(mu.logo, matrix(rnorm(prod(dim(net1$covariate))),nc=5)) )

net4 <- neuralnet(Y~X1+X2+X3+X4+X5+X1.1+X2.1+X3.1+X4.1+X5.1, data=mu.logo,
                  hidden=0,linear.output=TRUE, act.fct=function(x) {x})

plot(lm.fit$fitted.values, lm.fit$residuals)

plot(net4$net.result[[1]], net4$data$Y-net4$net.result[[1]])

plot(net4)


############## Keras ######################################################
library(keras)
mu.logo <- data.frame(read.csv("mu-logo.csv"))
y_train <- mu.logo[,1]
x_train <- as.matrix(mu.logo[,2:6])
model <- keras_model_sequential() 
model %>% 
  layer_dense(units = 1, input_shape = 5) 
  # layer_dense(units = 256, input_shape = 5) %>% 
  # layer_dropout(rate = 0.4) %>% 
  # layer_dense(units = 128) %>%
  # layer_dropout(rate = 0.3) %>%
  # layer_dense(units = 1)

### Step 3 : Compile Model ###
model %>% compile(
  loss = 'mse',
  optimizer = optimizer_rmsprop(),
  metrics = c('mean_absolute_error')
)

### Step 4 : Model Training ###
history <- model %>% fit(
  x_train, y_train, 
  batch_size = 256, 
  epoch = 100,
  validation_split = 0.2
)

plot(history)

y_hat <- model %>% predict(x_train)
plot(y_hat, y_train-y_hat)

