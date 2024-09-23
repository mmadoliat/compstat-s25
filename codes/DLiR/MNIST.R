library(keras)

### Step 1 : Data preprocessing ###
# Load MNIST (Modified National Institute of Standards and Technology) images datasets 
c(c(x_train, y_train), c(x_test, y_test)) %<-% dataset_mnist()

# Flatten images and transform RGB values into [0,1] range 
x_train <- array_reshape(x_train, c(nrow(x_train), 784))
x_test <- array_reshape(x_test, c(nrow(x_test), 784))
x_train <- x_train / 255
x_test <- x_test / 255

# Convert class vectors to binary class matrices
y_train <- to_categorical(y_train, 10)
y_test <- to_categorical(y_test, 10)

### Step 2 : Model definition ###
model <- keras_model_sequential() 
model %>% 
  layer_dense(units = 256, activation = 'relu', input_shape = c(784)) %>% 
  layer_dropout(rate = 0.4) %>% 
  layer_dense(units = 128, activation = 'relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 10, activation = 'softmax')

### Step 3 : Compile Model ###
model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)

### Step 4 : Model Training ###
history <- model %>% fit(
  x_train, y_train, 
  batch_size = 128, 
  epoch = 10,
  validation_split = 0.2
)

plot(history)

model %>% predict_classes(x_test[1:100,])

round(model %>% predict(x_test[1:9,]),5)

image(t(array_reshape(x_test,c(nrow(x_test), 28,28))[12,28:1,]), axes=FALSE, col = grey(seq(0, 1, length = 256)))
