diamonds <- ggplot2::diamonds
summary(diamonds$carat)
summary(diamonds$cut)

library(sloop)
otype(1:10)
otype(mtcars)
mle_obj <- stats4::mle(function(x = 1) (x - 2) ^ 2)
otype(mle_obj)

is.object(1:10)
sloop::otype(1:10)
is.object(mtcars)
sloop::otype(mtcars)

attr(1:10, "class")
attr(mtcars, "class")

typeof(1:10)
typeof(mtcars)
typeof(Matrix::Matrix(0,nr=2))
