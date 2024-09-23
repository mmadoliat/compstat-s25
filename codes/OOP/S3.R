library(sloop)

f <- factor(c("a", "b", "c"))
typeof(f)
attributes(f)

unclass(f)

print(f)
print(unclass(f))

ftype(t.test)
ftype(t.data.frame)

weighted.mean.Date
s3_get_method(weighted.mean.Date)
stats:::weighted.mean.Date

# Create and assign class in one step
x <- structure(list(), class = "my_class")
# Create, then set class
x <- list()
class(x) <- "my_class"

# Create a linear model
mod <- lm(log(mpg) ~ log(disp), data = mtcars)
class(mod)
print(mod)
# Turn it into a date (?!)
class(mod) <- "Date"
# Unsurprisingly this doesn't work very well
print(mod)

new_Date <- function(x = double()) {
  stopifnot(is.double(x))
  structure(x, class = "Date")
}
new_Date(c(-1, 0, 1))

x <- Sys.Date()
s3_dispatch(print(x))

x <- matrix(1:10, nrow = 2)
s3_dispatch(mean(x))
s3_dispatch(sum(Sys.time()))

s3_methods_generic("print")
s3_methods_class("ordered")

class(ordered("x"))
class(Sys.time())

s3_dispatch(print(ordered("x")))
s3_dispatch(print(Sys.time()))

new_secret <- function(x = double()) {
  stopifnot(is.double(x))
  structure(x, class = "secret")
}
print.secret <- function(x, ...) {
  print(strrep("x", nchar(x)))
  invisible(x)
}
x <- new_secret(c(15, 1, 456))
x

s3_dispatch(x[1])
x[1]

`[.secret` <- function(x, i) {
  print("error")
  new_secret(x[i])
}
x[1]

`[.secret` <- function(x, i) {
  x <- unclass(x)
  new_secret(x[i])
}
x[1]

`[.secret` <- function(x, i) {
  new_secret(NextMethod())
}
x[1]
s3_dispatch(x[1])

new_secret <- function(x, ..., class = character()) {
  stopifnot(is.double(x))
  structure(
    x,
    ...,
    class = c(class, "secret")
  )
}

new_supersecret <- function(x) {
  new_secret(x, class = "supersecret")
}
print.supersecret <- function(x, ...) {
  print(rep("xxxxx", length(x)))
  invisible(x)
}
x2 <- new_supersecret(c(15, 1, 456))
x2

`[.supersecret` <- function(x, ...) {
  new_supersecret(NextMethod())
}
x2[1:3]

########### Toy Example #####################
new_human <- function(name=NA, age=0, race="Black", gender="Female", ..., class = character()) {
  structure(
    list(name=name, age=age, race=race, gender=gender, ...),
    class = c(class, "human")
  )
}

print.human <- function(x) {
  print(paste0(x$name," is a ",x$age," years old ",x$race," ",x$gender,"."))
  invisible(x)
}

a <- new_human("Alice",18)
print(a)

new_human_Univ <- function(name=NA, age=0, race="Black", gender="Female", type="Student") {
  type <- match.arg(type, c("Student", "Staff", "Faculty"))
  new_human(name, age, race, gender, type=type, class="human_Univ")
}

print.human_Univ <- function(x) {
  print(x$type)
  NextMethod()
}

b <- new_human_Univ("Alex",25,"Asian","Male","Stu")
print(b)

s3_dispatch(print(a))
s3_dispatch(print(b))

plot.human <- function(x) {
  print("You can do some Human plot")
}

s3_dispatch(plot(a))
s3_dispatch(plot(b))

predict.human_Univ <- function(x) {
  print("You can do some Human Univ Prediction")
}

s3_dispatch(predict(a))
s3_dispatch(predict(b))

'+.human' <- function(x,y) {
  return(x$age+y$age)
}
a+b
