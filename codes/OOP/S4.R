library(methods)

setClass("Person", 
         slots = c(
           name = "character", 
           age = "numeric"
         )
)

john <- new("Person", name = "John Smith", age = NA_real_)

is(john)
john@name
slot(john, "age")

setGeneric("age", function(x) standardGeneric("age"))
setGeneric("age<-", function(x, value) standardGeneric("age<-"))

setMethod("age", "Person", function(x) x@age)
setMethod("age<-", "Person", function(x, value) {
  x@age <- value
  x
})
age(john) <- 50
age(john)

sloop::otype(john)
sloop::ftype(age)

setClass("Person", 
         slots = c(
           name = "character", 
           age = "numeric"
         ), 
         prototype = list(
           name = NA_character_,
           age = NA_real_
         )
)
me <- new("Person", name = "Hadley")
str(me)

setClass("Employee", 
         contains = "Person", 
         slots = c(
           boss = "Person"
         ),
         prototype = list(
           boss = new("Person")
         )
)
str(new("Employee"))

is(new("Person"))
is(new("Employee"))

is(john, "Person")

setClass("A", slots = c(x = "numeric"))
a <- new("A", x = 10)
setClass("A", slots = c(a_different_slot = "numeric"))
a

Person <- function(name, age = NA) {
  age <- as.double(age)
  new("Person", name = name, age = age)
}
Person("Hadley")

Person(mtcars)

Person("Hadley", age = c(30, 37))

setValidity("Person", function(object) {
  if (length(object@name) != length(object@age)) {
    "@name and @age must be same length"
  } else {
    TRUE
  }
})

alex <- Person("Alex", age = 30)
alex@age <- 1:10

validObject(alex)

args(getGeneric("show"))

setMethod("show", "Person", function(object) {
  cat(is(object)[[1]], "\n",
      "  Name: ", object@name, "\n",
      "  Age:  ", object@age, "\n",
      sep = ""
  )
})
john

setGeneric("name", function(x) standardGeneric("name"))
setMethod("name", "Person", function(x) x@name)
name(john)

setGeneric("name<-", function(x, value) standardGeneric("name<-"))
setMethod("name<-", "Person", function(x, value) {
  x@name <- value
  validObject(x)
  x
})
name(john) <- "Jon Smythe"
name(john)
name(john) <- letters
