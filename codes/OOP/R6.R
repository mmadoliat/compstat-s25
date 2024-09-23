library(R6)

Accumulator <- R6Class("Accumulator", list(
  sum = 0,
  add = function(x = 1) {
    self$sum <- self$sum + x 
    invisible(self)
  })
)

Accumulator

x <- Accumulator$new() 

x$add(4) 
x$sum

x$add(10)$add(10)$sum

Person <- R6Class("Person", list(
  name = NULL,
  age = NA,
  initialize = function(name, age = NA) {
    stopifnot(is.character(name), length(name) == 1)
    stopifnot(is.numeric(age), length(age) == 1)
    self$name <- name
    self$age <- age
  }
))
hadley <- Person$new("Hadley", age = "thirty-eight")
hadley <- Person$new("Hadley", age = 38)

Person <- R6Class("Person", list(
  name = NULL,
  age = NA,
  initialize = function(name, age = NA) {
    self$name <- name
    self$age <- age
  },
  print = function(...) {
    cat("Person: \n")
    cat("  Name: ", self$name, "\n", sep = "")
    cat("  Age:  ", self$age, "\n", sep = "")
    invisible(self)
  }
))
hadley2 <- Person$new("Hadley")
hadley2

AccumulatorChatty <- R6Class("AccumulatorChatty", 
                             inherit = Accumulator,
                             public = list(
                               add = function(x = 1) {
                                 cat("Adding ", x, "\n", sep = "")
                                 super$add(x = x)
                               }
                             )
)

x2 <- AccumulatorChatty$new()
x2$add(10)$add(1)$sum

class(hadley2)

Person <- R6Class("Person", 
                  public = list(
                    initialize = function(name, age = NA) {
                      private$name <- name
                      private$age <- age
                    },
                    print = function(...) {
                      cat("Person: \n")
                      cat("  Name: ", private$name, "\n", sep = "")
                      cat("  Age:  ", private$age, "\n", sep = "")
                    }
                  ),
                  private = list(
                    age = NA,
                    name = NULL
                  )
)

hadley3 <- Person$new("Hadley")
hadley3
hadley3$name

Rando <- R6::R6Class("Rando", active = list(
  random = function(value) {
    if (missing(value)) {
      runif(1)  
    } else {
      stop("Can't set `$random`", call. = FALSE)
    }
  }
))
x <- Rando$new()
x$random <- 5

Person <- R6Class("Person", 
                  private = list(
                    .age = NA,
                    .name = NULL
                  ),
                  active = list(
                    age = function(value) {
                      if (missing(value)) {
                        private$.age
                      } else {
                        stop("`$age` is read only", call. = FALSE)
                      }
                    },
                    name = function(value) {
                      if (missing(value)) {
                        private$.name
                      } else {
                        stopifnot(is.character(value), length(value) == 1)
                        private$.name <- value
                        self
                      }
                    }
                  ),
                  public = list(
                    initialize = function(name, age = NA) {
                      private$.name <- name
                      private$.age <- age
                    }
                  )
)

hadley4 <- Person$new("Hadley", age = 38)
hadley4$name
hadley4$name <- 10
hadley4$age <- 20

y1 <- Accumulator$new() 
y2 <- y1
y1$add(10)
c(y1 = y1$sum, y2 = y2$sum)

y1 <- Accumulator$new() 
y2 <- y1$clone()
y1$add(10)
c(y1 = y1$sum, y2 = y2$sum)
