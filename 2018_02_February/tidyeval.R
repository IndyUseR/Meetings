# Material: https://www.youtube.com/watch?v=nERXS3ssntw

# 1. R LANGUAGE IS A TREE

quote(f(x, "y", 1))
str(quote(f(x, "y", 1)))
expr <- parse(
  text = 'f(x, "y", 1)',
  keep.source = FALSE
)
str(expr)
str(expr[[1]])
View(expr)
View(quote(f(g(x, z), "y", 1)))
lobstr::ast(f(g(x, z), "y", 1))
# pryr::ast(f(g(x, z), "y", 1))
View(quote(y <- x * 10))

# Walk abstract syntax trees (ASTs):
# http://adv-r.had.co.nz/Expressions.html#ast-funs

# 2. DIFFERENT WAYS TO QUOTE

library(rlang)

f1 <- function(x) expr(x)
f2 <- function(x) enexpr(x)

f1(y + z)
f2(y + z)

# What about outside a function?
expr(x + y)
enexpr(x + y)

expr
enexpr

# Evaluated
mean(x + y) # requires x and y to have values

# Quoted
library(ggplot2)
ggplot(mtcars, aes(disp, mpg)) + # disp and mpg are quoted
  geom_point()

library(dplyr)
mtcars %>%
  filter(cyl > 2) %>%
  select(mpg:hp) %>%
  head(10)

# 3. UNQUOTING TO CHAIN TOGETHER LANGUAGE TREES

# install.packages("rlang_0.2.0.tar.gz", type = "source", repos = NULL)
xy <- expr(x + y)
expr(!!xy + z)
expr(1/!!xy)
expr(f(!!xy, y))

# 4. QUOTING AND UNQUOTING TOGETHER

my_summarize <- function(df, var){
  var <- enexpr(var)                    # Quote first.
  summarize(df, mean(!!var), sd(!!var)) # Then unquote.
}

my_summarize(mtcars, mpg)

#. 5. QUOSURES

my_mutate <- function(df, var){
  n <- 10 # `n` defined here
  var <- enexpr(var)
  mutate(df, !!var)
}

df <- tibble(x = 1)
n <- 100 # `n` also defined here
my_mutate(df, x + n) # Which `n` is used?

# Quosures capture the environment
# where the expression is DEFINED.
# Agrees with R's internals
# Related to lexical scoping.

my_mutate <- function(df, var){
  n <- 10
  var <- enquo(var) # Use enquo() instead of enexpr().
  mutate(df, !!var)
}

df <- tibble(x = 1)
n <- 100
my_mutate(df, x + n)
