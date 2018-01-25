library(drake)
load_basic_example()
plot_graph(my_plan)
outdated(my_plan)
max_useful_jobs(my_plan)
make(my_plan)
plot_graph(my_plan)
reg2 = function(d){ # Change one of your functions.
  d$x3 = d$x^3
  lm(y ~ x3, data = d)
}
outdated(my_plan) # Some targets depend on reg2().
plot_graph(my_plan)
