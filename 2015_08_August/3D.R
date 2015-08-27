setwd("~/Talks/3D Interactive Scatterplots")

# Install and load required packages

install.packages(c("rgl", "car"))
library("car")


# Prepare the data
# We'll use the iris data set in the following examples :

data(iris)
head(iris)
sep.l <- iris$Sepal.Length
sep.w <- iris$Sepal.Width
pet.l <- iris$Petal.Length

# iris data set gives the measurements of the variables sepal length and width, 
# petal length and width, respectively, for 50 flowers from each of 3 species of 
# iris. The species are Iris setosa, versicolor, and virginica.


# The function scatter3d
# The simplified formats are:
# scatter3d(formula, data)
# scatter3d(x, y, z)

# x, y, z are respectively the coordinates of points to be plotted. The 
# arguments y and z can be optional depending on the structure of x.

# formula: a model formula of form y ~ x + z. If you want to plot the points by 
# groups, you can use y ~ x + z | g where g is a factor dividing the data into 
# groups

# data: data frame within which to evaluate the formula


# Basic 3D scatter plots
# 3D plot with the regression plane
library(car)
scatter3d(x = sep.l, y = pet.l, z = sep.w)

# Note that, the plot can be manually rotated by holding down on the mouse or 
# touchpad. It can be also zoomed using the scroll wheel on a mouse or pressing 
# ctrl + using the touchpad on a PC or two fingers (up or down) on a mac.

# Change point colors and remove the regression surface:

scatter3d(x = sep.l, y = pet.l, z = sep.w,
          point.col = "blue", surface=FALSE)

# Plot the points by groups

scatter3d(x = sep.l, y = pet.l, z = sep.w, groups = iris$Species)

# Remove the surfaces

# To remove the grids only, the argument grid = FALSE can be used as follow:

scatter3d(x = sep.l, y = pet.l, z = sep.w, groups = iris$Species,
          grid = FALSE)


# Note that, the display of the surface(s) can be changed using the argument 
# fit. Possible values for fit are "linear", "quadratic", "smooth" and 
# "additive"

scatter3d(x = sep.l, y = pet.l, z = sep.w, groups = iris$Species,
          grid = FALSE, fit = "smooth")

# Remove surfaces. The argument surface = FALSE is used.

scatter3d(x = sep.l, y = pet.l, z = sep.w, groups = iris$Species,
          grid = FALSE, surface = FALSE)

# Add concentration ellipsoids

scatter3d(x = sep.l, y = pet.l, z = sep.w, groups = iris$Species,
          surface=FALSE, ellipsoid = TRUE)

# Remove the grids from the ellipsoids:

scatter3d(x = sep.l, y = pet.l, z = sep.w, groups = iris$Species,
          surface=FALSE, grid = FALSE, ellipsoid = TRUE)

# Change point colors by groups

# The argument surface.col is used. surface.col is a vector of colors for the 
# regression planes.
# For multi-group plots, the colors are used for the regression surfaces and for
# the points in the several groups.

scatter3d(x = sep.l, y = pet.l, z = sep.w, groups = iris$Species,
          surface=FALSE, grid = FALSE, ellipsoid = TRUE,
          surface.col = c("#999999", "#E69F00", "#56B4E9"))

# It's also possible to use color palettes from the RColorBrewer package:

library("RColorBrewer")
colors <- brewer.pal(n=3, name="Dark2")
scatter3d(x = sep.l, y = pet.l, z = sep.w, groups = iris$Species,
          surface=FALSE, grid = FALSE, ellipsoid = TRUE,
          surface.col = colors)


# Axes
# Change axis labels:
# The arguments xlab, ylab and zlab are used:

scatter3d(x = sep.l, y = pet.l, z = sep.w,
          point.col = "blue", surface=FALSE,
          xlab = "Sepal Length (cm)", ylab = "Petal Length (cm)",
          zlab = "Sepal Width (cm)")

# Remove axis scales
# axis.scales = FALSE

scatter3d(x = sep.l, y = pet.l, z = sep.w,
          point.col = "blue", surface=FALSE, 
          axis.scales = FALSE)

# By default, different colors are used for the 3 axes. The argument axis.col is
# used to specify colors for the 3 axes:

scatter3d(x = sep.l, y = pet.l, z = sep.w, groups = iris$Species,
          surface=FALSE, grid = FALSE, ellipsoid = TRUE,
          axis.col = c("black", "black", "black"))

# Add text labels for the points
# The arguments below are used:
# labels: text labels for the points, one for each point
# id.n: Number of relatively extreme points to identify automatically

scatter3d(x = sep.l, y = pet.l, z = sep.w, 
          surface=FALSE, labels = rownames(iris), id.n=nrow(iris))


# Export images
# The plot can be saved as png or pdf.
library("rgl")
rgl.snapshot(filename = "plot.png")

# The function rgl.postscript() is used to save the screenshot to a file in ps,
# eps, tex, pdf, svg or pgf format:
# rgl.postscript("plot.pdf",fmt="pdf") TAKES A FEW MINUTES - SHOW PRIOR



# Source: http://www.sthda.com/english/wiki/print.php?id=210
