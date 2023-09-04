### This script contains some examples for Principal Components###

# ref. http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

#####################################
####### PCA     ############# #######
#####################################

# Load required library
library(MASS)

# 1. Generate 2D data using bivariate normal
set.seed(123)  # For reproducibility
mu <- c(5, 5)  # mean for x and y
sigma <- matrix(c(4, 2, 2, 3), 2)  # covariance matrix

data <- mvrnorm(50, mu, sigma)
x <- data[,1]
y <- data[,2]

# Get the range for both x and y axes
plot_range <- range(c(x, y))

# 2. Plotting the Scatter Plot
plot(x, y, xlab="X", ylab="Y", pch=19, main="Scatter Plot with PC Directions", xlim=plot_range, ylim=plot_range)

# 3. Calculating Eigenvectors using PCA
pca_result <- prcomp(data, center=TRUE, scale.=TRUE)
eigenvectors <- pca_result$rotation[, 1:2]
# column j in eigenvectors is the j-th PC direction

# 4. Plotting the Eigenvectors with transparent colors
# Arrows for the positive direction of the PCs
# all vectors are scaled by 2 for better visualization
arrows(mean(x), mean(y), mean(x) + eigenvectors[1, 1]*2, mean(y) + eigenvectors[2, 1]*2, col=rgb(1,0,0,alpha=0.5), lwd=2)
arrows(mean(x), mean(y), mean(x) + eigenvectors[1, 2]*2, mean(y) + eigenvectors[2, 2]*2, col=rgb(0,0,1,alpha=0.5), lwd=2)
# Arrows for the negative direction of the PCs
arrows(mean(x), mean(y), mean(x) - eigenvectors[1, 1]*2, mean(y) - eigenvectors[2, 1]*2, col=rgb(1,0,0,alpha=0.5), lwd=2)
arrows(mean(x), mean(y), mean(x) - eigenvectors[1, 2]*2, mean(y) - eigenvectors[2, 2]*2, col=rgb(0,0,1,alpha=0.5), lwd=2)

legend("topleft", legend=c("1st PC", "2nd PC"), col=c(rgb(1,0,0,alpha=0.5), rgb(0,0,1,alpha=0.5)), lwd=2)

# 
selected_point <- data[25,]

# Plot the selected point in green
points(selected_point[1], selected_point[2], col="green", pch=19, cex=1.5)

# Calculate the projections of the selected point onto the PCA directions
proj_PC1 <- sum((selected_point - c(mean(x), mean(y))) * eigenvectors[,1]) * eigenvectors[,1] + c(mean(x), mean(y))
proj_PC2 <- sum((selected_point - c(mean(x), mean(y))) * eigenvectors[,2]) * eigenvectors[,2] + c(mean(x), mean(y))

# Plot the projection on PC1 using a red dashed line and mark the end point with a red circle
segments(selected_point[1], selected_point[2], proj_PC1[1], proj_PC1[2], col="red", lty=2)
points(proj_PC1[1], proj_PC1[2], col="red", pch=19, cex=1.2)

# Plot the projection on PC2 using a blue dashed line and mark the end point with a blue circle
segments(selected_point[1], selected_point[2], proj_PC2[1], proj_PC2[2], col="blue", lty=2)
points(proj_PC2[1], proj_PC2[2], col="blue", pch=19, cex=1.2)


#####################################
####### Dataset ############# #######
#####################################

# The functions prcomp() uses the singular value decomposition (SVD).

# Load necessary libraries
library(ggplot2)

# Load the mtcars dataset
data(mtcars)
# Write mtcars to a .txt file
# write.table(mtcars, file = "mtcars.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

mtcars <- mtcars[, -1] # remove column "mpg"
summary(mtcars)

# 1. Data Preparation
# In principal component analysis, variables are often scaled (i.e. standardized). 
# Scale the data to have mean 0 and unit variance
scaled_mtcars <- scale(mtcars)

# 2. PCA
# Compute principal components
# The functions prcomp() uses the singular value decomposition (SVD).
pca_result <- prcomp(scaled_mtcars, center = FALSE, scale. = FALSE)

# Display the summary of the PCA results
summary(pca_result)
# you have Standard deviation, Proportion of Varianc and Cumulative Proportion

# 3. Visualization
# Extract the first two principal components
# pca_result$x: n by num_cps 
PC1 <- pca_result$x[,1] # a vector of PC scores
PC2 <- pca_result$x[,2]

# Create a data frame with the principal components
df_pca <- data.frame(Car = rownames(mtcars), PC1, PC2)

# Plot the first two principal components with dots
ggplot(df_pca, aes(x=PC1, y=PC2, label=Car)) +
  geom_point(aes(color=Car), size=3) +  # this adds the dots
  geom_text(aes(label=Car), vjust=2, hjust=0.5, size=4) +  # adjusted for clarity with dots
  ggtitle("First two Principal Components of mtcars Dataset") +
  xlab("First Principal Component") +
  ylab("Second Principal Component") +
  theme_minimal() +
  theme(legend.position="none")  
# Individuals with a similar profile are grouped together.
# The colors of individuals points are pure for aesthetic purpose and do not have any meaning.

#####################################
# spree plot and plot of variances explained
#####################################

# Eigenvalues can be used to determine the number of principal components to retain after PCA 

# (1). Scree Plot
eigenvalues <- pca_result$sdev^2

# Plot eigenvalues
plot(eigenvalues, type="b", main="Scree Plot", 
     xlab="Principal Component", 
     ylab="Eigenvalue",
     pch=19, col="blue")

# (2). Variance Explained Plot

prop_varex <- eigenvalues / sum(eigenvalues)
cum_varex <- cumsum(prop_varex)

# Plot variance explained
plot(cum_varex, type="b", main="Variance Explained Plot", 
     xlab="Principal Component", 
     ylab="Cumulative Proportion of Variance Explained",
     pch=19, col="red", ylim=c(0,1))
abline(h=0.95, col="blue", lty=2) # line at 95% variance explained

# In the scree plot, search for a noticeable "bend" where the eigenvalues 
# begin to stabilize. The principal components situated before this bend capture 
# the bulk of the essential data, while those after it can frequently be disregarded 
# without sacrificing significant information.
# 
# For the variance explained graph, the goal is typically to identify the number 
# of principal components required to account for a specific percentage of the variance, 
# like 95%. A horizontal line drawn at 0.95 will assist in visually pinpointing this.


#####################################
#### individual case' contribution ##
#####################################

# Compute the scores
scores <- pca_result$x

# Compute the squared scores
squared_scores <- scores^2

# Compute the contribution of each observation to PC1 and PC2
contribution_PC1 <- squared_scores[,1] / sum(squared_scores[,1])
contribution_PC2 <- squared_scores[,2] / sum(squared_scores[,2])

# Print the contributions
print(data.frame(
                 Contribution_to_PC1 = contribution_PC1, 
                 Contribution_to_PC2 = contribution_PC2))

#####################################
####. variables correlation plot ####
#####################################. 

# 1. Extract Loadings
# pca_result$rotation: the matrix V, num_variables x num_pcs
loadings <- pca_result$rotation[, 1:2]

# 2. Visualize Loadings
# a visual representation of the contribution of each variable to 
# the first two principal components.

# Original plot with lines indicating axes
plot(loadings, type="n", xlim=c(-1,1), ylim=c(-1,1), 
     xlab="PC1", ylab="PC2", main="Loadings Plot")
abline(h=0, v=0, col="gray", lty=2)

# Add vectors
arrows(0, 0, loadings[,1], loadings[,2], angle=10, col="blue")

# Add text labels
text(loadings, labels=rownames(loadings), cex=0.7, pos=3)

# Add a circle with radius 1
symbols(x=0, y=0, circles=1, add=TRUE, inches=FALSE, lty=2)

# Variables that are close to each other are positively correlated.
# Variables that are orthogonal (perpendicular) to each other are not correlated.
# Variables that are on opposite ends of a plot (e.g., opposite sides of the origin along the same line) are negatively correlated.
# The closer a variable's loading is to the edge of the plot (a circle of radius 1), the more influence it has on that particular principal component.
# Variables that are closed to the center of the plot are less important for the first few components.

#####################################
###### contribution of variables ####
#####################################

## (1) contribution of the variable ##

# Compute squared loadings
squared_loadings <- pca_result$rotation^2 # num_case * num_pcs

# Compute the contribution of each variable to PC1 and PC2
contribution_PC1 <- squared_loadings[,1] / sum(squared_loadings[,1])
contribution_PC2 <- squared_loadings[,2] / sum(squared_loadings[,2])

# Print the contributions
print(data.frame( Contribution_to_PC1 = contribution_PC1, 
                 Contribution_to_PC2 = contribution_PC2))

## (2) contribution of the variable in variance ##

# Get eigenvalues (variances explained by each PC)
eigenvalues <- pca_result$sdev^2

# Compute the weighted total contribution of each variable across all PCs
weighted_total_contribution <- colSums(squared_loadings %*% diag(eigenvalues))

# Compute proportional total contribution
proportional_weighted_total_contribution <- weighted_total_contribution / sum(eigenvalues)

# Print the weighted total contributions
print(data.frame(Variable = rownames(squared_loadings), 
                 Weighted_Total_Contribution = weighted_total_contribution, 
                 Proportional_Weighted_Total_Contribution = proportional_weighted_total_contribution))
