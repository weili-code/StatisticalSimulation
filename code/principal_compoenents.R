### This script contains some examples for Principal Components###


set.seed(123)
n <- 10
mean <- c(0, 0)
sigma <- matrix(c(3, 1.5, 1.5, 2), 2)  # Covariance matrix
X <- MASS::mvrnorm(n, mean, sigma)

plot(X[,1], X[,2], xlab="X1", ylab="X2", main="Simulated Data with Principal Components")

cov_X <- cov(X)
eigenvalues_eigenvectors <- eigen(cov_X)
eigenvectors <- eigenvalues_eigenvectors$vectors

arrows(0, 0, eigenvectors[1,1]*3, eigenvectors[2,1]*3, col="red", lwd=2)
arrows(0, 0, eigenvectors[1,2]*3, eigenvectors[2,2]*3, col="blue", lwd=2)

pca <- prcomp(X, center = TRUE, scale. = FALSE)
Z <- pca$x

projected1 <- X %*% matrix(eigenvectors[,1], ncol = 1)
projected2 <- X %*% matrix(eigenvectors[,2], ncol = 1)

segments(X[,1], X[,2], projected1 * eigenvectors[1,1], projected1 * eigenvectors[2,1], col = "red", lty=2)
segments(X[,1], X[,2], projected2 * eigenvectors[1,2], projected2 * eigenvectors[2,2], col = "blue", lty=2)



##### Dataset ############# 

# Load necessary libraries
install.packages("ggplot2")
library(ggplot2)

# Load the mtcars dataset
data(mtcars)
summary(mtcars)

# 1. Data Preparation
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

# spree plot and plot of variances explained

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

# In the scree plot, you should look for an "elbow" where the 
# eigenvalues start to level off. The principal components to the 
# left of this elbow retain most of the useful information, and those to the right 
# can often be omitted without losing much information.
# 
# In the variance explained plot, the aim is often to determine how many 
# principal components are needed to explain a certain proportion of the 
# variance (e.g., 95%). The horizontal line at 0.95 will help you visually identify this.

####
####

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
arrows(0, 0, loadings[,1], loadings[,2], length=0.1, angle=10, col="blue")

# Add text labels
text(loadings, labels=rownames(loadings), cex=0.7, pos=3)

# Add a circle with radius 1
symbols(x=0, y=0, circles=1, add=TRUE, inches=FALSE, lty=2)
