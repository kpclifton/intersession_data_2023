
#### Lesson 1 (01-25-23): general preliminary investigation of data and data visualization ####

data <- read.csv('data/charmander.csv.gz')
dim(data) #get dimensions of data (rows and columns)
head(data) #see the first n rows and all the column names
data[1:5,1:5] #if too many columns can index first n rows and first n columns

summary(data[1:5,1:5]) #see summary statistics of each column, stats provided differ a bit based on data class
class(data$x_centroid) #numeric
class(data$X) #integer

library(ggplot2)
ggplot(data = data) +
  geom_point(aes(x = x_centroid, y = y_centroid), size=0.1) +
  theme_classic() 

ggplot(data = data) +
  geom_histogram(aes(x = area)) +
  theme_classic()

#### Lesson 2 (01-27-23): exploring spatial transcriptomics data ####

## storing gene expession data in separate data frame for convenience
gexp <- data[,5:ncol(data)]
gexp[1:5,1:5]

## each row is a different cell and each column is a different gene
dim(gexp) # return number of cells and number of genes

## rows can be summed to total number of genes expressed per cell 
length(rowSums(gexp)) # length of vector is the number of cells 

## columns can be summed to total number of copies of each gene detected across all cells 
length(colSums(gexp)) # length of vector is the number of genes

par(mfrow=c(1,2)) #set graphical parameters so subsequent figures will be drawn in an 1 row x 2 column array 
hist(rowSums(gexp), main = 'total genes per cell')
hist(colSums(gexp), main = 'total copies of each gene across all cells')

## log transform gene expression 
mat <- log10(gexp+1)

par(mfrow=c(1,2)) #set graphical parameters so subsequent figures will be drawn in an 1 row x 2 column array 
hist(colSums(gexp), main = 'total copies of each gene across all cells')
hist(colSums(mat), main = 'log total copies of each gene across all cells')

# ## determining gene prevalence by totaling how many cells have non-zero expression for each gene
# gprev <- colSums(gexp != 0)
# 
# par(mfrow=c(1,1)) #set graphical parameters so subsequent figures will be drawn in an 1 row x 1 column array 
# hist(gprev, main = 'prevalence of genes across cells')
# 
# ## ordering genes in decreasing order of prevalence
# prevalent.genes <- sort(gprev, decreasing = TRUE)
# prevalent.genes[1] #most prevalent gene
# prevalent.genes[length(prevalent.genes)] #least prevalent gene
# 
# ## plotting cells positions colored by which cells have non-zero expression of ERBB2
# p1 <- ggplot(data = data) +
#   geom_point(aes(x = x_centroid, y = y_centroid,
#                  col = ERBB2 > 0), size=0.1) +
#   theme_classic()
# 
# p2 <- ggplot(data = data) +
#   geom_point(aes(x = x_centroid, y = y_centroid,
#                  col = CRHBP > 0), size=0.1) +
#   theme_classic()
# 
# grid.arrange(p1, p2)

## determining variance of gene expression by applying the variance function on each column of gexp
gvar <- sapply(gexp, var)

par(mfrow=c(1,1)) #set graphical parameters so subsequent figures will be drawn in an 1 row x 1 column array 
hist(gvar, main = 'variance of genes')

## ordering genes in decreasing order of prevalence
variable.genes <- sort(gvar, decreasing = TRUE)
variable.genes[1:5] #most variable genes
variable.genes[(length(variable.genes)-4):length(variable.genes)] #least variable genes

par(mfrow=c(1,4))
hist(gexp$ERBB2)
hist(gexp$LUM)
hist(gexp$MPO)
hist(gexp$CRHBP)

## scaled gene expression 
mat <- scale(gexp) #center and divide the centered columns of gexp (genes) by their standard deviations
class(mat)
class(gexp)
mat <- data.frame(mat)

## determining variance of scaled gene expression by applying the variance function on each column of mat
gvar2 <- sapply(mat, var)

par(mfrow=c(1,2)) #set graphical parameters so subsequent figures will be drawn in an 1 row x 1 column array 
hist(gvar2, main = 'variance of scaled genes')
hist(gvar, main = 'variance of genes')

## ordering genes in decreasing order of prevalence
variable.genes <- sort(gvar2, decreasing = TRUE)
variable.genes[1:5] #most variable genes
variable.genes[(length(variable.genes)-4):length(variable.genes)] #least variable genes

par(mfrow=c(1,4))
hist(mat$MYBPC1)
hist(mat$CLDN4)
hist(mat$MMP1)
hist(mat$CD1C)

# ## relationship between ERBB2 vs. LUM
# ggplot(data = data) +
#   geom_point(aes(x = ERBB2, y = LUM), size=0.1) +
#   theme_classic() 
# 
# ## plotting cells positions colored by expression of ERBB2 and LUM
# p1 <- ggplot(data = data) +
#   geom_point(aes(x = x_centroid, y = y_centroid, col = ERBB2), size=0.1) +
#   scale_color_gradient(low = 'lightgrey', high='red') + 
#   theme_classic()
# 
# p2 <- ggplot(data = data) +
#   geom_point(aes(x = x_centroid, y = y_centroid, col = LUM), size=0.1) +
#   scale_color_gradient(low = 'lightgrey', high='red') +
#   theme_classic()
# 
# grid.arrange(p1, p2)

#### Lesson 4 (02-01-23): PCA ####

mat <- scale(log10(gexp+1)) ##scale and log-transformed counts
mat <- mat*1000 #how many genes would be in a cell if each cell expressed 1000 genes

##PCA
pcs <- prcomp(mat) ##centered and log-transformed counts

#should be as many PCs as genes
dim(mat) #cells x genes
dim(pcs$x) #cells x PCs

#scree plot: look at how much variance is explained by the first n PCs, n = 10
par(mfrow=c(1,1))
plot(1:10, pcs$sdev[1:10], type = "l")

## look at loadings (coefficients of the linear combination of genes for each PC)
head(sort(abs(pcs$rotation[,1]), decreasing=TRUE)) #highest absolute loadings, genes that contribute the most to PC1, 
head(sort(abs(pcs$rotation[,2]), decreasing=TRUE)) #genes that contribute the most to PC2,
head(sort(abs(pcs$rotation[,3]), decreasing=TRUE)) #genes that contribute the most to PC3,
head(sort(abs(pcs$rotation[,4]), decreasing=TRUE)) #genes that contribute the most to PC4,

## make a data visualization to explore our first two PCs
df <- data.frame(pcs$x[,1:2], gene = mat[,'KRT7']) 
head(df)
ggplot(data = df, aes(x = PC1, y = PC2, col=gene)) + 
  geom_point() + 
  theme_classic() +
  scale_color_gradient(low = 'lightgrey', high='red')

#### Lesson 5 (02-03-23): tSNE ####

library(Rtsne)
set.seed(0) ## for reproducibility
emb <- Rtsne(mat)
dim(mat)
dim(emb$Y)

df2 <- data.frame(emb$Y, gene1=mat[,'KRT7'], gene2=mat[,'LUM'])
head(df)
p1 <- ggplot(df2, aes(x = X1, y = X2, col=gene1)) + geom_point(size = 0.1) + 
  scale_color_gradient(low = 'lightgrey', high='red') + ggtitle('KRT7') + theme_classic()

