###
### Script for running a clustering model
### (the variational Bayesian Dirichlet process Gaussian mixture model, VB-DP-GMM)
###
### Contents:
###  - Load the required functions
###  - Set up the model settings
###  - Optimize the model hyper-parameters
###  - Fit the model with optimal hyper-parameters and select the best model
###  - Inspect the cluster assignments
###  - Inspect the observations in a cluster
###
### Input:
### X: Numeric matrix with samples (e.g., subjects) as rows and variables (e.g., metabolomic peaks) as columns. 'colnames( X )' should correspond to the names of the variables.
### covariate: Factor of sample group annotations (e.g., "Control" and "Case"). 'length( covariate )' should match 'nrow( X )'. Used only for making the PDF figure of clusters.
###
### Version history:
###
### 05.11.15 Added: Functionality for making a PDF figure of all the clusters.
### 16.10.15 File created.
###
###
### Author:
###
### Tommi Suvitaival
### tsvv@steno.dk
###



### INPUT


# LOAD YOUR DATA X HERE.


# REPLACE OR COMMENT OUT! Only creates a toy data set.

# X <- array( data=rnorm( n=10*50 ), dim=c( 10, 50 ) )
# X[ 1:5, 1:10 ] <- X[ 1:5, 1:10 ] + 2
# X[ 6:10, 11:20 ] <- X[ 6:10, 11:20 ] + 2
# X[ c( 1, 3, 4, 7, 9 ), 21:30 ] <- X[ c( 1, 3, 4, 7, 9 ), 21:30 ] - 2
# colnames( X ) <- c( paste( "PC(", 1:10, ":0)", sep="" ), paste( "TG(", 1:10, ":0)", sep="" ), paste( "Unknown (", 1:30, ")", sep="" ) )
X<-(as.matrix(CIMT_for_clustering_discovery))
# covariate <- factor( x=rep( x=c( "Control", "Case" ), each=ceiling( nrow( X ) / 2 ), length.out=nrow( X ) ), levels=c( "Control", "Case" ) )
#covariate<-factor(Medim_new$Country,levels = as.character(unique(Medim_new$Country)))
X.non.scaled <- t(X)
X <- scale( x = X.non.scaled )
#rownames(X)<-Medim_new$Compound.Name
# END OF REPLACE OR COMMENT OUT!


opts.clustering <- list()

## File paths for saving results.

# File path for saving the optimization result (optional)
opts.clustering$optimization$file.path.result <- file.path( "Z:", "AQAI - Ashfaq Ali", "CIMT_trial","Clustering","result_clustering_optimization_clust.RData" )


# File path for saving the optimal model (optional)
opts.clustering$fit$file.path.result <- file.path( "Z:", "AQAI - Ashfaq Ali", "CIMT_trial","Clustering", "result_clustering_fit.RData" )



### LIBRARIES


# The clustering algorithm comes from the "netresponse" bioconductor package.
# Install the required package if necessary.
# See the package documentation for details of the clustering algorithm.

if ( !( "netresponse" %in% .packages( all.available=TRUE ) ) ) {
    
    source( file="http://bioconductor.org/biocLite.R" )
    biocLite()
    biocLite( pkgs="netresponse" )
    
}

require( package="netresponse" )

# If necessary, install the package by running the lines below.

# Load the wrapper functions for optimizing and fitting the cluster model.
source( file=file.path( "Z:", "__Scripts", "R", "TSVV", "clustering", "vb_dp_gmm", "vb_dp_gmm.R" ) )

# Load the function for analyzing the cluster enrichments.
source( file=file.path( "Z:", "__Scripts", "R", "TSVV", "clustering", "compute_cluster_lipid_composition", "compute_cluster_lipid_composition.R" ) )

# Load the function for creating a table of the cluster assignments.
source( file=file.path( "Z:", "__Scripts", "R", "TSVV", "clustering", "cluster_map_assignments", "cluster_map_assignments.R" ) )

# Load the function for making a figure of the observations in the clusters.
source( file=file.path( "Z:", "__Scripts", "R", "TSVV", "clustering", "plot_observations_in_clusters", "plot_observations_in_clusters.R" ) )



### SETTINGS



# The pre-set values below are likely reasonable.


## Optimization of the hyper-parameters


# Number of repeats to run for the grid of hyper-parameter optimization
opts.clustering$optimization$N.repeats <- 10

# Tested values for the hyper-parameters
opts.clustering$optimization$prior.alpha <- seq( from=0.1, to=3, length.out=10 )
opts.clustering$optimization$prior.alphaKsi <- seq( from=0.1, to=3, length.out=10 )
opts.clustering$optimization$prior.betaKsi <- seq( from=0.1, to=3, length.out=10 )


## Selection of the optimal model


# Number of repeated runs of the hyper-parameter-optimized model
opts.clustering$fit$N.repeats <- 1000



### MODEL



## Optimization of the hyper-parameters


# Run the optimization
# This may take several hours, depending on the size of the optimization grid and the size of the data.
result.clustering.optimization <- vb_dp_gmm_optimize_hyperparameters(
    x = X,
    prior.alpha = opts.clustering$optimization$prior.alpha,
    prior.alphaKsi = opts.clustering$optimization$prior.alphaKsi,
    prior.betaKsi = opts.clustering$optimization$prior.betaKsi,
    N.repeats = opts.clustering$optimization$N.repeats,
    cluster.columns = TRUE,
    c.max = NULL,
    do.sort = TRUE,
    initial.K = 1,
    implicit.noise = 0,
    ite = Inf,
    min.size = 1,
    speedup = TRUE,
    threshold = 1e-05,
    verbose = TRUE
)


## Save/load the optimization result (optional)


# Save the result of the optimization (commented out; not run).
# If used, the given folder needs to exist.
save( result.clustering.optimization, file=opts.clustering$optimization$file.path.result, compress=TRUE )

# Load an earlier result of the optimization (commented out; not run).
# If used, the given file needs to exist.
 load( file=opts.clustering$optimization$file.path.result )


## Selection of the optimal model


# Fit the cluster model repeatedly with the optimal hyper-parameter setting
# and select the model with the best fit.
# This may take few hours, depending on the number of repeats and the size of the data.

result.clustering.fit <- vb_dp_gmm_fit(
    x = X,
    prior.alpha = result.clustering.optimization$optimal$param$prior.alpha,
    prior.alphaKsi = result.clustering.optimization$optimal$param$prior.alphaKsi,
    prior.betaKsi = result.clustering.optimization$optimal$param$prior.betaKsi,
    N.repeats = opts.clustering$fit$N.repeats,
    cluster.columns = TRUE,
    c.max = NULL,
    do.sort = TRUE,
    initial.K = 1,
    implicit.noise = 0,
    ite = Inf,
    min.size = 1,
    speedup = TRUE,
    threshold = 1e-05,
    verbose = TRUE
)


## Save/load the optimal model


# Save the result of the selection of the optimal model (commented out; not run).
# If used, the given folder needs to exist.
 save( result.clustering.fit, file=opts.clustering$fit$file.path.result, compress=TRUE )

# Load an earlier result of the selection of the optimal model (commented out; not run).
# If used, the given file needs to exist.
load( file=opts.clustering$fit$file.path.result )



### ANALYSIS



## Extract the cluster assignments


# Find the most likely cluster assignment for each item (i.e., column of the matrix 'X').
# (Technically, find the index of the maximum value at each row of the items-by-clusters assignment matrix.)
# Returns a vector with the length equal to the number of columns in the matrix X.
cluster.assignments <- apply( X = result.clustering.fit$optimal$posterior$qOFz,
                              MAR = 1,
                              FUN = which.max
)

# View the count of items in the clusters (commented out).
# table( cluster.assignments )


## The following work only if 'colnames( X )' correspond to the names of the variables.


## Create a table of the cluster assignments


cluster.assignments.table <- cluster_map_assignments(
    x = result.clustering.fit$optimal$posterior$qOFz,
    item.information = NULL,
    cluster.index.decreasing = TRUE
)

# View the cluster assignments table (commented out).
# View( cluster.assignments.table )


## Inspect the enrichment of functional groups/compound classes in the clusters.


# (i.e., multinomial test for each cluster-category pair)
# With LC-MS data, the function extracts the categories automatically from the names of the items.
# With GCxGC data, the categories will have to be supplied in 'names( cluster.assignments )' instead of the original item names (with 'category.as.name=TRUE').

cluster.enrichments <- compute_cluster_lipid_composition(
    clusters = cluster.assignments,
    p.adjustment = "BH",
    category.as.name = FALSE
)

# View the q-values of the enrichment test (commented out).
# View( cluster.enrichments$p.value.adjusted )


## Compute the cluster averages.


# Initialize the samples-by-clusters matrix.
X.cluster.mean <- array( dim=c( nrow( X ), max( cluster.assignments ) ) )
colnames( X.cluster.mean ) <- paste( "Cluster", 1:ncol( X.cluster.mean ) )
rownames( X.cluster.mean ) <- rownames( X )

X.cluster.sd <- X.cluster.mean

# Compute the sample-specific mean for each cluster.
for ( k in 1:max( cluster.assignments ) ) { # Go through the clusters.
    
    # Find the items assigned to the cluster 'k'.
    items.k <- which( cluster.assignments == k )
    
    if ( length( items.k ) > 0 ) {
        
        # Compute the sample-specific mean over the items assigned to cluster 'k'.
        X.cluster.mean[ , k ] <- rowMeans( X[ , items.k, drop=FALSE ] )
        X.cluster.sd[ , k ] <- apply( X=X[ , items.k, drop=FALSE ], MAR=1, FUN=sd )
        
    }
}


## Plot the observations in a cluster


# (Samples are on the x-axis, sorted by the sample-specific cluster mean;
# each observation in the cluster shown with a diagonal cross)

opts.plot <- list()

# Select the cluster to be plotted.
# Change this to plot one of the clusters from 1 to max( cluster.assignments ).
opts.plot$k <- 1

opts.plot$col.mean <- "#009FDAFF" # light blue
opts.plot$col.obs <- "#0019653F" # dark blue with 50 % opacity
opts.plot$col.ci <- "#E0DED8FF" # light gray
opts.plot$col.bl <- "#82786FFF" # dark gray
opts.plot$lwd.mean <- 1
opts.plot$pch.obs <- 4
opts.plot$cex.obs <- 1
opts.plot$lty.bl <- 3 # dotted

# Find the items that are assigned to cluster 'k'.
items.k <- which( cluster.assignments == opts.plot$k )
rank.samples.by.mean.k <- rank( x=result.clustering.fit$optimal$posterior$centroids[ opts.plot$k, ] )

# Initialize the plot.
plot( x = NA,
      xlim = c( 1, nrow( X ) ),
      ylim = range( X[ , items.k ], min( result.clustering.fit$optimal$posterior$centroids[ opts.plot$k, ]-2*result.clustering.fit$optimal$posterior$sds[ opts.plot$k, ]), max( result.clustering.fit$optimal$posterior$centroids[ opts.plot$k, ]+2*result.clustering.fit$optimal$posterior$sds[ opts.plot$k, ] ) ),
      main = paste( "Cluster ", opts.plot$k, " - ", length( items.k ), " Peaks", sep = "" ),
      xlab = "Samples",
      ylab = "Observation"
)

# Draw a line that corresponds to the average of the cluster over all samples.
# (i.e., the "baseline")
abline( h = mean( X[ , items.k ] ),
        col = opts.plot$col.bl,
        lwd = opts.plot$lwd.mean,
        lty = opts.plot$lty.bl
)

# Draw the confidence interval, the observations and the mean of each sample.

for ( i in 1:nrow( X ) ) { # Go through the samples.
    
    # Draw the confidence interval for sample 'i'.
    rect( xleft = rank.samples.by.mean.k[ i ]-0.25,
          xright = rank.samples.by.mean.k[ i ]+0.25,
          ybottom = result.clustering.fit$optimal$posterior$centroids[ opts.plot$k, i ]-2*result.clustering.fit$optimal$posterior$sds[ opts.plot$k, i ],
          ytop = result.clustering.fit$optimal$posterior$centroids[ opts.plot$k, i ]+2*result.clustering.fit$optimal$posterior$sds[ opts.plot$k, i ],
          col = opts.plot$col.ci,
          border = NA
    )
    
    # Draw the observations for sample 'i'.
    points( x = rep( x=rank.samples.by.mean.k[ i ], times=length( items.k ) ),
            y = X[ i, items.k ],
            col = opts.plot$col.obs,
            pch = opts.plot$pch.obs,
            cex = opts.plot$cex.obs
    )
    
    # Draw the mean for sample 'i'.
    segments( x0 = rank.samples.by.mean.k[ i ]-1/3,
              x1 = rank.samples.by.mean.k[ i ]+1/3,
              y0 = result.clustering.fit$optimal$posterior$centroids[ opts.plot$k, i ],
              col = opts.plot$col.mean,
              lwd = opts.plot$lwd.mean
    )
    
}

# Create a legend for the figure.
legend( x="bottomright",
        legend=c( "Mean (Sample)", "Confidence Interval", "Observation", "Mean (Global)" ),
        col=c( opts.plot$col.mean, opts.plot$col.ci, opts.plot$col.obs, opts.plot$col.bl ),
        pch=c( NA, NA, opts.plot$pch.obs, NA ),
        lty=c( 1, 1, NA, opts.plot$lty.bl ),
        lwd=c( opts.plot$lwd.mean, opts.plot$lwd.mean*3, NA, opts.plot$lwd.mean )
)

# Make a PDF figure of all the clusters.

if ( FALSE ) { # Not run.
    
    plot_observations_in_clusters( x = X,
                                   model = result.clustering.fit$optimal,
                                   file = file.path( "Z:", "AQAI - Ashfaq Ali", "MEDIM", "Medim_new_data", "cluster_profiles_3.pdf" ),
                                   covariate = covariate,
                                   sort.samples = TRUE,
                                   order.clusters.inverted = TRUE,
                                   ncol.figure = 4,
                                   col.mean = c( "#001965FF", "#C81F49FF", "#F58220FF", "#00AF41FF" ), # dark blue, red, orange, green; full opacity
                                   col.obs = c( "#0019653F", "#C81F493F", "#F582203F", "#00AF413F" ), # dark blue, red, orange, green; 50 % opacity
                                   col.ci = "#E0DED8FF", # light gray
                                   col.bl = "#82786FFF", # dark gray
                                   lwd.mean = 1,
                                   pch.obs = 4, # cross
                                   cex.obs = 1,
                                   lty.bl = 3 # dotted
    )
    
}

### TECHNICAL ANALYSIS



## NB: This appendix includes technical analysis that is not fully commented.


## Plot the fit of the model as a function of the number of clusters.

if ( FALSE ) { # Not run.
    
    opts.plot <- list()
    
    opts.plot$pch.optimal <- 8 # asterisk
    opts.plot$pch.other <- 4 # diagonal cross
    opts.plot$col.optimal <- "#009FDAFF" # light blue
    opts.plot$col.other <- "#0019653F" # dark blue with 50 % opacity
    
    plot( x = result.clustering.fit$repeats$K,
          y = -result.clustering.fit$repeats$energy,
          yaxt = "n",
          pch = opts.plot$pch.other,
          col = opts.plot$col.other,
          xlab = "Number of clusters",
          ylab = "",
          main = "Model fit v. the number of clusters"
    )
    
    points( x = result.clustering.fit$repeats$K[ which.min( result.clustering.fit$repeats$energy ) ],
            y = -result.clustering.fit$repeats$energy[ which.min( result.clustering.fit$repeats$energy ) ],
            pch = opts.plot$pch.optimal,
            col = opts.plot$col.optimal
    )
    
    title( ylab = "Variational lower bound for\nthe marginal likelihood of the observed data",
           line = 0.5
    )
    
    legend( x = "bottomright",
            legend = c( "Best model", "Other learned models" ),
            pch = c( opts.plot$pch.optimal, opts.plot$pch.other ),
            col = c( opts.plot$col.optimal, opts.plot$col.other )
    )
    
}
