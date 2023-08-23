# ##############################################################################
# Plots for provided genes lists
# ##############################################################################

#### MODULE SCORE

## @knitr heterogeneity_modules_scoring

MODULES_GENES = GENES_LIST

# Compute the score of the cells according to lists of genes
for( moduleName in names( MODULES_GENES))
{
  message(moduleName);
  if( length( MODULES_GENES[[moduleName]]) == 0)
  {
    warning( paste0( "List of genes in module '", moduleName, "' is empty, ignoring..."));
  } else
  {
    sc10x <- AddModuleScore( object = sc10x,
                             features = MODULES_GENES[ moduleName], # Must be a list
                             ctrl = MODULES_CONTROL_SIZE,           #length(MODULES_GENES[[ moduleName]]),
                             name = moduleName,
                             seed = SEED);
  }
}




## @knitr heterogeneity_modules_heatmap

# DoHeatmap replaced by use of iHeatmapr after testing several options
#htmltools::tagList( ggplotly( DoHeatmap(object = sc10x, features = topMarkersDF[["gene"]])));

# Get the matrix of module scores (not expression) for each cell and associated clusters from Seurat object
modulesScoreMat = t( as.matrix( sc10x[[ paste0(names( MODULES_GENES),1)]]));
clusterID = Idents( sc10x);

# Remove the extra numeric character added to modules names by Seurat
rownames( modulesScoreMat) = substr( rownames( modulesScoreMat), 1, nchar( rownames( modulesScoreMat))-1);

# Select genes in modules and reorder cells to group clusters together
clusterOrdering = order( clusterID);

modulesScoreMat = modulesScoreMat[, clusterOrdering];
clusterID = clusterID[ clusterOrdering];

# Prepare rows and columns annotation bars (module and cluster respectively)
#rowsAnnot = data.frame( Module = names( MODULES_GENES));
colsAnnot = data.frame( Cluster = clusterID);



# # Plot the 'non-interactive' heatmap
# pheatmap( modulesScoreMat,
#           cluster_rows = FALSE,
#           cluster_cols = FALSE,
#           #          annotation_row = rowsAnnot,
#           annotation_col = colsAnnot,
#           labels_row = originalRowNames,
#           annotation_colors = list( Cluster = clustersColor),
#           show_colnames = FALSE);

scaleRowsHeatmap = TRUE

Heatmap( if(scaleRowsHeatmap) t(scale(t(modulesScoreMat))) else modulesScoreMat, 
         name = "Module Score", 
         col = RColorBrewer::brewer.pal(name = "Reds", n = 4),
         rect_gp = gpar(col = "white", lwd = 0), # Border of cells
         column_names_rot = 0, #75
         show_column_names = FALSE,
         cluster_rows = FALSE,
         cluster_row_slices = FALSE, # Do not cluster rows split categories (based on logFC) to keep order of slices under control
         cluster_columns = TRUE,
         left_annotation = NULL,
         top_annotation = HeatmapAnnotation( df = colsAnnot,
                                             col = list("Cluster" = clustersColor)),
         # Make groups of columns based on categories in dataframe
         column_split = colsAnnot[["Cluster"]],
         column_title = character(0), # Default: character(0), Disable: NULL
         column_title_rot = 0,
         # Make groups of rows
         row_split = NULL,
         row_title = character(0), # Default: character(0), Disable: NULL
         row_title_rot = 0)



## @knitr heterogeneity_modules_expression_projection
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Plot scores of modules on dimensionality reduction figures
invisible( lapply( names(MODULES_GENES), function(featureName)
{
  print( FeaturePlot( sc10x, features = paste0(featureName, "1"), reduction = ifelse( exists("useReduction"), useReduction, "umap"), order = TRUE) +
           ggtitle( label = featureName) +
           theme( axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  legend.position = "none"));
}));

cat(" \n \n"); # Required for '.tabset'



## @knitr heterogeneity_modules_expression_violin

# Violinplot for module values by cluster
invisible( lapply( paste0(names(MODULES_GENES), 1), violinFeatureByCluster, seuratObject = sc10x, clustersColor = clustersColor, yLabel = "Score", addStats = FALSE, trimTitle = 1));

cat(" \n \n"); # Required for '.tabset'





## @knitr heterogeneity_modules_expression_violin_byHTO

# Violinplot for module values by cluster
invisible( lapply( paste0(names(MODULES_GENES), 1), violinFeatureByMetadata, seuratObject = sc10x, yLabel = "Score", addStats = FALSE, trimTitle = 1));

cat(" \n \n"); # Required for '.tabset'





## @knitr heterogeneity_modules_expression_projection_pngFile
# Plot expression values of individual module genes as png files (not in report)
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Compute 'pixel' dimensions to match figures in report (based on chunk params)
defaultDim = c( 800, 800);  # Fallback values if chunk params not available
defaultDpi = 72;            # (e.g. executing directly from R)

chunkDim    = knitr::opts_current$get( "fig.dim");
chunkDpi    = knitr::opts_current$get( "dpi");

figDim = if( is.null( chunkDim) || is.null( chunkDpi) ) defaultDim else chunkDim * chunkDpi;
figDpi = if( is.null( chunkDpi) ) defaultDpi else chunkDpi;


invisible( lapply( names( MODULES_GENES), function(moduleName)
{
  message(moduleName); # Just for tracking progress in console
  
  # Create subfolder for current module to store all png files
  pathCurrentModule = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "modulesExpressionIndividual_", ifelse(exists("useReduction"), useReduction, "umap")), moduleName);
  dir.create(pathCurrentModule, showWarnings = FALSE, recursive = TRUE);
  
  # Plots expression on projected cells (or error message if feature not found)
  invisible( lapply( MODULES_GENES[[moduleName]], function(featureName)
  {
    # Create png file
    png( file.path( pathCurrentModule, paste0( featureName, ".png")), 
         width = figDim[1],
         height = figDim[2],
         res = figDpi);
    
    print(FeaturePlot( sc10x, features = featureName, reduction = ifelse(exists("useReduction"), useReduction, "umap"), order = TRUE)  +
            theme( axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   legend.position = "none"))
    
    dev.off(); # Close file descriptor
  }));
  
}));




## @knitr heterogeneity_modules_expression_violin_pngFile
# Plot expression values of monitored genes as violinplot in a png file for each cluster (TODO: message if list empty)

# Compute 'pixel' dimensions to match figures in report (based on chunk params)
defaultDim = c( 800, 800);  # Fallback values if chunk params not available
defaultDpi = 72;            # (e.g. executing directly from R)

chunkDim    = knitr::opts_current$get( "fig.dim");
chunkDpi    = knitr::opts_current$get( "dpi");

figDim = if( is.null( chunkDim) || is.null( chunkDpi) ) defaultDim else chunkDim * chunkDpi;
figDpi = if( is.null( chunkDpi) ) defaultDpi else chunkDpi;


invisible( lapply( names( MODULES_GENES), function(moduleName)
{
  message(moduleName); # Just for tracking progress in console
  
  # Create subfolder for current module to store all png files
  pathCurrentModule = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "modulesExpressionIndividual_violin"), moduleName);
  dir.create(pathCurrentModule, showWarnings = FALSE, recursive = TRUE);
  
  # Violinplot for expression value of monitored genes by cluster (+ number of 'zero' and 'not zero' cells)
  invisible( lapply( MODULES_GENES[[moduleName]], function(featureName)
  {
    # Create png file
    png( file.path( pathCurrentModule, paste0( featureName, ".png")), 
         width = figDim[1],
         height = figDim[2],
         res = figDpi);
    
    violinFeatureByCluster(featureName, seuratObject = sc10x, clustersColor = clustersColor);
    
    dev.off(); # Close file descriptor
  }))
  
}));





## @knitr heterogeneity_aucell

BiocManager::install("AUCell")
library("AUCell")

MODULES_GENES = GENES_LIST

expMat = as.matrix(GetAssayData(sc10x, slot = "count"))

cellsByHTO = split( colnames(expMat), sc10x[["factorHTO"]])
expMatByHTO = lapply( cellsByHTO, function(x){ expMat[, x] })


# Calculate enrichment scores
auc_rankings_HTOlist = lapply(expMatByHTO, AUCell_buildRankings, GENES_LIST )#, aucMaxRank=nrow(cells_rankings)*0.05)
auc_HTOlist = lapply(auc_rankings_HTOlist,function(x) {AUCell_calcAUC(GENES_LIST, x)})


AUCell_calcAUC(GENES_LIST, auc_rankings_list[[1]])

# Optional: Set the assignment thresholds
par(mfrow=c(3,3))
set.seed(123)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)


