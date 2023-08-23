# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################




# PCA
# ---

## @knitr heterogeneity_pca

# Compute PCA on selected variable genes
nbPC=PCA_NPC
if(PCA_NPC>length(Cells(sc10x)))
{
  warning( paste0( "Number of cells in object (", length(Cells(sc10x)), ") smaller than requested number of PCs (", PCA_NPC,"), setting lower PC number..." ))
  nbPC = length(Cells(sc10x))
}           
sc10x <- RunPCA( object   = sc10x,
                 features = VariableFeatures( sc10x),
                 npcs     = nbPC,
                 verbose  = .VERBOSE);

# Identify clusters of cells by graph approach
nbPC_findclusters=FINDCLUSTERS_USE_PCA_NBDIMS
if(FINDCLUSTERS_USE_PCA_NBDIMS>nbPC)
{
  warning( paste0( "Number of computed PCs  (", nbPC, ") smaller than requested PCs for 'findclusters' (", FINDCLUSTERS_USE_PCA_NBDIMS,"), setting lower PC number (", nbPC, ")..." ))
  nbPC_findclusters = nbPC
}  
sc10x <- FindNeighbors(object    = sc10x,
                       reduction = "pca",
                       dims      = 1:nbPC_findclusters,
                       verbose   = .VERBOSE);

# Plot PCA, highlighting seurat clusters for combination of dimensions
invisible( apply( combn( 1:PCA_PLOTS_NBDIMS, 2), 2, function(dims)
{
  print( DimPlot( object=sc10x, reduction="pca", dims = dims, group.by = "orig.ident") +
           theme( legend.position = "none"));                               # Remove legend in this case...
}));




## @knitr heterogeneity_pca_umisCounts

# Same PCA plots but highlight UMI counts
invisible( apply( combn( 1:PCA_PLOTS_NBDIMS, 2), 2, function(dims)
{
  print( FeaturePlot( sc10x, feature = "nCount_RNA", reduction = "pca", dims = dims) +
           theme( legend.position = "none",
                  plot.margin = margin(0, 0, 0, 0, "cm"),
                  plot.title = element_blank()));
}));




## @knitr heterogeneity_pca_genesCounts

# Same PCA plots but highlight feature counts
invisible( apply( combn( 1:PCA_PLOTS_NBDIMS, 2), 2, function(dims)
{
  print( FeaturePlot( sc10x, feature = "nFeature_RNA", reduction = "pca", dims = dims) +
           theme( legend.position = "none",
                  plot.margin = margin(0, 0, 0, 0, "cm"),
                  plot.title = element_blank()));
}));




## @knitr heterogeneity_pca_correlations

# Isolate the PCA location of cells in first dimensions, together with UMI and genes counts for correlation analysis (makes use of cbind recycling to repeat values for each stacked PC)
relationToPC = suppressWarnings( cbind( stack( as.data.frame( Embeddings( sc10x, reduction = "pca")[ rownames( sc10x[[]]), paste0( "PC_", 1:PCA_PLOTS_NBDIMS) ])),
                                        sc10x[[ c( "nCount_RNA", "nFeature_RNA") ]],
                                        Cluster = Idents( sc10x)));

# Plot relationship of UMIs and genes counts with PCs (combine plots using '/' from 'patchwork' lib)
print( (ggplot( data = relationToPC, aes( x = values, y = nCount_RNA)) +
          facet_wrap( ~ind) +
          stat_binhex( bins = 60) +
          #geom_point( aes(col = Cluster), alpha = 0.5) +
          geom_smooth( method = 'lm') +
          stat_cor( method = "spearman") +
          ylab( "# UMIs") +
          theme( axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 plot.margin = unit( c( 1, 1, -0.5, 0.5), "lines")))
       /
         (ggplot( data = relationToPC, aes( x = values, y = nFeature_RNA)) +
            facet_wrap( ~ind) +
            stat_binhex( bins = 60) +
            #geom_point( aes(col = Cluster), alpha = 0.5) +
            geom_smooth( method = 'lm') +
            stat_cor( method = "spearman") +
            xlab( "PC values") +
            ylab( "# Genes")));




## @knitr heterogeneity_pca_loadings

# Plot PCA loadings
invisible( apply( combn( 1:PCA_PLOTS_NBDIMS, 2), 2, function(dims)
{
  namesPC=paste0( "PC_", dims);
  # Get the loading values for concerned PCs
  loadingsMatrix = Loadings( sc10x, reduction = "pca")[ , namesPC ];
  # Sort features by average absolute value and convert as DF with features names as column
  loadingsMatrix = head( loadingsMatrix[ order( apply( loadingsMatrix, 1, function(x){ mean( abs( x)) }), decreasing = TRUE), ], PCA_PLOTS_NBFEATURES);
  loadingsDF = data.frame( loadingsMatrix, features = rownames( loadingsMatrix));
  
  # Define symmetric and consistent axes for group of plots
  axesLimit = max( abs( loadingsMatrix));
  
  # Plot arrows and features name
  print( ggplot( data = loadingsDF, aes_( x = as.name( namesPC[1]), y = as.name( namesPC[2]))) +
           coord_cartesian( xlim = c( -axesLimit, axesLimit), ylim = c( -axesLimit, axesLimit)) +
           geom_text_repel( aes( label = features), max.iter = 10000) +
           geom_segment( x = 0 , y = 0, aes_( xend = as.name(namesPC[1]), yend = as.name(namesPC[2])), col = "#00000044", arrow = arrow( length = unit( 2, "mm"))));
}));




# DIMENSIONAL REDUCTION (TSNE/UMAP)
###################################

## @knitr heterogeneity_dimReduc
nbPC_dimreduc=DIMREDUC_USE_PCA_NBDIMS
if(DIMREDUC_USE_PCA_NBDIMS>nbPC)
{
  warning( paste0( "Number of computed PCs  (", nbPC, ") smaller than requested PCs for 'dimreduc' (", DIMREDUC_USE_PCA_NBDIMS,"), setting lower PC number (", nbPC, ")..." ))
  nbPC_dimreduc = nbPC
}

sc10x = RunUMAP( sc10x, dims = 1:nbPC_dimreduc);
sc10x = RunTSNE( sc10x, dims = 1:nbPC_dimreduc);

# Save resulting coordinates for all cells as 'tsv' files
write.table( Embeddings(sc10x, reduction = "umap"), 
             file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "cellsCoordinates_umap.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");

write.table( Embeddings(sc10x, reduction = "tsne"), 
             file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "cellsCoordinates_tsne.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");





## @knitr dimreduc_ggplotly_overlayBatches
# Interactive plot using plotly
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)


# Gather data to be visualized together (cell name + numID + metrics)
cellsData = cbind( "Cell" = colnames( sc10x), # cell names from rownames conflict with search field in datatable
                   sc10x[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                   "percent.mito" = if(length( mito.genes)) as.numeric(format(sc10x[["percent.mito", drop = TRUE]], digits = 5)) else NULL,
                   "percent.ribo" = if(length( ribo.genes)) as.numeric(format(sc10x[["percent.ribo", drop = TRUE]], digits = 5)) else NULL,
                   "Batch" = sc10x[["orig.ident", drop = TRUE]],
                   "HTO" = sc10x[["HTO_classification", drop = TRUE]]);

# Create text to show under cursor for each cell
hoverText = do.call(paste, c(Map( paste,
                                  c( "",
                                     "Cell ID: ",
                                     "# UMIs: ",
                                     "# Genes: ",
                                     if(length( mito.genes)) "% Mito: ",
                                     if(length( ribo.genes)) "% Ribo: "),
                                  cellsData[-ncol( cellsData)],
                                  sep = ""),
                             sep = "\n"));

## Do it in pure plotly
dimReducData = cbind(as.data.frame(Embeddings(sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"))),
                     Cluster = Idents( sc10x),
                     Batch   = sc10x[["orig.ident", drop = TRUE]],
                     HTO     = sc10x[["HTO_classification", drop = TRUE]]);

# Reusing 'hoverText' summarizing cell info (generated for split violin plots)
plotDimReduc = plot_ly(dimReducData,
                       x = as.formula( paste( "~", colnames( dimReducData)[1])),
                       y = as.formula( paste( "~", colnames( dimReducData)[2])),
                       #                       name = ~paste("Cluster", Cluster),
                       color = ~Batch,
                       height = 800) %>%
  add_trace( type = "scattergl",
             mode = "markers",
             marker = list( size = 5,
                            opacity = 0.6),
             text = hoverText,
             hoverinfo = "text+name") %>%
  #  hide_colorbar() %>%
  #  hide_legend() %>%
  config( displaylogo = FALSE,
          toImageButtonOptions = list(format='svg'),
          modeBarButtons = list(list('toImage'),
                                list('zoom2d', 'zoomIn2d', 'zoomOut2d', 'pan2d', 'resetScale2d')));

# Control layout using flex
# 'browsable' required in console, not in script/document
#browsable(
div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
     div(plotDimReduc, style = paste("flex : 0 0 800px;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;")))
#)





## @knitr dimreduc_ggplotly_overlayHTOs
# Interactive plot using plotly
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)


# Gather data to be visualized together (cell name + numID + metrics)
cellsData = cbind( "Cell" = colnames( sc10x), # cell names from rownames conflict with search field in datatable
                   sc10x[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                   "percent.mito" = if(length( mito.genes)) as.numeric(format(sc10x[["percent.mito", drop = TRUE]], digits = 5)) else NULL,
                   "percent.ribo" = if(length( ribo.genes)) as.numeric(format(sc10x[["percent.ribo", drop = TRUE]], digits = 5)) else NULL,
                   "Batch" = sc10x[["orig.ident", drop = TRUE]],
                   "HTO" = sc10x[["HTO_classification", drop = TRUE]]);

# Create text to show under cursor for each cell
hoverText = do.call(paste, c(Map( paste,
                                  c( "",
                                     "Cell ID: ",
                                     "# UMIs: ",
                                     "# Genes: ",
                                     if(length( mito.genes)) "% Mito: ",
                                     if(length( ribo.genes)) "% Ribo: "),
                                  cellsData[-ncol( cellsData)],
                                  sep = ""),
                             sep = "\n"));

## Do it in pure plotly
dimReducData = cbind(as.data.frame(Embeddings(sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"))),
                     Cluster = Idents( sc10x),
                     Batch   = sc10x[["orig.ident", drop = TRUE]],
                     HTO     = sc10x[["HTO_classification", drop = TRUE]]);

# Reusing 'hoverText' summarizing cell info (generated for split violin plots)
plotDimReduc = plot_ly(dimReducData,
                       x = as.formula( paste( "~", colnames( dimReducData)[1])),
                       y = as.formula( paste( "~", colnames( dimReducData)[2])),
                       #                       name = ~paste("Cluster", Cluster),
                       color = ~HTO,
                       height = 800) %>%
  add_trace( type = "scattergl",
             mode = "markers",
             marker = list( size = 5,
                            opacity = 0.6),
             text = hoverText,
             hoverinfo = "text+name") %>%
  #  hide_colorbar() %>%
  #  hide_legend() %>%
  config( displaylogo = FALSE,
          toImageButtonOptions = list(format='svg'),
          modeBarButtons = list(list('toImage'),
                                list('zoom2d', 'zoomIn2d', 'zoomOut2d', 'pan2d', 'resetScale2d')));

# Control layout using flex
# 'browsable' required in console, not in script/document
#browsable(
div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
     div(plotDimReduc, style = paste("flex : 0 0 800px;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;")))
#)






## @knitr compareBatchHTO_table

kable(table(sc10x[[c("orig.ident", "HTO_classification")]]))



# IDENTIFY CLUSTERS
###################

## @knitr heterogeneity_identifyClusters

sc10x <- FindClusters(object             = sc10x,
                      resolution         = FINDCLUSTERS_RESOLUTION,
                      algorithm          = FINDCLUSTERS_ALGORITHM,
                      temp.file.location = "/tmp/",
                      verbose            = .VERBOSE);

#DimPlot(sc10x)

# Show number of cells in each cluster
clustersCount = as.data.frame( table( Cluster = sc10x[[ "seurat_clusters" ]]), responseName = "CellCount");

# Define a set of colors for clusters (based on ggplot default)
clustersColor = hue_pal()( nlevels( Idents( sc10x)));
names( clustersColor) = levels( Idents( sc10x));


# Save cells cluster identity as determined with 'FindClusters'
write.table( data.frame(sc10x[["numID"]], identity = Idents(sc10x)), 
             file = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "cellsClusterIdentity.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");

# Also save cluster color attribution for reference
# Save cells cluster identity as determined with 'FindClusters'
write.table( clustersColor, 
             file = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "clustersColor.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = FALSE,
             sep="\t");


# Create datatable
datatable( clustersCount,
           class = "compact",
           rownames = FALSE,
           colnames = c("Cluster", "Nb. Cells"),
           options = list(dom = "<'row'rt>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                          autoWidth = FALSE,
                          columnDefs = list( # Center all columns
                                            list( targets = 0:(ncol(clustersCount)-1),
                                            className = 'dt-center')),
                          orderClasses = FALSE, # Disable flag for CSS to highlight columns used for ordering (for performance)
                          paging = FALSE, # Disable pagination (show all)
                          processing = TRUE,
                          scrollCollapse = TRUE,
                          scroller = TRUE,  # Only load visible data
                          scrollX = TRUE,
                          scrollY = "525px",
                          stateSave = TRUE)) %>%
  # Add color from cluster
  formatStyle( columns = "Cluster",
               color = styleEqual( names(clustersColor), clustersColor),
               fontWeight = 'bold');


# Compute a matrix of average expression value by cluster (for each gene)
geneExpByCluster = do.call( rbind, 
                            apply( as.matrix( GetAssayData( sc10x)), # expression values
                                   1,                                # by rows
                                   tapply,                           # apply by group
                                   INDEX = Idents( sc10x),           # clusters IDs
                                   mean,                             # summary function
                                   simplify = FALSE));               # do not auto coerce

# Save it as 'tsv' file
write.table( geneExpByCluster, 
             file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "normExpressionByCluster.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");




## @knitr heterogeneity_identifyClusters_splitStats
# Show UMIs Genes Mitochondrial and Ribosomal content split by cluster

# Gather data to be visualized together (cell name + numID + metrics)
cellsData = cbind( "Cell" = colnames( sc10x), # cell names from rownames conflict with search field in datatable
                   sc10x[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                   "percent.mito" = if(length( mito.genes)) as.numeric(format(sc10x[["percent.mito", drop = TRUE]], digits = 5)) else NULL,
                   "percent.ribo" = if(length( ribo.genes)) as.numeric(format(sc10x[["percent.ribo", drop = TRUE]], digits = 5)) else NULL,
                   "Batch" = sc10x[["orig.ident"]],
                   "HTO" = sc10x[["HTO_classification"]],
                   "Cluster" = Idents( sc10x));

# Create text to show under cursor for each cell
hoverText = do.call(paste, c(Map( paste,
                                  c( "",
                                     "Cell ID: ",
                                     "# UMIs: ",
                                     "# Genes: ",
                                     if(length( mito.genes)) "% Mito: ",
                                     if(length( ribo.genes)) "% Ribo: "),
                                  cellsData[-ncol( cellsData)],
                                  sep = ""),
                    sep = "\n"));

# Define size for panels (or assembled figure when using subplot)
panelWidth = 90 * nlevels(cellsData[["Cluster"]]);
panelHeight = 800;

# Generate plotly violin/jitter panels for #umis, #genes, and %mitochondrial stats
lypanel_umis  = plotViolinJitter( cellsData,
                                  xAxisFormula = ~as.numeric( Cluster),
                                  yAxisFormula = ~nCount_RNA,
                                  colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                  yAxisTitle = "# UMIs",
                                  hoverText = hoverText,
                                  traceName = ~paste( "Cluster", Cluster),
                                  xTicklabels = levels(cellsData[["Cluster"]]),
                                  panelWidth = panelWidth,
                                  panelHeight = panelHeight);

lypanel_genes = plotViolinJitter( cellsData,
                                  xAxisFormula = ~as.numeric( Cluster),
                                  yAxisFormula = ~nFeature_RNA,
                                  colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                  yAxisTitle = "# Genes",
                                  hoverText = hoverText,
                                  traceName = ~paste( "Cluster", Cluster),
                                  xTicklabels = levels(cellsData[["Cluster"]]),
                                  panelWidth = panelWidth,
                                  panelHeight = panelHeight);

lypanel_mitos = if(length( mito.genes)) plotViolinJitter( cellsData,
                                                          xAxisFormula = ~as.numeric( Cluster),
                                                          yAxisFormula = ~percent.mito,
                                                          colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                                          yAxisTitle = "% Mito",
                                                          hoverText = hoverText,
                                                          traceName = ~paste( "Cluster", Cluster),
                                                          xTicklabels = levels(cellsData[["Cluster"]]),
                                                          panelWidth = panelWidth,
                                                          panelHeight = panelHeight) else NULL;

lypanel_ribos = if(length( ribo.genes)) plotViolinJitter( cellsData,
                                                          xAxisFormula = ~as.numeric( Cluster),
                                                          yAxisFormula = ~percent.ribo,
                                                          colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                                          yAxisTitle = "% Ribo",
                                                          hoverText = hoverText,
                                                          traceName = ~paste( "Cluster", Cluster),
                                                          xTicklabels = levels(cellsData[["Cluster"]]),
                                                          panelWidth = panelWidth,
                                                          panelHeight = panelHeight) else NULL;

# Set panels as a list, define plotly config, and remove legend
panelsList = list( lypanel_umis, lypanel_genes, lypanel_mitos, lypanel_ribos);
panelsList = lapply( panelsList, config, displaylogo = FALSE,
                     toImageButtonOptions = list( format='svg'),
                     modeBarButtons = list( list('toImage'),
                                            list( 'zoom2d', 'pan2d', 'resetScale2d')));

# Group plotly violin/jitter panels so we can synchronise axes and use highlight on the full set
plotPanels = layout( subplot( panelsList,
                              nrows = 4,
                              shareX = TRUE,
                              titleY = TRUE),
                     xaxis = list(title = "Seurat Cluster",
                                  showgrid = TRUE,
                                  tickvals = seq(nlevels(cellsData[["Cluster"]])),
                                  ticktext = levels(cellsData[["Cluster"]])),
                     showlegend = FALSE, # Remove eventual legends (does not mix well with subplot and highlight)
                     autosize = TRUE);

# Control layout using flex
# 'browsable' required in console, not in script/document
#browsable(
div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: green;"),
     div( style = paste("display: flex; flex-flow: row nowrap; overflow-x: auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
          div(plotPanels, style = paste("flex : 0 0 auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;"))))
#)




## @knitr heterogeneity_dimReduc_interactivePlot
# Interactive plot using plotly
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# # Using Seurat::Dimplot does not allow to control tooltip properly (neither using HoverLocator on ggplot, or plotly tooltip + added aestetic in ggplot). Mostly because 'label' and 'combine' change the type of returned object)
# # Compute tSNE dimensional reductionsand plot projection
# plotDimReduc = ggplotly( suppressWarnings( DimPlot( sc10x,
#                                                     reduction = ,
#                                                     label=TRUE,
#                                                     label.size = 6)) +
#                                              theme( legend.position = "none",
#                                                     plot.margin = margin(0, 0, 0, 0, "cm")),
#                          height = 600) %>%
#   config( displaylogo = FALSE,
#           toImageButtonOptions = list(format='svg'),
#           modeBarButtons = list(list('toImage'),
#                                 list('zoom2d', 'pan2d', 'resetScale2d')));

## Do it in pure plotly
dimReducData = cbind(as.data.frame(Embeddings(sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"))),
                     Cluster = Idents( sc10x));

# Reusing 'hoverText' summarizing cell info (generated for split violin plots)
plotDimReduc = plot_ly(dimReducData,
                       x = as.formula( paste( "~", colnames( dimReducData)[1])),
                       y = as.formula( paste( "~", colnames( dimReducData)[2])),
                       #                       name = ~paste("Cluster", Cluster),
                       color = ~Cluster,
                       colors = clustersColor,
                       height = 800) %>%
  add_trace( type = "scattergl",
             mode = "markers",
             marker = list( size = 5,
                            opacity = 0.6),
             text = hoverText,
             hoverinfo = "text+name") %>%
  # Compute median of coordinates for each group and use it as new data for Cluster text annotations (type='scatter', mode='text')
  add_trace( data = as.data.frame( do.call( rbind, by( dimReducData, dimReducData[["Cluster"]], function(x){ return( data.frame( lapply( x[1:2], median), "Cluster" = x[1, "Cluster"]));}))),
             x = as.formula( paste( "~", colnames( dimReducData)[1])),
             y = as.formula( paste( "~", colnames( dimReducData)[2])),
             type = "scattergl",
             mode = "text",
             text = ~Cluster,
             textfont = list( color = '#000000', size = 20),
             hoverinfo = 'skip',
             showlegend = FALSE) %>%
  #  hide_colorbar() %>%
  #  hide_legend() %>%
  config( displaylogo = FALSE,
          toImageButtonOptions = list(format='svg'),
          modeBarButtons = list(list('toImage'),
                                list('zoom2d', 'zoomIn2d', 'zoomOut2d', 'pan2d', 'resetScale2d')));

# Control layout using flex
# 'browsable' required in console, not in script/document
#browsable(
div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
     div(plotDimReduc, style = paste("flex : 0 0 800px;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;")))
#)





## @knitr rmd_heterogeneity_dimReduc_plot_byHTO

#### Dimensionality reduction 2D (ggplot + plotly conversion)
dimReducWidth  = 800;
dimReducHeight = 800;
initialPointSize = 5;


# Combine Dimreduc data and cells informations (keep coordinates in first cols)
# + add batch and HTO information for facetting
dimReducData = cbind(as.data.frame(Embeddings(sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"))),
                     sc10x[[c("orig.ident", "HTO_classification", "seurat_clusters")]]);


# Plot a thumbnail highlighting cluster cells
ggFigure = ggplot(  dimReducData,
                    aes_string( x = colnames( dimReducData)[1],
                                y = colnames( dimReducData)[2],
                                color = "seurat_clusters")) +
  geom_point( stroke = 0) + # Just ignore point size and alpha here as it will be overriden after plotly conversion
  facet_wrap(~HTO_classification)

ggFigure = ggFigure +
  theme_classic() +
  theme( #axis.text  = element_blank(),
    #axis.title = element_blank(),
    #axis.ticks  = element_blank(),
    #axis.line = element_blank(),
    #legend.position = "none",
    plot.margin = margin( 0, 0, 0, 0, "cm"),
    plot.title = element_text( face = "bold",
                               size = rel( 1.8), #rel( 16/14),
                               hjust = 0.5,
                               vjust = 1,
                               margin = margin( b = 7)));

plotDimReducGG =  ggplotly( ggFigure, 
                            width = dimReducWidth, 
                            height = dimReducHeight) %>% 
  style( type = "scattergl",   # trick to use webGL for rendering the scatterplot since toWebGL does not have width/height arguments, and use of width/height in layout is deprecated
         marker.size = initialPointSize, # Set point size and opacity with plotly so we can synchronize them with slider (ggplot and plotly units don't match)
         marker.opacity = 0.3) %>% 
  config( displaylogo = FALSE,
          toImageButtonOptions = list( format='svg'),
          modeBarButtons = list(list( 'toImage'),
                                list( 'zoom2d', 
                                      'zoomIn2d', 
                                      'zoomOut2d', 
                                      'pan2d', 
                                      'resetScale2d')))

cat("\n")
plotDimReducGG
cat("\n")



## @knitr compareClustersHTO_table

kable(table(sc10x[[c("seurat_clusters", "HTO_classification")]]))



## @knitr heterogeneity_dimReduc_thumbnail
# Non-interactive with large labels for generating thumbnails when analyzing monitored and marker genes expression
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)
# If contains several values, loop on them so we don't have to call chunk separately (which creates a new <p>aragraph in html)

# Replot the projection with colored clusters and large labels (add title if several plotted in the loop)
reductionVector = if(exists("useReduction")) useReduction else "umap";
for(currentReduction in reductionVector)
{
  ggFigure = suppressWarnings( DimPlot( sc10x, reduction = currentReduction, label=TRUE, label.size = 10)) +
                                 scale_color_manual( name = "Cluster", values = clustersColor) +  # Fix legend title + set colors
                                 theme( axis.title.x = element_blank(),
                                        axis.title.y = element_blank(),
                                        legend.position = "none",
                                        plot.margin = margin( 0, 0, 0, 0, "cm"),
                                        plot.title = element_text( face = "bold",
                                                                   size = rel( 16/14),
                                                                   hjust = 0.5,
                                                                   vjust = 1,
                                                                   margin = margin( b = 7)));

  if( length( reductionVector) > 1) ggFigure = ggFigure + ggtitle( label = currentReduction);
  print( ggFigure);
}




# MARKER GENES
##############

## @knitr heterogeneity_markerGenes

# Identify marker genes
markers = FindAllMarkers( object          = sc10x,
                          test.use        = FINDMARKERS_METHOD,
                          only.pos        = FINDMARKERS_ONLYPOS,
                          min.pct         = FINDMARKERS_MINPCT,
                          logfc.threshold = FINDMARKERS_LOGFC_THR,
                          return.thresh   = FINDMARKERS_PVAL_THR,
                          random.seed     = SEED,
                          verbose         = .VERBOSE);

# Save markers list as 'tsv' table
write.table( markers,
             file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "MarkerGenes.tsv")),
             quote = FALSE,
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");

# Filter markers by cluster (TODO: check if downstream code works when no markers found)
topMarkers = by( markers, markers[["cluster"]], function(x)
{
  # Filter markers based on adjusted PValue
  x = x[ x[["p_val_adj"]] < FINDMARKERS_PVAL_THR, , drop = FALSE];
  # Sort by decreasing logFC
  x = x[ order(abs(x[["avg_log2FC"]]), decreasing = TRUE), , drop = FALSE ]
  # Return top ones
  return( if(is.null( FINDMARKERS_SHOWTOP)) x else head( x, n = FINDMARKERS_SHOWTOP));
});

# Merge marker genes in a single data.frame and render it as datatable
topMarkersDF = do.call( rbind, topMarkers);
# Select and order columns to be shown in datatable
topMarkersDT = topMarkersDF[c("gene", "cluster", "avg_log2FC", "p_val_adj")]




## @knitr heterogeneity_markerGenes_table

# Create datatable
datatable( topMarkersDT,
           class = "compact",
           filter="top",
           rownames = FALSE,
           colnames = c("Gene", "Cluster", "Avg. LogFC", "Adj. Pvalue"),
           caption = paste(ifelse( is.null( FINDMARKERS_SHOWTOP), "All", paste("Top", FINDMARKERS_SHOWTOP)), "marker genes for each cluster"),
           extensions = c('Buttons', 'Select'),
           options = list(dom = "<'row'<'col-sm-8'B><'col-sm-4'f>> <'row'<'col-sm-12'l>> <'row'<'col-sm-12'rt>> <'row'<'col-sm-12'ip>>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                          autoWidth = FALSE,
                          buttons = exportButtonsListDT,
                          columnDefs = list(
                            list( # Center all columns except first one
                              targets = 1:(ncol( topMarkersDT)-1),
                              className = 'dt-center'),
                            list( # Set renderer function for 'float' type columns (LogFC)
                              targets = ncol( topMarkersDT)-2,
                              render = htmlwidgets::JS("function ( data, type, row ) {return type === 'export' ? data : data.toFixed(4);}")),
                            list( # Set renderer function for 'scientific' type columns (PValue)
                              targets = ncol( topMarkersDT)-1,
                              render = htmlwidgets::JS( "function ( data, type, row ) {return type === 'export' ? data : data.toExponential(4);}"))),
                          #fixedColumns = TRUE, # Does not work well with filter on this column
                          #fixedHeader = TRUE, # Does not work well with 'scrollX'
                          lengthMenu = list(c( 10, 50, 100, -1),
                                            c( 10, 50, 100, "All")),
                          orderClasses = FALSE, # Disable flag for CSS to highlight columns used for ordering (for performance)
                          processing = TRUE,
                          #search.regex= TRUE, # Does not work well with 'search.smart'
                          search.smart = TRUE,
                          select = TRUE, # Enable ability to select rows
                          scrollCollapse = TRUE,
                          scroller = TRUE,  # Only load visible data
                          scrollX = TRUE,
                          scrollY = "525px",
                          stateSave = TRUE)) %>%
  # Add bar relative to logFC
  formatStyle( columns = "avg_log2FC",
               background = styleColorBar( data = range( topMarkersDT[["avg_log2FC"]]), 'lightblue', angle = -90),
               backgroundSize = '95% 50%',      # Set horizontal and vertical span in cell
               backgroundRepeat = 'no-repeat',
               backgroundPosition = 'center') %>%
  # Add color from cluster
  formatStyle( columns = "cluster",
               backgroundColor = styleEqual( names(clustersColor),
                                             scales::alpha(clustersColor, 0.3)));




## @knitr heterogeneity_markerGenes_heatmap

# DoHeatmap replaced by use of iHeatmapr after testing several options
#htmltools::tagList( ggplotly( DoHeatmap(object = sc10x, features = topMarkersDF[["gene"]])));

# Get the matrix of expression and associated clusters from Seurat object
expMat = as.matrix( GetAssayData( sc10x));
clusterID = Idents( sc10x);

# Select marker genes and reorder cells to group clusters together
topMarkersGenes = topMarkersDF[["gene"]];
clusterOrdering = order( clusterID);

expMat = expMat[topMarkersGenes, clusterOrdering];
clusterID = clusterID[clusterOrdering];

# Prepare rows and columns annotation bars (cluster groups)
rowsAnnot = data.frame(Markers = topMarkersDF[["cluster"]]);
colsAnnot = data.frame(Cluster = clusterID);


# Only create interactive heatmap if reasonable number of cells (prevents crashing when rendering and slow html result)
if(length( Cells( sc10x)) < PLOT_RASTER_NBCELLS_THRESHOLD)
{
  # Create a temporary file that will store the code for interactive heatmap and
  # knitted as separate doc to prevent inflating report with demanding figures
  tempFileName = tempfile( fileext = ".Rmd");
  writeLines( con = tempFileName, text = '
---
title: Marker genes
---

```{r heterogeneity_markerGenes_heatmapInteractive, echo = FALSE}
    # Create a matrix for text on mouse over events (customize "text" that is supposed to show value for adding info)
    hoverTextMarkers = expMat;
    #hoverTextMarkers[] = paste( colnames( expMat)[col( expMat)], rownames( expMat)[row( expMat)], sep = "<br><br><br>");
    hoverTextMarkers[] = paste( paste( "Cell ID:", sc10x[["numID", drop = FALSE]][colnames( expMat)[col( expMat)],]),
                                paste( "Value:", expMat),
                                paste( "Cluster:", sort( clusterID)[col( expMat)]),
                                paste( "Markers cl.:", topMarkersDF[["cluster"]][row( expMat)]),
                                sep = "<br>");


    # Create heatmap
    heatmapPlot = iheatmap(expMat,
                           row_labels = TRUE,
                           col_labels = FALSE,
                           text = hoverTextMarkers,
                           tooltip = setup_tooltip_options(row = TRUE,
                                                           col = TRUE,
                                                           value = TRUE,
                                                           prepend_row = "Gene: ",
                                                           prepend_col = "Cell: ",
                                                           prepend_value = ""),
                           row_annotation = rowsAnnot,
                           col_annotation = colsAnnot,
                           row_annotation_colors = list(Markers = clustersColor),
                           col_annotation_colors = list(Cluster = clustersColor),
                           layout = list(width = 900, height = 900));

    # Hack to customize modebar buttons as for plotly native objects
    # see: https://github.com/ropensci/iheatmapr/issues/38#issuecomment-445428166
    heatmapPlot_widget = iheatmapr::to_widget(heatmapPlot);
    # Edit the htmlwidget object itself
    heatmapPlot_widget[["x"]][["config"]][["displaylogo"]] = FALSE;
    heatmapPlot_widget[["x"]][["config"]][["toImageButtonOptions"]] = list( format="svg");  # This one does not seem to work... TODO: Check issue on GitHub
    heatmapPlot_widget[["x"]][["config"]][["modeBarButtons"]]       = list( list( "toImage"),
                                                                            list( "zoom2d", "zoomIn2d", "zoomOut2d", "pan2d", "resetScale2d"));

    # Render the heatmap
    heatmapPlot_widget;
```
');

  # Create a file name for the document knitted separately
  interactiveHeatmap_filename = paste0( outputFilesPrefix, "markerGenes_interactiveHeatmap.html");

  # knit the interactive heatmap in a separate document using current environment
  resultFile = rmarkdown::render( tempFileName,
                                  output_file = I( interactiveHeatmap_filename),
                                  output_dir = PATH_ANALYSIS_OUTPUT,
                                  envir=environment(),
                                  clean = TRUE,
                                  quiet = TRUE);

  cat( '\n<a target="_blank" href="', interactiveHeatmap_filename, '">Click here for interactive heatmap in a separate window</a>\n', sep = "");
}


# Prepare unique rows and cols names for pheatmap (annotation rows) and match with rowAnnots and colAnnots row names
originalRowNames = rownames( expMat);
originalColNames = colnames( expMat);
rownames( expMat) = make.unique( originalRowNames);
colnames( expMat) = make.unique( originalColNames);
rownames( rowsAnnot) = rownames( expMat);
rownames( colsAnnot) = colnames( expMat);

pheatmap( expMat,
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          annotation_row = rowsAnnot,
          annotation_col = colsAnnot,
          labels_row = originalRowNames,
          annotation_colors = list( Markers = clustersColor,
                                    Cluster = clustersColor),
          show_colnames = FALSE);




## @knitr heterogeneity_markerGenes_expression_projection
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Plot expression values of marker genes on dimreduc figures for each cluster (TODO: message if list empty)
invisible( lapply( names( topMarkers), function(clusterName)
{
  cat("##### Cl. <span style='border-radius: 3px; border: 3px solid ", clustersColor[clusterName], "; padding:0px 2px'>", clusterName, "</span>\n");

  # Highlight cells of current cluster on a dimreduc plot
  highlightClusterPlot(clusterName, seuratObject = sc10x, reduction = ifelse( exists("useReduction"), useReduction, "umap"));

  # Plots expression on projected cells
  invisible( lapply( topMarkers[[clusterName]][["gene"]], function(featureName)
    {
      print( FeaturePlot( sc10x, features = featureName, reduction = ifelse( exists("useReduction"), useReduction, "umap"), order = TRUE) +
               theme( axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      legend.position = "none"));
    }));

  cat(" \n \n"); # Required for '.tabset'
}));




## @knitr heterogeneity_markerGenes_expression_violin

# Plot expression values of marker genes as violinplot for each cluster (TODO: message if list empty)
invisible( lapply( names( topMarkers), function(clusterName)
{
  cat("##### Cl. <span style='border-radius: 3px; border: 3px solid ", clustersColor[clusterName], "; padding:0px 2px'>", clusterName, "</span>\n");

  # Remind cluster name in an empty figure to keep consistent alignment of panels between tabs
  plot( c( 0, 1), c( 0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n');
  text( x = 0.5, y = 0.5, paste( "Cluster", clusterName), cex = 2, col = clustersColor[clusterName]);

  # Violinplot for expression value of marker genes by cluster (+ number of 'zero' and 'not zero' cells)
  invisible( lapply( topMarkers[[clusterName]][["gene"]], violinFeatureByCluster, seuratObject = sc10x, clustersColor = clustersColor));

  cat(" \n \n"); # Required for '.tabset'
}));




# MONITORED GENES
#################

## @knitr heterogeneity_monitoredGenes

# Just remind the warning for genes names not in object
if(any( is.na( matchMonitoredGenes))) warning( paste( "Following gene(s) to be monitored could not be found in experimental data and will be ignored:", paste( monitoredGenesNotFound, collapse=" - ")));




## @knitr heterogeneity_monitoredGenes_heatmap

# DoHeatmap replaced by use of iHeatmapr after testing several options
#htmltools::tagList( ggplotly( DoHeatmap(object = sc10x, features = unlist(MONITORED_GENES))));

# Get the matrix of expression and associated clusters from Seurat object
expMat = as.matrix( GetAssayData( sc10x));
clusterID = Idents( sc10x);

# Select monitored genes and reorder cells to group clusters together
monitoredGenes = unlist( MONITORED_GENES);
clusterOrdering = order( clusterID);

expMat = expMat[monitoredGenes, clusterOrdering];
clusterID = clusterID[ clusterOrdering];

# Prepare rows and columns annotation bars (monitored group and cluster respectively)
rowsAnnot = data.frame( Monitored = fct_rev(fct_inorder(rep( names( MONITORED_GENES), sapply( MONITORED_GENES, length)))));
colsAnnot = data.frame( Cluster = clusterID);

# Only create interactive heatmap if reasonable number of cells (prevents crashing when rendering and slow html result)
if(length( Cells( sc10x)) < PLOT_RASTER_NBCELLS_THRESHOLD)
{
  # Create a temporary file that will store the code for interactive heatmap and
  # knitted as separate doc to prevent inflating report with demanding figures
  tempFileName = tempfile( fileext = ".Rmd");
  writeLines( con = tempFileName, text = '
---
title: Monitored genes
---

```{r heterogeneity_monitoredGenes_heatmapInteractive, echo = FALSE}
    # Create a matrix for text on mouse over events (customize "text" that is supposed to show value for adding info)
    hoverTextMarkers = expMat;
    hoverTextMarkers[] = paste( paste( "Cell ID:", sc10x[["numID", drop = FALSE]][colnames( expMat)[col( expMat)],]),
                                paste( "Value:", expMat),
                                paste( "Cluster:", sort( clusterID)[col( expMat)]),
                                paste( "Monitored Grp.:", rep(names(MONITORED_GENES), sapply(MONITORED_GENES, length))[row( expMat)]),
                                sep = "<br>");

    # Create heatmap
    heatmapPlot = iheatmap( expMat,
                            colorbar_position = 1,
                            row_labels = TRUE,
                            col_labels = FALSE,
                            text = hoverTextMarkers,
                            tooltip = setup_tooltip_options(row = TRUE,
                                                            col = TRUE,
                                                            value = TRUE,
                                                            prepend_row = "Gene: ",
                                                            prepend_col = "Cell: ",
                                                            prepend_value = ""),
                            row_annotation = rowsAnnot,
                            row_annotation_colors = list( Monitored = rainbow( nlevels( rowsAnnot[["Monitored"]]), s = 0.8)),
                            col_annotation = colsAnnot,
                            col_annotation_colors = list( Cluster = clustersColor),
                            layout = list(width = 900, height = 900));

    # Hack to customize modebar buttons as for plotly native objects
    # see: https://github.com/ropensci/iheatmapr/issues/38#issuecomment-445428166
    heatmapPlot_widget = iheatmapr::to_widget(heatmapPlot);
    # Edit the htmlwidget object itself
    heatmapPlot_widget[["x"]][["config"]][["displaylogo"]] = FALSE;
    heatmapPlot_widget[["x"]][["config"]][["toImageButtonOptions"]] = list( format="svg");  # This one does not seem to work... TODO: Check issue on GitHub
    heatmapPlot_widget[["x"]][["config"]][["modeBarButtons"]]       = list( list( "toImage"),
                                                                            list( "zoom2d", "zoomIn2d", "zoomOut2d", "pan2d", "resetScale2d"));

    # Render the heatmap
    heatmapPlot_widget;
```
');

  # Create a file name for the document knitted separately
  interactiveHeatmap_filename = paste0( outputFilesPrefix, "monitoredGenes_interactiveHeatmap.html");

  # knit the interactive heatmap in a separate document using current environment
  resultFile = rmarkdown::render( tempFileName,
                                  output_file = I( interactiveHeatmap_filename),
                                  output_dir = PATH_ANALYSIS_OUTPUT,
                                  envir=environment(),
                                  clean = TRUE,
                                  quiet = TRUE);

  cat( '\n<a target="_blank" href="', interactiveHeatmap_filename, '">Click here for interactive heatmap in a separate window</a>\n', sep = "");
}


# Prepare unique rows and cols names for pheatmap (annotation rows) and match with rowAnnots and colAnnots row names
originalRowNames = rownames( expMat);
originalColNames = colnames( expMat);
rownames( expMat) = make.unique( originalRowNames);
colnames( expMat) = make.unique( originalColNames);
rownames( rowsAnnot) = rownames( expMat);
rownames( colsAnnot) = colnames( expMat);
# Prepare colors of monitored genes groups for pheatmap (requires named vector matching data factor levels)
monitoredColors = rainbow( nlevels( rowsAnnot[["Monitored"]]), s = 0.8);
names( monitoredColors) = levels( rowsAnnot[["Monitored"]]);

pheatmap( expMat,
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          annotation_row = rowsAnnot,
          annotation_col = colsAnnot,
          labels_row = originalRowNames,
          annotation_colors = list( Monitored = monitoredColors,
                                    Cluster = clustersColor),
          show_colnames = FALSE);




## @knitr heterogeneity_monitoredGenes_expression_projection
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Plot expression values of monitored genes (TODO: message if list empty)
invisible( lapply( names( MONITORED_GENES), function(monitoredGroup)
{
  cat("#####", monitoredGroup, "\n");

  # Plots expression on projected cells (or error message if feature not found)
  invisible( lapply( MONITORED_GENES[[monitoredGroup]], function(featureName)
  {
    print(
      tryCatch( suppressWarnings( # Raise a warning # A warning is raised when all values are 0 (e.g. contamination genes after filtering)
                FeaturePlot( sc10x, features = featureName, reduction = ifelse(exists("useReduction"), useReduction, "umap"), order = TRUE)  +
                  theme( axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         legend.position = "none")),
                error = function(e){return(conditionMessage(e))})); # If gene not found in object (should have been removed earlier anyway)
  }));

  cat(" \n \n"); # Required for '.tabset'
}));




## @knitr heterogeneity_monitoredGenes_expression_violin

# Plot expression values of monitored genes as violinplot for each cluster (TODO: message if list empty)
invisible( lapply( names( MONITORED_GENES), function(monitoredGroup)
{
  cat("#####", monitoredGroup, "\n");

  # Violinplot for expression value of monitored genes by cluster (+ number of 'zero' and 'not zero' cells)
  invisible( lapply( MONITORED_GENES[[monitoredGroup]], violinFeatureByCluster, seuratObject = sc10x, clustersColor = clustersColor));

  cat(" \n \n"); # Required for '.tabset'
}));




# MODULES ANALYSIS
##################

## @knitr heterogeneity_modules

# Just remind the warning for genes names not in object, or modules that were transfered to individual monitoring of genes
if(any( is.na( matchModulesGenes)))
{
    warning( paste( "Following gene(s) from modules list could not be found in experimental data and will be ignored:",
                    paste( modulesGenesNotFound, collapse=" - ")));
}

if(any(modulesToTransfer))
{
    warning( paste0( "Following Module(s) contained very few genes (<",
                     MONITORED_GENES_SMALL_MODULE,
                     "). These genes were transfered to 'Monitored genes' to be analyzed individually: ",
                     paste( names(modulesToTransfer)[modulesToTransfer], collapse = " - ")));
}




## @knitr heterogeneity_modules_scoring

# Compute the score of the cells according to group of monitored genes
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

# Only create interactive heatmap if reasonable number of cells (prevents crashing when rendering and slow html result)
if(length( Cells( sc10x)) < PLOT_RASTER_NBCELLS_THRESHOLD)
{
  # Create a temporary file that will store the code for interactive heatmap and
  # knitted as separate doc to prevent inflating report with demanding figures
  tempFileName = tempfile( fileext = ".Rmd");
  writeLines( con = tempFileName, text = '
---
title: Modules
---

```{r heterogeneity_modules_heatmapInteractive, echo = FALSE}
    # Create a matrix for text on mouse over events (customize "text" that is supposed to show value for adding info)
    hoverTextMarkers = modulesScoreMat;
    hoverTextMarkers[] = paste( paste( "Cell ID:", sc10x[["numID", drop = FALSE]][colnames( modulesScoreMat)[col( modulesScoreMat)],]),
                                paste( "Value:", modulesScoreMat),
                                paste( "Cluster:", sort( clusterID)[col( modulesScoreMat)]),
                                sep = "<br>");

    # Create heatmap
    heatmapPlot = iheatmap( modulesScoreMat,
                            colorbar_position = 7, # Put it on third column to keep space for rownames
                            row_labels = TRUE,
                            col_labels = FALSE,
                            text = hoverTextMarkers,
                            tooltip = setup_tooltip_options(row = TRUE,
                                                            col = TRUE,
                                                            value = TRUE,
                                                            prepend_row = "Module: ",
                                                            prepend_col = "Cell: ",
                                                            prepend_value = ""),
                            layout = list(width = 900, height = 900)) %>%
                  add_col_annotation( annotation = colsAnnot,
                                      colors = list( Cluster = clustersColor),
                                      side = "bottom",
                                      show_colorbar = FALSE);

    # Hack to customize modebar buttons as for plotly native objects
    # see: https://github.com/ropensci/iheatmapr/issues/38#issuecomment-445428166
    heatmapPlot_widget = iheatmapr::to_widget(heatmapPlot);
    # Edit the htmlwidget object itself
    heatmapPlot_widget[["x"]][["config"]][["displaylogo"]] = FALSE;
    heatmapPlot_widget[["x"]][["config"]][["toImageButtonOptions"]] = list( format="svg");  # This one does not seem to work... TODO: Check issue on GitHub
    heatmapPlot_widget[["x"]][["config"]][["modeBarButtons"]]       = list( list( "toImage"),
                                                                            list( "zoom2d", "zoomIn2d", "zoomOut2d", "pan2d", "resetScale2d"));

    # Render the heatmap
    heatmapPlot_widget;
```
');

  # Create a file name for the document knitted separately
  interactiveHeatmap_filename = paste0( outputFilesPrefix, "modules_interactiveHeatmap.html");

  # knit the interactive heatmap in a separate document using current environment
  resultFile = rmarkdown::render( tempFileName,
                                  output_file = I( interactiveHeatmap_filename),
                                  output_dir = PATH_ANALYSIS_OUTPUT,
                                  envir=environment(),
                                  clean = TRUE,
                                  quiet = TRUE);

  cat( '\n<a target="_blank" href="', interactiveHeatmap_filename, '">Click here for interactive heatmap in a separate window</a>\n', sep = "");
}


# Prepare unique rows and cols names for pheatmap (annotation rows) and match with rowAnnots and colAnnots row names
originalRowNames = rownames( modulesScoreMat);
originalColNames = colnames( modulesScoreMat);
rownames( modulesScoreMat) = make.unique( originalRowNames);
colnames( modulesScoreMat) = make.unique( originalColNames);
#rownames( rowsAnnot) = rownames( modulesScoreMat);
rownames( colsAnnot) = colnames( modulesScoreMat);

# Plot the 'non-interactive' heatmap
pheatmap( modulesScoreMat,
          cluster_rows = FALSE,
          cluster_cols = FALSE,
#          annotation_row = rowsAnnot,
          annotation_col = colsAnnot,
          labels_row = originalRowNames,
          annotation_colors = list( Cluster = clustersColor),
          show_colnames = FALSE);




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


