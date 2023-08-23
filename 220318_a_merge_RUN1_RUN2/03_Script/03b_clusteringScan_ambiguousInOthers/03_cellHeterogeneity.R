# ##############################################################################
# This script analyses the heterogeneity of cells using dimension reduction 
# techniques. As opposed to other regular versions, this one loops on several
# resolution values, so all processing and plots must be generated in a single
# chunk. Some plots of regular reports are skipped (PCA, ...) and analysis stops
# after identification of marker genes.
#
# Important: once clusters are identified, the script restarts current analysis
# for each cluster, performing a recursive subclustering (until target recursion
# level reached). An output folder and an html report is created for each
# recursion level (report optionnaly included as iframe in upper level report).
#
# Results are stored in a recursive folder hierarchy:
# ResultsFolder
# |- Report.html (level 1 report, no subclustering)
# |- R0.1
# |  |- objectSeurat_R0.1_allClusters.RDS
# |  |- C1
# |  |  |- objectSeurat_R0.1C1.RDS (contains only cluster1 from level 1)
# |  |  |- Report_R0.1C1.html       (level 2 reports, 1 subclustering)
# |  |  |- R0.1
# |  |  |  |- objectSeurat_R0.1C1_R0.1_allClusters.RDS
# |  |  |  |- C1
# |  |  |  |  |- objectSeurat_R0.1C1_R0.1C1.RDS (contains only cluster1 from level 2)
# |  |  |  |  |- Report_R0.1C1_R0.1C1.html       (level 3 reports, 2 subclusterings)
# |  |  |  |
# |  |  |  |- ... (other clusters, level 2)
# |  |  |
# |  |  |- ... (other resolutions, level 2)
# |  |
# |  |- ... (other clusters, level 1)
# |
# |- ... (other resolutions, level 1)
# |
# 
# ##############################################################################




## @knitr heterogeneity
# From this chunk, all subsequent analyses/plots will be made in a large loop
# sweeping resolution values. As a consequence all titles (+tabsets) and layout
# are handled in this unique chunk instead of Rmd file. It means that we are
# limited to a single set of chunk params (specially for figures dimensions).

### Compute modules score (not to be done in loop)

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



### PCA
message( paste( "Computing PCA"));
nbPC=PCA_NPC
if(PCA_NPC>=length(Cells(sc10x)))
{
  warning( paste0( "Number of cells in object (", length(Cells(sc10x)), ") smaller or equal to requested number of PCs (", PCA_NPC,"), setting lower PC number (", length(Cells(sc10x))-1, ")..." ))
  nbPC = length(Cells(sc10x))-1
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
message( paste( "Computing Neighbors graph"));
sc10x <- FindNeighbors(object    = sc10x,
                       reduction = "pca",
                       dims      = 1:nbPC_findclusters,
                       verbose   = .VERBOSE);


### Compute the dimensionality reductions (plotted later with colored clusters)
# UMAP
umapFilename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "cellsCoordinates_umap.tsv"))
# Check if an output has already been computed earlier, eventually use it if available (see USE_PRECOMPUTED_DIMREDUC_CLUSTERING from analysis params)
if(USE_PRECOMPUTED_DIMREDUC_CLUSTERING && file.exists(umapFilename))
{
  message("Found previous output with 'umap' coordinates, loading it...");
  # Create dim reduc object from 'tsv' file and integrate it in seurat object
  sc10x[["umap"]] = CreateDimReducObject(embeddings = as.matrix( read.table( umapFilename,
                                                                             header = TRUE,
                                                                             row.names = 1,
                                                                             quote = "",
                                                                             sep = "\t")), 
                                         assay = "RNA", 
                                         key = "UMAP_")
} else { 
  
  message("Computing umap...");
  nbPC_dimreduc=DIMREDUC_USE_PCA_NBDIMS
  if(DIMREDUC_USE_PCA_NBDIMS>nbPC)
  {
    warning( paste0( "Number of computed PCs  (", nbPC, ") smaller than requested PCs for 'dimreduc' (", DIMREDUC_USE_PCA_NBDIMS,"), setting lower PC number (", nbPC, ")..." ))
    nbPC_dimreduc = nbPC
  }
  options('Seurat.warn.umap.uwot' = FALSE); # Prevent warning on changed UMAP method
  sc10x = RunUMAP( sc10x, dims = 1:nbPC_dimreduc, seed.use = SEED);

  # Save resulting coordinates for all cells as 'tsv' files
  write.table( Embeddings(sc10x, reduction = "umap"), 
               file= umapFilename, 
               quote = FALSE, 
               row.names = TRUE, 
               col.names = NA, # Add a blank column name for row names (CSV convention)
               sep="\t");
}

# Removed because not really used and crashes script early when number gets low
# # TSNE
# tsneFilename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "cellsCoordinates_tsne.tsv"))
# # Check if an output has already been computed earlier, eventually use it if available (see USE_PRECOMPUTED_DIMREDUC_CLUSTERING from analysis params)
# if(USE_PRECOMPUTED_DIMREDUC_CLUSTERING && file.exists(tsneFilename))
# {
#   message("Found previous output with 'tsne' coordinates, loading it...");
#   # Create dim reduc object from 'tsv' file and integrate it in seurat object
#   sc10x[["tsne"]] = CreateDimReducObject(embeddings = as.matrix( read.table( tsneFilename,
#                                                                              header = TRUE,
#                                                                              row.names = 1,
#                                                                              quote = "",
#                                                                              sep = "\t")), 
#                                          assay = "RNA", 
#                                          key = "tSNE_")
# } else {
#   message("Computing tsne...");
#   nbPC_dimreduc=DIMREDUC_USE_PCA_NBDIMS
#   if(DIMREDUC_USE_PCA_NBDIMS>nbPC)
#   {
#     warning( paste0( "Number of computed PCs  (", nbPC, ") smaller than requested PCs for 'dimreduc' (", DIMREDUC_USE_PCA_NBDIMS,"), setting lower PC number (", nbPC, ")..." ))
#     nbPC_dimreduc = nbPC
#   }
#   sc10x = RunTSNE( sc10x, dims = 1:nbPC_dimreduc, k.seed = SEED);
# 
#   write.table( Embeddings(sc10x, reduction = "tsne"), 
#                file= tsneFilename, 
#                quote = FALSE, 
#                row.names = TRUE, 
#                col.names = NA, # Add a blank column name for row names (CSV convention)
#                sep="\t");
# }
# 

# Prepare a list that will store for each resolution, a matrix of each gene 
# average expression value by cluster (also saved individually in loop)
geneExpByCluster_resolutionList = list();

#message(class(initialClustering))
#message(FINDCLUSTERS_RESOLUTION_VALUES);
#message(paste(ls(environment()), collaspe = " -- "));


# Flag special case when clustering is given as external file (only for level 1)
if(!is.null(initialClustering)) FINDCLUSTERS_RESOLUTION_VALUES = NA;


### Here starts the global loop iterating on FINDCLUSTERS_RESOLUTION values
for(FINDCLUSTERS_RESOLUTION in FINDCLUSTERS_RESOLUTION_VALUES)
{

  # Prepare a name for current resolution (for filenames, columns, etc...)
  currentResolutionName = if(is.na(FINDCLUSTERS_RESOLUTION)) "initialClustering" else paste0( "R", format( FINDCLUSTERS_RESOLUTION, scientific = FALSE));
  
  # Make a title (/tabset) for each resolution
  cat( "\n\n## ",currentResolutionName, "\n\n", sep="");

  message(paste("Resolution:", currentResolutionName));


  # Create output subdirectory for current resolution and update corresponding variables 
  PATH_ANALYSIS_OUTPUT_CURRENTRES = file.path( PATH_ANALYSIS_OUTPUT, currentResolutionName);
  dir.create( PATH_ANALYSIS_OUTPUT_CURRENTRES, recursive = TRUE, showWarnings = FALSE);
  outputFilesPrefix_CURRENTRES = paste0( outputFilesPrefix, currentResolutionName, "_");

  
  # Search clusters (or use initialClustering from file)
  clusteringResultFilename = file.path( PATH_ANALYSIS_OUTPUT_CURRENTRES, paste0( outputFilesPrefix_CURRENTRES, "cellsClusterIdentities.tsv"))
  if(is.na(FINDCLUSTERS_RESOLUTION))
  {
    # For level 1, one might use an eternal file for clusters definition (see INITIAL_CLUSTERING_PATH)
    Idents(sc10x) = initialClustering[,"identity", drop = FALSE]; # Already sorted as in object on loading
  } else
  {
    # Check if an output has already been computed earlier, eventually use it if available (see USE_PRECOMPUTED_DIMREDUC_CLUSTERING from analysis params)
    if(USE_PRECOMPUTED_DIMREDUC_CLUSTERING && file.exists(clusteringResultFilename))
    {
      message("Found previous clustering result, loading it...");
      Idents(sc10x) = factor(read.table( clusteringResultFilename,
                                  header = TRUE,
                                  row.names = 1,
                                  quote = "",
                                  sep = "\t")[["identity"]])
    } else {
      message("Computing clustering...");
      sc10x <- FindClusters(object             = sc10x,
                            resolution         = FINDCLUSTERS_RESOLUTION,
                            algorithm          = FINDCLUSTERS_ALGORITHM,
                            temp.file.location = "/tmp/",
                            verbose            = .VERBOSE);
    }
  }

  clustersCount = as.data.frame( table( Cluster = Idents(sc10x)), responseName = "CellCount");

  # Define a set of colors for clusters (based on ggplot default)
  clustersColor = hue_pal()( nlevels( Idents( sc10x)));
  names( clustersColor) = levels( Idents( sc10x));

  # Save cells cluster identity as determined with 'FindClusters'
  write.table( data.frame(sc10x[["numID"]], identity = Idents(sc10x)), 
              file = clusteringResultFilename, 
              quote = FALSE, 
              row.names = TRUE, 
              col.names = NA, # Add a blank column name for row names (CSV convention)
              sep="\t");

  # Also save cluster color attribution for reference
  # Save cells cluster identity as determined with 'FindClusters'
  write.table( clustersColor, 
              file = file.path( PATH_ANALYSIS_OUTPUT_CURRENTRES, paste0( outputFilesPrefix_CURRENTRES, "clusterColors.tsv")), 
              quote = FALSE, 
              row.names = TRUE, 
              col.names = FALSE,
              sep="\t");




  # Gather data to be visualized together (cell name + numID + metrics) for plotly hover
  cellsData = cbind( "Cell" = colnames( sc10x), # cell names from rownames conflict with search field in datatable
                     sc10x[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                     "percent.mito" = if("percent.mito" %in% colnames(sc10x[[]])) as.numeric(format(sc10x[["percent.mito", drop = TRUE]], digits = 5)) else NULL,
                     "percent.ribo" = if("percent.ribo" %in% colnames(sc10x[[]])) as.numeric(format(sc10x[["percent.ribo", drop = TRUE]], digits = 5)) else NULL,
                     "Cluster" = Idents( sc10x));

  # Create text to show under cursor for each cell for plotly hover
  hoverText = do.call(paste, c(Map( paste,
                                    c( "",
                                       "Cell ID: ",
                                       "# UMIs: ",
                                       "# Genes: ",
                                       if("percent.mito" %in% colnames(sc10x[[]])) "% Mito: ",
                                       if("percent.ribo" %in% colnames(sc10x[[]])) "% Ribo: "),
                                    cellsData[-ncol( cellsData)], # Do not include cluster in hover text
                                    sep = ""),
                      sep = "\n"));

  cat( "\n<br>\n");




  #### PLOT DIMENSIONALITY REDUCTION (TSNE/UMAP)

  cat( "\n\n### Dimensionality reduction {.tabset .tabset-fade .tabset-pills}\n\n");
  
  #### Plot dimensionality reduction results
  # Interactive plot using plotly
  # In regular reports, tSNE & UMAP are requested sequentially in the Rmd file.
  # Here we need to add an internal loop in R code (but still using previously 
  # defined variable 'useReduction' to specify which reduction method to use).
  cat( "\n\n<br>\n\n");

  for(useReduction in c( if(!is.null( originaCellsCoordinates)) "externalRef", "umap")) # Removed "tsne"
  {
    cat( "\n\n####", useReduction, "\n\n");

    message(paste("Dimensionality reduction:", useReduction));


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

    # Get coordinates corresponding to requested  representation (umap/tsne/external)
    plotCoordinates = if(useReduction == "externalRef") originaCellsCoordinates else Embeddings(sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"));

    dimReducData = data.frame( plotCoordinates,
                               Cluster = Idents( sc10x)[rownames(plotCoordinates)]);
    
    # Compute the median coordinates for each cluster to plot corresponding name 
    dimReducLabels = as.data.frame( do.call( rbind, 
                                             by( dimReducData, 
                                                 dimReducData[["Cluster"]], 
                                                 function(x)
                                                 { 
                                                   return( data.frame( lapply( x[colnames( plotCoordinates)], median), 
                                                                       "Cluster" = x[1, "Cluster"]));
                                                 })));
    

    # ## Do it in pure plotly
    # 
    # # Set figure dimensions
    # dimReducWidth  = 800;
    # dimReducHeight = 800;
    # # Set point size in a variable that will be used for the plot and the slider
    # initialPointSize = 5;
    # 
    # # Reusing 'hoverText' summarizing cell info (generated for split violin plots)
    # plotDimReduc = plot_ly(dimReducData,
    #                        x = as.formula( paste( "~", colnames( dimReducData)[1])),
    #                        y = as.formula( paste( "~", colnames( dimReducData)[2])),
    # #                       name   = ~paste("Cluster", Cluster),
    #                        color  = ~Cluster,
    #                        colors = clustersColor,
    #                        width  = dimReducWidth,
    #                        height = dimReducHeight) %>%
    #   add_trace( type = "scattergl",
    #              mode = "markers",
    #              marker = list( size = initialPointSize,
    #                             opacity = 0.6),
    #              text = hoverText,
    #              hoverinfo = "text+name") %>%
    # #  hide_colorbar() %>%
    # #  hide_legend() %>%
    #   config( displaylogo = FALSE,
    #           toImageButtonOptions = list(format='svg'),
    #           modeBarButtons = list(list('toImage'),
    #                                 list('zoom2d', 'zoomIn2d', 'zoomOut2d', 'pan2d', 'resetScale2d')));
    # 
    # 
    # # Set a unique html ID to access plot with html/javascript for external slider
    # figID = if(exists("figID")) figID+1 else 0;
    # plotDimReduc[["elementId"]] = paste0("plotDimReduc_withSlider", figID);
    # 
    # 
    # print(tagList(plotDimReduc))
    # 
    # cat( "\n<br>\n");
    # 
    # # Create a html slider using javascript for controlling point size
    # cat( paste0( '\n',
    #             '<p>Point size: <span id="valueTag_', figID, '">', initialPointSize, '</span></p>',
    #             '<input type="range" min="2" max="50" value="', initialPointSize, '" id="slider_', figID, '" style="width: ', dimReducWidth, 'px;">',
    #             '<script> document.getElementById("slider_', figID, '").oninput = function() { Plotly.restyle("', plotDimReduc[["elementId"]], '", {"marker.size": this.value}); document.getElementById("valueTag_', figID, '").innerHTML = this.value; } </script>',
    #             '\n'));
    # 
    # # Create a html/js block to change color of clusters and get selected values
    # cat( paste0( '\n', 
    #             'Color:&nbsp;',
    #             # Create a selector, each option corresponding to a plotly trace (defined here by color)
    #             '<select id="colorChange_traceSelector_plotly', figID, '" style="width: 200px;">', 
    #               paste('<option value="', c( 'undefined', 0:(nlevels( cellsData[["Cluster"]])-1)), '">', c('All', levels( cellsData[["Cluster"]])), '</option>', sep ="", collapse = " " ), 
    #             '</select>', 
    #             '&nbsp;&nbsp;',
    #             # Create HTML anchor for color picker
    #             '<input type="text" id="colorChange_colorPicker_plotly', figID, '" />', 
    #             # Call javascript that transforms it in an actual color picker (updates plotly color of selected trace on change, and resets its own value to 'clear')
    #             '<script> $("#colorChange_colorPicker_plotly', figID, '").spectrum({ allowEmpty: true, showInput: true, preferredFormat: "rgb", clickoutFiresChange: false, showPalette: true, showSelectionPalette: true, palette: [ ', paste('"', clustersColor, '"', collapse = ", ", sep = ""), ' ], localStorageKey: "myLocalPalette", change: function(color) { if(color) { Plotly.restyle("', plotDimReduc[["elementId"]], '", {"marker.color": color.toHexString()}, $("#colorChange_traceSelector_plotly', figID, '").val()); $("#colorChange_colorPicker_plotly', figID, '").spectrum("set", ""); } } }); </script>', 
    #             '&nbsp;&nbsp;Get current colors:&nbsp;',
    #             # Add button to extract current colors from plotly figure (hexadecimal format)
    #             '<button id="colorGetter_alertButton_plotly', figID, '" onclick="getColorTab(', plotDimReduc[["elementId"]], ', \'scattergl\', \'hex\')">Hexadecimal</button>',
    #             # Add button to extract current colors from plotly figure (RGB format)
    #             '<button id="colorGetter_alertButton_plotly', figID, '" onclick="getColorTab(', plotDimReduc[["elementId"]], ', \'scattergl\', \'rgb\')">RGB</button>',
    #             '\n'));
    
    # Actually use a simple raster plot as a bug shows an extra figure and the number of opengl panels is very limited...
    



    #### Dimensionality reduction 2D (ggplot)
    
    # Plot a thumbnail highlighting cluster cells
    ggFigure = ggplot(  dimReducData, 
                        aes_string( x = colnames( dimReducData)[1], 
                                    y = colnames( dimReducData)[2],
                                    color = "Cluster")) +
      geom_point( alpha = PLOT_DIMREDUC_GROUPS_ALPHA,
                  size = if(PLOT_DIMREDUC_GROUPS_POINTSIZE==0) Seurat:::AutoPointSize(cellsData) else PLOT_DIMREDUC_GROUPS_POINTSIZE)
    
    if(if(exists( "showDimReducLegend")) showDimReducLegend else TRUE)
    {
      ggFigure = ggFigure +
        geom_text( data = dimReducLabels,
                   aes( label = Cluster),
                   color = "black",
                   size = 7)
    }
    
    ggFigure = ggFigure +
      scale_color_manual( name = "Cluster", values = clustersColor) +
      theme_classic() +
      theme( #axis.text  = element_blank(),
        axis.title = element_blank(),
        #axis.ticks = element_blank(),
        #axis.line  = element_blank(),
        legend.position = "none",
        plot.margin = margin( 0, 0, 0, 0, "cm"),
        plot.title = element_text( face = "bold",
                                   size = rel( 1.8), #rel( 16/14),
                                   hjust = 0.5,
                                   vjust = 1,
                                   margin = margin( b = 7)));
    
    # Plot figure in a sub chunk with its own specific dimensions
    #subchunkify( ggFigure, options = "fig.dim=c(9,9), out.width='100%'");
    print(ggFigure)
    
  }
  
  cat(" \n \n"); # Required for '.tabset'
  


  #### Show number of cells in each cluster
  cat( "\n\n### Cells distribution\n\n");

  message("Cells cluster distribution...");


  print( knitr::kable( clustersCount,
                       align = "c")) #%>% kable_styling(bootstrap_options = c("condensed"), full_width = FALSE);

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
               file= file.path( PATH_ANALYSIS_OUTPUT_CURRENTRES, paste0( outputFilesPrefix_CURRENTRES, "normExpressionByCluster.tsv")),
               quote = FALSE,
               row.names = TRUE,
               col.names = NA, # Add a blank column name for row names (CSV convention)
               sep="\t");

  # Also save it in external list (all resolutions)
  geneExpByCluster_resolutionList[[currentResolutionName]] = geneExpByCluster;




  # #### Show cells stats (UMIs, Genes, Mito, Ribo) distribution split by cluster
  cat( "\n\n### Statistics by cluster\n\n");

  message("Statistics by cluster...");


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

  lypanel_mitos = if("percent.mito" %in% colnames(sc10x[[]])) plotViolinJitter( cellsData,
                                                                                xAxisFormula = ~as.numeric( Cluster),
                                                                                yAxisFormula = ~percent.mito,
                                                                                colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                                                                yAxisTitle = "% Mito",
                                                                                hoverText = hoverText,
                                                                                traceName = ~paste( "Cluster", Cluster),
                                                                                xTicklabels = levels(cellsData[["Cluster"]]),
                                                                                panelWidth = panelWidth,
                                                                                panelHeight = panelHeight) else NULL;

  lypanel_ribos = if("percent.ribo" %in% colnames(sc10x[[]])) plotViolinJitter( cellsData,
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

  print(tagList(plotPanels))

  cat( "\n<br>\n");


  
  
  #### Show the distribution of selected categorical metadata in clusters
  # Mostly used to check uniform repartition of batches in clusters
  
  cat( "\n\n### Categorical distributions\n\n");
  
  for( currentCategoryName in names(CATEGORICAL_IN_CLUSTERS))
  {
    cat( "\n\n#### ", currentCategoryName, "\n\n", sep = "");
    message(paste0("Categorical distribution: '", currentCategoryName, "'"));
    
    # Extract information from Seurat object metadata and attribute name
    categoryData = data.frame( sc10x[[CATEGORICAL_IN_CLUSTERS[currentCategoryName], drop = TRUE]], 
                               Idents(sc10x))
    names(categoryData) = c( currentCategoryName, "Cluster")
    
    # Plot proportional barplot (add stats as text)
    categoryPlot = ggplot( categoryData , aes_string(x = "Cluster", fill = currentCategoryName)) +
                     geom_bar( position = "fill") + 
                     geom_text( aes( label=paste0(..count.., "\n(", 100*signif( ..count.. / tapply( ..count.., ..x.., sum)[as.character( ..x..)], digits=3), "%)")), 
                                stat = 'count', 
                                position = position_fill( vjust = 0.5),
                                size = 3) +
                     geom_text( aes( label=after_stat(count), group = 1), stat='count', position = "fill", vjust = -0.5) +
                     ylab( "") +
                     theme( legend.text = element_text( size = 8), legend.position = "top", legend.box = "horizontal", legend.title = element_blank())
    

    print(categoryPlot)
      
    cat( "\n\n");
  }
  


  # MARKER GENES
  ##############

  #if(FALSE)
  if(nlevels( Idents( sc10x)) > 1)
  {
    #### Identification of marker genes for current clusters
    cat( "\n\n### Marker genes table\n\n");

    message("Marker genes table...");


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
                file= file.path( PATH_ANALYSIS_OUTPUT_CURRENTRES, paste0( outputFilesPrefix_CURRENTRES, "MarkerGenes.tsv")),
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

    # Create datatable
    print( htmltools::tagList( datatable( topMarkersDT,
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
                                                scales::alpha(clustersColor, 0.3)))));

    cat( "\n<br>\n");

    
    
    ## For external 'png' figure files
    # Match with figures dimensions in report (based on chunk params)
    defaultDim = c( 6, 6);  # Fallback values if chunk params not available
    defaultDpi = 72;        # (e.g. executing directly from R)
    
    chunkDim    = knitr::opts_current$get( "fig.dim");
    chunkDpi    = knitr::opts_current$get( "dpi");
    
    figDim = if( is.null( chunkDim) ) defaultDim else chunkDim;
    figDpi = if( is.null( chunkDpi) ) defaultDpi else chunkDpi;
    
    

    #### Markers expression on dimensionality reduction
    cat( "\n\n### Markers expression {.tabset .tabset-fade}\n\n");
    
    for(useReduction in c(if(!is.null(originaCellsCoordinates)) "externalRef", "umap")) # Removed "tsne"
    {
      cat( "\n\n####", useReduction, "{.tabset .tabset-fade}\n\n");
      
      message( paste( "Markers expression:", useReduction));


      # Get coordinates corresponding to requested representation (umap/tsne/external)
      plotCoordinates = if(useReduction == "externalRef") originaCellsCoordinates else Embeddings(sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"));

      dimReducData = data.frame( plotCoordinates,
                                 Cluster = Idents( sc10x)[rownames(plotCoordinates)]);

      # Compute the position of labels for groups using median of group coordinates
      dimReducLabels = cbind( as.data.frame( do.call( rbind, by( dimReducData[1:2], dimReducData[["Cluster"]], apply, 2, median))), Cluster = levels(dimReducData[["Cluster"]]))


      # Plot a thumbnail highlighting cluster cells
      ggFigure = ggplot( dimReducData,
                         aes_string( x = colnames( dimReducData)[1],
                                     y = colnames( dimReducData)[2],
                                     color = "Cluster")) +
                  geom_point( alpha = PLOT_DIMREDUC_GROUPS_ALPHA,
                              size = if(PLOT_DIMREDUC_GROUPS_POINTSIZE==0) Seurat:::AutoPointSize(dimReducData) else PLOT_DIMREDUC_GROUPS_POINTSIZE) +
                  geom_text( data = dimReducLabels,
                             aes( label = Cluster),
                             color = "black",
                             size = 7) +
                  scale_color_manual( name = "Cluster", values = clustersColor) +
                  theme_classic() +
                  theme( #axis.text  = element_blank(),
                         axis.title = element_blank(),
                         #axis.ticks  = element_blank(),
                         #axis.line = element_blank(),
                         legend.position = "none",
                         plot.margin = margin( 0, 0, 0, 0, "cm"),
                         plot.title = element_text( face = "bold",
                                                    size = rel( 16/14),
                                                    hjust = 0.5,
                                                    vjust = 1,
                                                    margin = margin( b = 7)));
      print( ggFigure);

      cat(" \n \n");


      # Plot expression values of marker genes on dimreduc figures for each cluster (TODO: message if list empty)
      invisible( lapply( names( topMarkers), function(clusterName)
      {        
        cat("\n\n##### Cl. <span style='border-radius: 3px; border: 3px solid ", clustersColor[clusterName], "; padding:0px 2px'>", clusterName, "</span>\n\n");

        ### Highlight cells of current cluster on a dimreduc plot

        # Separate layers to plot red points on top of gray ones
        ggFigure = ggplot( dimReducData[ dimReducData[["Cluster"]] != clusterName,],
                           aes_string( x = colnames( dimReducData)[1],
                                       y = colnames( dimReducData)[2])) +
                    geom_point( color = "#44444422",
                                alpha = PLOT_DIMREDUC_GROUPS_ALPHA,
                                size = if(PLOT_DIMREDUC_GROUPS_POINTSIZE==0) Seurat:::AutoPointSize(dimReducData) else PLOT_DIMREDUC_GROUPS_POINTSIZE) +
                    geom_point( data = dimReducData[ dimReducData[["Cluster"]] == clusterName,],
                                color = "#FF000088",
                                alpha = PLOT_DIMREDUC_GROUPS_ALPHA,
                                size = 1.5 * if(PLOT_DIMREDUC_GROUPS_POINTSIZE==0) Seurat:::AutoPointSize(dimReducData) else PLOT_DIMREDUC_GROUPS_POINTSIZE) +
                    labs( title = paste("Cluster", clusterName)) +
                    theme_classic() +
                    theme( #axis.text  = element_blank(),
                            axis.title = element_blank(),
                            #axis.ticks  = element_blank(),
                            #axis.line = element_blank(),
                            legend.position = "none",
                            plot.margin = margin( 0, 0, 0, 0, "cm"),
                            plot.title = element_text( face = "bold",
                                                      size = rel( 16/14),
                                                      hjust = 0.5,
                                                      vjust = 1,
                                                      margin = margin( b = 7)));
        print( ggFigure);
        
        ### Plot individual gene expression on dimensionality reduction
        
        # Plots expression on projected cells
        invisible( lapply( head( topMarkers[[clusterName]][["gene"]], FINDMARKERS_SHOWTOP_FIGS) , function(featureName)
        {
          # Get the requested feature values from Seurat object
          currentFeatureExpression = GetAssayData(sc10x, slot="data")[featureName, rownames(dimReducData)];
          
          # Set how data will be reorganized to control which points are drawn first
          cellsOrdering = if( is.null( PLOT_DIMREDUC_EXPRESSION_MAXONTOP)) 1:nrow(dimReducData) else order(currentFeatureExpression, decreasing= !PLOT_DIMREDUC_EXPRESSION_MAXONTOP);
          
          ggFigure = ggplot( cbind( dimReducData, currentFeatureExpression)[cellsOrdering, ],
                             aes_string( x = colnames( dimReducData)[1],
                                         y = colnames( dimReducData)[2],
                                         color = "currentFeatureExpression")) +
            geom_point( alpha = PLOT_DIMREDUC_EXPRESSION_ALPHA,
                        size = if(PLOT_DIMREDUC_EXPRESSION_POINTSIZE==0) Seurat:::AutoPointSize(dimReducData) else PLOT_DIMREDUC_EXPRESSION_POINTSIZE) +
            scale_color_gradient( low = "lightgrey", high = "blue") +
            labs( title = featureName) +
            theme_classic() +
            theme( #axis.text  = element_blank(),
              axis.title = element_blank(),
              #axis.ticks  = element_blank(),
              #axis.line = element_blank(),
              legend.position = "none",
              plot.margin = margin( 0, 0, 0, 0, "cm"),
              plot.title = element_text( face = "bold",
                                         size = rel( 16/14),
                                         hjust = 0.5,
                                         vjust = 1,
                                         margin = margin( b = 7)));

          print( ggFigure);
          
        }));
        
        ### Same but as external png files (with eventually more figures)
        
        if(RENDEREXTERNALFIGURES_MARKERS)
        {
          # Create subfolder for current resolution to store figure files
          pathCurrentFigures = file.path( PATH_ANALYSIS_OUTPUT_CURRENTRES, "exportedFigures", "markers", useReduction, clusterName);
          dir.create(pathCurrentFigures, showWarnings = FALSE, recursive = TRUE);

          # Plots expression on projected cells (loop on features/genes)
          invisible( mclapply( head( topMarkers[[clusterName]][["gene"]], FINDMARKERS_SHOWTOP_FIGS_EXT) , function(featureName)
          {
            figFilename = file.path( pathCurrentFigures, paste0( featureName, ".png"))
            ## Create png file
            #invisible( png( figFilename, 
            #                width = figDim[1],
            #                height = figDim[2],
            #                res = figDpi));
          
            # Get the requested feature values from Seurat object
            currentFeatureExpression = GetAssayData(sc10x, slot="data")[featureName, rownames(dimReducData)];
          
            # Set how data will be reorganized to control which points are drawn first
            cellsOrdering = if( is.null( PLOT_DIMREDUC_EXPRESSION_MAXONTOP)) 1:nrow(dimReducData) else order(currentFeatureExpression, decreasing= !PLOT_DIMREDUC_EXPRESSION_MAXONTOP);
            
            ggFigure = ggplot( cbind( dimReducData, currentFeatureExpression)[cellsOrdering, ],
                               aes_string( x = colnames( dimReducData)[1],
                                           y = colnames( dimReducData)[2],
                                           color = "currentFeatureExpression")) +
              geom_point( alpha = PLOT_DIMREDUC_EXPRESSION_ALPHA,
                          size = if(PLOT_DIMREDUC_EXPRESSION_POINTSIZE==0) Seurat:::AutoPointSize(dimReducData) else PLOT_DIMREDUC_EXPRESSION_POINTSIZE) +
              scale_color_gradient( low = "lightgrey", high = "blue") +
              labs( title = featureName) +
              theme_classic() +
              theme( #axis.text  = element_blank(),
                axis.title = element_blank(),
                #axis.ticks  = element_blank(),
                #axis.line = element_blank(),
                legend.position = "none",
                plot.margin = margin( 0, 0, 0, 0, "cm"),
                plot.title = element_text( face = "bold",
                                           size = rel( 16/14),
                                           hjust = 0.5,
                                           vjust = 1,
                                           margin = margin( b = 7)));
          
            #print( ggFigure);
            resNULL = ggsave( filename = figFilename,
                              ggFigure,
                              device = "png",
                              width = figDim[1],
                              height = figDim[2],
                              units = "in",
                              dpi = figDpi);
            
            return( NULL);
            
            ## Close png descriptor (twice to compensate for a bug when using png inside a rmarkdown rendering...)
            #tryCatch(invisible(dev.off()), error=function(e){});
            #tryCatch(invisible(dev.off()), error=function(e){});
            
            #include_graphics(figFilename);
            
          }, mc.cores = NBCORES_PLOT_EXTERNAL_FIGURES));

        } # /RENDEREXTERNALFIGURES_MARKERS
          
        cat(" \n \n"); # Required for '.tabset'
      }));

      
      
      cat(" \n \n"); # Required for '.tabset'
      
    }


    cat( "\n\n####  Violin {.tabset .tabset-fade}\n\n");
    
    message( "Markers expression: violin");

    
    ### Plot expression values of marker genes as violinplot for each cluster (TODO: message if list empty)
    invisible( lapply( names( topMarkers), function(clusterName)
    {
      cat("\n\n##### Cl. <span style='border-radius: 3px; border: 3px solid ", clustersColor[clusterName], "; padding:0px 2px'>", clusterName, "</span>\n\n");
      
      # Remind cluster name in an empty figure to keep consistent alignment of panels between tabs
      plot( c( 0, 1), c( 0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n');
      text( x = 0.5, y = 0.5, paste( "Cluster", clusterName), cex = 2, col = clustersColor[clusterName]);
      
      # Violinplot for expression value of marker genes by cluster (+ number of 'zero' and 'not zero' cells)
      invisible( lapply(  head( topMarkers[[clusterName]][["gene"]], FINDMARKERS_SHOWTOP_FIGS), violinFeatureByCluster, seuratObject = sc10x, clustersColor = clustersColor));
      

      if(RENDEREXTERNALFIGURES_MARKERS)
      {
        ### Same violinplots but as external png files (with eventually more figures)
        # Create subfolder for current resolution to store figure files
        pathCurrentFigures = file.path( PATH_ANALYSIS_OUTPUT_CURRENTRES, "exportedFigures", "markers", "violinplots", clusterName);
        dir.create(pathCurrentFigures, showWarnings = FALSE, recursive = TRUE);
        
        invisible( mclapply(  head( topMarkers[[clusterName]][["gene"]], FINDMARKERS_SHOWTOP_FIGS_EXT), function(featureName)
                    {
                      ## Create png file
                      #invisible( png( file.path( pathCurrentFigures, paste0( featureName, ".png")), 
                      #           width = figDim[1],
                      #           height = figDim[2],
                      #           res = figDpi));
        
                      resNULL = ggsave( filename = file.path( pathCurrentFigures, paste0( featureName, ".png")),
                                        violinFeatureByCluster(featureName, seuratObject = sc10x, clustersColor = clustersColor),
                                        device = "png",
                                        width = figDim[1],
                                        height = figDim[2],
                                        units = "in",
                                        dpi = figDpi);
        
                      return( NULL);
                      
                      ## Close png descriptor (twice to compensate for a bug when using png inside a rmarkdown rendering...)
                      #tryCatch(invisible(dev.off()), error=function(e){});
                      #tryCatch(invisible(dev.off()), error=function(e){});
                    }, mc.cores = NBCORES_PLOT_EXTERNAL_FIGURES));
      }
      
      cat(" \n \n"); # Required for '.tabset'
      
    }));
    
    cat(" \n \n"); # Required for '.tabset'
    
  }

  
  
  
  # MONITORED GENES (external figures only)
  #################
  
  ## For external 'png' figure files
  # Match with figures dimensions in report (based on chunk params)
  defaultDim = c( 6, 6);  # Fallback values if chunk params not available
  defaultDpi = 72;        # (e.g. executing directly from R)
  
  chunkDim    = knitr::opts_current$get( "fig.dim");
  chunkDpi    = knitr::opts_current$get( "dpi");
  
  figDim = if( is.null( chunkDim) ) defaultDim else chunkDim;
  figDpi = if( is.null( chunkDpi) ) defaultDpi else chunkDpi;
  
  if(RENDEREXTERNALFIGURES_MONITORED)
  {
    # Plot expression values of marker genes on dimreduc figures for each cluster (TODO: message if list empty)
    invisible( lapply( names( MONITORED_GENES), function(groupName)
    {
      
      # Dimension reduction (scatterplot)
      for(useReduction in c(if(!is.null(originaCellsCoordinates)) "externalRef", "umap")) # Removed "tsne"
      {
        
        message( paste0( "Monitored genes '", groupName , "': ", useReduction));
        
        # Get coordinates corresponding to requested representation (umap/tsne/external)
        plotCoordinates = if(useReduction == "externalRef") originaCellsCoordinates else Embeddings(sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"));
        
        dimReducData = data.frame( plotCoordinates,
                                   Cluster = Idents( sc10x)[rownames(plotCoordinates)]);
        
        # Compute the position of labels for groups using median of group coordinates
        dimReducLabels = cbind( as.data.frame( do.call( rbind, by( dimReducData[1:2], dimReducData[["Cluster"]], apply, 2, median))), Cluster = levels(dimReducData[["Cluster"]]))
        
        
        ### Plot individual gene expression on dimensionality reduction        
        
        # Create subfolder for current resolution to store figure files
        pathCurrentFigures = file.path( PATH_ANALYSIS_OUTPUT_CURRENTRES, "exportedFigures", "monitored", useReduction, groupName);
        dir.create(pathCurrentFigures, showWarnings = FALSE, recursive = TRUE);
        
        # Plots expression on projected cells (loop on features/genes)
        invisible( mclapply( MONITORED_GENES[[groupName]] , function(featureName)
        {
          figFilename = file.path( pathCurrentFigures, paste0( featureName, ".png"))
          ## Create png file
          #invisible( png( figFilename, 
          #                width = figDim[1],
          #                height = figDim[2],
          #                res = figDpi));
          
          # Get the requested feature values from Seurat object
          currentFeatureExpression = GetAssayData(sc10x, slot="data")[featureName, rownames(dimReducData)];
          
          # Set how data will be reorganized to control which points are drawn first
          cellsOrdering = if( is.null( PLOT_DIMREDUC_EXPRESSION_MAXONTOP)) 1:nrow(dimReducData) else order(currentFeatureExpression, decreasing= !PLOT_DIMREDUC_EXPRESSION_MAXONTOP);
          
          ggFigure = ggplot( cbind( dimReducData, currentFeatureExpression)[cellsOrdering, ],
                             aes_string( x = colnames( dimReducData)[1],
                                         y = colnames( dimReducData)[2],
                                         color = "currentFeatureExpression")) +
            geom_point( alpha = PLOT_DIMREDUC_EXPRESSION_ALPHA,
                        size = if(PLOT_DIMREDUC_EXPRESSION_POINTSIZE==0) Seurat:::AutoPointSize(dimReducData) else PLOT_DIMREDUC_EXPRESSION_POINTSIZE) +
            scale_color_gradient( low = "lightgrey", high = "blue") +
            labs( title = featureName) +
            theme_classic() +
            theme( #axis.text  = element_blank(),
              axis.title = element_blank(),
              #axis.ticks  = element_blank(),
              #axis.line = element_blank(),
              legend.position = "none",
              plot.margin = margin( 0, 0, 0, 0, "cm"),
              plot.title = element_text( face = "bold",
                                         size = rel( 16/14),
                                         hjust = 0.5,
                                         vjust = 1,
                                         margin = margin( b = 7)));
          
          #print( ggFigure);
          resNULL = ggsave( filename = figFilename,
                            ggFigure,
                            device = "png",
                            width = figDim[1],
                            height = figDim[2],
                            units = "in",
                            dpi = figDpi);
          
          return( NULL);
          
          ## Close png descriptor (twice to compensate for a bug when using png inside a rmarkdown rendering...)
          #tryCatch(invisible(dev.off()), error=function(e){});
          #tryCatch(invisible(dev.off()), error=function(e){});
          
          #include_graphics(figFilename);
          
        }, mc.cores = NBCORES_PLOT_EXTERNAL_FIGURES))
      }
      
      
      # Violinplots (same loop as dimreduc as we don't need to separate them in different paragraphs)
      # Create subfolder for current resolution to store figure files
      message( paste0( "Monitored genes '", groupName , "': violinplot"));
      
      # Create subfolder for current resolution to store figure files
      pathCurrentFigures = file.path( PATH_ANALYSIS_OUTPUT_CURRENTRES, "exportedFigures", "monitored", "violinplots", groupName);
      dir.create(pathCurrentFigures, showWarnings = FALSE, recursive = TRUE);
      
      invisible( mclapply( MONITORED_GENES[[groupName]], function(featureName)
      {
        ## Create png file
        #invisible( png( file.path( pathCurrentFigures, paste0( featureName, ".png")), 
        #           width = figDim[1],
        #           height = figDim[2],
        #           res = figDpi));
        
        resNULL = ggsave( filename = file.path( pathCurrentFigures, paste0( featureName, ".png")),
                          violinFeatureByCluster(featureName, seuratObject = sc10x, clustersColor = clustersColor),
                          device = "png",
                          width = figDim[1],
                          height = figDim[2],
                          units = "in",
                          dpi = figDpi);
        
        return( NULL);
        
        ## Close png descriptor (twice to compensate for a bug when using png inside a rmarkdown rendering...)
        #tryCatch(invisible(dev.off()), error=function(e){});
        #tryCatch(invisible(dev.off()), error=function(e){});
      }, mc.cores = NBCORES_PLOT_EXTERNAL_FIGURES));
      
    }))
  }  # /RENDEREXTERNALFIGURES_MONITORED
  
  
  
  
  
  
  # MODULES GENES (external figures only)
  #################
  
  ## For external 'png' figure files
  # Match with figures dimensions in report (based on chunk params)
  defaultDim = c( 6, 6);  # Fallback values if chunk params not available
  defaultDpi = 72;        # (e.g. executing directly from R)
  
  chunkDim    = knitr::opts_current$get( "fig.dim");
  chunkDpi    = knitr::opts_current$get( "dpi");
  
  figDim = if( is.null( chunkDim) ) defaultDim else chunkDim;
  figDpi = if( is.null( chunkDpi) ) defaultDpi else chunkDpi;
  
  if(RENDEREXTERNALFIGURES_MODULES)
  {
    # Plot expression values of marker genes on dimreduc figures for each cluster (TODO: message if list empty)
    invisible( lapply( names( MODULES_GENES), function(groupName)
    {
      
      # Create the module name as Seurat does (to get values from Seurat metadata)
      featureName = paste0(groupName, "1");
      
      
      # Dimension reduction (scatterplot)
      for(useReduction in c(if(!is.null(originaCellsCoordinates)) "externalRef", "umap")) # Removed "tsne"
      {
        
        message( paste0( "Modules genes '", groupName , "': ", useReduction));
        
        # Get coordinates corresponding to requested representation (umap/tsne/external)
        plotCoordinates = if(useReduction == "externalRef") originaCellsCoordinates else Embeddings(sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"));
        
        dimReducData = data.frame( plotCoordinates,
                                   Cluster = Idents( sc10x)[rownames(plotCoordinates)]);
        
        # Compute the position of labels for groups using median of group coordinates
        dimReducLabels = cbind( as.data.frame( do.call( rbind, by( dimReducData[1:2], dimReducData[["Cluster"]], apply, 2, median))), Cluster = levels(dimReducData[["Cluster"]]))
        
        
        ### Plot individual gene expression on dimensionality reduction        
        
        # Create subfolder for current resolution to store figure files
        pathCurrentFigures = file.path( PATH_ANALYSIS_OUTPUT_CURRENTRES, "exportedFigures", "modules", useReduction);
        dir.create(pathCurrentFigures, showWarnings = FALSE, recursive = TRUE);
        

        figFilename = file.path( pathCurrentFigures, paste0( groupName, ".png"))
        ## Create png file
        #invisible( png( figFilename, 
        #                width = figDim[1],
        #                height = figDim[2],
        #                res = figDpi));
        
        # Get the requested feature values from Seurat object metadata
        currentFeatureExpression = sc10x[[featureName, drop = TRUE]][rownames(dimReducData)];
        
        # Set how data will be reorganized to control which points are drawn first
        cellsOrdering = if( is.null( PLOT_DIMREDUC_EXPRESSION_MAXONTOP)) 1:nrow(dimReducData) else order(currentFeatureExpression, decreasing= !PLOT_DIMREDUC_EXPRESSION_MAXONTOP);
        
        ggFigure = ggplot( cbind( dimReducData, currentFeatureExpression)[cellsOrdering, ],
                           aes_string( x = colnames( dimReducData)[1],
                                       y = colnames( dimReducData)[2],
                                       color = "currentFeatureExpression")) +
          geom_point( alpha = PLOT_DIMREDUC_EXPRESSION_ALPHA,
                      size = if(PLOT_DIMREDUC_EXPRESSION_POINTSIZE==0) Seurat:::AutoPointSize(dimReducData) else PLOT_DIMREDUC_EXPRESSION_POINTSIZE) +
          scale_color_gradient( low = "lightgrey", high = "blue") +
          labs( title = featureName) +
          theme_classic() +
          theme( #axis.text  = element_blank(),
            axis.title = element_blank(),
            #axis.ticks  = element_blank(),
            #axis.line = element_blank(),
            legend.position = "none",
            plot.margin = margin( 0, 0, 0, 0, "cm"),
            plot.title = element_text( face = "bold",
                                       size = rel( 16/14),
                                       hjust = 0.5,
                                       vjust = 1,
                                       margin = margin( b = 7)));
        
        #print( ggFigure);
        resNULL = ggsave( filename = figFilename,
                          ggFigure,
                          device = "png",
                          width = figDim[1],
                          height = figDim[2],
                          units = "in",
                          dpi = figDpi);
        
        ## Close png descriptor (twice to compensate for a bug when using png inside a rmarkdown rendering...)
        #tryCatch(invisible(dev.off()), error=function(e){});
        #tryCatch(invisible(dev.off()), error=function(e){});
        
        #include_graphics(figFilename);
        
      }
      
      
      # Violinplots (same loop as dimreduc as we don't need to separate them in different paragraphs)
      # Create subfolder for current resolution to store figure files
      message( paste0( "Modules genes '", groupName , "': violinplot"));
      
      # Create subfolder for current resolution to store figure files
      pathCurrentFigures = file.path( PATH_ANALYSIS_OUTPUT_CURRENTRES, "exportedFigures", "modules", "violinplots");
      dir.create(pathCurrentFigures, showWarnings = FALSE, recursive = TRUE);

      ## Create png file
      #invisible( png( file.path( pathCurrentFigures, paste0( featureName, ".png")), 
      #           width = figDim[1],
      #           height = figDim[2],
      #           res = figDpi));
      
      resNULL = ggsave( filename = file.path( pathCurrentFigures, paste0( groupName, ".png")),
                        violinFeatureByCluster(featureName, seuratObject = sc10x, clustersColor = clustersColor),
                        device = "png",
                        width = figDim[1],
                        height = figDim[2],
                        units = "in",
                        dpi = figDpi);
      
      ## Close png descriptor (twice to compensate for a bug when using png inside a rmarkdown rendering...)
      #tryCatch(invisible(dev.off()), error=function(e){});
      #tryCatch(invisible(dev.off()), error=function(e){});

    }))
  }  # /RENDEREXTERNALFIGURES_MODULES
  
  
  
  
  #### Finalize (save sc10x object as RDS)
  message("Saving Seurat object...");

  # Save binary file of Seurat object with clusters fur current resolution
  seuratObjectPath = file.path( PATH_ANALYSIS_OUTPUT_CURRENTRES, paste0( outputFilesPrefix_CURRENTRES, "seuratObject_allClusters.RDS"));
  saveRDS( object = sc10x, file = seuratObjectPath);



  
  #### Recursive subclustering analysis
  # Subset Seurat object by cluster and restart analysis for each of them (use
  # call to rmarkdown::render as in launcher after updating specific variables
  # used for recursion).

  # Check if recursion should be started depending on current and max levels
  if(SUBCLUSTERING_RECURSION_CURRENTLEVEL < SUBCLUSTERING_RECURSION_NBLEVELS)
  {

    # Start a recursive subclustering for each identified cluster
    invisible(mclapply(levels( Idents( sc10x)), function(currentCluster)
    {

      ### Prepare subfolders, names, and subset seurat object by current cluster 

      # Prepare a name for current cluster reference (for filenames, columns, etc...)
      currentClusterName = if(is.na(FINDCLUSTERS_RESOLUTION)) currentCluster else paste0( "C", format( currentCluster, scientific = FALSE));

      message( paste( "Preparing subclustering for:", currentClusterName));
      
      nbCellsInCluster = sum(Idents( sc10x) == currentCluster);
      message( paste( "Number of cells:", nbCellsInCluster));
      
      if( nbCellsInCluster < SUBCLUSTERING_RECURSION_MINCELLS )
      {
        message( paste0( "Stopping subclustering recursion... Check argument 'SUBCLUSTERING_RECURSION_MINCELLS' (", SUBCLUSTERING_RECURSION_MINCELLS,")."));
      } else
      {
  
        # Create a subfolder to store recursive analyses for current cluster
        PATH_ANALYSIS_OUTPUT_CURRENTRES_CURRENTCLUSTER = file.path(PATH_ANALYSIS_OUTPUT_CURRENTRES, paste0("subGroup_", currentClusterName) );
        dir.create( PATH_ANALYSIS_OUTPUT_CURRENTRES_CURRENTCLUSTER, recursive = TRUE, showWarnings = FALSE);
        
        # Update output file names
        outputFilesPrefix_CURRENTRES_CURRENTCLUSTER = paste0( outputFilesPrefix_CURRENTRES, currentClusterName, "_");
  
        # Subset seurat object (backup saved before processing) with current cluster and save it
        sc10xBackupOriginal_subsetCurrentCluster = sc10xBackupOriginal[, Idents( sc10x) == currentCluster];
        seuratObjectPath_subsetCurrentCluster = file.path( PATH_ANALYSIS_OUTPUT_CURRENTRES_CURRENTCLUSTER, paste0( outputFilesPrefix_CURRENTRES_CURRENTCLUSTER, "seuratObject_subsetCluster.RDS"));
        saveRDS( object = sc10xBackupOriginal_subsetCurrentCluster, file = seuratObjectPath_subsetCurrentCluster);
  
  
        ### Render same (recursive) report for subset only (see 'launch_reports_compilation.R')
  
        # Make a local COPY (not duplication) of environment from initial one.
        # Required for parallel processing to prevent instances modifying content
        # of other ones (see environments storage properties)
        localRenderEnvironment = new.env();
        for(n in ls(initialRenderEnvironment, all.names=TRUE)) assign(n, get(n, initialRenderEnvironment), localRenderEnvironment)
  
        # Update variables of environment (saved at start, initially built by launcher) for next level
        # -- recursive suffix to be used for report name (eventually titles, etc...)
        localRenderEnvironment[["recursiveSuffix"]] = paste0( if(exists("recursiveSuffix")) recursiveSuffix, "_", currentResolutionName, "_", currentClusterName);
        # -- output files folders
        localRenderEnvironment[["PATH_ANALYSIS_OUTPUT"]] = PATH_ANALYSIS_OUTPUT_CURRENTRES_CURRENTCLUSTER;
        # -- output file prefix (initialy created in 'globalParams.R')
        localRenderEnvironment[["outputFilesPrefix"]] = outputFilesPrefix_CURRENTRES_CURRENTCLUSTER;
        # -- increment recursion level (set to 1 if notExist/firstLevel at begin of 'Rmd' file)
        localRenderEnvironment[["SUBCLUSTERING_RECURSION_CURRENTLEVEL"]] = SUBCLUSTERING_RECURSION_CURRENTLEVEL+1;
        # -- use Seurat object subset as base fr next level 
        localRenderEnvironment[["PATH_RDS_SEURAT_OBJECT"]] = seuratObjectPath_subsetCurrentCluster;
  
        # Create an output filename for report as in launcher + level suffix (first level fixed by launcher)
        reportOutputFilename = paste0( SCIENTIFIC_PROJECT_NAME, "_",
                                       EXPERIMENT_PROJECT_NAME, "_",
                                       ANALYSIS_STEP_NAME,
                                       if(exists("recursiveSuffix")) recursiveSuffix else "", "_",
                                       currentResolutionName, "_",
                                       currentClusterName,
                                       ".html");
  
        # Create temporary copy of Rmd for safe parallel-rendering (Rmd available 
        # since execution is done in Rmd folder by default). Use clustename in as
        # base pattern because parallel execution would give same name otherwise.
        rmdCopyName = tempfile( pattern = paste0("tempRmdCopy_", currentClusterName, '_'), tmpdir = WORKING_DIR, fileext = ".Rmd");
        stopifnot( all( file.copy( from = file.path( WORKING_DIR, paste0( "Report_", ANALYSIS_STEP_NAME, ".Rmd")),
                                   to   = rmdCopyName)));
  
  
        options(knitr.duplicate.label = "allow");
  
        message( paste( "Launching recursive call to 'rmarkdown::render'"));
        
        # Generate next level report, using same environment as for current one
        # (saved at start), with updated variable for managing recursion levels
        resultPath = rmarkdown::render( input = rmdCopyName,
                                        output_dir = PATH_ANALYSIS_OUTPUT_CURRENTRES_CURRENTCLUSTER,
                                        output_file  = I( reportOutputFilename),
                                        envir = localRenderEnvironment,
                                        quiet = TRUE);
  
        # Remove temporary 'Rmd' copy
        if(! file.remove( rmdCopyName)) warning( paste0( "Temporary file '", rmdCopyName, "' could not be removed..."));
  
        # TODO: Add resulting html file as an iframe in current report ?
      }
      
    }, mc.cores = NBCORES));
  }
  
  cat(" \n \n"); # Required for '.tabset'

} # End of global loop on FINDCLUSTERS_RESOLUTION values 




cat( "\n\n<br>\n\n");


