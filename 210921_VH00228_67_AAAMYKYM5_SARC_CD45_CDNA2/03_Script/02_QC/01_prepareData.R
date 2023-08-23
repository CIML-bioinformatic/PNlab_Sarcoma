# #########################################
# This script reads and filters sc10x  data
# #########################################




# READ DATA
# ---------

## @knitr loadData

# Read input files from named vector of folders (multimodal dataset gets read as a list)
scMatrixList = Read10X( PATH_10XFILES)

# Create the seurat object (RNA only)
sc10x = CreateSeuratObject( counts = scMatrixList[["Gene Expression"]],
                            min.cells = LOAD_MIN_CELLS,
                            min.features = LOAD_MIN_FEATURES,
                            project = "");
# TO REMOVE !!!!!!
#sc10x = sc10x[,c(1:500, 13000:13500)];

# Add HTO (reindex on sc10x cell IDs to get a perfect match) to seurat object
hto_assay = CreateAssayObject( scMatrixList[["Antibody Capture"]][, colnames(sc10x)])
sc10x[["HTO"]] = hto_assay

# Attribute a numeric ID to each cell (easier than barcode)
sc10x[["numID"]] = 1:length(Cells(sc10x));


### Remove eventual NULL (empty) list elements from list of genes to be monitored
monitoredGroupEmpty = sapply( MONITORED_GENES, is.null);
if(any( monitoredGroupEmpty)) warning( paste("Following group(s) of genes to be monitored will be ignored because empty:", paste( names(monitoredGroupEmpty)[monitoredGroupEmpty], collapse=" - ")));
MONITORED_GENES = MONITORED_GENES[! monitoredGroupEmpty];

# Check whether genes in MONITORED_GENES are actually found in assay object (case mismatch is not corrected anymore to avoid genes names confusion between species)
#matchMonitoredGenes = match( toupper( unlist( MONITORED_GENES)), toupper( rownames( GetAssayData( sc10x))));
matchMonitoredGenes = match( ( unlist( MONITORED_GENES)), ( rownames( GetAssayData( sc10x))));
monitoredGenesNotFound = unique( unlist( MONITORED_GENES)[is.na( matchMonitoredGenes)]);
if(any( is.na( matchMonitoredGenes))) warning( paste( "Following gene(s) to be monitored could not be found in experimental data and will be ignored:", paste( paste0("'", monitoredGenesNotFound, "'"), collapse=" - ")));
# Replace names by names found in assay (correcting eventual case mismatch, NA for not found)
MONITORED_GENES = relist( rownames( GetAssayData( sc10x))[ matchMonitoredGenes ], skeleton = MONITORED_GENES); # Does not work with NULL list elements (removed earlier)
# Finally remove names that did not match
MONITORED_GENES = lapply( MONITORED_GENES, na.omit);


### Remove eventual NULL (empty) list elements from list of genes in modules
modulesGroupEmpty = sapply( MODULES_GENES, is.null);
if(any( modulesGroupEmpty)) warning( paste("Following module(s) of genes will be ignored because empty:", paste( names(modulesGroupEmpty)[modulesGroupEmpty], collapse=" - ")));
MODULES_GENES = MODULES_GENES[! modulesGroupEmpty];

# Check whether genes in MODULES_GENES are actually found in assay object (case mismatch is not corrected anymore to avoid genes names confusion between species)
#matchModulesGenes = match( toupper( unlist( MODULES_GENES)), toupper( rownames( GetAssayData( sc10x))));
matchModulesGenes = match( ( unlist( MODULES_GENES)), ( rownames( GetAssayData( sc10x))));
modulesGenesNotFound = unique( unlist( MODULES_GENES)[is.na( matchModulesGenes)]);
if(any( is.na( matchModulesGenes))) warning( paste( "Following gene(s) from modules list could not be found in experimental data and will be ignored:", paste( paste0("'", modulesGenesNotFound, "'"), collapse=" - ")));
# Replace names by names found in assay (correcting eventual case mismatch, NA for not found)
MODULES_GENES = relist( rownames( GetAssayData( sc10x))[ matchModulesGenes ], skeleton = MODULES_GENES); # Does not work with NULL list elements (removed earlier)
# Finally remove names that did not match
MODULES_GENES = lapply( MODULES_GENES, na.omit);


### Transfer genes in very small modules (MODULES_GENES) to be analyzed individually (MONITORED_GENES)
modulesToTransfer = sapply(MODULES_GENES, length) < MONITORED_GENES_SMALL_MODULE
if(any(modulesToTransfer))
{
    warning( paste0( "Following Module(s) contained very few genes (<", MONITORED_GENES_SMALL_MODULE, "). These genes were transfered to 'Monitored genes' to be analyzed individually: ", paste( names(modulesToTransfer)[modulesToTransfer], collapse = " - ")));
    # Alter name so they can be recognized in list of Monitored genes
    names( MODULES_GENES)[modulesToTransfer] = paste0( "MOD_", names( MODULES_GENES)[modulesToTransfer]);
    MONITORED_GENES = c( MONITORED_GENES, MODULES_GENES[modulesToTransfer]);
    MODULES_GENES = MODULES_GENES[!modulesToTransfer];
}





# HTOs
# ----

## @knitr demuxHTO

# Normalize HTO counts
sc10x = NormalizeData(sc10x, assay = "HTO", normalization.method = "CLR")

# Demultiplex cells based on HTOs
sc10x = HTODemux(sc10x, assay = "HTO", positive.quantile = 0.99)


# Group cells based on the max HTO signal and show ridgeplot
Idents(sc10x) = "HTO_maxID"
RidgePlot(sc10x, assay = "HTO", features = rownames(sc10x[["HTO"]]), ncol = 2)

kable(table(sc10x[["HTO_classification.global"]]), col.names = NULL)

# Plot UMIs/genes distribution for singlets doublets and negative cells
#Idents(sc10x) <- "HTO_classification.global"
#VlnPlot(sc10x, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
#VlnPlot(sc10x, features = "nFeature_RNA", pt.size = 0.1, log = TRUE)

#Idents(sc10x) = "HTO_maxID"
#VlnPlot(sc10x, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
#Idents(sc10x) = "HTO_classification"
#VlnPlot(sc10x, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# Calculate tSNE embedding of the HTO data
DefaultAssay(sc10x) = "HTO"
sc10x = ScaleData(sc10x, features = rownames(sc10x), verbose = FALSE)
sc10x = RunPCA(sc10x, features = rownames(sc10x), approx = FALSE, verbose = FALSE)
sc10x = RunTSNE(sc10x, dims = 1:20, check_duplicates = FALSE)


## @knitr demuxHTO_tsne

# Plot showing singlet/doublet/negative cells
Idents(sc10x) <- "HTO_classification.global"
DimPlot(sc10x, reduction = "tsne")

# Show affectation (with all doublets combination)
Idents(sc10x) <- "HTO_classification"
DimPlot(sc10x, reduction = "tsne")


## @knitr demuxHTO_tsne_final

# Remove negative and doublet "cells" from dataset
Idents(sc10x) <- "HTO_classification.global"
sc10x <- subset(sc10x, idents = c("Negative", "Doublet"), invert = TRUE)

# Replot final HTO affectation
Idents(sc10x) <- "HTO_classification"
DimPlot(sc10x, reduction = "tsne")

# Reset default assay to RNA
DefaultAssay(sc10x) = "RNA"




# FILTER DATA
# -----------

## @knitr filterData_selection

### Identify mitocondrial genes in matrix
mito.genes = grep( pattern = "^mt-", x = rownames(GetAssayData(object = sc10x, slot = "counts")), value = TRUE, ignore.case = TRUE)
if(length(mito.genes)==0) 
{
  warning( "No mitochondrial genes could be identified in this dataset.");
} else 
{
  # Compute percentage of mitochondrial transcripts in each cell
  percent.mito <- PercentageFeatureSet(sc10x, features=mito.genes)
  # Add the mitocondrial gene percentage as meta information in the Seurat object
  sc10x[["percent.mito"]] <- percent.mito
}

### Identify ribosomal genes in matrix
ribo.genes = grep(pattern = "^rps|^rpl",  x = rownames(GetAssayData(object = sc10x, slot = "counts")), value = TRUE, ignore.case=TRUE)
if(length(ribo.genes)==0) 
{
  warning( "No ribosomal genes could be identified in this dataset.");
} else 
{
  # Compute percentage of ribosomal transcripts in each cell
  percent.ribo <- PercentageFeatureSet(sc10x, features=ribo.genes)
  # Add the ribosomal gene percentage as meta information in the Seurat object
  sc10x[["percent.ribo"]] <- percent.ribo
}


### Identify cells that will be rejected based on specified thresholds

# Reject cells based on UMI numbers
nUMI.drop = logical( length(Cells(sc10x)));
if( ! is.null( FILTER_UMI_MIN))
{
  nUMI.drop = nUMI.drop | (sc10x[["nCount_RNA", drop=TRUE]] < FILTER_UMI_MIN);
}
if( ! is.null( FILTER_UMI_MAX))
{
  nUMI.drop = nUMI.drop | (sc10x[["nCount_RNA", drop=TRUE]] > FILTER_UMI_MAX);
}

# Reject cells based on number of expressed genes
nGene.drop = logical( length(Cells(sc10x)));
if( ! is.null( FILTER_FEATURE_MIN))
{
  nGene.drop = nGene.drop | (sc10x[["nFeature_RNA", drop=TRUE]] < FILTER_FEATURE_MIN);
}
if( ! is.null( FILTER_FEATURE_MAX))
{
  nGene.drop = nGene.drop | (sc10x[["nFeature_RNA", drop=TRUE]] > FILTER_FEATURE_MAX);
}

# Identify cells with high percentage of mitocondrial genes
mito.drop = logical( length(Cells(sc10x)));
if( length(mito.genes) && (! is.null(FILTER_MITOPCT_MAX)))
{
  mito.drop = (sc10x[["percent.mito", drop=TRUE]] > FILTER_MITOPCT_MAX);
}

# Identify cells with low percentage of ribosomal genes
ribo.drop = logical( length(Cells(sc10x)));
if( length(ribo.genes) && (! is.null(FILTER_RIBOPCT_MIN)))
{
  ribo.drop = (sc10x[["percent.ribo", drop=TRUE]] < FILTER_RIBOPCT_MIN);
}


### Plot distributions of #UMIs, #Genes, %Mito, and %Ribo among cells

# The metrics used for filtering will be plotted using ggplot rasterized graphs when having to show a lot of points (instead of plotly interactive graphs)
if(length( Cells( sc10x)) < PLOT_RASTER_NBCELLS_THRESHOLD)
{
  # Do interactive plot using plotly (can cause overhead when viewing many point)

  # Gather data to be visualized together (cell name + numID + metrics)
  cellsData = cbind( "Cell" = colnames( sc10x), # cell names from rownames conflict with search field in datatable
                     sc10x[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                     "percent.mito" = if(length( mito.genes)) as.numeric(format(sc10x[["percent.mito", drop = TRUE]], digits = 5)) else NULL,
                     "percent.ribo" = if(length( ribo.genes)) as.numeric(format(sc10x[["percent.ribo", drop = TRUE]], digits = 5)) else NULL);

  # Create text to show under cursor for each cell
  hoverText = do.call(paste, c(Map( paste, 
                                    c( "", 
                                       "Cell ID: ", 
                                       "# UMIs: ", 
                                       "# Genes: ", 
                                       if(length( mito.genes)) "% Mito: ", 
                                       if(length( ribo.genes)) "% Ribo: "), 
                                    cellsData, 
                                    sep=""), 
                      sep="\n"));

  panelWidth =  190; # 800 when using subplot
  panelHeight = 400;
  # Generate plotly violin/jitter panels for #umis, #genes, %mitochondrial, and %ribosomal stats
  lypanel_umis  = plotViolinJitter( cbind( cellsData, outliers = nUMI.drop), 
                                    xAxisFormula = ~as.numeric(1), 
                                    yAxisFormula = ~nCount_RNA, 
                                    pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), 
                                    hoverText = hoverText, 
                                    xTicklabels = "# UMIs", 
                                    thresholdHigh = FILTER_UMI_MAX, 
                                    thresholdLow = FILTER_UMI_MIN, 
                                    panelWidth = panelWidth, 
                                    panelHeight = panelHeight);

  lypanel_genes = plotViolinJitter( cbind( cellsData, outliers = nGene.drop), 
                                    xAxisFormula = ~as.numeric(1), 
                                    yAxisFormula = ~nFeature_RNA, 
                                    pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), 
                                    hoverText = hoverText, 
                                    xTicklabels = "# Genes", 
                                    thresholdHigh = FILTER_FEATURE_MAX, 
                                    thresholdLow = FILTER_FEATURE_MIN, 
                                    panelWidth = panelWidth, 
                                    panelHeight = panelHeight);

  lypanel_mitos = if(length(mito.genes)) plotViolinJitter( cbind( cellsData, outliers = mito.drop), 
                                                           xAxisFormula = ~as.numeric(1), 
                                                           yAxisFormula = ~percent.mito, 
                                                           pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), 
                                                           hoverText = hoverText, 
                                                           xTicklabels = "% Mito", 
                                                           thresholdHigh = FILTER_MITOPCT_MAX, 
                                                           panelWidth = panelWidth, 
                                                           panelHeight = panelHeight) else NULL;


  lypanel_ribos = if(length(ribo.genes)) plotViolinJitter( cbind( cellsData, outliers = ribo.drop), 
                                                           xAxisFormula = ~as.numeric(1), 
                                                           yAxisFormula = ~percent.ribo, 
                                                           pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), 
                                                           hoverText = hoverText, 
                                                           xTicklabels = "% Ribo", 
                                                           thresholdLow = FILTER_RIBOPCT_MIN, 
                                                           panelWidth = panelWidth, 
                                                           panelHeight = panelHeight) else NULL;

  # Set panels as a list and define plotly config
  panelsList = list( lypanel_umis, lypanel_genes, lypanel_mitos, lypanel_ribos);
  panelsList = lapply(panelsList, config, displaylogo = FALSE, 
                      toImageButtonOptions = list(format='svg'), 
                      modeBarButtons = list(list('toImage'), 
                                            list('zoom2d', 'pan2d', 'resetScale2d')));

  # # Group plotly violin/jitter panels (for sizing, it uses layout of one of the plot)
  # plotPanels = layout( subplot( lypanel_umis,
  #                               lypanel_genes,
  #                               lypanel_mitos),
  #                      showlegend = FALSE) # Remove eventual legends (does not mix well with subplot)

  # Control layout using flex because subplot is limited regarding plot sizing and alignment
  # 'browsable' required in console, not in script/document
  #browsable(
  div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: green;"),
           div( style = paste("display: flex; flex-flow: row nowrap; overflow-x: auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
                lapply(panelsList, div, style = paste("flex : 0 0 auto; margin: 5px;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;"))))
  #)

} else
{
  # Do the plots as simple images with ggplot when having a lot of points

  # Generate plotly violin/jitter panels for #umis, #genes, %mitochondrial, and %ribosomal stats
  ggpanel_umis  = ggplot( cbind(sc10x[["nCount_RNA"]], outliers = nUMI.drop), aes( y = nCount_RNA)) + 
                    geom_jitter( aes(x = 1.5, col = outliers), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
                    geom_violin( aes(x = 0.5), width = 0.9, fill = "darkblue", col = "#AAAAAA", alpha = 0.3) + 
                    geom_hline( yintercept = c( FILTER_UMI_MIN, FILTER_UMI_MAX), 
                                color = c( "blue", "red")[!sapply( list( FILTER_UMI_MIN, FILTER_UMI_MAX), is.null)], 
                                alpha = 0.5,
                                size = 1) +
                    labs(x = "# UMIs", y = "") +
                    theme( axis.text.x = element_blank(), 
                           axis.ticks.x = element_blank(),
                           axis.title.y = element_blank(),
                           legend.position = "none") +
                    scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#444444"));

  ggpanel_genes = ggplot( cbind(sc10x[["nFeature_RNA"]], outliers = nGene.drop), aes( y = nFeature_RNA)) + 
                    geom_jitter( aes(x = 1.5, col = outliers), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
                    geom_violin( aes(x = 0.5), width = 0.9, fill = "darkblue", col = "#AAAAAA", alpha = 0.3) + 
                    geom_hline( yintercept = c( FILTER_FEATURE_MIN, FILTER_FEATURE_MAX), 
                                color = c( "blue", "red")[!sapply( list( FILTER_FEATURE_MIN, FILTER_FEATURE_MAX), is.null)], 
                                alpha = 0.5,
                                size = 1) +
                    labs(x = "# Genes", y = "") +
                    theme( axis.text.x = element_blank(), 
                           axis.ticks.x = element_blank(),
                           axis.title.y = element_blank(),
                           legend.position = "none") +
                    scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#444444"));

  ggpanel_mitos = if(length(mito.genes)) ggplot( cbind(sc10x[["percent.mito"]], outliers = mito.drop), aes( y = percent.mito)) + 
                                           geom_jitter( aes(x = 1.5, col = outliers), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
                                           geom_violin( aes(x = 0.5), width = 0.9, fill = "darkblue", col = "#AAAAAA", alpha = 0.3) + 
                                           geom_hline( yintercept = FILTER_MITOPCT_MAX, 
                                                       color = "red", 
                                                       alpha = 0.5,
                                                       size = 1) +
                                           labs(x = "% Mito", y = "") +
                                           theme( axis.text.x = element_blank(), 
                                                  axis.ticks.x = element_blank(),
                                                  axis.title.y = element_blank(),
                                                  legend.position = "none") +
                                           scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#444444")) else NULL;


  ggpanel_ribos = if(length(ribo.genes)) ggplot( cbind(sc10x[["percent.ribo"]], outliers = ribo.drop), aes( y = percent.ribo)) + 
                                           geom_jitter( aes(x = 1.5, col = outliers), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
                                           geom_violin( aes(x = 0.5), width = 0.9, fill = "darkblue", col = "#AAAAAA", alpha = 0.3) + 
                                           geom_hline( yintercept = FILTER_RIBOPCT_MIN, 
                                                       color = "blue", 
                                                       alpha = 0.5,
                                                       size = 1) +
                                           labs(x = "% Ribo", y = "") +
                                           theme( axis.text.x = element_blank(), 
                                                  axis.ticks.x = element_blank(),
                                                  axis.title.y = element_blank(),
                                                  legend.position = "none") +
                                           scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#444444")) else NULL;

  # Use patchwork library to assemble panels
  print( ggpanel_umis + ggpanel_genes + ggpanel_mitos + ggpanel_ribos + plot_layout( nrow = 1));
  
}


cat( "<br>Number of cells removed based on number of UMIs:", sum( nUMI.drop));
cat( "<br>Number of cells removed based on number of genes:", sum( nGene.drop));
if(exists( "mito.drop")) cat( "<br>Number of cells removed based on high percentage of mitochondrial transcripts:", sum( mito.drop));
if(exists( "ribo.drop")) cat( "<br>Number of cells removed based on low percentage of ribosomal transcripts:", sum( ribo.drop));
cat( "\n<br>\n");

# Identify cells to exclude as union of cells with low nb UMI, low nb expressed genes, high pct mito genes, low pct ribo genes
sc10x[["outlier"]] = nUMI.drop  | 
                     nGene.drop | 
                     ( if(exists( "mito.drop")) mito.drop else FALSE ) | 
                     ( if(exists( "ribo.drop")) ribo.drop else FALSE );

cat("<br><br>Removed cells after filters:", sum( unlist(sc10x[["outlier"]] )));
cat("<br>Remaining cells after filters:", sum( ! unlist(sc10x[["outlier"]] )));
cat("\n<br>\n");

### Record which cells got rejected

# Export the excluded cells to file
write.table( data.frame( cellid = names(which(nUMI.drop))), 
             file= file.path( PATH_ANALYSIS_OUTPUT, 
                              paste0( outputFilesPrefix, "excluded_cells_nbUMI.txt")), 
             quote = FALSE, 
             row.names = FALSE, 
             col.names = TRUE, 
             sep="\t");

write.table( data.frame( cellid = names(which(nGene.drop))), 
             file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "excluded_cells_nbGene.txt")), 
             quote = FALSE, 
             row.names = FALSE, 
             col.names = TRUE, 
             sep="\t");

if(exists( "mito.drop"))
{
  write.table( data.frame( cellid = names(which(mito.drop))), 
               file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "excluded_cells_pctMito.txt")), 
               quote = FALSE, 
               row.names = FALSE, 
               col.names = TRUE, 
               sep="\t");
}

if(exists( "ribo.drop"))
{
  write.table( data.frame( cellid = names(which(ribo.drop))), 
               file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "excluded_cells_pctRibo.txt")), 
               quote = FALSE, 
               row.names = FALSE, 
               col.names = TRUE, 
               sep="\t");
}




## @knitr filterData_summaryPlot

### Plot dispersion of excluded and non-excluded cells

# number of genes and number of UMIs
ggplot( sc10x[[]][order( sc10x[["outlier"]]),], # Plot FALSE first and TRUE after
        aes( x = nFeature_RNA, 
             y = nCount_RNA, 
             color = outlier)) + 
  geom_point( size = 0.5) +
  geom_vline( xintercept = c( FILTER_FEATURE_MIN, FILTER_FEATURE_MAX), 
              linetype = 2, 
              color = c( "blue", "red")[!sapply( list( FILTER_FEATURE_MIN, FILTER_FEATURE_MAX), is.null)], 
              alpha = 0.5) +
  geom_hline( yintercept = c( FILTER_UMI_MIN, FILTER_UMI_MAX), 
              linetype = 2, 
              color = c( "blue", "red")[!sapply( list( FILTER_UMI_MIN, FILTER_UMI_MAX), is.null)], 
              alpha = 0.5) +
  scale_color_manual( values = c( "TRUE" = "#FF0000AA", "FALSE" = "#44444444")) + 
  labs( x = "# Genes", y = "# UMIs") +
  theme( legend.position = "none")

# Mitochondrial vs ribosomal distributions
if(exists( "mito.drop") && exists( "ribo.drop"))
{
  ggplot( sc10x[[]][order( sc10x[["outlier"]]),], # Plot FALSE first and TRUE after
          aes( x = percent.ribo, 
               y = percent.mito, 
               color = outlier)) + 
    geom_point( size = 0.5) +
    geom_vline( xintercept = FILTER_RIBOPCT_MIN, linetype = 2, color = "blue", alpha = 0.5) +
    geom_hline( yintercept = FILTER_MITOPCT_MAX, linetype = 2, color = "red", alpha = 0.5) +
    scale_color_manual( values = c( "TRUE" = "#FF0000AA", "FALSE" = "#44444444")) + 
    labs( x = "% Ribosomal genes", y = "% Mitochondrial genes") +
    theme( legend.position = "none")
}




## @knitr filterData_filterObject

# Filter the excluded cells in the Seurat object
sc10x_nonFiltered = sc10x;
sc10x = sc10x[ , ! sc10x[[ "outlier", drop=TRUE ]] ];

# Save the list of remaining cells after selection during loading, and #Genes, #UMIs, pctMito, pctRibo
write.table( data.frame( cellid = Cells(sc10x)), 
             file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "selected_cells.txt")), 
             quote = FALSE, 
             row.names = FALSE, 
             col.names = TRUE, 
             sep="\t");




# NORMALIZE DATA
# --------------

## @knitr normalizeData

sc10x = NormalizeData( object = sc10x,
                       normalization.method = DATA_NORM_METHOD,
                       scale.factor = DATA_NORM_SCALEFACTOR,
                       verbose = .VERBOSE);

sc10x = ScaleData( object    = sc10x,
                   do.center = DATA_CENTER,
                   do.scale  = DATA_SCALE,
                   vars.to.regress = DATA_VARS_REGRESS,
                   verbose = .VERBOSE);


