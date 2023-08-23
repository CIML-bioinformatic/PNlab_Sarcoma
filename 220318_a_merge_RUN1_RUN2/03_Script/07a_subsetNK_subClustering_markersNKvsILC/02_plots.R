# ##############################################################################
# This script represents cells on a 2D plot (e.g. dimensionality reduction), and
# uses javascript link between report objects to select cells in plots and show
# their characteristics on attached graphs or tables. 
# 
# ##############################################################################


## @knitr reClustering
#### Redo clustering on this selected population (replace eventual previous identity)

nbPC_findclusters=FINDCLUSTERS_USE_PCA_NBDIMS
if(FINDCLUSTERS_USE_PCA_NBDIMS>nbPC)
{
  warning( paste0( "Number of computed PCs  (", nbPC, ") smaller than requested PCs for 'findclusters' (", FINDCLUSTERS_USE_PCA_NBDIMS,"), setting lower PC number (", nbPC, ")..." ))
  nbPC_findclusters = nbPC
}  
sc10x <- FindNeighbors(object    = sc10x,
                       reduction = "pca",
                       dims      = 1:nbPC_findclusters);

sc10x <- FindClusters(object             = sc10x,
                      resolution         = 0.5,
                      algorithm          = FINDCLUSTERS_ALGORITHM,
                      temp.file.location = "/tmp/",
                      verbose            = .VERBOSE);

# Define a set of colors for clusters (based on ggplot default)
clustersColor = hue_pal()( nlevels( Idents( sc10x)));
names( clustersColor) = levels( Idents( sc10x));




## @knitr plotsPreparation
#### Prepare SharedData object and variables for plots 

# For each cell, gather all data to be used for plots generation
cellsData = cbind( "Cell" = colnames( sc10x), # cell names from rownames conflict with search field in datatable
                    sc10x[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                    "percent.mito" = if("percent.mito" %in% colnames(sc10x[[]])) as.numeric(format(sc10x[["percent.mito", drop = TRUE]], digits = 5)) else NULL,
                    "percent.ribo" = if("percent.ribo" %in% colnames(sc10x[[]])) as.numeric(format(sc10x[["percent.ribo", drop = TRUE]], digits = 5)) else NULL,
                    "Cluster" = Idents( sc10x),
                    "HTO" = sc10x[["factorHTO", drop = TRUE]],
                    cellsCoordinates,
                    jitterCoords = jitter(rep(0, ncol( sc10x)), amount = 0.15) + 0.17); # Precompute jitter for eventual dotplots so coordinates are not recomputed at each refresh

# Create text to show under cursor for each cell for plotly hover
hoverText = do.call(paste, c(Map( paste,
                                  c( "",
                                      "Cell ID: ",
                                      "# UMIs: ",
                                      "# Genes: ",
                                      if("percent.mito" %in% colnames(sc10x[[]])) "% Mito: ",
                                      if("percent.ribo" %in% colnames(sc10x[[]])) "% Ribo: ",
                                      "Cluster: ",
                                      "HTO: "),
                                  cellsData[1:(ncol(cellsData)-3)], # Do not include coordinates in hover text
                                  sep = ""),
                    sep = "\n"));

# Create SharedData object (using plotly accessor) for use by all plots/widgets
dataShare = highlight_key(cellsData);

# Compute the median coordinates for each cluster to plot corresponding name 
dimReducLabels = as.data.frame( do.call( rbind, 
                                         by( cellsData, 
                                             cellsData[["Cluster"]], 
                                             function(x)
                   { 
                     return( data.frame( lapply( x[colnames( cellsCoordinates)], median), 
                                         "Cluster" = x[1, "Cluster"]));
                   })));




## @knitr plotClustersSize
#### Show number of cells in each cluster

# Compute the number of cells by cluster
#clustersCount = as.data.frame( table( Cluster = Idents(sc10x)), responseName = "Total");

# Compute the cluster counts sliced by HTO condition
clustersCountByHTO = table(cbind(Idents( sc10x), sc10x[["factorHTO"]]))

# Add it to the report
print( knitr::kable( addmargins(clustersCountByHTO),
                    align   = "c",
                    caption = "Number of cells by population (Cluster) and condition (HTO)")) #%>% kable_styling(bootstrap_options = c("condensed"), full_width = FALSE);

# Get a table of proportion (for selected margin)
#print( knitr::kable( addmargins( prop.table( clustersCountByHTO, margin = 1) * 100,
#                                 margin = 2),
#                     align   = "c",
#                     caption = "Proportion (pct) of cells in conditions (HTO) for each population (Cluster)")) #%>% kable_styling(bootstrap_options = c("condensed"), full_width = FALSE);




# Table of counts sliced by sample/replicate and HTO
#kable(table(sc10x[[c("orig.ident", "factorHTO")]]))

# Compute a matrix of average expression value by cluster (for each gene)
geneExpByCluster = do.call( rbind, 
                            apply( as.matrix( GetAssayData( sc10x, 
                                                            slot = "data", 
                                                            assay="RNA")),  # normalized expression values
                                   1,                                       # by rows
                                   tapply,                                  # apply by group
                                   INDEX = Idents( sc10x),                  # clusters IDs
                                   mean,                                    # summary function
                                   simplify = FALSE));                      # do not auto coerce

# Save it as 'tsv' file
write.table( geneExpByCluster, 
             file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "normExpressionByCluster.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");

# Compute a matrix of average expression value by cluster and HTOs (for each gene)
geneExpByClusterAndHTO = do.call( rbind, 
                                  apply( as.matrix( GetAssayData( sc10x, 
                                                                  slot = "data", 
                                                                  assay="RNA")),  # normalized expression values
                                         1,                                       # by rows
                                         tapply,                                  # apply by group
                                         INDEX = paste( Idents( sc10x),
                                                        sc10x[["factorHTO", drop = TRUE]], 
                                                        sep = "_"),               # clusters + HTOs IDs
                                         mean,                                    # summary function
                                         simplify = FALSE));                      # do not auto coerce

# Save it as 'tsv' file
write.table( geneExpByClusterAndHTO, 
             file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "normExpressionByClusterAndHTO.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");




## @knitr chisquare_test_clusters_HTO
# Test independance between rows and columns (clusters and sample condition/HTO)
# to detect counts over expected values under random condition.
chisq = chisq.test(clustersCountByHTO)

# Show the expected counts considering margins statistics, to contrast with
# observed counts shown above
print( knitr::kable( addmargins(round(chisq[["expected"]],2)),
                     align   = "c",
                     caption = "Expected number of cells by population (Cluster) and condition (HTO) from X².\nTo compare with observed counts.")) #%>% kable_styling(bootstrap_options = c("condensed"), full_width = FALSE);

# Show result as simple table
pander(chisq)


## @knitr chisquare_test_clusters_HTO_plotResidualsAndContributions
# Set a layout to combine pearson residuals and contribution plots
def.par <- par(no.readonly = TRUE) # save default, for resetting...
layout(matrix(c(1,2), nrow=1, byrow = TRUE)) # Two next figures on same row

# Plot Pearson residuals
corrplot(chisq$residuals, 
         is.cor=FALSE, 
         addCoef.col="darkgrey", 
         tl.srt=60, 
         cl.ratio=0.35, 
         #col = rev(colorRampPalette(RColorBrewer::brewer.pal(name="RdBu", 11))(100)),
         mar=c(0,0,2,0), 
         title = "Pearson residuals")

# Plot contributions (see: http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r)
contrib = round(100*(chisq$residuals^2)/chisq$statistic,3)
corrplot(contrib, 
         is.cor=FALSE, 
         addCoef.col="darkgrey", 
         tl.srt=60, 
         cl.ratio=0.35, 
         col = colorRampPalette(c("white", "yellow", "darkorange", "darkred"))(100),
         mar=c(0,0,2,0), 
         title = "Contribution to X² statistic (%)")

# Reset layout to previously saved default
par(def.par)


## @knitr chisquare_test_clusters_HTO_plotMosaic
# Plot as mosaic
mosaicplot(clustersCountByHTO, shade = TRUE, las = 2, main = "X² residuals", xlab="", ylab="")




## @knitr plotDimReducInteractive_colorClusters
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables

#### Dimensionality reduction 2D (plotly)
# Parameters defined in chunk
# dimReducWidth  = 800;
# dimReducHeight = 800;
# initialPointSize = 5;

plotDimReduc = plot_ly( dataShare, # 'SharedData' object for plots interactions
                        x = as.formula( paste( "~", colnames( cellsCoordinates)[1])),
                        y = as.formula( paste( "~", colnames( cellsCoordinates)[2])),
                        #                       name = ~paste("Cluster", Cluster),
                        color = ~Cluster,
                        colors = clustersColor,
                        #symbol = ~HTO,
                        width = dimReducWidth,
                        height = dimReducHeight) %>%
  add_trace( type = "scattergl",
             mode = "markers",
             marker = list( size = initialPointSize,
                            opacity = 0.6),
             text = hoverText,
             hoverinfo = "text+name", 
             showlegend = if(exists( "showLegend")) showLegend else TRUE)

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  plotDimReduc = plotDimReduc %>% 
    # Use precomputed median of coordinates for each group as new data for Cluster text annotations (type='scatter', mode='text')
    add_trace( data = dimReducLabels,
               x = as.formula( paste( "~", colnames( cellsCoordinates)[1])),
               y = as.formula( paste( "~", colnames( cellsCoordinates)[2])),
               type = "scattergl",
               mode = "text",
               symbol = NULL,
               text = ~Cluster,
               textfont = list( color = '#000000', size = 20),
               hoverinfo = 'skip',
               showlegend = FALSE)
}

plotDimReduc = plotDimReduc %>% 
  #  hide_colorbar() %>%
  #  hide_legend() %>%
  config( displaylogo = FALSE,
          toImageButtonOptions = list( format='svg'),
          modeBarButtons = list(list( 'toImage'),
                                list( 'zoom2d', 
                                      'zoomIn2d', 
                                      'zoomOut2d', 
                                      'pan2d', 
                                      'resetScale2d')))

# Set a unique html ID to access plot with html/javascript for external slider
figID = if(exists("figID")) figID+1 else 0;
plotDimReduc[["elementId"]] = paste0("plotDimReduc_withSlider", figID);

# Render plot
(div(htmltools::tagList(plotDimReduc)));

# Create a html slider using javascript for controlling point size
cat( paste0( '\n',
             '<p>Point size: <span id="valueTag_', figID, '">', initialPointSize, '</span></p>',
             '<input type="range" min="2" max="50" value="', initialPointSize, '" id="slider_', figID, '" style="width: ', dimReducWidth, 'px;">',
             '<script> document.getElementById("slider_', figID, '").oninput = function() { Plotly.restyle("', plotDimReduc[["elementId"]], '", {"marker.size": this.value}); document.getElementById("valueTag_', figID, '").innerHTML = this.value; } </script>',
             '\n'));

# Create a html/js block to change color of clusters and get selected values
cat( paste0( '\n', 
             'Color:&nbsp;',
             # Create a selector, each option corresponding to a plotly trace (defined here by color)
             '<select id="colorChange_traceSelector_plotly', figID, '" style="width: 200px;">', 
             paste('<option value="', c( 'undefined', 0:(nlevels( cellsData[["Cluster"]])-1)), '">', c('All', levels( cellsData[["Cluster"]])), '</option>', sep ="", collapse = " " ), 
             '</select>', 
             '&nbsp;&nbsp;',
             # Create HTML anchor for color picker
             '<input type="text" id="colorChange_colorPicker_plotly', figID, '" />', 
             # Call javascript that transforms it in an actual color picker (updates plotly color of selected trace on change, and resets its own value to 'clear')
             '<script> $("#colorChange_colorPicker_plotly', figID, '").spectrum({ allowEmpty: true, showInput: true, preferredFormat: "rgb", clickoutFiresChange: false, showPalette: true, showSelectionPalette: true, palette: [ ', paste('"', clustersColor, '"', collapse = ", ", sep = ""), ' ], localStorageKey: "myLocalPalette", change: function(color) { if(color) { Plotly.restyle("', plotDimReduc[["elementId"]], '", {"marker.color": color.toHexString()}, $("#colorChange_traceSelector_plotly', figID, '").val()); $("#colorChange_colorPicker_plotly', figID, '").spectrum("set", ""); } } }); </script>', 
             '&nbsp;&nbsp;Get current colors:&nbsp;',
             # Add button to extract current colors from plotly figure (hexadecimal format)
             '<button id="colorGetter_alertButton_plotly', figID, '" onclick="getColorTab(', plotDimReduc[["elementId"]], ', \'scattergl\', \'hex\')">Hexadecimal</button>',
             # Add button to extract current colors from plotly figure (RGB format)
             '<button id="colorGetter_alertButton_plotly', figID, '" onclick="getColorTab(', plotDimReduc[["elementId"]], ', \'scattergl\', \'rgb\')">RGB</button>',
             '\n'));
cat('\n \n')


## @knitr plotDimReducRaster_colorClusters
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables

# Prepare dimreduc figure as ggplot
ggFigure = ggplot(  cellsData, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2],
                                color = "Cluster")) +
  geom_point( alpha = PLOT_DIMREDUC_GROUPS_ALPHA, 
              size  = PLOT_DIMREDUC_GROUPS_POINTSIZE)

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  ggFigure = ggFigure +
    geom_text( data = dimReducLabels,
               aes( label = Cluster),
               color = "black",
               size = 5)
}

ggFigure = ggFigure +
  scale_color_manual( name = "Cluster", values = clustersColor) +
  theme_classic() +
  theme( #axis.text  = element_blank(),
    axis.title = element_blank(),
    #axis.ticks  = element_blank(),
    #axis.line = element_blank(),
    legend.position = if(if(exists( "showLegend")) showLegend else TRUE) NULL else "none",
    plot.margin = margin( 0, 0, 0, 0, "cm"),
    plot.title = element_text( face = "bold",
                               size = rel( 1.8), #rel( 16/14),
                               hjust = 0.5,
                               vjust = 1,
                               margin = margin( b = 7)));
print( ggFigure);
cat('\n \n')




## @knitr plotDimReducInteractive_colorClusters_symbolHTOs
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables

#### Dimensionality reduction 2D (plotly)
# Parameters defined in chunk
# dimReducWidth  = 800;
# dimReducHeight = 800;
# initialPointSize = 5;

plotDimReduc = plot_ly( dataShare, # 'SharedData' object for plots interactions
                        x = as.formula( paste( "~", colnames( cellsCoordinates)[1])),
                        y = as.formula( paste( "~", colnames( cellsCoordinates)[2])),
                        #                       name = ~paste("Cluster", Cluster),
                        color = ~Cluster,
                        colors = clustersColor,
                        symbol = ~HTO,
                        width = dimReducWidth,
                        height = dimReducHeight) %>%
  add_trace( type = "scattergl",
             mode = "markers",
             marker = list( size = initialPointSize,
                            opacity = 0.6),
             text = hoverText,
             hoverinfo = "text+name", 
             showlegend = if(exists( "showLegend")) showLegend else TRUE)

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  plotDimReduc = plotDimReduc %>% 
    # Use precomputed median of coordinates for each group as new data for Cluster text annotations (type='scatter', mode='text')
    add_trace( data = dimReducLabels,
               x = as.formula( paste( "~", colnames( cellsCoordinates)[1])),
               y = as.formula( paste( "~", colnames( cellsCoordinates)[2])),
               type = "scattergl",
               mode = "text",
               symbol = NULL,
               text = ~Cluster,
               textfont = list( color = '#000000', size = 20),
               hoverinfo = 'skip',
               showlegend = FALSE)
}

plotDimReduc = plotDimReduc %>% 
  #  hide_colorbar() %>%
  #  hide_legend() %>%
  config( displaylogo = FALSE,
          toImageButtonOptions = list( format='svg'),
          modeBarButtons = list(list( 'toImage'),
                                list( 'zoom2d', 
                                      'zoomIn2d', 
                                      'zoomOut2d', 
                                      'pan2d', 
                                      'resetScale2d')))

# Set a unique html ID to access plot with html/javascript for external slider
figID = if(exists("figID")) figID+1 else 0;
plotDimReduc[["elementId"]] = paste0("plotDimReduc_withSlider", figID);

# Render plot
(div(htmltools::tagList(plotDimReduc)));

# Create a html slider using javascript for controlling point size
cat( paste0( '\n',
             '<p>Point size: <span id="valueTag_', figID, '">', initialPointSize, '</span></p>',
             '<input type="range" min="2" max="50" value="', initialPointSize, '" id="slider_', figID, '" style="width: ', dimReducWidth, 'px;">',
             '<script> document.getElementById("slider_', figID, '").oninput = function() { Plotly.restyle("', plotDimReduc[["elementId"]], '", {"marker.size": this.value}); document.getElementById("valueTag_', figID, '").innerHTML = this.value; } </script>',
             '\n'));
cat('\n \n')


## @knitr plotDimReducRaster_colorClusters_symbolHTOs
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables

# Prepare dimreduc figure as ggplot
ggFigure = ggplot(  cellsData, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2],
                                color = "Cluster",
                                pch = "HTO")) +
  geom_point( alpha = PLOT_DIMREDUC_GROUPS_ALPHA, 
              size  = PLOT_DIMREDUC_GROUPS_POINTSIZE)

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  ggFigure = ggFigure +
    geom_text( data = dimReducLabels,
               aes( label = Cluster, pch = NULL),
               color = "black",
               size = 3.5)
}

ggFigure = ggFigure +
  scale_color_manual( name = "Cluster", values = clustersColor) +
  theme_classic() +
  theme( #axis.text  = element_blank(),
    axis.title = element_blank(),
    #axis.ticks  = element_blank(),
    #axis.line = element_blank(),
    legend.position = if(if(exists( "showLegend")) showLegend else TRUE) NULL else "none",
    plot.margin = margin( 0, 0, 0, 0, "cm"),
    plot.title = element_text( face = "bold",
                               size = rel( 1.8), #rel( 16/14),
                               hjust = 0.5,
                               vjust = 1,
                               margin = margin( b = 7)));
print( ggFigure);
cat('\n \n')




## @knitr plotDimReducInteractiveGG_colorClusters_facetHTOs
# Using 'ggplotly' as facetting is not trivial with plotly 
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables

#### Facetted dimensionality reduction 2D (ggplot + plotly conversion)
# Parameters defined in chunk
# dimReducWidth  = 800;
# dimReducHeight = 800;
# initialPointSize = 5;

# Prepare dimreduc figure as ggplot
ggFigure = ggplot(  cellsData, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2],
                                color = "Cluster")) +
  geom_point( stroke = 0) + # Just ignore point size and alpha here as it will be overriden after plotly conversion
  facet_wrap(facets = vars(HTO))

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  ggFigure = ggFigure +
    geom_text( data = dimReducLabels, # 'repel' not implemented by ggplotly yet
               aes( label = Cluster),
               color = "black",
               size = 3.5)
}

ggFigure = ggFigure +
  scale_color_manual( name = "Cluster", values = clustersColor) +
  theme_classic() +
  theme( #axis.text  = element_blank(),
    axis.title = element_blank(),
    #axis.ticks  = element_blank(),
    #axis.line = element_blank(),
    legend.position = if(if(exists( "showLegend")) showLegend else TRUE) NULL else "none",
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
         marker.opacity = 0.6) %>% 
  config( displaylogo = FALSE,
          toImageButtonOptions = list( format='svg'),
          modeBarButtons = list(list( 'toImage'),
                                list( 'zoom2d', 
                                      'zoomIn2d', 
                                      'zoomOut2d', 
                                      'pan2d', 
                                      'resetScale2d')))

# Set a unique html ID to access plot with html/javascript for external slider
figID = if(exists("figID")) figID+1 else 0;
plotDimReducGG[["elementId"]] = paste0("plotDimReducGG_withSlider", figID);

# Render plot. ggplot conversion generates warnings about incompatibility between 'scatter' and 'hoveron'
div(htmltools::tagList(plotDimReducGG));

# Create a html slider using javascript for controlling point size
cat( paste0( '\n',
             '<p>Point size: <span id="valueTag_', figID, '">', initialPointSize, '</span></p>',
             '<input type="range" min="2" max="50" value="', initialPointSize, '" id="slider_', figID, '" style="width: ', dimReducWidth, 'px;">',
             '<script> document.getElementById("slider_', figID, '").oninput = function() { Plotly.restyle("', plotDimReducGG[["elementId"]], '", {"marker.size": this.value}); document.getElementById("valueTag_', figID, '").innerHTML = this.value; } </script>',
             '\n'));
cat('\n \n')


## @knitr plotDimReducRaster_colorClusters_facetHTOs
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables

# Prepare dimreduc figure as ggplot
ggFigure = ggplot(  cellsData, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2],
                                color = "Cluster"))  +
  geom_point( alpha = PLOT_DIMREDUC_GROUPS_ALPHA, 
              size  = PLOT_DIMREDUC_GROUPS_POINTSIZE) +
  facet_wrap(facets = vars(HTO))

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  ggFigure = ggFigure +
    geom_text( data = dimReducLabels,
               aes( label = Cluster),
               color = "black",
               size = 2.5)
}

ggFigure = ggFigure +
  scale_color_manual( name = "Cluster", values = clustersColor) +
  theme_classic() +
  theme( #axis.text  = element_blank(),
    axis.title = element_blank(),
    #axis.ticks  = element_blank(),
    #axis.line = element_blank(),
    legend.position = if(if(exists( "showLegend")) showLegend else TRUE) NULL else "none",
    plot.margin = margin( 0, 0, 0, 0, "cm"),
    plot.title = element_text( face = "bold",
                               size = rel( 1.8), #rel( 16/14),
                               hjust = 0.5,
                               vjust = 1,
                               margin = margin( b = 7)));
print( ggFigure);
cat('\n \n')




## @knitr plotDimReducInteractiveGG_colorHTOs_facetClusters
# Using 'ggplotly' as facetting is not trivial with plotly 
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables

#### Facetted dimensionality reduction 2D (ggplot + plotly conversion)
# Parameters defined in chunk
# dimReducWidth  = 800;
# dimReducHeight = 800;
# initialPointSize = 5;

# Prepare dimreduc figure as ggplot
ggFigure = ggplot(  cellsData, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2],
                                color = "HTO")) +
  geom_point( stroke = 0) + # Just ignore point size and alpha here as it will be overriden after plotly conversion
  facet_wrap(facets = vars(Cluster), scales = "free") # Wrap and rescale around cluster coordinates

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  ggFigure = ggFigure +
    geom_text( data = dimReducLabels, # 'repel' not implemented by ggplotly yet
               aes( label = Cluster),
               color = "black",
               size = 3.5)
}

ggFigure = ggFigure +
  scale_color_manual( name = "HTO", values = HTOsColor) +
  theme_classic() +
  theme( #axis.text  = element_blank(),
    axis.title = element_blank(),
    #axis.ticks  = element_blank(),
    #axis.line = element_blank(),
    legend.position = if(if(exists( "showLegend")) showLegend else TRUE) NULL else "none",
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
         marker.opacity = 0.6) %>% 
  config( displaylogo = FALSE,
          toImageButtonOptions = list( format='svg'),
          modeBarButtons = list(list( 'toImage'),
                                list( 'zoom2d', 
                                      'zoomIn2d', 
                                      'zoomOut2d', 
                                      'pan2d', 
                                      'resetScale2d')))

# Set a unique html ID to access plot with html/javascript for external slider
figID = if(exists("figID")) figID+1 else 0;
plotDimReducGG[["elementId"]] = paste0("plotDimReducGG_withSlider", figID);

# Render plot
div(htmltools::tagList(plotDimReducGG));

# Create a html slider using javascript for controlling point size
cat( paste0( '\n',
             '<p>Point size: <span id="valueTag_', figID, '">', initialPointSize, '</span></p>',
             '<input type="range" min="2" max="50" value="', initialPointSize, '" id="slider_', figID, '" style="width: ', dimReducWidth, 'px;">',
             '<script> document.getElementById("slider_', figID, '").oninput = function() { Plotly.restyle("', plotDimReducGG[["elementId"]], '", {"marker.size": this.value}); document.getElementById("valueTag_', figID, '").innerHTML = this.value; } </script>',
             '\n'));
cat('\n \n')


## @knitr plotDimReducRaster_colorHTOs_facetClusters
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables

# Prepare dimreduc figure as ggplot
ggFigure = ggplot(  cellsData, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2],
                                color = "HTO"))  +
  geom_point( alpha = PLOT_DIMREDUC_GROUPS_ALPHA, 
              size  = PLOT_DIMREDUC_GROUPS_POINTSIZE) +
  facet_wrap(facets = vars(Cluster), scales = "free") # Wrap and rescale around cluster coordinates

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  ggFigure = ggFigure +
    geom_text( data = dimReducLabels,
               aes( label = Cluster),
               color = "black",
               size = 2.5)
}

ggFigure = ggFigure +
  scale_color_manual( name = "HTO", values = HTOsColor) +
  theme_classic() +
  theme( #axis.text  = element_blank(),
    axis.title = element_blank(),
    #axis.ticks  = element_blank(),
    #axis.line = element_blank(),
    legend.position = if(if(exists( "showLegend")) showLegend else TRUE) NULL else "none",
    plot.margin = margin( 0, 0, 0, 0, "cm"),
    plot.title = element_text( face = "bold",
                               size = rel( 1.8), #rel( 16/14),
                               hjust = 0.5,
                               vjust = 1,
                               margin = margin( b = 7)));
print( ggFigure);
cat('\n \n')




## @knitr plotDimReducRaster
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables

#### Dimensionality reduction 2D (ggplot)

# Plot a thumbnail highlighting cluster cells
ggFigure = ggplot(  cellsData, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2],
                                color = "Cluster")) +
             geom_point( alpha = PLOT_DIMREDUC_GROUPS_ALPHA,
                         size = if(PLOT_DIMREDUC_GROUPS_POINTSIZE==0) Seurat:::AutoPointSize(cellsData) else PLOT_DIMREDUC_GROUPS_POINTSIZE)

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  ggFigure = ggFigure +
               geom_text_repel( data = dimReducLabels,
                               aes( label = Cluster),
                               color = "black",
                               size = 7)
}

ggFigure = ggFigure +
             scale_color_manual( name = "Cluster", values = clustersColor) +
             theme_classic() +
             theme( #axis.text  = element_blank(),
                    axis.title = element_blank(),
                    #axis.ticks  = element_blank(),
                    #axis.line = element_blank(),
                    legend.position = if(if(exists( "showLegend")) showLegend else TRUE) NULL else "none",
                    plot.margin = margin( 0, 0, 0, 0, "cm"),
                    plot.title = element_text( face = "bold",
                                               size = rel( 1.8), #rel( 16/14),
                                               hjust = 0.5,
                                               vjust = 1,
                                               margin = margin( b = 7)));
print( ggFigure);




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




# MODULES ANALYSIS
##################

## @knitr heterogeneity_modules
# Just remind the warning(s) for genes names not in object, or modules that were transfered to individual monitoring of genes
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





# MONITORED GENES
#################

## @knitr heterogeneity_monitoredGenes
# Just remind the warning for genes names not in object
if(any( is.na( matchMonitoredGenes))) warning( paste( "Following gene(s) to be monitored could not be found in experimental data and will be ignored:", paste( monitoredGenesNotFound, collapse=" - ")));


## @knitr heterogeneity_monitoredGenes_expression_projection
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Plot expression values of monitored genes (TODO: message if list empty)
invisible( lapply( names( MONITORED_GENES), function(monitoredGroup)
{
  cat("\n#####", monitoredGroup, "\n");
  
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
  cat("\n#####", monitoredGroup, "\n");
  
  # Violinplot for expression value of monitored genes by cluster (+ number of 'zero' and 'not zero' cells)
  invisible( lapply( MONITORED_GENES[[monitoredGroup]], violinFeatureByCluster, seuratObject = sc10x, clustersColor = clustersColor));
  
  cat(" \n \n"); # Required for '.tabset'
}));



