# ##############################################################################
# This script represents cells on a 2D plot (e.g. dimensionality reduction), and
# uses javascript link between report objects to select cells in plots and show
# their characteristics on attached graphs or tables. 
# 
# ##############################################################################




## @knitr plotsPreparation
#### Prepare SharedData object and variables for plots 

# For each cell, gather all data to be used for plots generation
cellsData = cbind( "Cell" = colnames( sc10x), # cell names from rownames conflict with search field in datatable
                    sc10x[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                    "percent.mito" = if("percent.mito" %in% colnames(sc10x[[]])) as.numeric(format(sc10x[["percent.mito", drop = TRUE]], digits = 5)) else NULL,
                    "percent.ribo" = if("percent.ribo" %in% colnames(sc10x[[]])) as.numeric(format(sc10x[["percent.ribo", drop = TRUE]], digits = 5)) else NULL,
                    "Cluster" = Idents( sc10x),
                    cellsCoordinates,
                    jitterCoords = jitter(rep(0, ncol( sc10x)), amount = 0.15) + 0.17);

# Create text to show under cursor for each cell for plotly hover
hoverText = do.call(paste, c(Map( paste,
                                  c( "",
                                      "Cell ID: ",
                                      "# UMIs: ",
                                      "# Genes: ",
                                      if("percent.mito" %in% colnames(sc10x[[]])) "% Mito: ",
                                      if("percent.ribo" %in% colnames(sc10x[[]])) "% Ribo: "),
                                  cellsData[1:(ncol(cellsData)-4)], # Do not include cluster and coordinates in hover text
                                  sep = ""),
                    sep = "\n"));

# Create SharedData object (using plotly accessor) for use by all plots/widgets
dataShare = highlight_key(cellsData);




## @knitr plotClustersSize
#### Show number of cells in each cluster

clustersCount = as.data.frame( table( Cluster = Idents(sc10x)), responseName = "CellCount");

print( knitr::kable( clustersCount,
                      align = "c")) #%>% kable_styling(bootstrap_options = c("condensed"), full_width = FALSE);




## @knitr plotsPanelsRender

#### Define size of plot panels to be assembled using flex CSS

# Dimensionality reduction 2D (plotly)
dimReducWidth  = 800;
dimReducHeight = 800;

# Violinplots assembled figure (using subplot)
violinWidth = 100 * nlevels(cellsData[["Cluster"]]);
violinHeight = 800;

# Datatable width
datatableWidth = dimReducWidth + violinWidth;




#### Plot cells according to selected coordinates (dimensionality reduction) 

plotDimReduc = plot_ly( dataShare, # 'SharedData' object for plots interactions
                        x = as.formula( paste( "~", colnames( cellsCoordinates)[1])),
                        y = as.formula( paste( "~", colnames( cellsCoordinates)[2])),
#                       name = ~paste("Cluster", Cluster),
                        color = ~Cluster,
                        colors = clustersColor,
                        width = dimReducWidth,
                        height = dimReducHeight) %>%
  add_trace( type = "scattergl",
              mode = "markers",
              marker = list( size = 5,
                            opacity = 0.6),
              text = hoverText,
              hoverinfo = "text+name") %>%
  # Compute median of coordinates for each group and use it as new data for Cluster text annotations (type='scatter', mode='text')
  add_trace( data = as.data.frame( do.call( rbind, by( cellsData, cellsData[["Cluster"]], function(x){ return( data.frame( lapply( x[colnames( cellsCoordinates)], median), "Cluster" = x[1, "Cluster"]));}))),
              x = as.formula( paste( "~", colnames( cellsCoordinates)[1])),
              y = as.formula( paste( "~", colnames( cellsCoordinates)[2])),
              type = "scattergl",
              mode = "text",
              text = ~Cluster,
              textfont = list( color = '#000000', size = 20),
              hoverinfo = 'skip',
              showlegend = FALSE) %>%
#  hide_colorbar() %>%
#  hide_legend() %>%
  config( displaylogo = FALSE,
          toImageButtonOptions = list( format='svg'),
          modeBarButtons = list(list( 'toImage'),
                                list( 'select2d',
                                      'lasso2d'),
                                list( 'zoom2d', 
                                      'zoomIn2d', 
                                      'zoomOut2d', 
                                      'pan2d', 
                                      'resetScale2d'))) %>% 
  highlight( on = 'plotly_selected');




#### Show cells stats (UMIs, Genes, Mito, Ribo) distribution split by cluster

# Generate plotly violin/jitter panels for #umis, #genes, and %mitochondrial stats
lypanel_umis  = plotViolinJitter( dataShare,
                                  xAxisFormula = ~as.numeric( Cluster),
                                  yAxisFormula = ~nCount_RNA,
                                  colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                  yAxisTitle = "# UMIs",
                                  hoverText = hoverText,
                                  traceName = ~paste( "Cluster", Cluster),
                                  xTicklabels = levels(cellsData[["Cluster"]]),
                                  panelWidth = violinWidth,
                                  panelHeight = violinHeight);

lypanel_genes = plotViolinJitter( dataShare,
                                  xAxisFormula = ~as.numeric( Cluster),
                                  yAxisFormula = ~nFeature_RNA,
                                  colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                  yAxisTitle = "# Genes",
                                  hoverText = hoverText,
                                  traceName = ~paste( "Cluster", Cluster),
                                  xTicklabels = levels(cellsData[["Cluster"]]),
                                  panelWidth = violinWidth,
                                  panelHeight = violinHeight);

lypanel_mitos = if("percent.mito" %in% colnames(sc10x[[]])) plotViolinJitter( dataShare,
                                                                              xAxisFormula = ~as.numeric( Cluster),
                                                                              yAxisFormula = ~percent.mito,
                                                                              colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                                                              yAxisTitle = "% Mito",
                                                                              hoverText = hoverText,
                                                                              traceName = ~paste( "Cluster", Cluster),
                                                                              xTicklabels = levels(cellsData[["Cluster"]]),
                                                                              panelWidth = violinWidth,
                                                                              panelHeight = violinHeight) else NULL;

lypanel_ribos = if("percent.ribo" %in% colnames(sc10x[[]])) plotViolinJitter( dataShare,
                                                                              xAxisFormula = ~as.numeric( Cluster),
                                                                              yAxisFormula = ~percent.ribo,
                                                                              colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                                                              yAxisTitle = "% Ribo",
                                                                              hoverText = hoverText,
                                                                              traceName = ~paste( "Cluster", Cluster),
                                                                              xTicklabels = levels(cellsData[["Cluster"]]),
                                                                              panelWidth = violinWidth,
                                                                              panelHeight = violinHeight) else NULL;

# Set panels as a list, define plotly config, and remove legend
panelsList = list( lypanel_umis, lypanel_genes, lypanel_mitos, lypanel_ribos);
panelsList = lapply( panelsList, config, displaylogo = FALSE,
                      toImageButtonOptions = list( format='svg'),
                      modeBarButtons = list( list('toImage'),
                                             list( 'select2d',
                                                   'lasso2d'),
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
                      autosize = FALSE) %>% 
             highlight( on = 'plotly_selected');




#### Show table with cells informations

# Create datatable
tableCells = datatable( dataShare,
                        width = datatableWidth,
                        class = "compact",
                        filter="top",
                        rownames = FALSE,
                        colnames = c("Cell", "Num. ID", "# UMIs", "# Features", "% Mito", "% Ribo", "Identity", "Coord_1", "Coord_2", "GraphicJitter"),
                        extensions = c('Buttons', 'Select'),
                        options = list(dom = "<'row'<'col-sm-8'B><'col-sm-4'f>> <'row'<'col-sm-12'l>> <'row'<'col-sm-12'rt>> <'row'<'col-sm-12'ip>>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                                        autoWidth = TRUE,
                                        buttons = exportButtonsListDT,
                                        columnDefs = list(
                                          list( # Center all columns except first one
                                            targets = 1:(ncol( cellsData)-1),
                                            className = 'dt-center'),
                                          list( # Set renderer function for 'float' type columns (LogFC)
                                            targets = 4:5, # %Mito and %Ribo
                                            render = htmlwidgets::JS("function ( data, type, row ) {return type === 'export' ? data : data.toFixed(4);}")),
                                          list( # Hide unnecessary columns (coordinates)
                                            targets = 7:9,
                                            visible = FALSE,
                                            searchable = FALSE)),
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
  # Add bars relative to numeric values
  formatStyle( columns = "nCount_RNA",
               background = styleColorBar( data = range( cellsData[["nCount_RNA"]]), 'lightblue', angle = -90),
               backgroundSize = '95% 50%',      # Set horizontal and vertical span in cell
               backgroundRepeat = 'no-repeat',
               backgroundPosition = 'center') %>%
  formatStyle( columns = "nFeature_RNA",
               background = styleColorBar( data = range( cellsData[["nFeature_RNA"]]), 'lightblue', angle = -90),
               backgroundSize = '95% 50%',      # Set horizontal and vertical span in cell
               backgroundRepeat = 'no-repeat',
               backgroundPosition = 'center') %>%
  formatStyle( columns = "percent.mito",
               background = styleColorBar( data = c(0, 100), 'lightblue', angle = -90),
               backgroundSize = '95% 50%',      # Set horizontal and vertical span in cell
               backgroundRepeat = 'no-repeat',
               backgroundPosition = 'center') %>%
  formatStyle( columns = "percent.ribo",
               background = styleColorBar( data = c(0, 100), 'lightblue', angle = -90),
               backgroundSize = '95% 50%',      # Set horizontal and vertical span in cell
               backgroundRepeat = 'no-repeat',
               backgroundPosition = 'center') %>%
  # Add color from cluster
  formatStyle( columns = "Cluster",
               backgroundColor = styleEqual( names(clustersColor),
                                             scales::alpha(clustersColor, 0.3)));



#### Assemble plot panels using flex and render

# 'browsable' required in console, not in script/document
#browsable(
  div( style = paste("display: flex; flex-wrap: wrap; width: 100%; height: 100%;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),  
        div( style = paste("flex: 0 0 auto; display: flex; flex-wrap: nowrap; flex-direction: row;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: yellow;"),  
          div( tableCells, style = paste("flex : 0 0 auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;"))),
        div( style = paste("flex: 0 0 auto; display: flex; flex-wrap: nowrap; flex-direction: row;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: green;"),  
          div( plotDimReduc, style = paste("flex : 0 0 auto;",  if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;")),
          div( plotPanels, style = paste("flex : 0 0 auto",  if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;"))))

#)

