# ##################################################
# The aim of this script is to find the most
# variable genes among the cell population
# ##################################################


## FIND THE MOST VARIABLE GENES USING THE SEURAT PACKAGE METHOD

## @knitr findVariableGenes_seuratMethod

# Find most variable features
sc10x = FindVariableFeatures( object = sc10x, 
                              selection.method = "vst",
                              loess.span = 0.3,
                              clip.max = "auto",
                              mean.function = ExpMean, 
                              dispersion.function = LogVMR,
                              nfeatures = VARIABLE_FEATURES_MAXNB,
                              verbose = .VERBOSE);

variablesGenesStats = paste(length( VariableFeatures( sc10x)), "/", nrow( sc10x));


## @knitr findVariableGenes_summaryPlot

# Prepare a Variance/Expression plot highlighting variable genes and add names of most variable genes
suppressMessages( suppressWarnings( LabelPoints( plot = VariableFeaturePlot( sc10x) + theme(legend.position = "none"), 
                                                 points = head( VariableFeatures( sc10x), 10), 
                                                 repel = TRUE)));




## @knitr findVariableGenes_summaryTable

# Extract variable genes info as data.frame
variableAnnotationsDT = head( HVFInfo( object = sc10x, assay = "RNA", selection.method = 'vst')[VariableFeatures( sc10x),], VARIABLE_FEATURES_SHOWTOP);
variableAnnotationsDT = cbind("Gene" = rownames(variableAnnotationsDT), variableAnnotationsDT);

# Create a table in report containing information about top variable genes
datatable( variableAnnotationsDT, # Set annotation names as column instead of rownames so datatable handles column search properly
           class = "compact",
           filter="top",
           rownames = FALSE,
           colnames = c("Gene", "Avg. Expression", "Variance", "Std. Variance"),
           extensions = c('Buttons', 'Select'),
           options = list(dom = "<'row'<'col-sm-8'B><'col-sm-4'f>> <'row'<'col-sm-12'l>> <'row'<'col-sm-12'rt>> <'row'<'col-sm-12'ip>>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                          autoWidth = FALSE,
                          buttons = exportButtonsListDT,
                          columnDefs = list( 
                            list( # Center all columns except first one
                              targets = 1:(ncol( variableAnnotationsDT)-1),
                              className = 'dt-center'),
                            list( # Set renderer function for 'float' type columns
                              targets = 1:(ncol( variableAnnotationsDT)-1),
                              render = htmlwidgets::JS( "function ( data, type, row ) {return type === 'export' ? data : data.toFixed(4);}"))), 
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
  formatStyle( columns = "mean",
               background = styleColorBar( data = range( variableAnnotationsDT[["mean"]]), 'lightblue', angle = -90),
               backgroundSize = '95% 50%',      # Set horizontal and vertical span in cell
               backgroundRepeat = 'no-repeat',
               backgroundPosition = 'center') %>%
  formatStyle( columns = "variance",
               background = styleColorBar( data = range( variableAnnotationsDT[["variance"]]), 'lightblue', angle = -90),
               backgroundSize = '95% 50%',      # Set horizontal and vertical span in cell
               backgroundRepeat = 'no-repeat',
               backgroundPosition = 'center') %>%
  formatStyle( columns = "variance.standardized",
               background = styleColorBar( data = range( variableAnnotationsDT[["variance.standardized"]]), 'lightblue', angle = -90),
               backgroundSize = '95% 50%',      # Set horizontal and vertical span in cell
               backgroundRepeat = 'no-repeat',
               backgroundPosition = 'center');



