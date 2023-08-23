# ##############################################################################
# Plots for provided genes lists
# ##############################################################################

#### HEATMAPS

## @knitr heatmap_function_definition

# Function plotting heatmaps of selected matrix in tabs with variations on rows organisation and scaling
# Creates external png in first run, and includes it in report in secind run to allow dynamic image size
# without altering final layout with 'png' and 'dev.off' calls during rendering... See 'FLAG_secondRun'.
# Used below to show DE genes for all comparisons
# DEBUG: avgExpressionMatrix=expressionValuesMatrix; basePath = NULL; scaleArgs = c("none", "row"); baseMain = ""; columnsClustering = FALSE; columnsCategories = NULL; columnsCategoriesColor = list(); splitColsByCatIndex = if(!is.null(columnsCategories)) 1 else NULL; rotateColsCat = FALSE; rowsClustering = FALSE; rowsCategories = NULL; rowsCategoriesColor = list(); slpitRowsByFC = FALSE; rotateRowsCat = FALSE; width = 900; height = 400+(15*nrow(avgExpressionMatrix))
plotHeatmaps = function(avgExpressionMatrix, 
                        basePath = NULL, 
                        scaleArgs = c("none", "row"), 
                        baseMain = "", 
                        columnsClustering = FALSE, 
                        columnsCategories = NULL, 
                        columnsCategoriesColor = list(), 
                        splitColsByCatIndex = if(!is.null(columnsCategories)) 1 else NULL, # Column index in 'columnsCategories' to use for splitting matrix columns
                        rotateColsCat = 0, # 0, 90, or NULL to hide columns categories names
                        #rowsClustering = FALSE, # Fixed, depending on the representation decided in function
                        rowsCategories = NULL,   # Mostly for logFC (col name "Log2 FC")
                        rowsCategoriesColor = list(), 
                        splitRowsByFC = FALSE, # Split Pos and Neg LogFC. Only works if a rows category is named "Log2 FC".
                        rotateRowsCat = 90, # 0, 90, or NULL to hide rows categories names
                        width = 900, 
                        height = 400+(15*nrow(avgExpressionMatrix)),
                        headerLevel = 6) # Header level to start from for rmarkdown titles (typically one more than current title level, which defines a tabset itself)
{
  
  # Prepare rows and columns annotation if corresponding dataframe provided
  rowsAnnot = NULL;
  if(!is.null(rowsCategories))
  {
    rowsAnnot = rowAnnotation( df = rowsCategories,
                               col = rowsCategoriesColor, # Random colors if NULL 
                               annotation_name_rot = 90) # 75
  }
  columnsAnnot = NULL;
  if(!is.null(columnsCategories))
  {
    columnsAnnot = HeatmapAnnotation( df = columnsCategories,
                                      col = columnsCategoriesColor) # Random colors if NULL 
  }
  
  message("## 01: Rows ordered by LogFC (as given)")
  for(scaleArg in (if(nrow(avgExpressionMatrix) > 1) scaleArgs else "none")) # No need to rescale if single row of data
  {
    cat(paste0("\n\n", strrep("#", headerLevel), " Rows ordering: none ==== Rows scaling: ",scaleArg, " {.tabset .tabset-fade .tabset-pills} \n\n"))
    
    main = paste0(baseMain, " - order:none - scaling:", scaleArg)
    
    # Render the plot with ComplexHeatmap
    hm = Heatmap( if(scaleArg == "rows") t(scale(t(avgExpressionMatrix))) else avgExpressionMatrix, 
                  name = "Average\nExpression", 
                  #col = RColorBrewer::brewer.pal(name = "Reds", n = 4),
                  col = rev(rocket(20)),
                  rect_gp = gpar(col = "white", lwd = 1), # Border of cells
                  column_names_rot = 90, #75
                  cluster_rows = FALSE,
                  cluster_row_slices = FALSE, # Do not cluster rows split categories (based on logFC) to keep order of slices under control
                  cluster_columns = columnsClustering,
                  left_annotation = rowsAnnot,
                  top_annotation = columnsAnnot,
                  # Make groups of columns based on categories in dataframe
                  column_split = if(!is.null(splitColsByCatIndex)) columnsCategories[[splitColsByCatIndex]],
                  column_title = if(is.null(rotateColsCat)) NULL else character(0), # Default: character(0), Disable: NULL
                  column_title_rot = if(is.null(rotateColsCat)) 0 else rotateColsCat, # NULL to disable name of columns categories (see line above)
                  # Make groups of rows based on sign of logFC
                  row_split = if(splitRowsByFC) factor(ifelse(rowsCategories[["Log2 FC"]] >0, "Pos", "Neg"), levels = c( "Pos", "Neg")) else NULL,
                  row_title = if(is.null(rotateRowsCat)) NULL else character(0), # Default: character(0), Disable: NULL
                  row_title_rot = if(is.null(rotateRowsCat)) 0 else rotateRowsCat) # NULL to disable name of rows categories (see line above)
    
    
    # Only render to file if a base path is specified
    if(!is.null(basePath))
    {
      filenamePNG = paste0( basePath, "_ORDER_none_SCALE_", scaleArg, ".png")
      
      if(FLAG_secondRun) # Check if this is a second run for rendering only (see analysisParams)
      {
        # Don't mess with RMD layout using 'png' and 'dev.off'...
        # Integrate the figure that was generated during first run
        cat( paste0( "\n![Heatmap ORDER:none SCALE:", scaleArg, "](", filenamePNG,")\n" ))
      } else # If it's the first run, create the figure as external file
      {
        # First run, generate figure as external file (messes up layout: figures in wring tabset)
        # Don't include them in report, will be done in secind run.
        png( filenamePNG, 
             width = width, 
             height = height)
        draw( hm)
        dev.off()
        cat("\n\nFirst (computing) run. Run report generation a second time to get all tables/figures and proper layout...")
      }
    }
    else
    {
      draw( hm)
    }
    
    cat("\n\n")
  }
  
  if(nrow(avgExpressionMatrix) > 1) # Clustering requires at least 2 values, and no need for different ordering anyway when only one row...
  {
    message("## 02: Rows ordered by hierarchical clustering")
    for(scaleArg in scaleArgs)
    {
      cat(paste0("\n\n", strrep("#", headerLevel), " Rows ordering: hClust ==== Rows scaling: ",scaleArg, " {.tabset .tabset-fade .tabset-pills} \n\n"))
      
      main = paste0(baseMain, " - order:hClust - scaling:", scaleArg)
      
      # Render the plot with ComplexHeatmap
      hm = Heatmap( if(scaleArg == "rows") t(scale(t(avgExpressionMatrix))) else avgExpressionMatrix, 
                    name = "Average\nExpression", 
                    #col = RColorBrewer::brewer.pal(name = "Reds", n = 4),
                    col = rev(rocket(20)),
                    rect_gp = gpar(col = "white", lwd = 1), # Border of cells
                    column_names_rot = 90, #75
                    cluster_rows = TRUE,
                    cluster_row_slices = FALSE, # Do not cluster rows split categories (based on logFC) to keep order of slices under control
                    clustering_distance_rows = "euclidean",
                    cluster_columns = columnsClustering,
                    left_annotation = rowsAnnot,
                    top_annotation = columnsAnnot,
                    # Make groups of columns based on categories in dataframe
                    column_split = if(!is.null(splitColsByCatIndex)) columnsCategories[[splitColsByCatIndex]],
                    column_title = if(is.null(rotateColsCat)) NULL else character(0), # Default: character(0), Disable: NULL
                    column_title_rot = if(is.null(rotateColsCat)) 0 else rotateColsCat, # NULL to disable name of columns categories (see line above)
                    # Make groups of rows based on sign of logFC
                    row_split = if(splitRowsByFC) factor(ifelse(rowsCategories[["Log2 FC"]] >0, "Pos", "Neg"), levels = c( "Pos", "Neg")) else NULL,
                    row_title = if(is.null(rotateRowsCat)) NULL else character(0), # Default: character(0), Disable: NULL
                    row_title_rot = if(is.null(rotateRowsCat)) 0 else rotateRowsCat) # NULL to disable name of rows categories (see line above)
      
      # Only render to file if a base path is specified
      if(!is.null(basePath))
      {
        filenamePNG = paste0( basePath, "_ORDER_hClust_SCALE_", scaleArg, ".png")
        
        if(FLAG_secondRun) # Check if this is a second run for rendering only (see analysisParams)
        {
          # Don't mess with RMD layout using 'png' and 'dev.off'...
          # Integrate the figure that was generated during first run
          cat( paste0( "\n![Heatmap ORDER:hClust SCALE:", scaleArg, "](", filenamePNG,")\n" ))
        } else # If it's the first run, create the figure as external file
        {
          # First run, generate figure as external file (messes up layout: figures in wring tabset)
          # Don't include them in report, will be done in secind run.
          png( filenamePNG, 
               width = width+50, # Add 50 pixels to accomodate for dendrogram 
               height = height)
          draw( hm)
          dev.off()
          cat("\n\nFirst (computing) run. Run report generation a second time to get all tables/figures and proper layout...")
        }
      }
      else
      {
        draw( hm)
      }
      
      cat("\n\n")
    }
    
    
    message("## 03: Rows ordered by average expression")
    for(scaleArg in "none") # No need for rows scaling when ordering by average expression... (?)
    {
      cat(paste0("\n\n", strrep("#", headerLevel), " Rows ordering: avgExp ==== Rows scaling: ",scaleArg, " {.tabset .tabset-fade .tabset-pills} \n\n"))
      
      # Here we recompute a rows annotation object with order from sorted matrix
      avgExpressionMatrix_orderRowsByAvgExp = avgExpressionMatrix[ order(rowMeans(avgExpressionMatrix)), , drop = FALSE]
      rowsAnnot_orderRowsByAvgExp = NULL;
      if(!is.null(rowsCategories))
      {
        rowsAnnot_orderRowsByAvgExp = rowAnnotation( df = rowsCategories[ rownames(avgExpressionMatrix_orderRowsByAvgExp), , drop = FALSE],
                                                     col = rowsCategoriesColor, # Random colors if NULL 
                                                     annotation_name_rot = 75)
      }
      
      main = paste0(baseMain, " - order:avgExp - scaling:", scaleArg)
      
      # Render the plot with ComplexHeatmap
      hm = Heatmap( if(scaleArg == "rows") t(scale(t(avgExpressionMatrix_orderRowsByAvgExp))) else avgExpressionMatrix_orderRowsByAvgExp, 
                    name = "Average\nExpression", 
                    #col = RColorBrewer::brewer.pal(name = "Reds", n = 4),
                    col = rev(rocket(20)),
                    rect_gp = gpar(col = "white", lwd = 1), # Border of cells
                    column_names_rot = 90, #75
                    cluster_rows = FALSE,
                    cluster_row_slices = FALSE, # Do not cluster rows split categories (based on logFC) to keep order of slices under control
                    cluster_columns = columnsClustering,
                    left_annotation = rowsAnnot_orderRowsByAvgExp,
                    top_annotation = columnsAnnot,
                    # Make groups of columns based on categories in dataframe
                    column_split = if(!is.null(splitColsByCatIndex)) columnsCategories[[splitColsByCatIndex]],
                    column_title = if(is.null(rotateColsCat)) NULL else character(0), # Default: character(0), Disable: NULL
                    column_title_rot = if(is.null(rotateColsCat)) 0 else rotateColsCat, # NULL to disable name of columns categories (see line above)
                    # Make groups of rows based on sign of logFC
                    row_split = if(splitRowsByFC) factor(ifelse(rowsCategories[rownames(avgExpressionMatrix_orderRowsByAvgExp), "Log2 FC"] >0, "Pos", "Neg"), levels = c( "Pos", "Neg")) else NULL,
                    row_title = if(is.null(rotateRowsCat)) NULL else character(0), # Default: character(0), Disable: NULL
                    row_title_rot = if(is.null(rotateRowsCat)) 0 else rotateRowsCat) # NULL to disable name of rows categories (see line above)
      
      # Only render to file if a base path is specified
      if(!is.null(basePath))
      {
        filenamePNG = paste0( basePath, "_ORDER_avgExp_SCALE_", scaleArg, ".png")
        
        if(FLAG_secondRun) # Check if this is a second run for rendering only (see analysisParams)
        {
          # Don't mess with RMD layout using 'png' and 'dev.off'...
          # Integrate the figure that was generated during first run
          cat( paste0( "\n![Heatmap ORDER:avgExp SCALE:", scaleArg, "](", filenamePNG,")\n" ))
        } else # If it's the first run, create the figure as external file
        {
          # First run, generate figure as external file (messes up layout: figures in wring tabset)
          # Don't include them in report, will be done in secind run.
          png( filenamePNG, 
               width = width, 
               height = height)
          draw( hm)
          dev.off()
          cat("\n\nFirst (computing) run. Run report generation a second time to get all tables/figures and proper layout...")
        }
      }
      else
      {
        draw( hm)
      }
      
      cat("\n\n")
    }
  }
}




## @knitr heatmaps_genes
# On first run create heatmap figures as external png files (variable dimensions 
# depending on number of rows and columns). Put them into report in second run.

# Loop on list of genes
for(currentListName in names(GENES_LIST))
{
  message( paste0( "Heatmaps for list: ", currentListName))
  
  # Create dir where png figures will be stored (and searched for in second run)
  pathHeatmapsOutput = file.path( PATH_ANALYSIS_OUTPUT,
                                  "heatmaps_expression",
                                  path_sanitize(currentListName, replacement = "_"))
  dir.create(pathHeatmapsOutput, recursive = TRUE, showWarnings = FALSE)
  
  cat( paste0( "\n\n#### Genes list: ", currentListName, " {.tabset .tabset-fade}\n\n"))
  
  ## 02 & 03
  # Get the precomputed average expression of genes for all clusters/HTOs combinations
  expressionValuesMatrix0203 = geneExpByClusterAndHTO[GENES_LIST[[currentListName]], , drop = FALSE]
  # Extract groups from colnames (created from factors HTO and Idents) for column annotation
  columnsCategories0203 = data.frame(do.call(rbind, strsplit(x = colnames(expressionValuesMatrix0203), split = "_")))
  # Convert to factor and reorder columns according to known factor levels
  colnames(columnsCategories0203) = c("Cluster", "HTO")
  columnsCategories0203[["Cluster"]] = factor(columnsCategories0203[["Cluster"]])
  columnsCategories0203[["HTO"]] = factor(columnsCategories0203[["HTO"]], levels = HTO_FACTOR_LEVELS)
  orderColumns = order(columnsCategories0203[["Cluster"]], columnsCategories0203[["HTO"]])
  expressionValuesMatrix0203 = expressionValuesMatrix0203[, orderColumns, drop = FALSE]
  columnsCategories0203 = columnsCategories0203[orderColumns, , drop = FALSE]
  
  # Filter eventual rows with all identical values (incompatible with hclust)
  sdZero = apply(expressionValuesMatrix0203, 1, sd)==0
  if(any(sdZero))
  {
    warning("Removed at least one line with all identical values (incompatible with hclust)")
    expressionValuesMatrix0203 = expressionValuesMatrix0203[ !sdZero , , drop = FALSE]
  }
  
  
  message("# Plots 02: All clusters & HTOs (by cluster)")
  cat("\n\n##### All clusters & HTOs (by cluster) {.tabset .tabset-fade .tabset-pills}\n\n")
  # Here make columns groups from category "Cluster" before clustering
  plotHeatmaps( avgExpressionMatrix = expressionValuesMatrix0203, 
                basePath = file.path(pathHeatmapsOutput, "AllClustersAndHTOs_CatByCluster"),
                scaleArgs = c("none", "rows"), # Generate figures with and without rows scaling 
                baseMain = currentListName, 
                columnsClustering = TRUE, 
                columnsCategories = columnsCategories0203, 
                columnsCategoriesColor = list("Cluster" = clustersColor, "HTO" = HTOsColor), 
                splitColsByCatIndex = 1, # Index of category for which to split the columns. NULL to ignore.
                rotateColsCat = NULL, # 0 or 90, NULL to hide columns categories names
                #rowsClustering = FALSE, # Fixed, depends on plots made in function
                rowsCategories = NULL,
                rowsCategoriesColor = list(), 
                splitRowsByFC = FALSE, # Split Pos and Neg LogFC. Only works if a rows category is named "Log2 FC". Cannot choose which category as for columns because a boolean operation is made on continuous log FC.
                rotateRowsCat = 90, # 0 or 90, NULL to hide rows categories names
                width = 290+(20*ncol(expressionValuesMatrix0203)), # logFC(60)+legend(100)+legendClusters(130)+cells (adds 50 internally when doing clustering)
                height = 300+(20*nrow(expressionValuesMatrix0203))) # colnames(200)+colcategories(50)+dendrogram(50)+cells
  
  message("# Plots 03: All clusters & HTOs (by HTO)")
  cat("\n\n##### All clusters & HTOs (by HTO) {.tabset .tabset-fade .tabset-pills}\n\n")
  # Here make columns groups from category "HTO" before clustering
  plotHeatmaps( avgExpressionMatrix = expressionValuesMatrix0203, 
                basePath = file.path(pathHeatmapsOutput, "AllClustersAndHTOs_CatByHTO"),
                scaleArgs = c("none", "rows"), # Generate figures with and without rows scaling 
                baseMain = currentListName, 
                columnsClustering = TRUE, 
                columnsCategories = columnsCategories0203, 
                columnsCategoriesColor = list("Cluster" = clustersColor, "HTO" = HTOsColor), 
                splitColsByCatIndex = 2, # Index of category for which to split the columns. NULL to ignore.
                rotateColsCat = 0, # 0 or 90, NULL to hide columns categories names
                #rowsClustering = FALSE, # Fixed, depends on plots made in function
                rowsCategories = NULL, 
                rowsCategoriesColor = list(), 
                splitRowsByFC = FALSE, # Split Pos and Neg LogFC. Only works if a rows category is named "Log2 FC". Cannot choose which category as for columns because a boolean operation is made on continuous log FC.
                rotateRowsCat = 90, # 0 or 90, NULL to hide rows categories names
                width = 340+(20*ncol(expressionValuesMatrix0203)), # logFC(60)+legend(100)+geneames(50)+legendClusters(130)+cells (adds 50 internally when doing clustering)
                height = 340+(20*nrow(expressionValuesMatrix0203))) # colnames(200)+colCategories(50)+dendrogram(50)+colCategoriesNames(40)+cells
  
  
  ## 06 Each cluster separately (for rows scaling) (but plot for all HTOs)
  for( currentIdentity in levels(columnsCategories0203[["Cluster"]]))
  {
    expressionValuesMatrix0607 = expressionValuesMatrix0203[, columnsCategories0203[["Cluster"]] == currentIdentity, drop = FALSE]
    columnsCategories0607 = columnsCategories0203[columnsCategories0203[["Cluster"]] == currentIdentity, , drop = FALSE]
    
    # Filter eventual rows with all identical values (incompatible with hclust)
    sdZero = apply(expressionValuesMatrix0607, 1, sd)==0
    if(any(sdZero))
    {
      warning("Removed at least one line with all identical values (incompatible with hclust)")
      expressionValuesMatrix0607 = expressionValuesMatrix0607[ !sdZero , , drop = FALSE]
    }
    
    message( paste0("Identity: ", currentIdentity))
    message( "# Plots 06: Selected clusters (all HTOs)")
    cat( paste0( "\n\n##### Selected cluster: ",  currentIdentity, " (all HTOs) {.tabset .tabset-fade .tabset-pills}\n\n"))
    # Here make columns groups from category "Cluster" before clustering
    plotHeatmaps( avgExpressionMatrix = expressionValuesMatrix0607, 
                    basePath = file.path(pathHeatmapsOutput, paste0("SelectedCluster_", gsub("\\W", "", path_sanitize(currentIdentity, replacement = "")), "_AllHTOs")),
                    scaleArgs = c("none", "rows"), # Generate figures with and without rows scaling 
                    baseMain = paste0(currentListName, " (", currentIdentity, ")"), 
                    columnsClustering = FALSE, 
                    columnsCategories = columnsCategories0607, 
                    columnsCategoriesColor = list("Cluster" = clustersColor, "HTO" = HTOsColor), 
                    splitColsByCatIndex = 1, # Index of category for which to split the columns. NULL to ignore.
                    rotateColsCat = NULL, # 0 or 90, NULL to hide columns categories names
                    #rowsClustering = FALSE, # Fixed, depends on plots made in function
                    rowsCategories = NULL, 
                    rowsCategoriesColor = list(), 
                    splitRowsByFC = FALSE, # Split Pos and Neg LogFC. Only works if a rows category is named "Log2 FC". Cannot choose which category as for columns because a boolean operation is made on continuous log FC.
                    rotateRowsCat = 90, # 0 or 90, NULL to hide rows categories names
                    width = 290+(20*ncol(expressionValuesMatrix0607)), # logFC(60)+legend(100)+legendClusters(130)+cells (adds 50 internally when doing clustering)
                    height = 250+(20*nrow(expressionValuesMatrix0607))) # colnames(200)+colcategories(50)+cells+!dendrogram(removed,50)
    
    # message("# Plots 07: Selected clusters (all HTOs) (+clustCols)")
    # cat( paste0( "\n\n##### Selected cluster: ",  currentIdentity, " (all HTOs) (+clustCols) {.tabset .tabset-fade .tabset-pills}\n\n"))
    # # Here make columns groups from category "Cluster" before clustering
    # plotHeatmaps( avgExpressionMatrix = expressionValuesMatrix0607,
    #               basePath = file.path(pathHeatmapsOutput, "SelectedClusters_AllHTOs_clustCols"),
    #               scaleArgs = c("none", "rows"), # Generate figures with and without rows scaling
    #               baseMain = paste0(currentListName, " (", currentIdentity, ")"),
    #               columnsClustering = TRUE,
    #               columnsCategories = columnsCategories0607,
    #               columnsCategoriesColor = list("Cluster" = clustersColor, "HTO" = HTOsColor),
    #               splitColsByCatIndex = 1, # Index of category for which to split the columns. NULL to ignore.
    #               rotateColsCat = NULL, # 0 or 90, NULL to hide columns categories names
    #               #rowsClustering = FALSE, # Fixed, depends on plots made in function
    #               rowsCategories = NULL,
    #               rowsCategoriesColor = list(),
    #               splitRowsByFC = FALSE, # Split Pos and Neg LogFC. Only works if a rows category is named "Log2 FC". Cannot choose which category as for columns because a boolean operation is made on continuous log FC.
    #               rotateRowsCat = 90, # 0 or 90, NULL to hide rows categories names
    #               width = 290+(20*ncol(expressionValuesMatrix0607)), # logFC(60)+legend(100)+legendClusters(130)+cells (adds 50 internally when doing clustering)
    #               height = 300+(20*nrow(expressionValuesMatrix0607))) # colnames(200)+colcategories(50)+cells+dendrogram(50)
  }
  
  
  ## 08 All cells and all HTOs
  message("# Plots 08: Global by HTO")
  cat("\n\n##### Global by HTO {.tabset .tabset-fade .tabset-pills}\n\n")
  # Compute the average expression of detected genes within the cells selections used for testing
  expressionValuesMatrix08 = geneExpByHTO[GENES_LIST[[currentListName]], , drop = FALSE]
  
  # Filter eventual rows with all identical values (incompatible with hclust)
  sdZero = apply(expressionValuesMatrix08, 1, sd)==0
  if(any(sdZero))
  {
    warning("Removed at least one line with all identical values (incompatible with hclust)")
    expressionValuesMatrix08 = expressionValuesMatrix08[ !sdZero , , drop = FALSE]
  }
  
  # Render the plots with graphic variations in tabs (scaling, ordering)
  plotHeatmaps( avgExpressionMatrix = expressionValuesMatrix08, 
                basePath = file.path(pathHeatmapsOutput, "GlobalByHTO"),
                scaleArgs = c("none", "rows"), # Generate figures with and without rows scaling 
                baseMain = currentListName, 
                columnsClustering = FALSE, 
                columnsCategories = NULL, 
                columnsCategoriesColor = list(), 
                splitColsByCatIndex = NULL, # Index of category for which to split the columns. NULL to ignore.
                rotateColsCat = NULL, # 0 or 90, NULL to hide columns categories names
                #rowsClustering = FALSE, # Fixed, depends on plots made in function
                rowsCategories = NULL, 
                rowsCategoriesColor = list(), 
                splitRowsByFC = FALSE, # Split Pos and Neg LogFC. Only works if a rows category is named "Log2 FC". Cannot choose which category as for columns because a boolean operation is made on continuous log FC.
                rotateRowsCat = 90, # 0 or 90, NULL to hide rows categories names
                width = max(c(250, 160+(20*ncol(expressionValuesMatrix08)))), # logFC(60)+legend(100)+legendClusters(130)+cells (adds 50 internally when doing clustering)
                height = max(c(230, 120+(20*nrow(expressionValuesMatrix08)))))  # colnames(120)+cells+!dendrogram(removed,50)
  
}



