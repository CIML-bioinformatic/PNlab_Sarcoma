# ##############################################################################
# Plots for provided genes lists
# ##############################################################################

#### CELLS COUNTS (individual genes)

## @knitr cellsCounts

# Get counts for genes of interest from Seurat object

# Loop on list of genes
for(currentListName in names(GENES_LIST))
{
  cat(paste0("\n\n#### ",  currentListName, "\n\n"))
      
  message( paste0( "Heatmaps for list: ", currentListName))
 
  currentCounts = as.matrix(GetAssayData(sc10x, assay = "RNA", slot = "counts"))[ GENES_LIST[[currentListName]], ]
  
  cat("\n\n")
  print( kable( apply(currentCounts, 1, function(x){table(x>0)}),
                caption = "Total number of positive (>0) cells"))
  
  
  print( kable( apply(currentCounts, 1, function(x)
    {
      tapply(x, sc10x[["factorHTO", drop = TRUE]], function(y)(sum(y>0)))
    }), caption = "Positive (count>0) cells by HTO"))
  
  
  # Convert counts to long format for processing with ggplot
  currentCountsLong = reshape2::melt(currentCounts, 
                                     varnames=c("geneName", "cellBarcode"),
                                     value.name = "Count")
  
  # Inject metadata for groupings with ggplot
  currentCountsLong[["Cluster"]] = Idents(sc10x)[ currentCountsLong[["cellBarcode"]]]
  
  
  countDistribution =  ggplot( currentCountsLong[ currentCountsLong[["Count"]]>0, ], 
                               aes( x = Count,
                                    #y= after_stat(density),
                                    fill = Cluster)) +
    geom_histogram(position = "dodge") +
    #facet_grid(  geneName ~ Cluster, scales = 'fixed')
    facet_grid( rows = vars(geneName), scales = 'fixed')
    
  print(countDistribution)
  
}





## @knitr cellsCounts_genesCombinations

for(currentListName in names(GENES_LIST))
{
  
  cat(paste0("\n\n#### ",  currentListName, "\n\n"))
  
  message( paste0( "Combinatory expression for list: ", currentListName))

  # Get count of current genes list from Seurat object  
  currentCounts = as.matrix(GetAssayData(sc10x, assay = "RNA", slot = "counts"))[ GENES_LIST[[currentListName]], ]
  
# combinedCounts = as.data.frame(t((currentCounts>0)+as.integer(0))) %>% 
#                    group_by_all() %>% 
#                    count() %>%        # Count number of rows for each group
#                    ungroup() %>%      # Reset grouping (count result is by group already)
#                    mutate( nbGenesPos = rowSums(across(-n))) %>% # Count number of positive genes in each line/group
#                    arrange(desc(nbGenesPos)) # Sort by number of positive genes
# 
# cellsCounts = as.data.frame(t((currentCounts>0)+as.integer(0))) %>% 
#                 mutate( nbGenesPos = rowSums(across())) %>% 
#                 arrange(desc(nbGenesPos))
# 
# # Make groups from each combination
# groupedCounts = as.data.frame(t((currentCounts>0)+as.integer(0))) %>% 
#                   group_by_all() # Make groups from each combination
#
# # Groups-related informations
# groupedCounts %>% tally(sort = TRUE)
# groupedCounts %>% group_keys()
# groupedCounts %>% group_indices()
# groupedCounts %>% group_rows() %>% head()
# groupedCounts %>% group_vars()
# groupedCounts %>% summarise()
# 
# # combinedCounts
# combinedCounts = groupedCounts %>% 
#                    tally(sort = TRUE) %>% # Count number of rows for each group
#                    ungroup() %>%          # Reset grouping (count result is by group already)
#                    mutate( nbGenesPos = rowSums(across(-n)) ) %>% # Count number of positive genes in each line/group
#                    arrange(desc(nbGenesPos)) # Sort by number of positive genes

  # Flag positive genes in cells
  posGenesDF = as.data.frame(t((currentCounts>0)+as.integer(0)))
  # Count number of positive genes in each cell
  countPosGenes = rowSums(posGenesDF)
  
  # Add computed count to Seurat object metadata (for umap plots)
  nameCurrentMD = paste0(currentListName, "_nbGenesPos")
  sc10x = AddMetaData(sc10x, countPosGenes, col.name = nameCurrentMD)
  
  # Select 3 largest counts of positive genes (to be highlighted in plots)
  largestCounts = tail( sort (unique( countPosGenes)), 3)
  
  # Attribute a color the the 3 largest counts
  largestCountsColor = c("yellow2", "darkorange", "red")
  names(largestCountsColor) = as.character(largestCounts)
  
  
  ### Histogram of counts
  # Histogram for number of cells relative to number of positive genes in list
  histogramCounts = ggplot( data.frame( countPosGenes, 
                                        Cluster = Idents(sc10x)[names(countPosGenes)],
                                        Condition = sc10x[["factorHTO"]][ names(countPosGenes),]), 
                            aes(x = countPosGenes)) + 
                      geom_bar( aes( fill = as.character(countPosGenes)),
                                show.legend = FALSE) +
                      geom_text( aes(label = stat(count)), 
                                 stat = "count", 
                                 colour = "black", 
                                 size = 3,
                                 vjust = -0.5) +
                      ylab( paste0( "Nb. of cells (total: ", ncol(sc10x), ")")) +
                      xlab( paste0( "Nb. of positive genes in list (list size:", nrow(currentCounts), ")")) +
                      scale_fill_manual( values = largestCountsColor, 
                                         na.value = "grey", 
                                         name = NULL) +
                      ggtitle( currentListName) +
                      scale_x_continuous( labels = min(countPosGenes):max(countPosGenes), 
                                          breaks = min(countPosGenes):max(countPosGenes))
  
  print(histogramCounts)
  
  #print(histogramCounts + facet_wrap(~Cluster, ncol = 1))
  
  print(histogramCounts + facet_wrap(~Cluster, nrow = 1))
  
  print(histogramCounts + facet_wrap(~Condition, nrow = 1))
  
  print(histogramCounts + facet_grid( rows = vars( Cluster), 
                                      cols = vars( Condition)))
  
  
  

  ### Dimensionality reduction plot showing number of positive genes
  print( DimPlot( sc10x, 
                  group.by = nameCurrentMD, 
                  order = TRUE, 
                  split.by = "factorHTO") +
           scale_color_manual( values = largestCountsColor, 
                               na.value = "grey", 
                               name = NULL)) 
  
  
  
  
  ### Heatmaps of genes by cells group within the 3 largest counts
  # Create dir where png figures will be stored (and searched for in second run)
  pathHeatmapsOutput = file.path( PATH_ANALYSIS_OUTPUT,
                                  "heatmaps_binary",
                                  path_sanitize(currentListName, replacement = "_"))
  dir.create(pathHeatmapsOutput, recursive = TRUE, showWarnings = FALSE)
  
  # Prepare filname and requested dimensions for external png figure
  filenamePNG_byCluster           = file.path( pathHeatmapsOutput, "topRepresented_byCluster.png")
  filenamePNG_byCount             = file.path( pathHeatmapsOutput, "topRepresented_byCount.png")
  filenamePNG_byCondition         = file.path( pathHeatmapsOutput, "topRepresented_byCondition.png")
  filenamePNG_byConditionAndCount = file.path( pathHeatmapsOutput, "topRepresented_byConditionAndCount.png")
  width  = 250 + (15*ncol(posGenesDF))
  height = 150+(4*sum(countPosGenes %in% largestCounts))
  

  # Prepare row/cells annotations to be added to heatmap
  rowsAnnot = rowAnnotation( df = data.frame( Condition = sc10x[["factorHTO"]][ countPosGenes %in% largestCounts,],
                                              Pos.Genes = countPosGenes[countPosGenes %in% largestCounts],
                                              Cluster = Idents(sc10x)[countPosGenes %in% largestCounts]),
                             col = list( Cluster = clustersColor,
                                         Pos.Genes = largestCountsColor,
                                         Condition = HTOsColor),
                             annotation_name_rot = 90) # 75
  
  

  
  
  if(FLAG_secondRun) # Check if this is a second run for rendering only (see analysisParams)
  {
    # Don't mess with RMD layout using 'png' and 'dev.off'...
    # Integrate the figure that was generated during first run
    cat( paste0( "\n\n",
                 "![", currentListName, " - Top represented cells by cluster](", filenamePNG_byCluster,"){width=20%} ",
                 "&nbsp;&nbsp;&nbsp;",
                 "![", currentListName, " - Top represented cells by count](", filenamePNG_byCount,"){width=20%} ",
                 "&nbsp;&nbsp;&nbsp;",
                 "![", currentListName, " - Top represented cells by condition](", filenamePNG_byCondition,"){width=20%} ",
                 "&nbsp;&nbsp;&nbsp;",
                 "![", currentListName, " - Top represented cells by condition and count](", filenamePNG_byConditionAndCount,"){width=20%} \n\n" ))

    
  } else # If it's the first run, create the figure as external file
  {
    # First run, generate figure as external file (messes up layout: figures in wring tabset)
    # Don't include them in report, will be done in secind run.
    
    # Split groups of cells by identity/cluster
    png( filenamePNG_byCluster, 
         width = width,
         height= height)
    
    hm_byCluster = Heatmap( posGenesDF[countPosGenes %in% largestCounts,],
                            col = c("#DDDDDD", "red"),
                            name = "Positive",
                            column_title = "Split by cluster",
                            rect_gp = gpar(col = "white", lwd = 1),
                            show_row_names = FALSE, 
                            clustering_distance_rows = "binary",
                            clustering_distance_columns = "binary",
                            right_annotation = rowsAnnot,
                            split = Idents(sc10x)[countPosGenes %in% largestCounts],
                            cluster_row_slices = FALSE)
    draw(hm_byCluster)
    dev.off()
    
    # Split groups of cells by number of positive genes in list
    png( filenamePNG_byCount, 
         width = width,
         height= height)
    
    hm_byCount = Heatmap( posGenesDF[countPosGenes %in% largestCounts,],
                          col = c("#DDDDDD", "red"),
                          name = "Positive",
                          column_title = "Split by positive count",
                          rect_gp = gpar(col = "white", lwd = 1),
                          show_row_names = FALSE, 
                          clustering_distance_rows = "binary",
                          clustering_distance_columns = "binary",
                          right_annotation = rowsAnnot,
                          split = countPosGenes[countPosGenes %in% largestCounts],
                          cluster_row_slices = FALSE)
    draw(hm_byCount)
    dev.off()
    
    # Split groups of cells by HTO/condition
    png( filenamePNG_byCondition, 
         width = width,
         height= height)
    
    hm_byCondition = Heatmap( posGenesDF[countPosGenes %in% largestCounts,],
                              col = c("#DDDDDD", "red"),
                              name = "Positive",
                              column_title = "Split by condition",
                              rect_gp = gpar(col = "white", lwd = 1),
                              show_row_names = FALSE, 
                              clustering_distance_rows = "binary",
                              clustering_distance_columns = "binary",
                              right_annotation = rowsAnnot,
                              split = sc10x[["factorHTO"]][ countPosGenes %in% largestCounts, ],
                              cluster_row_slices = FALSE)
    draw(hm_byCondition)
    dev.off()
    
    # Split groups of cells by HTO/condition
    png( filenamePNG_byConditionAndCount, 
         width = width + 150, # Account for non-rotated categories names (split)
         height= height)
    
    # Prepare factor levels for combination of condition and count (properly ordered)
    factorLevels_condition_count_DF = expand.grid( sort(unique(countPosGenes)), levels(sc10x[["factorHTO", drop = TRUE]])) # First factor varies fastest
    factorLevels_condition_count = apply( factorLevels_condition_count_DF, 1, function(x) {paste(x[2], x[1], collapse = " ")})
    # Create the split factor with previously computed levels order
    splitFactor = factor( paste( sc10x[["factorHTO"]][ countPosGenes %in% largestCounts, ], countPosGenes[countPosGenes %in% largestCounts]), levels = factorLevels_condition_count)
    
    hm_byConditionAndCount = Heatmap( posGenesDF[countPosGenes %in% largestCounts,],
                                      col = c("#DDDDDD", "red"),
                                      name = "Positive",
                                      column_title = "Split by condition and counts",
                                      rect_gp = gpar(col = "white", lwd = 1),
                                      show_row_names = FALSE, 
                                      clustering_distance_rows = "binary",
                                      clustering_distance_columns = "binary",
                                      right_annotation = rowsAnnot,
                                      split = splitFactor,
                                      row_title_rot = 0,
                                      cluster_row_slices = FALSE)
    draw(hm_byConditionAndCount)
    dev.off()
    
    
    # Dont include anything in the report yet since calls to 'png' and 'dev.off'
    # interferes with the layout, will be done on second run.
    cat("\n\nFirst (computing) run. Run report generation a second time to get all tables/figures and proper layout...")
  }
  
  
  cat(" \n \n"); # Required for '.tabset'
}










