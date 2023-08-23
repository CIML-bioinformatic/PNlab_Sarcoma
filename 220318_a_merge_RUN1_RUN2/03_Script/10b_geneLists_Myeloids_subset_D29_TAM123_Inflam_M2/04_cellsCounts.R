# ##############################################################################
# Plots for provided genes lists
# ##############################################################################

#### CELLS COUNTS (individual genes)

## @knitr cellsCounts


# Loop on list of genes
for(currentListName in names(GENES_LIST))
{
  cat(paste0("\n\n#### ",  currentListName, "\n\n"))
      
  # Get counts for genes of interest from Seurat object
  currentCounts = as.matrix(GetAssayData(sc10x, assay = "RNA", slot = "counts"))[ GENES_LIST[[currentListName]], ]
  
  cat("\n\n")
  print( kable( apply(currentCounts, 1, function(x){table(x>0)}),
                caption = "Total number of positive (>0) cells"))
  
  countsByHTO = apply(currentCounts, 1, function(x)
  {
    tapply(x, sc10x[["factorHTO", drop = TRUE]], function(y)(sum(y>0)))
  })
  
  print( kable( countsByHTO, caption = "Positive (count>0) cells by HTO"))
  
  # Chi2 analysis
  # Test independance between rows and columns (clusters and sample condition/HTO)
  # to detect counts over expected values under random condition.
  chisq = chisq.test(countsByHTO)
  
  # Show the expected counts considering margins statistics, to contrast with
  # observed counts shown above
  print( knitr::kable( addmargins(round(chisq[["expected"]],2)),
                       align   = "c",
                       caption = "Expected number of cells by population (Cluster) and condition (HTO) from X².\nTo compare with observed counts.")) #%>% kable_styling(bootstrap_options = c("condensed"), full_width = FALSE);
  
  lowExpectedCountsGenes = apply( chisq[["expected"]], 2, function(x) {any(x<5)})
  
  if(any(lowExpectedCountsGenes))
  {
    warning(paste("Some features show a too low expected number of positive cells for chi², filtering them:", paste( colnames(chisq[["expected"]])[lowExpectedCountsGenes], collapse = " - " )))
    
    chisq = chisq.test(countsByHTO[, !lowExpectedCountsGenes])
    
    print( knitr::kable( addmargins(round(chisq[["expected"]],2)),
                         align   = "c",
                         caption = "Expected number of cells by population (Cluster) and condition (HTO) from X².\nTo compare with observed counts.")) #%>% kable_styling(bootstrap_options = c("condensed"), full_width = FALSE);
    
    # Show result as simple table
    pander(chisq)
  }
  if(sum(!lowExpectedCountsGenes)>=3) 
  {
  
    # Set a layout to combine pearson residuals and contribution plots
    def.par <- par(no.readonly = TRUE) # save default, for resetting...
    layout(matrix(c(1,2), nrow=1, byrow = TRUE)) # Two next figures on same row
    
    # Plot Pearson residuals
    corrplot(chisq$residuals, 
             is.cor=FALSE, 
             addCoef.col="darkgrey", 
             tl.srt=60, 
             cl.ratio=0.35, 
             col = rev(colorRampPalette(RColorBrewer::brewer.pal(name="RdBu", 11))(100)),
             #col = colorRampPalette(c("red", "white", "blue"))(100),
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
      geom_bar(position = "dodge") +
      #facet_grid(  geneName ~ Cluster, scales = 'fixed')
      facet_grid( rows = vars(geneName), scales = 'free')
      
    #print(countDistribution)
  } else
  {
    cat("\n\n<BR>Number of features too low to apply chi² test...<BR>")
  }
}




## @knitr cellsCounts_violin

# Same approach but using norm/scaled data and violin/jitter plots

# Loop on list of genes
for(currentListName in names(GENES_LIST))
{
  cat(paste0("\n\n#### ",  currentListName, "\n\n"))
  
  currentExpression = as.matrix(GetAssayData(sc10x, assay = "RNA", slot = "data"))[ GENES_LIST[[currentListName]], ]
  
  # Convert counts to long format for processing with ggplot
  currentExpressionLong = reshape2::melt(currentExpression, 
                                     varnames=c("geneName", "cellBarcode"),
                                     value.name = "Expression")
  
  # Convert factor to character to prevent confusion when using as index
  currentExpressionLong[["cellBarcode"]] = as.character(currentExpressionLong[["cellBarcode"]])
  
  # Inject metadata for groupings with ggplot
  currentExpressionLong[["Cluster"]] = Idents(sc10x)[ currentExpressionLong[["cellBarcode"]]]
  
  # Prepare counting of zero and non-zero values by category to be added as text
  countsByAnnotation_pos = tapply( currentExpressionLong[["Expression"]], 
                                   list(geneName = currentExpressionLong[["geneName"]],
                                        Cluster = currentExpressionLong[["Cluster"]]),
                                   function(x)
                                   {
                                     return( (sum(x>0)))
                                     #return( paste0(sum(x==0), ' | ', sum(x>0), " →"))
                                   })
  
  countsByAnnotation_neg = tapply( currentExpressionLong[["Expression"]], 
                                   list(geneName = currentExpressionLong[["geneName"]],
                                        Cluster = currentExpressionLong[["Cluster"]]),
                                   function(x)
                                   {
                                     return( (sum(x==0)))
                                   })
  
  # Make a subselection for violinplot:
  #   Only consider expression values >0
  #   Only select genes with more than 2 points >0 in each group  (otherwise group gets dropped and shifts violin)
  selectedViolin = (currentExpressionLong[["Expression"]]>0) & (currentExpressionLong[["geneName"]] %in% rownames(countsByAnnotation_pos)[apply(countsByAnnotation_pos, 1, function(x){all(x>1)})])
  #selectedViolin = (currentExpressionLong[["Expression"]]>0)
  
  countsBoxplots = ggplot( currentExpressionLong, 
                           aes( x = geneName,
                                #y= after_stat(density),
                                group = Cluster,
                                color = Cluster)) +
    geom_jitter( aes( y = Expression), 
                 position = position_jitterdodge( jitter.width = 0.5, 
                                                  jitter.height = 0.05, 
                                                  dodge.width = 0.8, 
                                                  seed = SEED), 
                 alpha = 0.3,
                 show.legend = TRUE) +
    geom_violin( data = currentExpressionLong[selectedViolin , ],
                 aes( y = Expression), 
                 alpha = 0.3, 
                 position = position_dodge( width = 0.8), 
                 color = "#777777",draw_quantiles = 0.5) +
    # geom_label(data = reshape2::melt(countsByAnnotation_neg),
    #           aes(label = value, y = 0, fill = NULL), 
    #           alpha = 0.75,
    #           position = position_dodge( width = 0.8),
    #           hjust = 0.5) +
    geom_text( data = reshape2::melt(countsByAnnotation_neg),
               aes( label = value, 
                    y = - (max(currentExpressionLong[["Expression"]])/10)), 
               position = position_dodge( width = 0.8),
               hjust = 0.5,
               show.legend = FALSE) +
    geom_text( data = reshape2::melt(countsByAnnotation_pos),
               aes( label = value, 
                    y = max(currentExpressionLong[["Expression"]]) + (max(currentExpressionLong[["Expression"]])/10)), 
               position = position_dodge( width = 0.8),
               hjust = 0.5,
               show.legend = FALSE) +
    coord_flip() +
    facet_grid( rows = vars(geneName), scales = 'free', switch="both") +
    xlab("") +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

  print(countsBoxplots)
}




## @knitr cellsCounts_violin_eachCluster_byHTO

# Same approach but using norm/scaled data and violin/jitter plots

# Loop on list of genes
for(currentListName in names(GENES_LIST))
{
  cat(paste0("\n\n#### ",  currentListName, " {.tabset .tabset-fade} \n\n"))
  
  
  for(currentCluster in levels( Idents( sc10x)))
  {
    cat(paste0("\n\n##### ",  currentCluster, "\n\n"))
    
    currentExpression = as.matrix(GetAssayData(sc10x, assay = "RNA", slot = "data"))[ GENES_LIST[[currentListName]], Idents( sc10x) == currentCluster]
    
    # Convert counts to long format for processing with ggplot
    currentExpressionLong = reshape2::melt(currentExpression, 
                                           varnames=c("geneName", "cellBarcode"),
                                           value.name = "Expression")
    
    # Convert factor to character to prevent confusion when using as index
    currentExpressionLong[["cellBarcode"]] = as.character(currentExpressionLong[["cellBarcode"]])
    
    # Inject metadata for groupings with ggplot
    currentExpressionLong[["Cluster"]] = Idents(sc10x)[ currentExpressionLong[["cellBarcode"]]]
    currentExpressionLong[["HTO"]] = sc10x[["factorHTO", drop = TRUE]][ currentExpressionLong[["cellBarcode"]]]
    
    # Prepare counting of zero and non-zero values by category to be added as text
    countsByAnnotation_pos = tapply( currentExpressionLong[["Expression"]], 
                                     list(geneName = currentExpressionLong[["geneName"]],
                                          HTO = currentExpressionLong[["HTO"]]),
                                     function(x)
                                     {
                                       return( (sum(x>0)))
                                       #return( paste0(sum(x==0), ' | ', sum(x>0), " →"))
                                     })
    
    countsByAnnotation_neg = tapply( currentExpressionLong[["Expression"]], 
                                     list(geneName = currentExpressionLong[["geneName"]],
                                          HTO = currentExpressionLong[["HTO"]]),
                                     function(x)
                                     {
                                       return( (sum(x==0)))
                                     })
    
    # Make a subselection for violinplot:
    #   Only consider expression values >0
    #   Only select genes with more than 2 points >0 in each group  (otherwise group gets dropped and shifts violin)
    selectedViolin = (currentExpressionLong[["Expression"]]>0) & (currentExpressionLong[["geneName"]] %in% rownames(countsByAnnotation_pos)[apply(countsByAnnotation_pos, 1, function(x){all(x>1)})])
    #selectedViolin = (currentExpressionLong[["Expression"]]>0)
    
    countsBoxplots = ggplot( currentExpressionLong, 
                             aes( x = geneName,
                                  #y= after_stat(density),
                                  group = HTO,
                                  color = HTO)) +
      geom_jitter( aes( y = Expression), 
                   position = position_jitterdodge( jitter.width = 0.5, 
                                                    jitter.height = 0.05, 
                                                    dodge.width = 0.8, 
                                                    seed = SEED), 
                   alpha = 0.3,
                   show.legend = TRUE) +
      geom_violin( data = currentExpressionLong[selectedViolin , ],
                   aes( y = Expression), 
                   alpha = 0.3, 
                   position = position_dodge( width = 0.8), 
                   color = "#777777",draw_quantiles = 0.5) +
      # geom_label(data = reshape2::melt(countsByAnnotation_neg),
      #           aes(label = value, y = 0, fill = NULL), 
      #           alpha = 0.75,
      #           position = position_dodge( width = 0.8),
      #           hjust = 0.5) +
      geom_text( data = reshape2::melt(countsByAnnotation_neg),
                 aes( label = value, 
                      y = - (max(currentExpressionLong[["Expression"]])/10)), 
                 position = position_dodge( width = 0.8),
                 hjust = 0.5,
                 show.legend = FALSE) +
      geom_text( data = reshape2::melt(countsByAnnotation_pos),
                 aes( label = value, 
                      y = max(currentExpressionLong[["Expression"]]) + (max(currentExpressionLong[["Expression"]])/10)), 
                 position = position_dodge( width = 0.8),
                 hjust = 0.5,
                 show.legend = FALSE) +
      coord_flip() +
      facet_grid( rows = vars(geneName), scales = 'free', switch="both") +
      xlab("") +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
    
    print(countsBoxplots)
  }
}




## @knitr expression_bubbles_eachCluster_byHTO

# Same approach but using norm/scaled data and violin/jitter plots

# Loop on list of genes
for(currentListName in names(GENES_LIST))
{
  cat(paste0("\n\n#### ",  currentListName, " {.tabset .tabset-fade} \n\n"))
  
  
  for(currentCluster in levels( Idents( sc10x)))
  {
    cat(paste0("\n\n##### ",  currentCluster, "\n\n"))
    
    currentExpression = as.matrix(GetAssayData(sc10x, assay = "RNA", slot = "data"))[ GENES_LIST[[currentListName]], Idents( sc10x) == currentCluster]
    
    
    # Convert counts to long format for processing with ggplot
    currentExpressionLong = reshape2::melt(currentExpression, 
                                           varnames=c("geneName", "cellBarcode"),
                                           value.name = "Expression")
    
    # Convert factor to character to prevent confusion when using as index
    currentExpressionLong[["cellBarcode"]] = as.character(currentExpressionLong[["cellBarcode"]])
    
    # Inject metadata for groupings with ggplot
    currentExpressionLong[["Cluster"]] = Idents(sc10x)[ currentExpressionLong[["cellBarcode"]]]
    currentExpressionLong[["HTO"]] = sc10x[["factorHTO", drop = TRUE]][ currentExpressionLong[["cellBarcode"]]]
    
    # Prepare counting of zero, non-zero cells and average expression values by category
    pos =  reshape2::melt( tapply( currentExpressionLong[["Expression"]], 
                                   data.frame( geneName = factor(currentExpressionLong[["geneName"]]),
                                               HTO = factor(currentExpressionLong[["HTO"]]),
                                               Cluster = factor(currentExpressionLong[["Cluster"]])),
                                   function(x)
                                   {
                                     return( sum(x>0))
                                   }),
                           value.name = "NbCellsPositive")
    
    neg =  reshape2::melt( tapply( currentExpressionLong[["Expression"]], 
                                   data.frame( geneName = factor(currentExpressionLong[["geneName"]]),
                                               HTO = factor(currentExpressionLong[["HTO"]]),
                                               Cluster = factor(currentExpressionLong[["Cluster"]])),
                                   function(x)
                                   {
                                     return( sum(x==0))
                                   }),
                           value.name = "NbCellsNegative")
    
    avg =  reshape2::melt( tapply( currentExpressionLong[["Expression"]], 
                                   data.frame( geneName = factor(currentExpressionLong[["geneName"]]),
                                               HTO = factor(currentExpressionLong[["HTO"]]),
                                               Cluster = factor(currentExpressionLong[["Cluster"]])),
                                   mean),
                           value.name = "AvgExpression")
    
    # Merge computed values as single dataframe for plotting
    expressionSummaryDF = merge(merge( pos, neg), avg)
    
    
    
    # Plot
    print( ggplot( expressionSummaryDF,
                   aes( x = HTO, 
                        y = geneName)) +
             geom_point( aes( size = (NbCellsPositive/(NbCellsPositive+NbCellsNegative)) * 100,
                              fill = AvgExpression,
                              color = AvgExpression), 
                         show.legend = TRUE,
                         shape = 21,
                         stroke = 0.2) +
             scale_fill_gradient( #colours = c("white", "blue", "red"),
               low = "blue", 
               high = "red",
               name = "Average expression") +
             scale_color_gradient( #colours = c("white", "blue", "red"),
               low = "blue", 
               high = "red",
               name = "Average expression") +
             #scale_size_continuous(name = "% Expression in cluster") +
             scale_size_continuous(name = "% Positive Cells", range = c(1,10)) +
             scale_x_discrete(position = "top") +
             scale_y_discrete(position = "left") +
             # scale_x_discrete(guide = guide_axis(n.dodge = 1)) + 
             theme(strip.text.y.right = element_text(angle = 90),
                   axis.text.x = element_text(angle = 45, hjust=0)) +
             labs( x = "Cluster", 
                   y = "") +
             facet_wrap( vars(Cluster), ncol =1, strip.position = "right") )
    
    
  }
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
  largestCountsColor = c("darkgrey", "darkgrey", "darkgrey") #c("yellow2", "darkorange", "red")
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








