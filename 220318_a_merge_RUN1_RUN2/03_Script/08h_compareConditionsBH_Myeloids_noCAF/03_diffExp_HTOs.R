# ##############################################################################
# This script performs the differential expression analysis between groups of
# cells defined by HTOs (see analysisParams.R).
# ##############################################################################




## @knitr diffExp_computeDE
# Perform the DE analyses with specified parameters/conditions and store results
# in a list 

# Define the function selecting the groups of cells on specified conditions/HTOs 
# and selected cluster (union of cells from specified values, from Idents(), "AllCells" = all)
selectCellsForComparisons = function(currentComparison, idents = "AllCells")
{
  # DEBUG: currentComparison = HTO_DIFFEXP_COMPARISONLIST[["Panth_D21_vs_D29"]]
  # DEBUG: idents = "AllCells"
  # DEBUG: boxplot(lapply(compareCells, function(cellsName){as.numeric(GetAssay(sc10x)["Gzmc", cellsName])}))
  
  # currentComparison should be a list with 2 elements, each being a vector of
  # condition(s) (HTO names)
  if(length(currentComparison)!=2) 
  {
    warning( "2 conditions should be provided for DE comparison. Skipping current comparison...")
    return( NULL);
  }
  if( !((length( idents) >= 1) && all(idents %in% c("AllCells", levels( Idents( sc10x))))) )
  {
    warning( "Values for 'idents' are not in 'AllCells' or levels of 'Idents(sc10x)': ", paste( paste0( "'", levels( Idents( sc10x)), "'"), collapse = " - "), ". Skipping current comparison..." )
    return( NULL);
  }
  
  # Build list of cells to be compared from conditions/HTOs names and eventually selected cluster
  compareCells = lapply(currentComparison, function(currentCondition)
  {
    selectedCellIdents = (Idents(sc10x) %in% idents) | (idents=="AllCells") # Use unary operator to 'or' all values against TRUE when idents=="AllCells"
    selectedCellHTOs   = sc10x[["factorHTO", drop = TRUE]] %in% currentCondition # 'factorHTO' built when loading object
    colnames(sc10x)[selectedCellIdents & selectedCellHTOs] # Return cells passing both selections
  })
  
  # Remove eventually overlapping cells (with a warning, should not happen with HTOs)
  selectedCellsName = unlist( compareCells)
  duplicatedCellsName = unique( selectedCellsName[duplicated( selectedCellsName)])
  if(length( duplicatedCellsName))
  {
    warning( paste( length( duplicatedCellsName), "cells were selected in multiple conditions to be compared, removing them"))
    compareCells = lapply( compareCells, function(cellsName) { cellsName[!cellsName %in% duplicatedCellsName]})
  }
  
  # Abort (return NULL with a warning) on eventual empty group
  emptyGroups = (sapply( compareCells, length) < 3) # Only works because previous filterings (selection and duplicated cells filter) don't remove the list element but gives it character(0), and because number of groups is tested on function start. If adding a filter which removes list element, this test must be adjusted (eventually without reporting empty group(s) name). 
  if( any( emptyGroups))
  {
    warning( paste0( "Selected group for DE analysis contains less than 3 cells: ", paste( paste0("'", names( compareCells)[emptyGroups], "'"), collapse = " - "),". Skipping current comparison..."))
    return( NULL)
  }
  
  return(compareCells)
}
  

# Function performing the DE analysis using previously selected set of cells
diffExp = function(compareCells)
{
  message( paste0( "nbCells: ", length(compareCells[[1]]), " - ", length(compareCells[[2]])))
 
  if(is.null(compareCells)) return(NULL)
  
  # Identify DE genes
  resultDE = FindMarkers( object          = sc10x,
                          ident.1         = compareCells[[1]],
                          ident.2         = compareCells[[2]],
                          test.use        = HTO_DIFFEXP_METHOD,
                          only.pos        = HTO_DIFFEXP_ONLYPOS,
                          min.pct         = HTO_DIFFEXP_MINPCT,
                          logfc.threshold = HTO_DIFFEXP_LOGFC_THR,
                          random.seed     = SEED,
                          verbose         = .VERBOSE)
  
  
  if(HTO_DIFFEXP_FORCE_ADJUST_PVAL_BH)
  {
    # Force recomputing adjusted PValues using BH instead of Bonferroni (default)
    resultDE[["p_val_adj"]] = p.adjust(resultDE[["p_val"]], method = "BH")
    resultDE = resultDE[ order( resultDE[["p_val_adj"]]), , drop = FALSE] # Update ordering
  }
  
  return( resultDE)
}


# Define cluster groups to test for each comparison (union of cells if several)
diffExpHTO_identsToTest        = c("AllCells", levels(Idents(sc10x)))
names(diffExpHTO_identsToTest) = diffExpHTO_identsToTest

# Perform cells selection on selected conditions (HTO) and clusters (idents).
# The result is a list of identities/clusters nested in a list of comparisons.
selectedCells_nestedList = lapply(HTO_DIFFEXP_COMPARISONLIST, function(currentComparison)
{
  lapply(diffExpHTO_identsToTest, function(currentIdent)
  {
    selectCellsForComparisons(currentComparison, idents=currentIdent) # Call the previously defined function for cell selection on each defined comparison (see analysisParams.R)
  })
})

# Now perform the actual DE analysis for each selection of cells.
# The result is a list of identities/clusters nested in a list of comparisons.
resultsDE_nestedList_notFiltered = mclapply( selectedCells_nestedList, 
                                             lapply, 
                                             diffExp,
                                             mc.cores = NBCORES) # Unfiltered, used for volcanoplot

# Remove NULL or empty 
# removeNullOrEmtpyFromList = function(currentList)
# {
#   return(currentList[sapply(currentList, length)!=0]) # Does not remove empty data.frame (0 rows)
# }

# Filter based on PValue threshold
resultsDE_nestedList = lapply( resultsDE_nestedList_notFiltered, 
                               lapply, 
                               function(x) 
                                 {
                                   if(is.null(x)) return(NULL)
                                   return(x[ x[["p_val_adj"]] <= HTO_DIFFEXP_PVAL_THR, , drop = FALSE]) 
                                 });




# Save DEG lists as 'tsv' tables in (created specific subfolders (only if table not empty)
for(currentComparisonName in names(resultsDE_nestedList))
{
  for(currentIdent in names(resultsDE_nestedList[[currentComparisonName]])) 
  {
    # Create the folder (even if no result)
    pathDiffExpOutput = file.path(PATH_ANALYSIS_OUTPUT, "DiffExpByHTO", currentComparisonName)
    dir.create(pathDiffExpOutput, showWarnings = FALSE, recursive = TRUE);
    
    # Write file only if result is not empty
    if(!is.null( resultsDE_nestedList[[currentComparisonName]][[currentIdent]]))
    {
      if(nrow( resultsDE_nestedList[[currentComparisonName]][[currentIdent]]))
      {
        write.table( resultsDE_nestedList[[currentComparisonName]][[currentIdent]],
                     file= file.path( pathDiffExpOutput, paste0( outputFilesPrefix, "DiffExp_" , currentComparisonName, "_", gsub(" ", "_", currentIdent), ".tsv")),
                     quote = FALSE,
                     row.names = TRUE, 
                     col.names = NA, # Add a blank column name for row names (CSV convention)
                     sep="\t");
      }
    }
    #else
    #{
    #  resultsDE_nestedList[[currentIdent]][[currentComparisonName]]=NULL;
    #}
  }
}

# Filter top differentially expressed genes and remove eventual empty elements
topDiffExp_nestedList = lapply(resultsDE_nestedList, function(resultsCurrentComparison)
{
  # Return top ones (sorted by PValue as 'FindMarkers' output), ignore eventual 
  # NULLs from invalid conditions specifications. Add a column 'gene'.
  for(currentIdentity in names(resultsCurrentComparison))
  {
    # Remove an eventual empty data.frame (or NULL value)
    if(is.null(resultsCurrentComparison[[currentIdentity]]) || (nrow(resultsCurrentComparison[[currentIdentity]])) == 0)
    {
      resultsCurrentComparison[[currentIdentity]] = NULL
    } else
    {
      # Add a "gene" column from rownames
      resultsCurrentComparison[[currentIdentity]] = data.frame(gene = rownames(resultsCurrentComparison[[currentIdentity]]), 
                                                               resultsCurrentComparison[[currentIdentity]])
      # Eventually select the top rows
      if(!is.null( HTO_DIFFEXP_SELECT_TOP)) 
      {
        resultsCurrentComparison[[currentIdentity]] = head( resultsCurrentComparison[[currentIdentity]], 
                                                            n = HTO_DIFFEXP_SELECT_TOP)      
      }
    }
  }
  return(resultsCurrentComparison)
});


# Filter eventual empty comparisons (with a warning). Nested empty identities
# filtered in previous block.
emptyComparisons = (sapply(topDiffExp_nestedList, length) == 0)
if(any(emptyComparisons)) 
{
  warning( paste0( "Some comparisons selected for differential expression analysis returned no result and will be ignored: ",
                   paste( paste0( "'", names(topDiffExp_nestedList)[emptyComparisons], "'"), collapse = " - "), "."))
  topDiffExp_nestedList[emptyComparisons] = NULL
}




# @knitr functional_enrichment_preparation
# In the original version of the script (no HTOs) this loop applied to a list of
# data.frames (one for each cluster).
# Here we use the same structure but nested in a first 'comparison' level list.
# This chunk extracts genes names from enrichment nested lists, converts them to 
# 'entrez' IDS, prepare a general conversion dictionnary, and remove eventual 
# empty list elements.

# Convert the list of marker genes to ENTREZ IDs for functional analyses
genesDF_nestedList = lapply(topDiffExp_nestedList, lapply, "[[", "gene")

# Prepare list of all unique marker gene names to be converted (here the union 
# of all DE genes across conditions).
allGenes = unique( unlist( genesDF_nestedList))

# Convert unique SYMBOLs to ENTREZ IDs
universeConversion = bitr( allGenes, 
                           fromType="SYMBOL", 
                           toType=c("ENTREZID"), 
                           OrgDb=FUNCTIONAL_ORGANISM_LIBRARY)

# Detect genes for which conversion did not succeed
conversionFailed = !(allGenes %in% universeConversion[["SYMBOL"]]) 
if(any( conversionFailed)) 
{
  warning( paste( "Some DE genes symbol could not be converted to ENTREZ ID for functional analysis :", paste(allGenes[conversionFailed], collapse = " - ")))
}

# Filter genes for which conversion failed
genesNames_nestedList = lapply( genesDF_nestedList, lapply, function(genesIDs)
{
  return(genesIDs[! genesIDs %in% allGenes[conversionFailed]])
})

# Remove eventual empty categories / columns
genesNames_nestedList = lapply(genesNames_nestedList, function(x)
{
  x[ sapply(x, length) > 0 ]
})

# Prepare a named vector (dictionnary) for efficient conversion (and to eventually be used as universe for enrichment analysis)
universe = universeConversion[["ENTREZID"]]
names(universe) = universeConversion[["SYMBOL"]]
# And a reverse dictionnary for converting back to genes names after analysis
universeRev = names(universe)
names(universeRev) = universe

# Make a dataframe of converted genes IDs
genesNames_nestedList_entrez = lapply( genesNames_nestedList, lapply, function(genesName)
{
  return( universe[genesName])
}) 




# @knitr functional_kegg_pathway
# Perform the Enrichment analysis for each comparison and store it in 
# corresponding element of nested list.
# Requires UNIVERSE_IS_UNIONGENES logical value

# Make the KEGG (pathway, not MODULE) enrichment analysis for lists of DE genes (by comparison then identity)
enrichmentList_KEGG = lapply(genesNames_nestedList_entrez, function(list_currentComparison)
{
  
  # Return a named list element so we can merge the lists later and maintain 
  # identification of enrichment type
  return(list("K-Pathways" = lapply( list_currentComparison,
                               function(entrezIDs_currentIdentity)
                               {
                                 # Prepare a data frame with expected columns from enrichment function as it can return various formats (DF, empty DF, NULL)
                                 enrichRes = data.frame( matrix( ncol = 9, nrow = 0, dimnames = list( NULL, c( "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"))))
                                 
                                 message( paste0( "KEGG Pathway enrichment analysis..."));
                                 message( paste(entrezIDs_currentIdentity, collapse = " - "));
                                 
                                 # Append enrichment result to created data.frame (helps to keep expected format/colnames regardless the returned value)
                                 enrichRes = rbind( enrichRes, as.data.frame( enrichKEGG( gene          = entrezIDs_currentIdentity,
                                                                                          organism      = FUNCTIONAL_ORGANISM_CODE,
                                                                                          pAdjustMethod = KEGG_PADJUST_METHOD,
                                                                                          pvalueCutoff  = KEGG_PVALUE_CUTOFF,
                                                                                          universe      = if(UNIVERSE_IS_UNIONGENES) universe)))
                                 
                                 if(nrow( enrichRes))
                                 {
                                   # Convert back ENTREZ ID returned by enrichment function to genes names
                                   enrichRes[["geneID"]] = sapply( strsplit( enrichRes[["geneID"]], "/"), 
                                                                   function(entrezIDs)
                                                                   {
                                                                     paste( universeRev[entrezIDs], collapse = ",") # Default expected separator for plotting function (see RFutils::GO_plotClusterAndGenes)
                                                                   })
                                 }
                                 
                                 return(enrichRes)
                               })))
  
})




# @knitr functional_kegg_module
# Perform the Enrichment analysis for each comparison and store it in 
# corresponding element of nested list.
# Requires UNIVERSE_IS_UNIONGENES logical value

# Make the KEGG MODULE enrichment analysis for lists of DE genes (by comparison then identity)
enrichmentList_MODULE = lapply(genesNames_nestedList_entrez, function(list_currentComparison)
{
  
  # Return a named list element so we can merge the lists later and maintain 
  # identification of enrichment type
  return(list("K-Modules" = lapply( list_currentComparison,
                                 function(entrezIDs_currentIdentity)
                                 {
                                   # Prepare a data frame with expected columns from enrichment function as it can return various formats (DF, empty DF, NULL)
                                   enrichRes = data.frame( matrix( ncol = 9, nrow = 0, dimnames = list( NULL, c( "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"))))
                                   
                                   message( paste0( "KEGG Module enrichment analysis..."));
                                   message( paste(entrezIDs_currentIdentity, collapse = " - "));
                                   
                                   # Append enrichment result to created data.frame (helps to keep expected format/colnames regardless the returned value)
                                   enrichRes = rbind( enrichRes, as.data.frame( enrichMKEGG( gene          = entrezIDs_currentIdentity,
                                                                                             organism      = FUNCTIONAL_ORGANISM_CODE,
                                                                                             pAdjustMethod = KEGG_PADJUST_METHOD,
                                                                                             pvalueCutoff  = KEGG_PVALUE_CUTOFF,
                                                                                             universe      = if(UNIVERSE_IS_UNIONGENES) universe)))
                                   
                                   if(nrow( enrichRes))
                                   {
                                     # Convert back ENTREZ ID returned by enrichment function to genes names
                                     enrichRes[["geneID"]] = sapply( strsplit( enrichRes[["geneID"]], "/"), 
                                                                     function(entrezIDs)
                                                                     {
                                                                       paste( universeRev[entrezIDs], collapse = ",") # Default expected separator for plotting function (see RFutils::GO_plotClusterAndGenes)
                                                                     })
                                   }
                                   
                                   return(enrichRes)
                                 })))
  
})



## @knitr functional_go

#### Make the GO analysis against selected GO categories (CC, MF, BP) for markers of each cluster, and store in a list.
enrichmentList_GO = lapply(genesNames_nestedList_entrez, 
                           function(list_currentComparison)
                           {
                             lapply( GO_CATEGORIES, function(currentCategory)
                             {
                               # Loop on samples/cluster (nested list)
                               lapply( list_currentComparison, 
                                       function(entrezIDs_currentIdentity) # Use ENTREZ ID previously computed for KEGG (SYMBOL could be ok but not recommended).
                                       {
                                         # Prepare a data frame with expected columns from enrichment function as it can return various formats (DF, empty DF, NULL)
                                         enrichResDF = data.frame( matrix( ncol = 9, nrow = 0, dimnames = list( NULL, c( "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"))))
                                         
                                         # Do the enrichment analysis (returns an object 'enrichResult')
                                         enrichRes = enrichGO( gene          = entrezIDs_currentIdentity,
                                                               OrgDb         = FUNCTIONAL_ORGANISM_LIBRARY,
                                                               ont           = currentCategory,
                                                               pAdjustMethod = GO_PADJUST_METHOD,
                                                               pvalueCutoff  = GO_PVALUE_CUTOFF,
                                                               universe      = if(UNIVERSE_IS_UNIONGENES) universe)
                                         
                                         # Eventually remove redundant terms using dedicated 'simplify' function
                                         if(!is.null( GO_SIMPLIFY_CUTOFF))
                                         {
                                           enrichRes = simplify( enrichRes,
                                                                 cutoff = GO_SIMPLIFY_CUTOFF)
                                         }
                                         
                                         # Append enrichment result to created data.frame (helps to keep expected format/colnames regardless the returned value)
                                         enrichResDF = rbind( enrichResDF, as.data.frame( enrichRes))
                                         
                                         if(nrow( enrichResDF))
                                         {
                                           # Convert back ENTREZ ID returned by enrichment function to genes names
                                           enrichResDF[["geneID"]] = sapply( strsplit( enrichResDF[["geneID"]], "/"), 
                                                                             function(entrezIDs)
                                                                             {
                                                                               paste( universeRev[entrezIDs], collapse = ",") # Default expected separator for plotting function (see RFutils::GO_plotClusterAndGenes)
                                                                             })
                                         }
                                         
                                         return(enrichResDF)
                                       })
                               
                             })
                           })



## @knitr diffExp_plots
# Merge results and plot differential expression table and enrichment analyses 
# results 

# Save original 'par' to restore it later (GO_plotClusterAndGenes alters it and does not restore it properly after plot)
initialPar = par()

# Merge the different categories of enrichment analyses in a single list
enrichmentList = Map(c, enrichmentList_KEGG, enrichmentList_MODULE, enrichmentList_GO)
# Clean the enrichment list (removes empty enrichment results and empty categories)
enrichmentList = lapply(enrichmentList, GO_checkArgument_tableResultList) # 'GO_*' functions initially made for topGO results (see RFutils)

# Main loop on comparisons (as defined in 'AnalysesParams.R')
for(currentComparisonName in names(topDiffExp_nestedList))
{
  
  cat( "\n\n### ", currentComparisonName, " {.tabset .tabset-fade} \n\n", sep="")
  
  
  #### Plotting  DEGS count summary (before/after filter)
  
  # Count DEGs passing PValue threshold
  countDEG = do.call( rbind, 
                      lapply( resultsDE_nestedList_notFiltered[[currentComparisonName]], # '_notFiltered' contains DE results for all tested genes
         function(x)
         {
           categoryDEG = factor( ifelse( x[["p_val_adj"]] < HTO_DIFFEXP_PVAL_THR, 
                                         ifelse( x[["avg_log2FC"]]>0, 
                                                 "Pos", 
                                                 "Neg"), 
                                         "Not DE"), 
                                 levels = c("Pos", "Neg", "Not DE"))
           table(categoryDEG)
         }))
  
  # Report as formatted table
  cat("\n\n");
  print( kable( addmargins( countDEG, margin = 2),
                align   = "c",
                caption = "Differential expression analysis statistics (before eventual filtering). See Volcano/MA plots."))
  

  # Count DEGs after filtering (+ eventual top selection)
  countUsedDEG = do.call( rbind, 
                          lapply( topDiffExp_nestedList[[currentComparisonName]], # Results already filtered based on PValue threshold (+ top selection), no need to recompute as above.
                                  function(x)
                                  {
                                    categoryDEG = factor( ifelse( x[["avg_log2FC"]]>0, 
                                                                  "Pos", 
                                                                  "Neg"),
                                                          levels = c("Pos", "Neg"))
                                    table(categoryDEG)
                                  }))
  
  # Report as formatted table
  cat("\n\n");
  print( kable( addmargins( countUsedDEG, margin = 2),
                align   = "c",
                caption = paste0("Summary of DEG genes selected for enrichments analysis and report table (", if(!is.null( HTO_DIFFEXP_SELECT_TOP)) paste0("top ", HTO_DIFFEXP_SELECT_TOP, ")") else "no 'top' filter applied)", ". Full list in external csv.")))
  
  
  
  
  #### Plotting DEG TABLE
  
  # Render the table of differentialy expressed genes
  cat("\n\n");
  
  # Create a copy of topDiffExp_nestedList that will eventually be filtered for
  # reporting as datatable.
  currentTopDiffExp_listForDatatable = topDiffExp_nestedList[[currentComparisonName]]
  if(!is.null(HTO_DIFFEXP_TOP_DATATABLE))
  {
    currentTopDiffExp_listForDatatable = lapply(currentTopDiffExp_listForDatatable, function(x)
      {
        if(is.null(x)) return( NULL) # Should not happen as topDiffExp_nestedList is filtered earlier
        
        # Sort by adjusted PValues before selecting top genes
        head( x[order( x[["p_val_adj"]]), , drop = FALSE], 
              HTO_DIFFEXP_TOP_DATATABLE)
      })
  }
  
  # Create a data frame for all clusters/Identities (stored as list elements),
  # adding a "Group" column from list names to trace origin.
  groupedDEG = data.frame( "Group" = factor(rep( names( currentTopDiffExp_listForDatatable), 
                                                 sapply( currentTopDiffExp_listForDatatable, nrow))),
                           do.call( rbind, currentTopDiffExp_listForDatatable));
  
  #message( paste( "NB DEGs:", nrow(groupedDEG)))
  
  # Transform gene name to HTML with genecards link
  groupedDEG[["gene"]] = paste0("<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=", groupedDEG[["gene"]], "' target='_blank' rel='noopener noreferrer'>", groupedDEG[["gene"]], "</a>")
  
  
  if(!FLAG_secondRun) # check if this is a second run for rendering only (see analysisParams)
  {
    cat("\n\nFirst (computing) run. Run report generation a second time to get all tables/figures and proper layout...")
  } else
  {
      
    # Create datatable
    summaryTable = datatable( groupedDEG,
                              class = "compact",
                              filter="top",
                              rownames = FALSE,
                              #colnames = ,
                              caption = if(! (is.null(HTO_DIFFEXP_TOP_DATATABLE) && is.null(HTO_DIFFEXP_SELECT_TOP)) ) paste0("Top ", min( HTO_DIFFEXP_TOP_DATATABLE, HTO_DIFFEXP_SELECT_TOP), " DEGs by cluster/identity"),
                              extensions = c('Buttons', 'Select'),
                              options = list(dom = "<'row'<'col-sm-8'B><'col-sm-4'f>> <'row'<'col-sm-12'l>> <'row'<'col-sm-12'rt>> <'row'<'col-sm-12'ip>>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                                             autoWidth = FALSE,
                                             buttons = exportButtonsListDT,
                                             columnDefs = list(
                                               list( # Center all columns except first one
                                                 targets = 1:(ncol( groupedDEG)-1),
                                                 className = 'dt-center'),
                                               list( # Set renderer function for 'float' type columns (LogFC)
                                                 targets = c(3,4,5),
                                                 render = htmlwidgets::JS( "function ( data, type, row, meta ) {return type === 'export' || (!Number.isFinite(data)) ? data : data.toPrecision(4);}")),                                              
                                               list( # Set renderer function for 'scientific' type columns (PValue)
                                                 targets = c(2,6), # pvalue, adj. pvalue
                                                 render = htmlwidgets::JS( "function ( data, type, row, meta ) {return type === 'export' || (!Number.isFinite(data)) ? data : data.toExponential(4);}"))),
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
                                             stateSave = TRUE),
                              escape = FALSE) %>%
      # Add  bar relative to logFC
      formatStyle( columns = "avg_log2FC",
                   backgroundColor = styleInterval(seq( -max( abs( groupedDEG[["avg_log2FC"]])), max( abs(groupedDEG[["avg_log2FC"]])), length.out = 99), colorRampPalette(c("blue", "white", "darkred"))(100)))
    
    
    print( htmltools::tagList(summaryTable))
    
    cat("\n\n");
  
    
    
    
    #### Plotting VOLCANO + MA
  
    for(currentIdentity in names(topDiffExp_nestedList[[currentComparisonName]]))
    {
      cat("\n\n#### <span style='border-radius: 3px; border: 3px solid #7777CC; padding:0px 2px'> DE stats:", currentIdentity, "</span> \n\n")

      # Get the non-filtered result of DE analysis  for current condition and identity
      DEdf_noFilter = resultsDE_nestedList_notFiltered[[currentComparisonName]][[currentIdentity]]
      
      #DEdf_noFilter[["p_val_BH"]] = p.adjust(DEdf_noFilter[["p_val"]], method = "BH")
      
      
      ## Volcano-plot
      
      # Set gene name as an actual column
      DEdf_noFilter[["name"]] = rownames(DEdf_noFilter)
      # Define categories Pos/Neg based on specified threshold
      DEdf_noFilter[["category"]] = factor( ifelse( DEdf_noFilter[["p_val_adj"]] < HTO_DIFFEXP_PVAL_THR, 
                                                    ifelse( DEdf_noFilter[["avg_log2FC"]]>0, 
                                                            "Pos", 
                                                            "Neg"), 
                                                    "Not DE"), 
                                            levels = c("Pos", "Neg", "Not DE"))
      
      # Report name as label for top5 (based on avg PValue) of each Pos/Neg category
      DEdf_noFilter[["labelDE"]] = ""
      for(currentCat in c("Pos", "Neg"))
      {
        # Select rows of current category
        currentSelection = DEdf_noFilter[ DEdf_noFilter[["category"]] == currentCat, , drop = FALSE]
        if(nrow(currentSelection) > 0)
        {
          selectedLabels = head( currentSelection[ order(currentSelection[["p_val_adj"]]), "name"], 10)
          DEdf_noFilter[selectedLabels, "labelDE"] = selectedLabels;
        }
      }
      
      # Prepare graphic adjustment values
      categoryColors = c("red", "blue", "black") #c("darkred", "darkblue", "black")
      names(categoryColors) = levels( DEdf_noFilter[["category"]])
      categorySizes = c(1, 1, 0.3)
      names(categorySizes) = levels( DEdf_noFilter[["category"]])
      categoryNames = paste0( levels( DEdf_noFilter[["category"]]), " (", table( DEdf_noFilter[["category"]]), ")")
      names(categoryNames) = levels( DEdf_noFilter[["category"]])
      
      # Create base for volcano plot
      volcano_base = ggplot( DEdf_noFilter, 
                             aes( x = -log10(p_val), # avg_log2FC
                                  y = avg_log2FC, # -log10(p_val)
                                  color = category,
                                  size  = category)) +
        geom_point( alpha = 0.6) +
        scale_color_manual( categoryNames, values = categoryColors, name = NULL) +
        scale_size_manual( values = categorySizes, guide = "none") +
        labs( y = "Average log2 FC", x = "-log10 P-Value") + # labs( x = "Average log2 FC", y = "-log10 P-Value") +
        ggtitle( "Volcano-plot (threshold on adjusted P-value)") +
        theme_classic()
      
      
      volcano = NULL
      # Add labels for top identified top 5 genes of each Pos/Neg category
      if(any( DEdf_noFilter[["labelDE"]] != "" ))
      {
        volcano = volcano_base +
          geom_text_repel( data = DEdf_noFilter[ DEdf_noFilter[["category"]] %in% c("Pos", "Neg"), ], # Select all DEGs to repel on points ("" label)
                           aes(label = labelDE),
                           size = 4,
                           #xlim = c(-0.2, 0.2),
                           #nudge_x = ifelse(DEdf_noFilter[ DEdf_noFilter[["category"]] %in% c("Pos", "Neg"), "avg_log2FC"]>0, 0.5, -0.5),
                           #nudge_x = -DEdf_noFilter[ DEdf_noFilter[["category"]] %in% c("Pos", "Neg"), "avg_log2FC"]*.25,
                           #nudge_x = DEdf_noFilter[ DEdf_noFilter[["category"]] %in% c("Pos", "Neg"), "avg_log2FC"]*1.25,
                           force = 15,
                           min.segment.length = 0,
                           show.legend = FALSE,
                           segment.size = 0.15)
      }else 
      {
        volcano = volcano_base
      }
      
      
      
      ## MA-plot
      
      # Compute the average expression of detected genes within the cells selections used for testing
      expressionValuesMatrix = do.call(cbind, lapply( selectedCells_nestedList[[currentComparisonName]][[currentIdentity]], 
                                                      function(cellsID)
                                                      {
                                                        apply( as.matrix( GetAssayData( sc10x, 
                                                                                        slot = "data", 
                                                                                        assay="RNA")[ DEdf_noFilter[["name"]], cellsID, drop = FALSE]), # normalized expression values
                                                               1, 
                                                               mean)
                                                      })) 
      # Add average expression value to initital df
      DEdf_noFilter[["avgExp"]] = rowMeans(expressionValuesMatrix)
      
      # Create base for MA plot
      MA_base = ggplot( DEdf_noFilter,
                        aes( x = avgExp,
                             y = avg_log2FC,
                             color = category,
                             size  = category)) +
        geom_point( alpha = 0.6) +
        scale_color_manual( categoryNames, values = categoryColors, name = NULL) +
        scale_size_manual( values = categorySizes, guide = "none") +
        labs( x = "Average expression", y = "Average log2 FC") +
        ggtitle( "MA-plot") +
        theme_classic()
      
      MA = NULL
      # Add labels for top identified top genes of each Pos/Neg category
      if(any( DEdf_noFilter[["labelDE"]] != "" ))
      {
        MA = MA_base +
          geom_text_repel( data = DEdf_noFilter[ DEdf_noFilter[["category"]] %in% c("Pos", "Neg"), ], # Select all DEGs to repel on points ("" label)
                           aes(label = labelDE),
                           size = 4,
                           force = 15,
                           min.segment.length = 0,
                           show.legend = FALSE,
                           segment.size = 0.15)
      } else 
      {
        MA = MA_base
      }
        
      # Combine volcano and MA plot in a single panel and plot
      print( (volcano / MA) + 
               plot_annotation( title = paste0(currentComparisonName, " / ",  currentIdentity)))
      cat("\n\n")
      
      
      ## Umap of cells used for DE
      
      # Get list of cells selected for current DE analysis
      cellsList = selectedCells_nestedList[[currentComparisonName]][[currentIdentity]] 
      
      categoryColorsUMAP = c("grey", "orange", "darkgreen")
      names(categoryColorsUMAP) = c("Unselected", names(cellsList))
      categoryNames = c("Unselected", paste0( names(cellsList), " (", sapply( cellsList, length), ")"))
      names(categoryNames) = c("Unselected", names(cellsList))
      
      umapCellsDE = DimPlot( sc10x, 
                             cells.highlight = selectedCells_nestedList[[currentComparisonName]][[currentIdentity]],
                             label = TRUE,
                             label.size = 3,
                             label.color = "black",
                             label.box = TRUE,
                             repel = TRUE) +
                      ggtitle(paste0("Cells selection for DE analysis: ", currentComparisonName, " / ", currentIdentity, " -> ", sum(sapply( cellsList, length)))) +
                      scale_color_manual( labels = categoryNames, values = categoryColorsUMAP)
      
      print( umapCellsDE)
      cat("\n\n")
    
      
      ## Volcano/MA with specific genes highlighted 
      
      if( exists("MONITORED_GENES") && (length(MONITORED_GENES)>0) )
      {
        
        # Make a copy of non-filtered DE data_frame and set label/category for lists of genes
        DEdf_noFilter_monitored = DEdf_noFilter
        # Add a column that will store the name of group and one for corresponding labels
        DEdf_noFilter_monitored[["monitoredGroup"]] = ""
        DEdf_noFilter_monitored[["labelMonitor"]] = ""
          
        # Set category and label to genes of interest
        for(currentCategory in names(MONITORED_GENES))
        {
          DEdf_noFilter_monitored[MONITORED_GENES[[currentCategory]], "monitoredGroup"] = currentCategory
          DEdf_noFilter_monitored[MONITORED_GENES[[currentCategory]], "labelMonitor"] = MONITORED_GENES[[currentCategory]]
        }
        
        # Replace black by darkgrey for text color
        categoryColors = c("red", "blue", "darkgrey") 
        names(categoryColors) = levels( DEdf_noFilter[["category"]])
        
        # Create new graphic adjustment values (use of 'new_scale_color' from ggnewscale)
        categoryColors_monitored = c("darkorange", "darkgreen", "magenta4", brewer.pal(n = length(MONITORED_GENES), name = "Set3"))[1:length(MONITORED_GENES)]#c("gold1", "springgreen2", brewer.pal(n = length(MONITORED_GENES), name = "Set3"))[1:length(MONITORED_GENES)]
        names(categoryColors_monitored) = names(MONITORED_GENES)
        
        # Reuse volcano base figure and add labels annotations
        volcano = volcano_base +
          new_scale_color() + # Reset color scale for points
          geom_label_repel(  data = DEdf_noFilter_monitored[DEdf_noFilter_monitored[["monitoredGroup"]] %in% names(MONITORED_GENES),],
                            aes(label = labelMonitor, color = monitoredGroup ),
                            fill = "white",
                            size = 4,
                            force = 5,
                            min.segment.length = 0,
                            show.legend = FALSE,
                            segment.size = 0.3) +
          geom_point(data = DEdf_noFilter_monitored[DEdf_noFilter_monitored[["monitoredGroup"]] %in% names(MONITORED_GENES),],
                     pch = 1, size = 2.5, aes(color = monitoredGroup), stroke=1.5) +
          scale_color_manual( names(categoryColors_monitored), values = categoryColors_monitored, name = NULL)
        # volcano
        
        
        # volcano = volcano_base +
        #   geom_label_repel( data = DEdf_noFilter_monitored[DEdf_noFilter_monitored[["monitoredGroup"]] %in% names(MONITORED_GENES),],
        #                     aes(label = labelMonitor, fill = monitoredGroup, color = category ),
        #                     size = 4,
        #                     force = 5,
        #                     min.segment.length = 0,
        #                     show.legend = FALSE,
        #                     segment.size = 0.3) +
        #   scale_fill_manual( names(categoryColors_monitored), values = categoryColors_monitored, name = NULL)+
        #   new_scale_color() + # Reset color scale for points
        #   geom_point(data = DEdf_noFilter_monitored[DEdf_noFilter_monitored[["monitoredGroup"]] %in% names(MONITORED_GENES),],
        #              pch = 1, size = 2.5, aes(color = monitoredGroup), stroke=1.5) +
        #   scale_color_manual( names(categoryColors_monitored), values = categoryColors_monitored, name = NULL) 
        
        
        
        # Reuse MA base figure and add labels annotations
        MA = MA_base +
          new_scale_color() + # Reset color scale for points
          geom_label_repel(  data = DEdf_noFilter_monitored[DEdf_noFilter_monitored[["monitoredGroup"]] %in% names(MONITORED_GENES),],
                             aes(label = labelMonitor, color = monitoredGroup ),
                             fill = "white",
                             size = 4,
                             force = 5,
                             min.segment.length = 0,
                             show.legend = FALSE,
                             segment.size = 0.3) +
          geom_point(data = DEdf_noFilter_monitored[DEdf_noFilter_monitored[["monitoredGroup"]] %in% names(MONITORED_GENES),],
                     pch = 1, size = 2.5, aes(color = monitoredGroup), stroke=1.5) +
          scale_color_manual( names(categoryColors_monitored), values = categoryColors_monitored, name = NULL)
        # MA
        
        # Combine volcano and MA plot in a single panel and plot
        print( (volcano / MA) + 
                 plot_annotation( title = paste0(currentComparisonName, " / ",  currentIdentity, " - Monitored Genes")))
        cat("\n\n")
      }
    }
  }
  
  
  #### Plotting HEATMAPS
    
  # Function plotting heatmaps of selected matrix in tabs with variations on rows organisation and scaling
  # Creates external png in first run, and includes it in report in secind run to allow dynamic image size
  # without altering final layout with 'png' and 'dev.off' calls during rendering...
  # Used below to show DE genes for all comparisons
  # DEBUG: avgExpressionMatrix=expressionValuesMatrix; basePath = NULL; scaleArgs = c("none", "row"); baseMain = ""; columnsClustering = FALSE; columnsCategories = NULL; columnsCategoriesColor = list(); splitColsByCatIndex = if(!is.null(columnsCategories)) 1 else NULL; rotateColsCat = FALSE; rowsClustering = FALSE; rowsCategories = NULL; rowsCategoriesColor = list(); slpitRowsByFC = FALSE; rotateRowsCat = FALSE; width = 900; height = 400+(15*nrow(avgExpressionMatrix))
  plotHeatmapsDE = function(avgExpressionMatrix, 
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
                            height = 400+(15*nrow(avgExpressionMatrix)))
  {
    
    # Prepare rows and columns annotation if corresponding dataframe provided
    rowsAnnot = NULL;
    if(!is.null(rowsCategories))
    {
      rowsAnnot = rowAnnotation( df = rowsCategories,
                                 col = rowsCategoriesColor, # Random colors if NULL 
                                 annotation_name_rot = 75)
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
      cat(paste0("\n\n###### Rows ordering: logFC ==== Rows scaling: ",scaleArg, " {.tabset .tabset-fade .tabset-pills} \n\n"))
      
      main = paste0(baseMain, " - order:logFC - scaling:", scaleArg)
      
      # Render the plot with ComplexHeatmap
      hm = Heatmap( if(scaleArg == "rows") t(scale(t(avgExpressionMatrix))) else avgExpressionMatrix, 
                    name = "Average\nExpression", 
                    #col = RColorBrewer::brewer.pal(name = "Reds", n = 4),
                    col = rev(rocket(20)),
                    rect_gp = gpar(col = "white", lwd = 1), # Border of cells
                    column_names_rot = 75,
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
        filenamePNG = paste0( basePath, "_ORDER_logFC_SCALE_", scaleArg, ".png")
        
        if(FLAG_secondRun) # Check if this is a second run for rendering only (see analysisParams)
        {
          # Don't mess with RMD layout using 'png' and 'dev.off'...
          # Integrate the figure that was generated during first run
          cat( paste0( "\n![Heatmap ORDER:logFC SCALE:", scaleArg, "](", filenamePNG,")\n" ))
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
        cat(paste0("\n\n###### Rows ordering: hClust ==== Rows scaling: ",scaleArg, " {.tabset .tabset-fade .tabset-pills} \n\n"))

        main = paste0(baseMain, " - order:hClust - scaling:", scaleArg)
        
        # Render the plot with ComplexHeatmap
        hm = Heatmap( if(scaleArg == "rows") t(scale(t(avgExpressionMatrix))) else avgExpressionMatrix, 
                      name = "Average\nExpression", 
                      #col = RColorBrewer::brewer.pal(name = "Reds", n = 4),
                      col = rev(rocket(20)),
                      rect_gp = gpar(col = "white", lwd = 1), # Border of cells
                      column_names_rot = 75,
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
        cat(paste0("\n\n###### Rows ordering: avgExp ==== Rows scaling: ",scaleArg, " {.tabset .tabset-fade .tabset-pills} \n\n"))
        
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
                      column_names_rot = 75,
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
  
  # Use above function for plotting heatmaps for each comparison and identity
  for(currentIdentity in names(topDiffExp_nestedList[[currentComparisonName]]))
  {
    # Create output folder for heatmaps as external files (png) 
    pathDiffExpOutput = file.path(PATH_ANALYSIS_OUTPUT, "DiffExpByHTO", currentComparisonName, paste0("heatmaps_", currentIdentity))
    dir.create(pathDiffExpOutput, recursive = TRUE, showWarnings = FALSE)
    
    # Get result of DEG analysis for current comparison/identity
    currentGenesDE_DF = topDiffExp_nestedList[[currentComparisonName]][[currentIdentity]]
    nbCurrentGenes = nrow(currentGenesDE_DF)
    
    message( paste0( "NB DEGs in ", currentComparisonName, "/", currentIdentity, ": ", nbCurrentGenes))
    
    # Eventually select the top (by adjusted PValue, see PLOT_HEATMAPS_MAXGENES)
    currentGenesDE_DF = currentGenesDE_DF[order(currentGenesDE_DF[["p_val_adj"]]), ]
    if(!is.null(PLOT_HEATMAPS_MAXGENES))
    {
      currentGenesDE_DF = head(currentGenesDE_DF, PLOT_HEATMAPS_MAXGENES)
    }
    
    # Order by "avg_log2FC" (for convenience of first plot)
    currentGenesDE_DF = currentGenesDE_DF[order(currentGenesDE_DF[["avg_log2FC"]], decreasing = TRUE), ]
    
    cat("\n\n#### <span style='border-radius: 3px; border: 3px solid #77CC77; padding:0px 2px'> Heatmap:", currentIdentity, if((!is.null(PLOT_HEATMAPS_MAXGENES)) && nbCurrentGenes>PLOT_HEATMAPS_MAXGENES) paste0("(top", PLOT_HEATMAPS_MAXGENES, "/", nbCurrentGenes, ")"),"</span> {.tabset .tabset-fade}\n\n")
    
    ## Prepare rows annotations (common for all plots)
    rowsCategories = currentGenesDE_DF[,"avg_log2FC", drop = FALSE] # Only LogFC for now
    colnames(rowsCategories) = c("Log2 FC")
    # Prepare corresponding color scales (see ComplexHeatmap annotations)
    colBreaksLogFC = rev(RColorBrewer::brewer.pal(5,"RdBu"))[c(1,3,5)] # Get values for a colorblind-friendly symmetric red/blue scale 
    maxAbsLogFC = max( abs( rowsCategories[["Log2 FC"]]))
    colfunLogFC = colorRamp2( c( -maxAbsLogFC, 0, maxAbsLogFC), colBreaksLogFC) # Create the interpolation function used by ComplexHeatmap
    rowsCategoriesColor = list("Log2 FC" = colfunLogFC)
    
    ## 01
    message("# Plots 01: DE cells groups")
    cat("\n\n##### DE cells groups {.tabset .tabset-fade .tabset-pills}\n\n")
    # Compute the average expression of detected genes within the cells selections used for testing
    expressionValuesMatrix01 = do.call(cbind, lapply( selectedCells_nestedList[[currentComparisonName]][[currentIdentity]], 
                                                    function(cellsID)
                                                    {
                                                      apply( as.matrix( GetAssayData( sc10x, 
                                                                                      slot = "data", 
                                                                                      assay="RNA")[currentGenesDE_DF[["gene"]], cellsID, drop = FALSE]), # normalized expression values
                                                             1, 
                                                             mean)
                                                    }))
    # Render the plots with graphic variations in tabs (scaling, ordering)
    # Creates external png and includes it in report to allow dynamic image size
    plotHeatmapsDE( avgExpressionMatrix = expressionValuesMatrix01, 
                    basePath = file.path(pathDiffExpOutput, "DE_CellsGroups"),
                    scaleArgs = c("none"), # Scaling on rows irrelevant for 2 columns only 
                    baseMain = paste0(currentComparisonName, " (", currentIdentity, ")"), 
                    columnsClustering = FALSE, # Only 2 columns...
                    columnsCategories = NULL, 
                    columnsCategoriesColor = list(), 
                    splitColsByCatIndex = NULL, # Index of category for which to split the columns. NULL to ignore.
                    rotateColsCat = NULL, # 0 or 90, NULL to hide columns categories names
                    #rowsClustering = FALSE, # Fixed, depends on plots made in function
                    rowsCategories = rowsCategories, # Rows annotations: mostly for logFC (col name "Log2 FC", continuous)
                    rowsCategoriesColor = rowsCategoriesColor, 
                    splitRowsByFC = TRUE, # Split pos and neg LogFC. Only works if a rows category is named "Log2 FC". Cannot choose which category as for columns because a boolean operation is made on continuous log FC.
                    rotateRowsCat = 90, # 0 or 90, NULL to hide rows categories names
                    width = max(c(250, 160+(20*ncol(expressionValuesMatrix01)))), # logFC(60)+legend(100)+cells (adds 50 internally when doing clustering)
                    height = max(c(230,130+(20*nrow(expressionValuesMatrix01))))) # colnames(130)+cells (min size 150 to allow vertical layout of legends)
    
    
    
    ## 02 & 03
    # Get the precomputed average expression of detected genes for all clusters/HTOs combinations
    expressionValuesMatrix0203 = geneExpByClusterAndHTO[currentGenesDE_DF[["gene"]], , drop = FALSE]
    # Extract groups from colnames (created from factors HTO and Idents) for column annotation
    columnsCategories0203 = data.frame(do.call(rbind, strsplit(x = colnames(expressionValuesMatrix0203), split = "_")))
    # Convert to factor and reorder columns according to known factor levels
    colnames(columnsCategories0203) = c("Cluster", "HTO")
    columnsCategories0203[["Cluster"]] = factor(columnsCategories0203[["Cluster"]])
    columnsCategories0203[["HTO"]] = factor(columnsCategories0203[["HTO"]], levels = HTO_FACTOR_LEVELS)
    orderColumns = order(columnsCategories0203[["Cluster"]], columnsCategories0203[["HTO"]])
    expressionValuesMatrix0203 = expressionValuesMatrix0203[, orderColumns, drop = FALSE]
    columnsCategories0203 = columnsCategories0203[orderColumns, , drop = FALSE]
    
    message("# Plots 02: All clusters & HTOs (by cluster)")
    cat("\n\n##### All clusters & HTOs (by cluster) {.tabset .tabset-fade .tabset-pills}\n\n")
    # Here make columns groups from category "Cluster" before clustering
    plotHeatmapsDE( avgExpressionMatrix = expressionValuesMatrix0203, 
                    basePath = file.path(pathDiffExpOutput, "AllClustersAndHTOs_CatByCluster"),
                    scaleArgs = c("none", "rows"), # Generate figures with and without rows scaling 
                    baseMain = paste0(currentComparisonName, " (", currentIdentity, ")"), 
                    columnsClustering = TRUE, 
                    columnsCategories = columnsCategories0203, 
                    columnsCategoriesColor = list("Cluster" = clustersColor, "HTO" = HTOsColor), 
                    splitColsByCatIndex = 1, # Index of category for which to split the columns. NULL to ignore.
                    rotateColsCat = NULL, # 0 or 90, NULL to hide columns categories names
                    #rowsClustering = FALSE, # Fixed, depends on plots made in function
                    rowsCategories = rowsCategories, # Rows annotations: mostly for logFC (col name "Log2 FC", continuous)
                    rowsCategoriesColor = rowsCategoriesColor, 
                    splitRowsByFC = TRUE, # Split Pos and Neg LogFC. Only works if a rows category is named "Log2 FC". Cannot choose which category as for columns because a boolean operation is made on continuous log FC.
                    rotateRowsCat = 90, # 0 or 90, NULL to hide rows categories names
                    width = 290+(20*ncol(expressionValuesMatrix0203)), # logFC(60)+legend(100)+legendClusters(130)+cells (adds 50 internally when doing clustering)
                    height = 300+(20*nrow(expressionValuesMatrix0203))) # colnames(200)+colcategories(50)+dendrogram(50)+cells
    
    message("# Plots 03: All clusters & HTOs (by HTO)")
    cat("\n\n##### All clusters & HTOs (by HTO) {.tabset .tabset-fade .tabset-pills}\n\n")
    # Here make columns groups from category "HTO" before clustering
    plotHeatmapsDE( avgExpressionMatrix = expressionValuesMatrix0203, 
                    basePath = file.path(pathDiffExpOutput, "AllClustersAndHTOs_CatByHTO"),
                    scaleArgs = c("none", "rows"), # Generate figures with and without rows scaling 
                    baseMain = paste0(currentComparisonName, " (", currentIdentity, ")"), 
                    columnsClustering = TRUE, 
                    columnsCategories = columnsCategories0203, 
                    columnsCategoriesColor = list("Cluster" = clustersColor, "HTO" = HTOsColor), 
                    splitColsByCatIndex = 2, # Index of category for which to split the columns. NULL to ignore.
                    rotateColsCat = 0, # 0 or 90, NULL to hide columns categories names
                    #rowsClustering = FALSE, # Fixed, depends on plots made in function
                    rowsCategories = rowsCategories, # Rows annotations: mostly for logFC (col name "Log2 FC", continuous)
                    rowsCategoriesColor = rowsCategoriesColor, 
                    splitRowsByFC = TRUE, # Split Pos and Neg LogFC. Only works if a rows category is named "Log2 FC". Cannot choose which category as for columns because a boolean operation is made on continuous log FC.
                    rotateRowsCat = 90, # 0 or 90, NULL to hide rows categories names
                    width = 340+(20*ncol(expressionValuesMatrix0203)), # logFC(60)+legend(100)+geneames(50)+legendClusters(130)+cells (adds 50 internally when doing clustering)
                    height = 340+(20*nrow(expressionValuesMatrix0203))) # colnames(200)+colCategories(50)+dendrogram(50)+colCategoriesNames(40)+cells
    
    
    ## 04 & 05 Selection of clusters and HTOs involved in DE analysis
    # Select columns relevant for current DE analysis from previously gathered expression matrix (clusters and HTOs)
    columnsSelection = columnsCategories0203[["HTO"]] %in% unlist(HTO_DIFFEXP_COMPARISONLIST[[currentComparisonName]]) & (if(currentIdentity == "AllCells") TRUE else (columnsCategories0203[["Cluster"]] == currentIdentity))
    # Skip if similar as previous (when comparison implies all clusters & HTOs) 
    if(!all(columnsSelection))
    {
      expressionValuesMatrix0405 = expressionValuesMatrix0203[, columnsSelection, drop = FALSE]
      columnsCategories0405 = columnsCategories0203[columnsSelection, , drop = FALSE]
      
      message("# Plots 04: Selected clusters & HTOs (by cluster)")
      cat("\n\n##### Selected clusters & HTOs (by cluster) {.tabset .tabset-fade .tabset-pills}\n\n")
      # Here make columns groups from category "Cluster" before clustering
      plotHeatmapsDE( avgExpressionMatrix = expressionValuesMatrix0405, 
                      basePath = file.path(pathDiffExpOutput, "SelectedClustersAndHTOs_CatByCluster"),
                      scaleArgs = c("none", "rows"), # Generate figures with and without rows scaling 
                      baseMain = paste0(currentComparisonName, " (", currentIdentity, ")"), 
                      columnsClustering = FALSE, 
                      columnsCategories = columnsCategories0405, 
                      columnsCategoriesColor = list("Cluster" = clustersColor, "HTO" = HTOsColor), 
                      splitColsByCatIndex = 1, # Index of category for which to split the columns. NULL to ignore.
                      rotateColsCat = NULL, # 0 or 90, NULL to hide columns categories names
                      #rowsClustering = FALSE, # Fixed, depends on plots made in function
                      rowsCategories = rowsCategories, # Rows annotations: mostly for logFC (col name "Log2 FC", continuous)
                      rowsCategoriesColor = rowsCategoriesColor, 
                      splitRowsByFC = TRUE, # Split Pos and Neg LogFC. Only works if a rows category is named "Log2 FC". Cannot choose which category as for columns because a boolean operation is made on continuous log FC.
                      rotateRowsCat = 90, # 0 or 90, NULL to hide rows categories names
                      width = 290+(20*ncol(expressionValuesMatrix0405)), # logFC(60)+legend(100)+legendClusters(130)+cells (adds 50 internally when doing clustering)
                      height = 250+(20*nrow(expressionValuesMatrix0405))) # colnames(200)+colcategories(50)+cells+!dendrogram(removed,50)
      
      message("# Plots 05: Selected clusters & HTOs (by HTO)")
      cat("\n\n##### Selected clusters & HTOs (by HTO) {.tabset .tabset-fade .tabset-pills}\n\n")
      # Here make columns groups from category "HTO" before clustering
      plotHeatmapsDE( avgExpressionMatrix = expressionValuesMatrix0405, 
                      basePath = file.path(pathDiffExpOutput, "SelectedClustersAndHTOs_CatByHTO"),
                      scaleArgs = c("none", "rows"), # Generate figures with and without rows scaling 
                      baseMain = paste0(currentComparisonName, " (", currentIdentity, ")"), 
                      columnsClustering = FALSE, 
                      columnsCategories = columnsCategories0405, 
                      columnsCategoriesColor = list("Cluster" = clustersColor, "HTO" = HTOsColor), 
                      splitColsByCatIndex = 2, # Index of category for which to split the columns. NULL to ignore.
                      rotateColsCat = 0, # 0 or 90, NULL to hide columns categories names
                      #rowsClustering = FALSE, # Fixed, depends on plots made in function
                      rowsCategories = rowsCategories, # Rows annotations: mostly for logFC (col name "Log2 FC", continuous)
                      rowsCategoriesColor = rowsCategoriesColor, 
                      splitRowsByFC = TRUE, # Split Pos and Neg LogFC. Only works if a rows category is named "Log2 FC". Cannot choose which category as for columns because a boolean operation is made on continuous log FC.
                      rotateRowsCat = 90, # 0 or 90, NULL to hide rows categories names
                      width = 340+(20*ncol(expressionValuesMatrix0405)), # logFC(60)+legend(100)+categoriesSpace(50)+legendClusters(130)+cells (adds 50 internally when doing clustering)
                      height = 290+(20*nrow(expressionValuesMatrix0405))) # colnames(200)+colCategories(50)+colCategoriesNames(40)+cells+!dendrogram(removed,50)
    }
    

    ## 06 Selection of clusters involved in DE analysis (but plot for all HTOs)
    # Skip if similar as previous (02 & 03: AllCells)
    if(currentIdentity != "AllCells")
    { 
      expressionValuesMatrix0607 = expressionValuesMatrix0203[, columnsCategories0203[["Cluster"]] == currentIdentity, drop = FALSE]
      columnsCategories0607 = columnsCategories0203[columnsCategories0203[["Cluster"]] == currentIdentity, , drop = FALSE]
      
      message("# Plots 06: Selected clusters (all HTOs)")
      cat("\n\n##### Selected clusters (all HTOs) {.tabset .tabset-fade .tabset-pills}\n\n")
      # Here make columns groups from category "Cluster" before clustering
      plotHeatmapsDE( avgExpressionMatrix = expressionValuesMatrix0607, 
                      basePath = file.path(pathDiffExpOutput, "SelectedClusters_AllHTOs"),
                      scaleArgs = c("none", "rows"), # Generate figures with and without rows scaling 
                      baseMain = paste0(currentComparisonName, " (", currentIdentity, ")"), 
                      columnsClustering = FALSE, 
                      columnsCategories = columnsCategories0607, 
                      columnsCategoriesColor = list("Cluster" = clustersColor, "HTO" = HTOsColor), 
                      splitColsByCatIndex = 1, # Index of category for which to split the columns. NULL to ignore.
                      rotateColsCat = NULL, # 0 or 90, NULL to hide columns categories names
                      #rowsClustering = FALSE, # Fixed, depends on plots made in function
                      rowsCategories = rowsCategories, # Rows annotations: mostly for logFC (col name "Log2 FC", continuous)
                      rowsCategoriesColor = rowsCategoriesColor, 
                      splitRowsByFC = TRUE, # Split Pos and Neg LogFC. Only works if a rows category is named "Log2 FC". Cannot choose which category as for columns because a boolean operation is made on continuous log FC.
                      rotateRowsCat = 90, # 0 or 90, NULL to hide rows categories names
                      width = 290+(20*ncol(expressionValuesMatrix0607)), # logFC(60)+legend(100)+legendClusters(130)+cells (adds 50 internally when doing clustering)
                      height = 250+(20*nrow(expressionValuesMatrix0607))) # colnames(200)+colcategories(50)+cells+!dendrogram(removed,50)
      
      message("# Plots 07: Selected clusters (all HTOs) (+clustCols)")
      cat("\n\n##### Selected clusters (all HTOs) (+clustCols) {.tabset .tabset-fade .tabset-pills}\n\n")
      # Here make columns groups from category "Cluster" before clustering
      plotHeatmapsDE( avgExpressionMatrix = expressionValuesMatrix0607, 
                      basePath = file.path(pathDiffExpOutput, "SelectedClusters_AllHTOs_clustCols"),
                      scaleArgs = c("none", "rows"), # Generate figures with and without rows scaling 
                      baseMain = paste0(currentComparisonName, " (", currentIdentity, ")"), 
                      columnsClustering = TRUE, 
                      columnsCategories = columnsCategories0607, 
                      columnsCategoriesColor = list("Cluster" = clustersColor, "HTO" = HTOsColor), 
                      splitColsByCatIndex = 1, # Index of category for which to split the columns. NULL to ignore.
                      rotateColsCat = NULL, # 0 or 90, NULL to hide columns categories names
                      #rowsClustering = FALSE, # Fixed, depends on plots made in function
                      rowsCategories = rowsCategories, # Rows annotations: mostly for logFC (col name "Log2 FC", continuous)
                      rowsCategoriesColor = rowsCategoriesColor, 
                      splitRowsByFC = TRUE, # Split Pos and Neg LogFC. Only works if a rows category is named "Log2 FC". Cannot choose which category as for columns because a boolean operation is made on continuous log FC.
                      rotateRowsCat = 90, # 0 or 90, NULL to hide rows categories names
                      width = 290+(20*ncol(expressionValuesMatrix0607)), # logFC(60)+legend(100)+legendClusters(130)+cells (adds 50 internally when doing clustering)
                      height = 300+(20*nrow(expressionValuesMatrix0607))) # colnames(200)+colcategories(50)+cells+dendrogram(50)
    }


    ## 08 All cells and all HTOs
    message("# Plots 08: Global by HTO")
    cat("\n\n##### Global by HTO {.tabset .tabset-fade .tabset-pills}\n\n")
    # Compute the average expression of detected genes within the cells selections used for testing
    expressionValuesMatrix08 = geneExpByHTO[currentGenesDE_DF[["gene"]], , drop = FALSE]
    # Render the plots with graphic variations in tabs (scaling, ordering)
    plotHeatmapsDE( avgExpressionMatrix = expressionValuesMatrix08, 
                    basePath = file.path(pathDiffExpOutput, "HTOs"),
                    scaleArgs = c("none", "rows"), # Generate figures with and without rows scaling 
                    baseMain = paste0(currentComparisonName, " (", currentIdentity, ")"), 
                    columnsClustering = FALSE, 
                    columnsCategories = NULL, 
                    columnsCategoriesColor = list(), 
                    splitColsByCatIndex = NULL, # Index of category for which to split the columns. NULL to ignore.
                    rotateColsCat = NULL, # 0 or 90, NULL to hide columns categories names
                    #rowsClustering = FALSE, # Fixed, depends on plots made in function
                    rowsCategories = rowsCategories, # Rows annotations: mostly for logFC (col name "Log2 FC", continuous)
                    rowsCategoriesColor = rowsCategoriesColor, 
                    splitRowsByFC = TRUE, # Split Pos and Neg LogFC. Only works if a rows category is named "Log2 FC". Cannot choose which category as for columns because a boolean operation is made on continuous log FC.
                    rotateRowsCat = 90, # 0 or 90, NULL to hide rows categories names
                    width = max(c(250, 160+(20*ncol(expressionValuesMatrix08)))), # logFC(60)+legend(100)+legendClusters(130)+cells (adds 50 internally when doing clustering)
                    height = max(c(230, 120+(20*nrow(expressionValuesMatrix08)))))  # colnames(120)+cells+!dendrogram(removed,50)
    
  }
  
  
  # End tabset
  #cat("\n\n### {.toc-ignore}\n\n")
  

  
  
  #### Plotting FUNCTIONAL ENRICHMENTS

  # Check whether we have some enrichments left before we save and plot them
  if(length( enrichmentList[[currentComparisonName]]))
  {
      message(paste( "FUNCTIONAL ENRICHMENTS -", currentComparisonName))

      # Format the full table of enrichments (+ add 'score') to a single dataframe (summary table) and save it to csv file
      resNULL = GO_writeEnrichmentTable( enrichmentList[[currentComparisonName]],
                                         fileName = file.path( PATH_ANALYSIS_OUTPUT,
                                                               paste0( outputFilesPrefix,
                                                                       "EnrichmentList_",
                                                                       ifelse(UNIVERSE_IS_UNIONGENES, "bg", "nobg"),
                                                                       ".tsv")),
                                         pValueColumnName = "p.adjust",
                                         pValueToScore = function(pVal){-log10(pVal)})
      
      # Filter top terms for not having massive datatable in report (full table generated as external file)
      topForDataTable = GO_filterEnrichment_firstRows( enrichmentList[[currentComparisonName]],
                                                       nbRows = TOPTERMS_DATATABLE,
                                                       sortBy = "p.adjust") 
      
      # Use this function again on filtered dataset, but just for formatting and score computation
      enrichmentDF = GO_writeEnrichmentTable( topForDataTable,
                                              fileName = NULL,
                                              pValueColumnName = "p.adjust",
                                              pValueToScore = function(pVal){-log10(pVal)})
      # Convert as factor so datatable provides a convenient filtering
      enrichmentDF[["sampleName"]] = factor(enrichmentDF[["sampleName"]])
      
      # Filter top enrichments for graphic summary representation
      topSummary = GO_filterEnrichment_firstRows( enrichmentList[[currentComparisonName]],
                                                  nbRows = TOPTERMS_SUMMARY,
                                                  sortBy = "p.adjust")
      
      # Prepare an eventual different filtering for 'AllCells' in case it must be plotted separately
      topSummary_allCells =  GO_filterEnrichment_firstRows( enrichmentList[[currentComparisonName]],
                                                            nbRows = SEPARATE_ALLCELLS_TOPTERMS,
                                                            sortBy = "p.adjust")
      
      # Plot each selected category (KEGG/Modules/BP/MF/CC) we found enrichments for
      for(currentCategoryName in names(enrichmentList[[currentComparisonName]]))
      {
        message(paste("Plotting enrichment category:", currentCategoryName))
        
        cat("\n\n#### <span style='border-radius: 3px; border: 3px solid #CC7777; padding:0px 2px'> Enrich:", currentCategoryName, "</span> {.tabset .tabset-fade .tabset-pills}\n\n")
  
        if(!FLAG_secondRun) # check if this is a second run for rendering only (see analysisParams)
        {
          cat("\n\nFirst (computing) run. Run report generation a second time to get all tables/figures and proper layout...")
        } else
        {
          # Plot the summary figure combining and organizing enrichments for all samples
          if(currentCategoryName %in% names(topSummary)) # Check if summary has something to show (empty enrichments filtered earlier)
          {
            # Eventually plot global enrichment (AllCells) separately from cluster-specific results
            if(SEPARATE_ALLCELLS_ENRICHMENTS && ("AllCells" %in% names(topSummary[[currentCategoryName]])))
            {
              # Create a temp nested list structure (required by plot function) to hold 'AllCells' data only (for currentCategoryName only)
              tempSummaryList = list()
              tempSummaryList[[currentCategoryName]] = topSummary_allCells[[currentCategoryName]]["AllCells"]
              
              message(paste("Enrichment summary ('AllCells' only) (top", SEPARATE_ALLCELLS_TOPTERMS, "terms)"))
              cat("\n\n##### Summary (AllCells) (top", SEPARATE_ALLCELLS_TOPTERMS, "terms)\n\n")
              
              # Truncate eventual very long term names to prevent error 'figure margins too large' 
              tempSummaryList = lapply( tempSummaryList, lapply, function(x)
              {
                x[["Description"]] = str_trunc(x[["Description"]], 75)
                return(x)
              })
              
              ## DEBUG: save enrichment list as RDS for external plot debug
              #enrichName = paste0(currentComparisonName, "__", currentCategoryName, "__AllCellsOnly")
              #outputPathTest = file.path(PATH_ANALYSIS_OUTPUT, paste0("tableGO__", enrichName, ".RDS"))
              #message(outputPathTest)
              #saveRDS(object = tempSummaryList, file = outputPathTest)
              #message(paste0("Now Plotting: ", enrichName))
              
              GO_plotClusterAndGenes( tempSummaryList, # Keep the nested list structure since plotting function can handle multiple categories
                                      col.names = c( featureID="ID",
                                                     featureName="Description",
                                                     pValue="p.adjust",
                                                     genesNames="geneID"),
                                      main = currentCategoryName,
                                      create.dev = NULL)
              # Restore original 'par' saved at chunk start (GO_plot function does not restore it properly)
              par(initialPar)
              
              cat(" \n \n");
            }
            
            # Create a temp nested list structure (eventually remove 'AllCells' data)
            tempSummaryList = list()
            tempSummaryList[[currentCategoryName]] = topSummary[[currentCategoryName]]
            if(SEPARATE_ALLCELLS_ENRICHMENTS) tempSummaryList[[currentCategoryName]][["AllCells"]] = NULL
            
            # Make sure we have something to plot after eventual removing of "AllCells"
            if(length(tempSummaryList[[currentCategoryName]])>0)
            {
              message(paste("Enrichment summary for all groups (top", TOPTERMS_SUMMARY, "terms)"))
              cat("\n\n##### Summary (top", TOPTERMS_SUMMARY, "terms /group)\n\n")
              
              # Truncate eventual very long term names to prevent error 'figure margins too large' 
              tempSummaryList = lapply( tempSummaryList, lapply, function(x)
              {
                x[["Description"]] = str_trunc(x[["Description"]], 75)
                return(x)
              })
              
              ## DEBUG: save enrichment list as RDS for external plot debug
              #enrichName = paste0(currentComparisonName, "__", currentCategoryName, "")
              #outputPathTest = file.path(PATH_ANALYSIS_OUTPUT, paste0("tableGO__", enrichName, ".RDS"))
              #message(outputPathTest)
              #saveRDS(object = tempSummaryList, file = outputPathTest)
              #message(paste0("Now Plotting: ", enrichName))
              
              GO_plotClusterAndGenes( tempSummaryList,
                                      col.names = c( featureID="ID",
                                                     featureName="Description",
                                                     pValue="p.adjust",
                                                     genesNames="geneID"),
                                      main = currentCategoryName,
                                      create.dev = NULL)
              # Restore original 'par' saved at chunk start (GO_plot function does not restore it properly)
              par(initialPar)
              
              cat(" \n \n");
            }
  
          
            # Draw a summary table for all terms for this category (from dataframe/csv with all categories merged)
            message(paste("Summary table (all terms)"))
            cat("\n\n##### Summary table (all terms)\n\n")
            # Select enrichment for current enrichment category
            currentEnrichmentDF = enrichmentDF[ enrichmentDF[["typeName"]] == currentCategoryName,]
            
            cat(" \n \n");
            
            # Create datatable
            summaryTable = datatable( currentEnrichmentDF,
                                      class = "compact",
                                      filter="top",
                                      rownames = FALSE,
                                      #colnames = ,
                                      #caption = ,
                                      extensions = c('Buttons', 'Select'),
                                      options = list(dom = "<'row'<'col-sm-8'B><'col-sm-4'f>> <'row'<'col-sm-12'l>> <'row'<'col-sm-12'rt>> <'row'<'col-sm-12'ip>>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                                                     autoWidth = FALSE,
                                                     buttons = exportButtonsListDT,
                                                     columnDefs = list(
                                                       list( # Center all columns except first one
                                                         targets = 1:(ncol( currentEnrichmentDF)-1),
                                                         className = 'dt-center'),
                                                       list( # Set renderer function for 'float' type columns (LogFC)
                                                         targets = 11, # score
                                                         render = htmlwidgets::JS( "function ( data, type, row, meta ) {return type === 'export' || (!Number.isFinite(data)) ? data : data.toPrecision(4);}")),
                                                       list( # Set renderer function for 'scientific' type columns (PValue)
                                                         targets = c(6, 7, 8), # pvalue, adj. pvalue, qvalue, score
                                                         render = htmlwidgets::JS( "function ( data, type, row, meta ) {return type === 'export' || (!Number.isFinite(data)) ? data : data.toExponential(4);}"))),
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
              # Add  bar relative to logFC
              formatStyle( columns = "score",
                           backgroundColor = styleInterval(seq( min( currentEnrichmentDF[["score"]]), max( currentEnrichmentDF[["score"]]), length.out = 99), colorRampPalette(c("blue", "darkred"))(100)))
            
            cat(" \n \n")
            
            print( htmltools::tagList(summaryTable))
            
            cat(" \n \n")
          }
        } # secondRun (rendering)

        
        # Plot enrichments for each identity/sample/cluster
        for( currentIdentity in names(enrichmentList[[currentComparisonName]][[currentCategoryName]]))
        {
          message(paste("Individual enrichments for identity:", currentIdentity))
          # Get the data frame of enrichments for current identity
          currentEnrichTop = enrichmentList[[currentComparisonName]][[currentCategoryName]][[currentIdentity]]
          # Filter top terms (done internally by dotplot but we do it for pathview too)
          if(!is.null(TOPTERMS_FIGURE))
          {
            currentEnrichTop = currentEnrichTop[order( currentEnrichTop[["p.adjust"]]),] # Sort it
            currentEnrichTop = head(currentEnrichTop, TOPTERMS_FIGURE) # Select first occurences
          }
          
          # Create the dotplot figure as plotly object
          dotplotFigure = dotPlotly( currentEnrichTop, topTerms = TOPTERMS_FIGURE)
  
          if(!is.null( dotplotFigure)) # Function returns NULL on empty enrichments
          {
            if(!FLAG_secondRun) # check if this is a second run for rendering only (see analysisParams)
            {
              cat("\n\nFirst (computing) run. Run report generation a second time to get all tables/figures and proper layout...")
            } else
            {
              cat("\n\n#####", currentIdentity, "(top", TOPTERMS_FIGURE, "terms)\n\n")
              
              message(paste("Dotplot"))
              cat("\n###### Dotplot")
              
              cat(" \n \n")
              print( htmltools::tagList(dotplotFigure))
              cat(" \n \n")
              
              # Add title for next section (KEGG graphs)
              if(currentCategoryName == "K-Pathways") cat("\n###### KEGG graphs:", currentComparisonName, "-", currentIdentity, "\n\n")
            }
            
            # Extract pathway graph from KEGG and show DEG genes
            if(currentCategoryName == "K-Pathways")
            {
              
              # Get DEGs logFC for current comparison/identity
              currentDEGenes_DF = topDiffExp_nestedList[[currentComparisonName]][[currentIdentity]]
              currentGenesLogFC_DE = currentDEGenes_DF[["avg_log2FC"]]
              names(currentGenesLogFC_DE) = rownames(currentDEGenes_DF)

              # Get logFC for all tested genes (not signif only)
              currentAllGenes_DF = resultsDE_nestedList_notFiltered[[currentComparisonName]][[currentIdentity]]
              currentGenesLogFC_All = currentAllGenes_DF[["avg_log2FC"]]
              names(currentGenesLogFC_All) = rownames(currentAllGenes_DF)

              # Loop on enriched pathways
              for(currentPathwayID in rownames(currentEnrichTop)[order(currentEnrichTop[["Description"]])]) # Sort by alphabetic order of description text
              {
                # Create destination folder
                pathEnrichPathwayOutput = file.path(PATH_ANALYSIS_OUTPUT, "DiffExpByHTO", currentComparisonName, paste0("pathways_", currentIdentity))
                
                # Create a name suffix four output png file created by pathview
                outSuffix = paste(currentComparisonName, currentIdentity, sep = "_")
                
                # Create external figures in first run only to avoid messing the layout
                if(!FLAG_secondRun) # check if this is a second run for rendering only (see analysisParams)
                {
                  cat("\n\nFirst (computing) run. Run report generation a second time to get all tables/figures and proper layout...")
                  
                  dir.create(pathEnrichPathwayOutput, recursive = TRUE, showWarnings = FALSE)
                  
                  # Switch to destination folder because pathview only writes in "current directory"
                  savedCurrentDir = getwd() # Save current path to restore later
                  setwd(pathEnrichPathwayOutput)
                  
                  # Create figure showing logFC for DEG genes only
                  message(paste("Pathview logFC DEG only"))
                  pathview( gene.data = currentGenesLogFC_DE,
                            gene.idtype = "SYMBOL",
                            pathway.id=gsub("^.{3}", "", currentPathwayID),  # Remove species prefix
                            species = "mmu", 
                            key.pos = "bottomright", 
                            both.dirs = TRUE, # LogFC is 'bidirectional'
                            kegg.native=TRUE, 
                            same.layer=FALSE, 
                            node.sum = "max.abs", # Node summary function when multiple genes mapped to it (using "max.abs", if value is larger in 'all genes' than 'DE', it means a gene with larger FC has not been declared significative, probably based on PValue)
                            low = 'blue', 
                            mid = 'lightgrey', 
                            high = 'red',
                            na.col = 'white', #"#DDFFDD", # Genes not in tested matrix (DE analysis), ('white' makes it similar to 'not-genes' rectangles too)
                            bins = 10,
                            cex = 0.25,
                            out.suffix = paste0(outSuffix, "_logFC_DEGenes"))

                  # Same figure showing logFC for all tested genes
                  message(paste("Pathview logFC all tested genes"))
                  pathview( gene.data = currentGenesLogFC_All,
                            gene.idtype = "SYMBOL",
                            pathway.id=gsub("^.{3}", "", currentPathwayID),  # Remove species prefix
                            species = "mmu", 
                            key.pos = "bottomright", 
                            both.dirs = TRUE, # LogFC is 'bidirectional'
                            kegg.native=TRUE, 
                            same.layer=FALSE, 
                            node.sum = "max.abs", # Node summary function when multiple genes mapped to it (using "max.abs", if value is larger in 'all genes' than 'DE', it means a gene with larger FC has not been declared significative, probably based on PValue)
                            low = 'blue', 
                            mid = 'lightgrey', 
                            high = 'red',
                            na.col = 'white', #"#DDFFDD", # Genes not in tested matrix (DE analysis), ('white' makes it similar to 'not-genes' rectangles too)
                            bins = 10,
                            cex = 0.25,
                            out.suffix = paste0(outSuffix, "_logFC_allGenes"))
                  
                  
                  # Repeat with no legend (could be obscuring information)
                  # Not included in report (just available as external file)
                  message(paste("Pathview logFC DEG only (no legend)"))
                  pathview( gene.data = currentGenesLogFC_DE,
                            gene.idtype = "SYMBOL",
                            pathway.id=gsub("^.{3}", "", currentPathwayID),  # Remove species prefix
                            species = "mmu", 
                            key.pos = "bottomright", 
                            plot.col.key = FALSE,
                            both.dirs = TRUE, # LogFC is 'bidirectional'
                            kegg.native=TRUE, 
                            same.layer=FALSE, 
                            node.sum = "max.abs", # Node summary function when multiple genes mapped to it (using "max.abs", if value is larger in 'all genes' than 'DE', it means a gene with larger FC has not been declared significative, probably based on PValue)
                            low = 'blue', 
                            mid = 'lightgrey', 
                            high = 'red',
                            na.col = 'white', #"#DDFFDD", # Genes not in tested matrix (DE analysis), ('white' makes it similar to 'not-genes' rectangles too)
                            bins = 10,
                            cex = 0.25,
                            out.suffix = paste0(outSuffix, "_logFC_DEGenes_noLegend"))
                  message(paste("Pathview logFC all tested genes (no legend)"))
                  pathview( gene.data = currentGenesLogFC_All,
                            gene.idtype = "SYMBOL",
                            pathway.id=gsub("^.{3}", "", currentPathwayID),  # Remove species prefix
                            species = "mmu", 
                            key.pos = "bottomright", 
                            plot.col.key = FALSE,
                            both.dirs = TRUE, # LogFC is 'bidirectional'
                            kegg.native=TRUE, 
                            same.layer=FALSE, 
                            node.sum = "max.abs", # Node summary function when multiple genes mapped to it (using "max.abs", if value is larger in 'all genes' than 'DE', it means a gene with larger FC has not been declared significative, probably based on PValue)
                            low = 'blue', 
                            mid = 'lightgrey', 
                            high = 'red',
                            na.col = 'white', #"#DDFFDD", # Genes not in tested matrix (DE analysis), ('white' makes it similar to 'not-genes' rectangles too)
                            bins = 10,
                            cex = 0.25,
                            out.suffix = paste0(outSuffix, "_logFC_allGenes_noLegend"))
                  
                  
                  # Restore previous "current directory"
                  setwd(savedCurrentDir)
                  
                } else # In second run, integrate figures generated in first one
                {
                  # Reconstruct filename as created by pathview...
                  outFilepath_all = file.path(pathEnrichPathwayOutput, paste0(currentEnrichTop[currentPathwayID,"ID"], ".", outSuffix, "_logFC_allGenes.png"))
                  outFilepath_DEG = file.path(pathEnrichPathwayOutput, paste0(currentEnrichTop[currentPathwayID,"ID"], ".", outSuffix, "_logFC_DEGenes.png"))
                  outFilepath_all_noLegend = file.path(pathEnrichPathwayOutput, paste0(currentEnrichTop[currentPathwayID,"ID"], ".", outSuffix, "_logFC_allGenes_noLegend.png"))
                  outFilepath_DEG_noLegend = file.path(pathEnrichPathwayOutput, paste0(currentEnrichTop[currentPathwayID,"ID"], ".", outSuffix, "_logFC_DEGenes_noLegend.png"))
                  # DEBUG: file.exists(outFilepath)
                  # Include generated file in report
                  cat( paste0( "\n- **", 
                               currentEnrichTop[currentPathwayID,"Description"]," (", currentPathwayID, "):**  \n",
                               "LogFC all => ![LogFC all Genes](", outFilepath_all,"){width=20%} ",
                               "&nbsp;&nbsp;&nbsp;",
                               "[LogFC DE Genes](", outFilepath_DEG,"){target='_blank'}",
                               "&nbsp;&nbsp;&nbsp;",
                               "[LogFC all (no legend)](", outFilepath_all_noLegend,"){target='_blank'}",
                               "&nbsp;&nbsp;&nbsp;",
                               "[LogFC DE Genes (no legend)](", outFilepath_DEG_noLegend,"){target='_blank'}",
                               "\n\n"))
                }
              }
            }
            
            cat(" \n \n")
          }
          
          cat(" \n \n")
        }
      }
      
  } else cat("\n\n<br>No enrichment found.<br>\n\n")
  
  
  
  
  #### Plotting GSEA MITOCARTA ENRICHMENTS
  
  for(currentIdentity in names(topDiffExp_nestedList[[currentComparisonName]]))
  {
    message(paste("GSEA Mitocarta for identity:", currentIdentity))
    cat("\n\n#### <span style='border-radius: 3px; border: 3px solid #F7BD58; padding:0px 2px'> GSEA Mitocarta:", currentIdentity, "</span> \n\n")
    
    # Create output folder for heatmaps as external files (png) 
    gseaOutputPath = file.path(PATH_ANALYSIS_OUTPUT, "DiffExpByHTO", currentComparisonName, paste0("gsea_mitocarta_", currentIdentity))
    
    # Get named vector of sorted logFC for current DE analysis (used for analysis AND PLOTTING)
    genesByLFC = resultsDE_nestedList_notFiltered[[currentComparisonName]][[currentIdentity]][["avg_log2FC"]]
    names(genesByLFC) = rownames(resultsDE_nestedList_notFiltered[[currentComparisonName]][[currentIdentity]])
    genesByLFC = sort(genesByLFC, decreasing = TRUE)
    
    # Prepare filename for gsea result as tsv table
    gseaResultFilepath = file.path( gseaOutputPath, paste0("gsea_mitocarta_", currentComparisonName, "_", currentIdentity, ".tsv"))
    
    gseaResult = data.frame()
    
    if(FLAG_secondRun) # Check if this is a second run for rendering only (see analysisParams)
    {
      if(file.exists( gseaResultFilepath))
      {
        gseaResult = read.table( gseaResultFilepath, 
                                 sep = "\t")

      }
    } else
    {
      # Make the GSEA enrichment analysis (cannot be done twice as results may vary between runs)
      gseaResult = as.data.frame( fgseaMultilevel( pathways = MITOCARTA_PATHWAYS, stats = genesByLFC))
      
      # Transform 'leadingEdge' column from list to csv
      gseaResult[["leadingEdge"]] = sapply( gseaResult[["leadingEdge"]], paste, collapse = ", ")
      
      dir.create(gseaOutputPath, recursive = TRUE, showWarnings = FALSE)
      write.table( gseaResult, 
                   file = gseaResultFilepath, 
                   sep = "\t",
                   quote = TRUE)
    }
    
    # Filter enrichment from total table and sort
    gseaResultFiltered = data.frame()
    if(nrow( gseaResult) > 0)
    {
      selectedResults = gseaResult[["padj"]] < GSEA_PVALUE_CUTOFF
      if(any(selectedResults))
      {
        gseaResultFiltered = gseaResult[ selectedResults, ]
        gseaResultFiltered = gseaResultFiltered[order( gseaResultFiltered[["padj"]]), ]
      }
    }

    if(nrow(gseaResultFiltered) > 0)
    {
      ## Plot table
      message(paste("GSEA datatable"))
      cat(" \n\n <br>  \n\n")
      cat("\n\n##### Summary table")
      cat(" \n \n");
      
      # Create datatable
      summaryTable = datatable( gseaResultFiltered,
                                class = "compact",
                                filter="top",
                                rownames = FALSE,
                                #colnames = ,
                                #caption = ,
                                extensions = c('Buttons', 'Select'),
                                options = list(dom = "<'row'<'col-sm-8'B><'col-sm-4'f>> <'row'<'col-sm-12'l>> <'row'<'col-sm-12'rt>> <'row'<'col-sm-12'ip>>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                                               autoWidth = FALSE,
                                               buttons = exportButtonsListDT,
                                               columnDefs = list(
                                                 list( # Center all columns except first one
                                                   targets = 1:(ncol( gseaResultFiltered)-1),
                                                   className = 'dt-center'),
                                                 list( # Set renderer function for 'scientific' type columns (PValue)
                                                   targets = 1:2, 
                                                   render = htmlwidgets::JS( "function ( data, type, row, meta ) {return type === 'export' || (!Number.isFinite(data)) ? data : data.toExponential(4);}"))),
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
                                               stateSave = TRUE))
      
      
      cat(" \n \n")
      print( htmltools::tagList(summaryTable))
      
      cat(" \n\n <br>  \n\n")
      
      cat("\n\n##### Summary figure")
      cat(" \n \n")
      
      ## Plot enrichment graph (summary for all enrichments)

      pngFilepath = file.path( gseaOutputPath, "Enrichments_summary.png")
      
      if(FLAG_secondRun) # Check if this is a second run for rendering only (see analysisParams)
      {
        # Don't mess with RMD layout using 'png' and 'dev.off'...
        # Integrate the figure that was generated during first run
        cat( paste0( "\n![GSEA Summary](", pngFilepath,")\n" ))
      } else # If it's the first run, create the figure as external file
      {
        # First run, generate figure as external file (messes up layout: figures in wring tabset)
        # Don't include them in report, will be done in secind run.
        png( pngFilepath, 
             width = 700, 
             height = 150+ ( nrow(gseaResultFiltered) *150 ))

        message(paste("GSEA summary figure"))
        plotGseaTable( pathways = MITOCARTA_PATHWAYS[gseaResultFiltered[["pathway"]]], 
                       stats = genesByLFC, 
                       gseaResultFiltered)
        
        dev.off()
        cat("\n\nFirst (computing) run. Run report generation a second time to get all tables/figures and proper layout...")
      }
      
      cat(" \n \n")
      
      
      
      ## Plot individual enriched pathway details
      for(currentPathwayName in gseaResultFiltered[["pathway"]])
      {
        cat(" \n\n <br>  \n\n")
        cat("\n\n##### Enrichment detail:", currentPathwayName)
        cat(" \n \n")
        message(paste("GSEA detail:", currentPathwayName))
        print( plotEnrichment( pathway = MITOCARTA_PATHWAYS[[currentPathwayName]],
                               stats   = genesByLFC))
        cat(" \n \n")
      }
      
    }
  }
  
  
  
  
  #### Plotting GSEA ENRICHMENTS FROM MISGDB IDs
  
  for(currentIdentity in names(topDiffExp_nestedList[[currentComparisonName]]))
  {
    message(paste("GSEA msigdb for identity:", currentIdentity))
    cat("\n\n#### <span style='border-radius: 3px; border: 3px solid #ad1cba; padding:0px 2px'> GSEA msigdb:", currentIdentity, "</span> \n\n")
    
    # Load gmt files
    MSIGDB_PATHWAYS_MOUSE = unlist( lapply( MSIGDB_GMT_PATHS_MOUSE, gmtPathways), recursive = FALSE)
    MSIGDB_PATHWAYS_HUMAN = unlist( lapply( MSIGDB_GMT_PATHS_HUMAN, gmtPathways), recursive = FALSE)
    
    # Convert human genes IDS to corresponding mouse IDs (library babelgene)
    MSIGDB_PATHWAYS_HUMAN2MOUSE = lapply( lapply( MSIGDB_PATHWAYS_HUMAN, 
                                                  orthologs, 
                                                  species = "mouse"),
                                          "[[",
                                          "symbol")
    names(MSIGDB_PATHWAYS_HUMAN2MOUSE) = paste0( names(MSIGDB_PATHWAYS_HUMAN2MOUSE),
                                                 "_IDsConvertedToMouse")
    # Merge original mouse and 'human to mouse' converted
    MSIGDB_PATHWAYS = c(MSIGDB_PATHWAYS_MOUSE, MSIGDB_PATHWAYS_HUMAN2MOUSE)
    
    # Create output folder for heatmaps as external files (png) 
    gseaOutputPath = file.path(PATH_ANALYSIS_OUTPUT, "DiffExpByHTO", currentComparisonName, paste0("gsea_msigdb_", currentIdentity))
    
    # Get named vector of sorted logFC for current DE analysis (used for analysis AND PLOTTING)
    genesByLFC = resultsDE_nestedList_notFiltered[[currentComparisonName]][[currentIdentity]][["avg_log2FC"]]
    names(genesByLFC) = rownames(resultsDE_nestedList_notFiltered[[currentComparisonName]][[currentIdentity]])
    genesByLFC = sort(genesByLFC, decreasing = TRUE)
    
    # Convert to ENTRE IDs ???
    #mapIds(org.Mm.eg.db, keys=rownames(resultsDE_nestedList_notFiltered[[1]][[1]]), column="ENTREZID", keytype="SYMBOL", multiVals="first")
    
    # Prepare filename for gsea result as tsv table
    gseaResultFilepath = file.path( gseaOutputPath, paste0("gsea_msigdb_", currentComparisonName, "_", currentIdentity, ".tsv"))
    
    gseaResult = data.frame()
    
    if(FLAG_secondRun) # Check if this is a second run for rendering only (see analysisParams)
    {
      if(file.exists( gseaResultFilepath))
      {
        gseaResult = read.table( gseaResultFilepath, 
                                 sep = "\t")
        
      }
    } else
    {
      
      # Make the GSEA enrichment analysis (cannot be done twice as results may vary between runs)
      gseaResult = as.data.frame( fgseaMultilevel( pathways = MSIGDB_PATHWAYS, stats = genesByLFC))
      
      # Transform 'leadingEdge' column from list to csv
      gseaResult[["leadingEdge"]] = sapply( gseaResult[["leadingEdge"]], paste, collapse = ", ")
      
      dir.create(gseaOutputPath, recursive = TRUE, showWarnings = FALSE)
      write.table( gseaResult, 
                   file = gseaResultFilepath, 
                   sep = "\t",
                   quote = TRUE)
    }
    
    # Filter enrichment from total table and sort
    gseaResultFiltered = data.frame()
    if(nrow( gseaResult) > 0)
    {
      selectedResults = gseaResult[["padj"]] < GSEA_PVALUE_CUTOFF
      if(any(selectedResults))
      {
        gseaResultFiltered = gseaResult[ selectedResults, ]
        gseaResultFiltered = gseaResultFiltered[order( gseaResultFiltered[["padj"]]), ]
      }
    }
    
    if(nrow(gseaResultFiltered) > 0)
    {
      ## Plot table
      message(paste("GSEA datatable"))
      cat(" \n\n <br>  \n\n")
      cat("\n\n##### Summary table")
      cat(" \n \n");
      
      # Create datatable
      summaryTable = datatable( gseaResultFiltered,
                                class = "compact",
                                filter="top",
                                rownames = FALSE,
                                #colnames = ,
                                #caption = ,
                                extensions = c('Buttons', 'Select'),
                                options = list(dom = "<'row'<'col-sm-8'B><'col-sm-4'f>> <'row'<'col-sm-12'l>> <'row'<'col-sm-12'rt>> <'row'<'col-sm-12'ip>>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                                               autoWidth = FALSE,
                                               buttons = exportButtonsListDT,
                                               columnDefs = list(
                                                 list( # Center all columns except first one
                                                   targets = 1:(ncol( gseaResultFiltered)-1),
                                                   className = 'dt-center'),
                                                 list( # Set renderer function for 'scientific' type columns (PValue)
                                                   targets = 1:2, 
                                                   render = htmlwidgets::JS( "function ( data, type, row, meta ) {return type === 'export' || (!Number.isFinite(data)) ? data : data.toExponential(4);}"))),
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
                                               stateSave = TRUE))
      
      
      cat(" \n \n")
      print( htmltools::tagList(summaryTable))
      
      cat(" \n\n <br>  \n\n")
      
      cat("\n\n##### Summary figure")
      cat(" \n \n")
      
      ## Plot enrichment graph (summary for all enrichments)
      
      pngFilepath = file.path( gseaOutputPath, "Enrichments_summary.png")
      
      if(FLAG_secondRun) # Check if this is a second run for rendering only (see analysisParams)
      {
        # Don't mess with RMD layout using 'png' and 'dev.off'...
        # Integrate the figure that was generated during first run
        cat( paste0( "\n![GSEA Summary](", pngFilepath,")\n" ))
      } else # If it's the first run, create the figure as external file
      {
        # First run, generate figure as external file (messes up layout: figures in wring tabset)
        # Don't include them in report, will be done in secind run.
        png( pngFilepath, 
             width = 700, 
             height = 150+ ( nrow(gseaResultFiltered) *150 ))
        
        message(paste("GSEA summary figure"))
        plotGseaTable( pathways = MSIGDB_PATHWAYS[gseaResultFiltered[["pathway"]]], 
                       stats = genesByLFC, 
                       gseaResultFiltered)
        
        dev.off()
        cat("\n\nFirst (computing) run. Run report generation a second time to get all tables/figures and proper layout...")
      }
      
      cat(" \n \n")
      
      
      
      ## Plot individual enriched pathway details
      for(currentPathwayName in gseaResultFiltered[["pathway"]])
      {
        cat(" \n\n <br>  \n\n")
        cat("\n\n##### Enrichment detail:", currentPathwayName)
        cat(" \n \n")
        message(paste("GSEA detail:", currentPathwayName))
        print( plotEnrichment( pathway = MSIGDB_PATHWAYS[[currentPathwayName]],
                               stats   = genesByLFC))
        cat(" \n \n")
      }
      
    }
  }
  
  
  
  
  ### PLOTTING SELECTED KEGG PATHWAYS
  
  # Extract pathway graph from KEGG and show DEG genes
  if( exists("KEGG_MANUAL_PATHWAY_IDS") && (length(KEGG_MANUAL_PATHWAY_IDS)>0) )
  {
    for(currentIdentity in names(topDiffExp_nestedList[[currentComparisonName]]))
    {
      message(paste("Manually selected KEGG pathways, projection of DE genes for identity:", currentIdentity))
      cat("\n\n#### <span style='border-radius: 3px; border: 3px solid #39e6e0; padding:0px 2px'> KEGG manual selection:", currentIdentity, "</span> \n\n")
    
      # Get DEGs logFC for current comparison/identity
      currentDEGenes_DF = topDiffExp_nestedList[[currentComparisonName]][[currentIdentity]]
      currentGenesLogFC_DE = currentDEGenes_DF[["avg_log2FC"]]
      names(currentGenesLogFC_DE) = rownames(currentDEGenes_DF)
      
      # Get logFC for all tested genes (not signif only)
      currentAllGenes_DF = resultsDE_nestedList_notFiltered[[currentComparisonName]][[currentIdentity]]
      currentGenesLogFC_All = currentAllGenes_DF[["avg_log2FC"]]
      names(currentGenesLogFC_All) = rownames(currentAllGenes_DF)
      
      # Loop on enriched pathways
      for(currentPathwayID in KEGG_MANUAL_PATHWAY_IDS)
      {
        # Create destination folder
        pathEnrichPathwayOutput = file.path(PATH_ANALYSIS_OUTPUT, "DiffExpByHTO", currentComparisonName, paste0("manual_selected_pathways_", currentIdentity))
        
        # Create a name suffix four output png file created by pathview
        outSuffix = paste(currentComparisonName, currentIdentity, sep = "_")
        
        # Create external figures in first run only to avoid messing the layout
        if(!FLAG_secondRun) # check if this is a second run for rendering only (see analysisParams)
        {
          cat("\n\nFirst (computing) run. Run report generation a second time to get all tables/figures and proper layout...")
          
          dir.create(pathEnrichPathwayOutput, recursive = TRUE, showWarnings = FALSE)
          
          # Switch to destination folder because pathview only writes in "current directory"
          savedCurrentDir = getwd() # Save current path to restore later
          setwd(pathEnrichPathwayOutput)
          
          # Create figure showing logFC for DEG genes only
          message(paste("Pathview logFC DEG only"))
          pathview( gene.data = currentGenesLogFC_DE,
                    gene.idtype = "SYMBOL",
                    pathway.id=currentPathwayID,
                    species = "mmu", 
                    key.pos = "bottomright", 
                    both.dirs = TRUE, # LogFC is 'bidirectional'
                    kegg.native=TRUE, 
                    same.layer=FALSE, 
                    node.sum = "max.abs", # Node summary function when multiple genes mapped to it (using "max.abs", if value is larger in 'all genes' than 'DE', it means a gene with larger FC has not been declared significative, probably based on PValue)
                    low = 'blue', 
                    mid = 'lightgrey', 
                    high = 'red',
                    na.col = 'white', #"#DDFFDD", # Genes not in tested matrix (DE analysis), ('white' makes it similar to 'not-genes' rectangles too)
                    bins = 10,
                    cex = 0.25,
                    out.suffix = paste0(outSuffix, "_logFC_DEGenes"))
          
          # Same figure showing logFC for all tested genes
          message(paste("Pathview logFC all tested genes"))
          pathview( gene.data = currentGenesLogFC_All,
                    gene.idtype = "SYMBOL",
                    pathway.id=currentPathwayID,
                    species = "mmu", 
                    key.pos = "bottomright", 
                    both.dirs = TRUE, # LogFC is 'bidirectional'
                    kegg.native=TRUE, 
                    same.layer=FALSE, 
                    node.sum = "max.abs", # Node summary function when multiple genes mapped to it (using "max.abs", if value is larger in 'all genes' than 'DE', it means a gene with larger FC has not been declared significative, probably based on PValue)
                    low = 'blue', 
                    mid = 'lightgrey', 
                    high = 'red',
                    na.col = 'white', #"#DDFFDD", # Genes not in tested matrix (DE analysis), ('white' makes it similar to 'not-genes' rectangles too)
                    bins = 10,
                    cex = 0.25,
                    out.suffix = paste0(outSuffix, "_logFC_allGenes"))
          
          
          # Repeat with no legend (could be obscuring information)
          # Not included in report (just available as external file)
          message(paste("Pathview logFC DEG only (no legend)"))
          pathview( gene.data = currentGenesLogFC_DE,
                    gene.idtype = "SYMBOL",
                    pathway.id=currentPathwayID,
                    species = "mmu", 
                    key.pos = "bottomright", 
                    plot.col.key = FALSE,
                    both.dirs = TRUE, # LogFC is 'bidirectional'
                    kegg.native=TRUE, 
                    same.layer=FALSE, 
                    node.sum = "max.abs", # Node summary function when multiple genes mapped to it (using "max.abs", if value is larger in 'all genes' than 'DE', it means a gene with larger FC has not been declared significative, probably based on PValue)
                    low = 'blue', 
                    mid = 'lightgrey', 
                    high = 'red',
                    na.col = 'white', #"#DDFFDD", # Genes not in tested matrix (DE analysis), ('white' makes it similar to 'not-genes' rectangles too)
                    bins = 10,
                    cex = 0.25,
                    out.suffix = paste0(outSuffix, "_logFC_DEGenes_noLegend"))
          
          message(paste("Pathview logFC all tested genes (no legend)"))
          pathview( gene.data = currentGenesLogFC_All,
                    gene.idtype = "SYMBOL",
                    pathway.id=currentPathwayID,
                    species = "mmu", 
                    key.pos = "bottomright", 
                    plot.col.key = FALSE,
                    both.dirs = TRUE, # LogFC is 'bidirectional'
                    kegg.native=TRUE, 
                    same.layer=FALSE, 
                    node.sum = "max.abs", # Node summary function when multiple genes mapped to it (using "max.abs", if value is larger in 'all genes' than 'DE', it means a gene with larger FC has not been declared significative, probably based on PValue)
                    low = 'blue', 
                    mid = 'lightgrey', 
                    high = 'red',
                    na.col = 'white', #"#DDFFDD", # Genes not in tested matrix (DE analysis), ('white' makes it similar to 'not-genes' rectangles too)
                    bins = 10,
                    cex = 0.25,
                    out.suffix = paste0(outSuffix, "_logFC_allGenes_noLegend"))
          
          
          # Restore previous "current directory"
          setwd(savedCurrentDir)
          
        } else # In second run, integrate figures generated in first one
        {
          # Reconstruct filename as created by pathview...
          outFilepath_all = file.path(pathEnrichPathwayOutput, paste0("mmu", currentPathwayID, ".", outSuffix, "_logFC_allGenes.png"))
          outFilepath_DEG = file.path(pathEnrichPathwayOutput, paste0("mmu", currentPathwayID, ".", outSuffix, "_logFC_DEGenes.png"))
          outFilepath_all_noLegend = file.path(pathEnrichPathwayOutput, paste0("mmu", currentPathwayID, ".", outSuffix, "_logFC_allGenes_noLegend.png"))
          outFilepath_DEG_noLegend = file.path(pathEnrichPathwayOutput, paste0("mmu", currentPathwayID, ".", outSuffix, "_logFC_DEGenes_noLegend.png"))
          if(all(file.exists(c(outFilepath_all, outFilepath_DEG, outFilepath_all_noLegend, outFilepath_DEG_noLegend))))
          {
            # DEBUG: file.exists(outFilepath)
            # Include generated file in report
            cat( paste0( "\n- **", 
                         currentPathwayID," (", currentPathwayID, "):**  \n",
                         "LogFC all => ![LogFC all Genes](", outFilepath_all,"){width=20%} ",
                         "&nbsp;&nbsp;&nbsp;",
                         "[LogFC DE Genes](", outFilepath_DEG,"){target='_blank'}",
                         "&nbsp;&nbsp;&nbsp;",
                         "[LogFC all (no legend)](", outFilepath_all_noLegend,"){target='_blank'}",
                         "&nbsp;&nbsp;&nbsp;",
                         "[LogFC DE Genes (no legend)](", outFilepath_DEG_noLegend,"){target='_blank'}",
                         "\n\n"))
          } else
          {
            cat(paste0("\n- Figures for Pathway '", currentPathwayID, "' could not be found. Make sure this reference exists for current species.\n\n"))
          }
        }
      }
    }
  }
}




