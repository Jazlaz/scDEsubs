scDEsubs <- function( org, mRNAexpr, mRNAnomenclature, pathways, 
                    DEtool, DEpar, CORtool, CORpar, subpathwayType,
                    rankedList=NULL, verbose=TRUE, grouping = NA)
{
    # Create the adjacency matrix
    dataNet <- .constructNetwork(org=org, 
                            mRNAexpr=mRNAexpr, 
                            mRNAnomenclature=mRNAnomenclature, 
                            pathways=pathways)
    mRNAexpr <- dataNet[['mRNAexpr']]
    edgeList <- dataNet[['edgeList']]

    # Filter the nodes and the edges of the organism's pathways network
    lens     <- suppressWarnings(split(1:ncol(mRNAexpr), 1:2))
    lens     <- unname(sapply(lens, function(x) { length(x) }))
    net      <- .pruneNetwork(  edgeList=edgeList,
                                mRNAexpr=mRNAexpr, 
                                DEGchoice=DEtool, 
                                DEGthresh=DEpar, 
                                classes=c(rep(1, lens[1]), rep(2, lens[2]) ),
                                corr_threshold=CORpar, org=org,
                                CORtool=CORtool,
                                rankedList=rankedList,
                                verbose=verbose, grouping=grouping,
                                mRNAnomenclature=mRNAnomenclature)

    edgeList <- net[['edgeList']]
    DEgenes  <- net[['DEgenes']]
    output   <- list('org'=org, 'mRNAnomenclature'=mRNAnomenclature,
                    'edgeList'=edgeList, 'DEgenes'=DEgenes)
    # Choose specific subpathways
    subTypes <- subpathwayTypes(grouping=subpathwayType)

    # Subpathway analysis
    for (subType in subTypes)
    {
        subs <- .subpathwayAnalysis( edgeList=edgeList, 
                                    method=subType,
                                    DEgenes=DEgenes,
                                    org=org,
                                    verbose=verbose)
        output <- c(output, subs)
    }

    return( output )
}


.DEanalysis <- function( count.matrix, DEGchoice, classes, grouping=NA, DEGthresh, ...)
{
    # 
    # Differential expression analysis using various DE analysis tools from
    # 
    # Soneson,C. and Delorenzi,M. (2013) A comparison of methods for 
    # differential expression analysis of RNA-seq data. BMC bioinformatics, 
    # 14(1), 1.
    #
    argList<-list(...)
    if ( missing(count.matrix) ) { message('Please supply a matrix.') }
    if ( missing(classes) )      { message('Please supply the classes.') }
    if ( missing(DEGchoice) )    { message('Please supply a type.') }

    supportedMethods <- c('All', 'Wilcoxon', 'edgeRLRT', 'limmatrend','MASTcpm')
    
    
    if ( DEGchoice == 'limmatrend')
    {
      message("limmatrend")
      
      dge <- DGEList(count.matrix, group = grouping)
      dge <- calcNormFactors(dge)
      design <- model.matrix(~grouping)
      y <- new("EList")
      y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
      fit <- lmFit(y, design = design)
      fit <- eBayes(fit, trend = TRUE, robust = TRUE)
      tt <- topTable(fit, n = Inf, adjust.method = "BH")
      
      res.pvalue2 <- tt$P.Value
      names(res.pvalue2) <- rownames(tt)
      saveRDS(res.pvalue2, file = "pvals/res_limmatrend.rds")
      return(res.pvalue2)
    }  
    if ( DEGchoice == 'MASTcpm'){
      message("MAST, CPM (including detection rate)")
      stopifnot(all(names(grouping) == colnames(count.matrix)))
      grp <- grouping
      cdr <- scale(colMeans(count.matrix > 0))
      dge <- DGEList(counts = count.matrix)
      dge <- edgeR::calcNormFactors(dge)
      cpms <- cpm(dge)
      sca <- FromMatrix(exprsArray = log2(cpms + 1), 
                        cData = data.frame(wellKey = names(grp), 
                                           grp = grp, cdr = cdr))
      zlmdata <- zlm(~cdr + grp, sca)
      mast <- lrTest(zlmdata, "grp")
      df1 = data.frame(pval = mast[, "hurdle", "Pr(>Chisq)"],
                      row.names = names(mast[, "hurdle", "Pr(>Chisq)"]))
      df_return<- df1[,1]
      names(df_return)<-rownames(df1)
      res_MAST<-df_return
      return(df_return)
    }
    
    if ( DEGchoice == 'Wilcoxon')
    {
      message("Wilcoxon")
      
      tmm <- edgeR::calcNormFactors(count.matrix)
      tpmtmm <- edgeR::cpm(count.matrix, lib.size = tmm * colSums(count.matrix))
      idx <- 1:nrow(tpmtmm)
      names(idx) <- rownames(tpmtmm)
      wilcox_p <- sapply(idx, function(i) {
        wilcox.test(tpmtmm[i, ] ~ grouping)$p.value
      })
      res.pvalue2<-wilcox_p
      return(res.pvalue2)
    }
    
    if ( DEGchoice == 'edgeRLRT')
    {
      message("edgeRLRT")
      dge <- DGEList(count.matrix, group = grouping)
      dge <- calcNormFactors(dge)
      design <- model.matrix(~grouping)
      dge <- estimateDisp(dge, design = design)
      fit <- edgeR::glmFit(dge, design = design)
      lrt <- glmLRT(fit)
      tt <- topTags(lrt, n = Inf)
      res.pvalue2<-tt[["table"]][["PValue"]]
      names(res.pvalue2)<-rownames(tt[,4])
      return(res.pvalue2)
    }
    
    if ( DEGchoice == 'All')
    {
      #limmatrend
      {
      message("limmatrend")
      
        dge <- DGEList(count.matrix, group = grouping)
        dge <- calcNormFactors(dge)
        design <- model.matrix(~grouping)
        y <- new("EList")
        y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
        fit <- lmFit(y, design = design)
        fit <- eBayes(fit, trend = TRUE, robust = TRUE)
        tt <- topTable(fit, n = Inf, adjust.method = "BH")
        
      res.pvalue2 <- tt$P.Value
      names(res.pvalue2) <- rownames(tt)
      res.limmatrend<- res.pvalue2
      }
      #Wilcoxon
      {
        message("Wilcoxon")

          tmm <- edgeR::calcNormFactors(count.matrix)
          tpmtmm <- edgeR::cpm(count.matrix,
                               lib.size = tmm * colSums(count.matrix))
          idx <- 1:nrow(tpmtmm)
          names(idx) <- rownames(tpmtmm)
          wilcox_p <- sapply(idx, function(i) {
            wilcox.test(tpmtmm[i, ] ~ grouping)$p.value
          })
      res.wilcox<-wilcox_p
      }
      #edgeRLRT
      {
      message("edgeRLRT")
      dge <- DGEList(count.matrix, group = grouping)
      dge <- calcNormFactors(dge)
      design <- model.matrix(~grouping)
      dge <- estimateDisp(dge, design = design)
      fit <- edgeR::glmFit(dge, design = design)
      lrt <- edgeR::glmLRT(fit)
      tt <- edgeR::topTags(lrt, n = Inf)
      res.edgeRLRT<-tt$table$PValue
      names(res.edgeRLRT)<-row.names(tt$table)
      }
      #MAST
      {
        message("MAST, CPM (including detection rate)")
        stopifnot(all(names(grouping) == colnames(count.matrix)))
        grp <- grouping
        cdr <- scale(colMeans(count.matrix > 0))
        dge <- DGEList(counts = count.matrix)
        dge <- edgeR::calcNormFactors(dge)
        cpms <- cpm(dge)
        sca <- FromMatrix(exprsArray = log2(cpms + 1), 
                          cData = data.frame(wellKey = names(grp), 
                                             grp = grp, cdr = cdr))
        zlmdata <- zlm(~cdr + grp, sca)
        mast <- lrTest(zlmdata, "grp")
        df1 = data.frame(pval = mast[, "hurdle", "Pr(>Chisq)"],
                         row.names = names(mast[, "hurdle", "Pr(>Chisq)"]))
        df_return<- df1[,1]
        names(df_return)<-rownames(df1)
        res_MAST<-df_return
      }
      
      ##################################
      
      A<-res.limmatrend[which(res.limmatrend < DEGthresh)]
      B<-res_MAST[which(res_MAST < DEGthresh)]
      C<-res.wilcox[which(res.wilcox < DEGthresh)]
      D<-res.edgeRLRT[which(res.edgeRLRT < DEGthresh)]
      X<-c(A,B,C,D)
      genes<- c(X[intersect(intersect(names(A),names(B)),names(C))], 
                X[intersect(intersect(names(B),names(C)),names(D))],
                X[intersect(intersect(names(A),names(C)),names(D))],
                X[intersect(intersect(names(A),names(B)),names(D))])

      return(genes)
    }
    
    
    if  ( !is.null(DEGchoice) )
    {
      unsupportedOptions <- DEGchoice[!DEGchoice %in% supportedMethods]
      if ( length(unsupportedOptions) > 0 )
      {
        message('Option ', unsupportedOptions, ' not supported.')
        message('Supported options are ', 
                paste0(supportedMethods, collapse=', '), '.')
        return( NULL ) 
      }
    }


    return( NULL )
}

