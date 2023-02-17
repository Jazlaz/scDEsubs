test_scDEsubs <- function() 
{
   message('Testing scDEsubs...', appendLF=FALSE)

	load(system.file('extdata', 'data.RData', package='scDEsubs'))

	scDEsubs.run <- scDEsubs(   org='hsa', 
	                        mRNAexpr=SingleCellmRNAexpr, 
	                        mRNAnomenclature='entrezgene', 
	                        pathways='All', 
	                        DEtool='edgeRLRT', 
	                        DEpar=0.05,
	                        CORtool='pearson', 
	                        CORpar=0.6, 
	                        subpathwayType=NULL,
	                        rankedList=NULL,
	                        grouping = grouping)

   all.equal(scDEsubs.run, scDEsubs.out)

   message('done.')
}
