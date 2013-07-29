# Get gene list.
# This function fetches GDS datasets from GEO for the IDs specified.
# It then runs Significance Analysis of Microarrays (SAM) on each
# dataset and returns a list of up- or down-regulated genes.

# NOTE: To select the delta threshold for choosing cutoffs, we choose
# the delta associated with the lowest estimated FDR. We also throw out
# any genes without gene IDs.

get_gene_list <- function(id){

	test_id <- paste("GDS", id, sep="")
	test_gds <- getGEO(test_id, AnnotGPL=T)
	test_eset <- GDS2eSet(test_gds)
	# The gene names are in fData(test_eset)[,3]
	#gene_ids <- fData(test_eset)[,c(1,3)]
	gene_ids <- toupper(as.vector(fData(test_eset)[,3]))
	# Rank normalize the data
	test_eset <- normalize.ExpressionSet.quantiles(test_eset)
	# Optional: get actual expression data
	#eset_out <- exprs(test_eset)

	# Run SAM
	sam_out <- sam(test_eset, pData(test_eset)[,2], method="wilc.stat")

	# Here, need to make sam_out@chip <- ""
	sam_out@chip <- ""
	# I'm taking the gene list provided by the Delta with the lowest
	# FDR that still produces significant genes...we should look into
	# this, since the FDR is not always 0.05
	fdr <- as.numeric(show(sam_out)$FDR)
	sum <- summary(sam_out, show(sam_out)$Delta[which.min(fdr[fdr>0])])

	cutlow <- sum@mat.fdr[6]
	#cutup <- sum@mat.fdr[7]
	d.val <- sum@mat.sig$d.value
	reg <- ifelse(d.val <= cutlow, "down", "up")
	
	# This gets the rows of the up- and down-regulated genes
	rows <- sum@mat.sig$Row
	#output <- cbind(gene_ids[rows,],"fold"=sam_out@fold[rows])
	output <- data.frame("gene"=gene_ids[rows],"fold"=sum@mat.sig$R.fold,"reg"=reg)
	# I'm throwing out any genes without gene names
	idx <- which(output$gene == "")
	output[-idx,]
	
}

# Calculate the enrichment score.
# This function takes a drug signature and gene signature and
# applies the process described by Butte et. al. (from Lamb et. al.)
# to calculate an enrichment score for the pair.
# This function is called by calc_dds, which uses the enrichment score
# to calculate an overall drug disease score.
calc_es <- function(drug,subgene){
	ord <- order(subgene$fold)
	ord_drug <- drug[order(drug)]
	dis_gene <- subgene$gene[ord]
	v_j <- match(dis_gene, names(ord_drug))
	# Tossing NAs for now
	idx <- which(is.na(v_j))
	v_j <- v_j[-idx]
	dis_gene <- dis_gene[-idx]

	a <- max((1:length(dis_gene))/length(dis_gene) - v_j/length(ord_drug))		
	b <- max(v_j/length(ord_drug) - (0:(length(dis_gene)-1))/length(dis_gene))		

	ifelse(a > b, a, -b)
}

# Calculate the drug disease score.
# This function takes a drug signature and disease signature, then
# splits the disease signature into up- and down-regulated gene sets.
# The enrichment score is calculated for each set by calc_es(), and
# the es_up and es_down scores are used to calculate the dds.
calc_dds <- function(drug,dat){

	uidx <- which(dat$reg == "up")
	ur <- dat[uidx,]
	dr <- dat[-uidx,]

	es_up <- calc_es(drug,ur)
	es_down <- calc_es(drug,dr)

	# Return the dds
	ifelse((sign(es_up) + sign(es_down)) == 0, (es_up - es_down),0)
}

# Calculate p-values.
# This function produces a p-value for the observed drug disease score
# by generating a "null" distribution of drug disease scores for each
# disease signature using random ranks.

# I guess we're assuming the null dist is symmetric...
# This definitely needs to be changed somehow...need to use Levene's Test...

calc_pval <- function(dat,dds){
	scores <- vector("numeric",100)
	for(i in 1:100){
		drug <- sample(1:length(dat),length(dat),replace=F)
		scores[i] <- calc_dds(drug,dat)
	}
	
	(length(which(scores > abs(dds))) + length(which(scores < -abs(dds))))/100			
}