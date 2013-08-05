# Get annotation.
# This function is called when it is not clear within a dataset which
# samples are cases and which are controls. The user is prompted to proivde
# an ordered vector indicating case/control in whataver format they choose.

get_annotation <- function(pd){
	cat("Here is the sample information for your dataset:\n")
	print(pd)
	cat("\n")
	
	ans <- readline(prompt=cat("Please provide a vector indiciating which samples are cases \nand which are controls,e.g. c('case','control','case')\n"))
	x <- eval(parse(text=ans))

	x
}


# Get gene list.
# This function fetches GDS datasets from GEO for the IDs specified.
# It then runs a straightforward linear model via limma that calculates
# the mean effect of case vs. control and uses the contrast to
# generate a list of significant genes. We report the genes that have
# non-empty gene IDs and p-values < 0.05 at an FDR of < 0.05. We also
# have a prompt system in place to ask the user to properly annotate the
# microrarray.

get_gene_list <- function(id){

	test_id <- paste("GDS", id, sep="")
	test_gds <- getGEO(test_id, AnnotGPL=T)
	test_eset <- GDS2eSet(test_gds)
	gene_ids <- toupper(as.vector(fData(test_eset)[,3]))
	# Rank normalize the data
	test_eset <- normalize.ExpressionSet.quantiles(test_eset)

	# Limma analysis. We're keeping it simple:
	# case and control groups with a direct contrast.

	pd <- pData(test_eset)

	if(length(levels(pd[,2])) != 2 & length(levels(pd[,3])) != 2){
		f <- get_annotation(pd)
	} else if(length(levels(pd[,2])) == 2 & length(levels(pd[,3])) != 2){
		ans <- ""
		while(!(ans %in% c("y","n"))){
			ans <- readline(prompt=paste("Using", colnames(pd)[2], "for annotation. OK? (y/n)",sep=" "))
			ans <- tolower(substr(ans,1,1))
		}

		if(ans == "y"){
			f <- pd[,2]
		} else {
			f <- get_annotation(pd)
		}
	} else if(length(levels(pd[,2])) != 2 & length(levels(pd[,3])) == 2){
		ans <- ""
		while(!(ans %in% c("y","n"))){
			ans <- readline(prompt=paste("Using", colnames(pd)[3], "for annotation. OK? (y/n)",sep=" "))
			ans <- tolower(substr(ans,1,1))
		}

		if(ans == "y"){
			f <- pd[,3]
		} else {
			f <- get_annotation(pd)
		}

	} else {
		ans <- ""
		while(!(ans %in% c(2,3))){
			ans <- readline(prompt=paste("Should I use 2:", colnames(pd)[2], "or 3:", colnames(pd)[3], " for annotation? (2/3) ",sep=" "))
		}	
		ans <- as.numeric(ans)
		f <- pd[,ans]
	}

	f <- factor(gsub(" ","", f))
	f <- factor(gsub("[[:punct:]]","",f))

	design <- model.matrix(~0+f)
	colnames(design) <- levels(f)

	cat("I'm here 1\n")

	fit <- lmFit(test_eset, design)

	cat("I'm here 2\n")
	ctrst <- as.character(paste(levels(f), collapse="-"))
	cat("ctrst is", ctrst, "with class", class(ctrst), "\n")
	cont.matrix <- makeContrasts(cc=ctrst,levels=design)

	cat("I'm here 3\n")
	fit2 <- contrasts.fit(fit,cont.matrix)
	fit2 <- eBayes(fit2)

	cat("I'm here 4\n")
	tt <- topTable(fit2,number=100000)
	tt <- tt[which(tt$adj.P.Val < 0.05),] # This is FDR!
	output <- tt[,c("Gene.symbol", "logFC")]
	names(output) <- c("gene", "fold")
	reg <- ifelse(output$fold < 0, "down", "up")
	output <- cbind(output, reg)
	# throw out entries with no gene names
	output <- output[-which(output$gene == ""),]
	output$gene <- toupper(output$gene)
	output
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
	if(length(idx) > 0){
		v_j <- v_j[-idx]
		dis_gene <- dis_gene[-idx]
	}

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

calc_pval <- function(dat,drugname,dds){
	scores <- vector("numeric",100)
	vec <- 1:length(drugname)
	for(i in 1:100){
		idx <- sample(vec, dim(dat)[1], replace=F)
		drug <- vec[idx]
		names(drug) <- drugname[idx]
		scores[i] <- calc_dds(drug,dat)
	}
	
	(length(which(scores > abs(dds))) + length(which(scores < -abs(dds))))/100			
}