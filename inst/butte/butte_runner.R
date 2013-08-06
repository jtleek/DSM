# Run the butte analysis
# This function expects as input a vector of numeric IDs,
# the type of dataset each ID refers to (currently, only GDS is allowed),
# and the number of drugs to test against (just the first n drugs from the
# connectivity map dataset).
butte_runner <- function(ids, id_types="GDS", ndrug=10){
	knit("butte_template.Rmd")
	time <- format(Sys.time(), "%Y_%m_%d_%H%M%S")
	fn <- paste("butte_analysis_",time,".html",sep="")
	cat(paste("Converting to html file: ",fn,"\n",sep=""))
	markdownToHTML("butte_template.md",output=fn)
	system("rm butte_template.md")	
}
