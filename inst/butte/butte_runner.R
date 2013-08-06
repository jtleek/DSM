# Run the butte analysis
# This function expects as input a vector of numeric IDs,
# the type of dataset each ID refers to (currently, only GDS is allowed),
# and the number of drugs to test against (just the first n drugs from the
# connectivity map dataset).
butte_runner <- function(ids, id_types="GDS", ndrug=10){
	knit("butte_template.Rmd")
	markdownToHTML("butte_template.md",output="butte_analysis.html")	
}