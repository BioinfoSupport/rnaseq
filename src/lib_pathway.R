
library(tidyverse)
library(patchwork)





#' require cond, ensembl_id , dir (a partition of the genes in cond)
#' A: a 2 column annotation table with ensembl_id,pathway_id
#' P: pathways metadata with: pathway_id,pathway_name
enrichment_analysis <- function(res,A,P) {
	
	# Compute universe sizes for each condition
	universe_sizes <- select(res,cond,ensembl_id) |>
		inner_join(A,by="ensembl_id",relationship="many-to-many") |>
		group_by(cond) |>
		summarize(universe_size=n_distinct(ensembl_id))
	
	# Compute geneset sizes for each condition/direction
	genesets_sizes <- select(res,cond,dir,ensembl_id) |>
		inner_join(A,by="ensembl_id",relationship="many-to-many") |>
		group_by(cond,dir) |>
		summarize(geneset_size = n_distinct(ensembl_id))
	
	# Compute pathway size for each condition/pathway
	pathway_sizes <- select(res,cond,ensembl_id) |>
		inner_join(A,by="ensembl_id",relationship="many-to-many") |>
		group_by(cond,pathway_id) |>
		summarize(pathway_size=n_distinct(ensembl_id))
	
	overlaps <- select(res,cond,dir,ensembl_id) |>
		inner_join(A,by="ensembl_id",relationship="many-to-many") |>
		group_by(cond,dir,pathway_id) |>
		summarize(overlap_size = n_distinct(ensembl_id))
	
	inner_join(pathway_sizes,genesets_sizes,by="cond",relationship = "many-to-many") |>
		inner_join(universe_sizes,by="cond",relationship = "many-to-one") |>
		left_join(overlaps,by=c("cond","dir","pathway_id"),relationship = "one-to-one") |>
		replace_na(list(overlap_size=0)) |>
		mutate(pval = phyper(pmax(0,overlap_size-1L),pathway_size,universe_size-pathway_size,geneset_size,lower.tail=FALSE)) |>
		group_by(cond) |>
		mutate(padj = p.adjust(pval,method="fdr")) |>
		ungroup() |>
		left_join(P,by="pathway_id") |>
		relocate(pathway_id,pathway_name,cond,dir) |>
		arrange(pval)
}



# Plot enrichment analysis results
ggenrichdot <- function(enrich,overlap_label=FALSE,max_padj=0.05) {
	p <- enrich |>
		ungroup() |>
		filter(pathway_id %in% enrich$pathway_id[enrich$padj<=max_padj]) |>
		mutate(x = sprintf("%s (%d)",cond,geneset_size)) |>
		mutate(y = sprintf("%s [%s](%d)",pathway_name,pathway_id,pathway_size)) |>
		mutate(y = reorder(y,pval,min)) |>
		ggplot(aes(y = y,x = x)) +
		facet_wrap(~dir,scales = "free_x",drop = TRUE,nrow=1) +
		geom_point(aes(size=2*pi*(sqrt(overlap_size/pathway_size/2/pi) + 0.1)^2),color="black",data=~filter(.,padj<=max_padj)) +
	  geom_point(aes(size=2*pi*(sqrt(overlap_size/pathway_size/2/pi) + 0.1)^2),color="grey",data=~filter(.,padj>max_padj)) +
		geom_point(aes(size=overlap_size/pathway_size,color=pmin(-log10(padj),10))) +
		scale_size_area(labels=scales::percent) + 
		scale_color_gradient(low="#EFD500",high="#0094CD",limits=c(0,6)) +
		labs(size="% set",color=expression(-log[10](padj)))
	if (overlap_label) p <- p + geom_text(aes(label=overlap_size),color="grey",size=3)
	p + xlab("") + ylab("") +
		theme_bw() +
		theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
}





reactome <- function(org = c("Homo sapiens","Mus musculus"),dir="data/pathway/reactome") {
  org <- match.arg(org)
  P = read.table(
    file.path(dir,"Ensembl2Reactome_All_Levels.txt.gz"),
    sep="\t",comment="",quote="",
    colClasses = c("NULL","character","NULL","character","NULL","character"),
    col.names = c(NA,"pathway_id",NA,"pathway_name",NA,"species")
  ) |>
    filter(species %in% org) |>
    distinct(pathway_id,pathway_name)
  
  A = read.table(
    file.path(dir,"Ensembl2Reactome_All_Levels.txt.gz"),
    sep="\t",comment="",quote="",
    colClasses = c("character","character","NULL","NULL","NULL","NULL"),
    col.names = c("ensembl_id","pathway_id",NA,NA,NA,NA)
  ) |>
    filter(pathway_id %in% P$pathway_id) |>
    arrange(pathway_id)
  
  H = read.table(
    file.path(dir,"ReactomePathwaysRelation.txt"),
    sep="\t",comment="",quote="",
    colClasses = c("character","character"),
    col.names = c("parent_id","child_id")
  ) |>
    filter(parent_id %in% P$pathway_id,child_id %in% P$pathway_id)
  
  list(P=P,A=A,H=H)
}





read_gmt <- function(con) {
  txt <- readLines(con)
  pat <- "^([^\t]*)\t([^\t]*)\t(.*)$"
  z <- tibble(
    pathway_id = str_extract(txt,pat,1),
    pathway_url = str_extract(txt,pat,2),
    members = str_split(str_extract(txt,pat,3),"\t")
  )
}

gsea <- function(org = c("Homo sapiens"),dir="data/pathway/msigdb_v2023.2.Hs_GMTs") {
  org <- match.arg(org)
  gmt <- read_gmt(file.path(dir,"msigdb.v2023.2.Hs.symbols.gmt"))
  
  P <- select(gmt,pathway_id = pathway_id,pathway_name = pathway_id,pathway_url=pathway_url)
  A <- gmt |> 
    select(pathway_id,members) |>
    unnest(members) |>
    dplyr::rename(symbol=members)
  list(P=P,A=A)
}
