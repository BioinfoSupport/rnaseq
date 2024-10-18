
library(tidyverse)
library(patchwork)




#' Perform multiple genes enrichment analyses
#' 
#' @param enrich a data.frame with 3 columns cond, ensembl_id, dir. 
#'   For each tested condition (cond), all the genes of the universe must appear 
#'   and are partitioned according to dir.
#' @param A a 2 column data.frame with genes annotations (must contain columns `ensembl_id` and `pathway_id`)
#' @param P a data.frame of pathways descriptions with: pathway_id,pathway_name
enrichment_analysis <- function(res,A,P,min_pathway_size=5,min_overlap_size=2) {
	
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
	
	# Compute overlaps between pathways and each condition
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
	  mutate(padj = case_when( # This correction is not ideal: better would be to run p.adjust only with tested pathways
	    pathway_size < min_pathway_size ~ 1,
	    overlap_size < min_overlap_size ~ 1,
	    TRUE ~ padj
	  )) |>
		ungroup() |>
		left_join(P,by="pathway_id") |>
		relocate(pathway_id,pathway_name,cond,dir) |>
		arrange(pval)
}



# Plot enrichment analysis results
ggenrichdot <- function(enrich,overlap_label=FALSE,max_padj=0.05,y_glue="{pathway_name} ({pathway_size})",x_glue="{cond}/{dir} ({geneset_size})") {
	p <- enrich |>
		group_by(pathway_id) |>
		filter(any(padj<=max_padj)) |>
	  ungroup() |>
	  mutate(x = as.factor(str_glue(x_glue))) |>
	  mutate(y = as.factor(str_glue(y_glue))) |>
		mutate(y = fct_reorder(y,pval,min)) |>
		ggplot(aes(y = y,x = x)) +
  		geom_point(aes(size=2*pi*(sqrt(overlap_size/pathway_size/2/pi) + 0.1)^2),color="black",data=~filter(.,padj<=max_padj)) +
  	  geom_point(aes(size=2*pi*(sqrt(overlap_size/pathway_size/2/pi) + 0.1)^2),color="grey",data=~filter(.,padj>max_padj)) +
  		geom_point(aes(size=overlap_size/pathway_size,color=-log10(padj))) +
  		scale_size_area(labels=scales::percent,limits=c(0,1),oob=scales::squish) + 
  		scale_color_gradient(low="#EFD500",high="#0094CD",limits=c(0,6),oob=scales::squish) +
  		labs(size="% set",color=expression(-log[10](padj)))
	if (overlap_label) {
	  p <- p + geom_text(aes(label=overlap_size),color="grey",size=3,hjust=0,position=position_nudge(x=0.2),data=~filter(.,overlap_size>0))
	}
	p + xlab("") + ylab("") +
		theme_bw() +
		theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1))
}




# Load reactome pathways dataset
reactome <- function(org = c("Homo sapiens","Mus musculus"),dir="data/pathway/reactome") {
  org <- match.arg(org)
  P = read.table(
    file.path(dir,"Ensembl2Reactome_All_Levels.txt.gz"),
    sep="\t",comment="",quote="",
    colClasses = c("NULL","character","NULL","character","NULL","character"),
    col.names = c(NA,"pathway_id",NA,"pathway_name",NA,"species")
  ) |>
    filter(species %in% org) |>
    distinct(pathway_id,pathway_name) |>
    mutate(pathway_name = make.unique(pathway_name))
  
  A = read.table(
    file.path(dir,"Ensembl2Reactome_All_Levels.txt.gz"),
    sep="\t",comment="",quote="",
    colClasses = c("character","character","NULL","NULL","NULL","NULL"),
    col.names = c("ensembl_id","pathway_id",NA,NA,NA,NA)
  ) |>
    filter(pathway_id %in% P$pathway_id) |>
    mutate(ensembl_id = str_replace(ensembl_id,"\\..*","")) |>
    distinct() |>
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




msigdb_generic <- function(db_file,namespace_id) {
	misgdb <- DBI::dbConnect(RSQLite::SQLite(),db_file)
	#DBI::dbListTables(misgdb)
	#DBI::dbListTables(misgdb) %>% enframe(name=NULL,value="tbl") %>% mutate(content=map(tbl,~DBI::dbGetQuery(misgdb,str_glue("SELECT * FROM {.} LIMIT 5")))) |> pull(content,name=tbl)
	P <- DBI::dbGetQuery(misgdb,"SELECT id as pathway_id,standard_name as pathway_name,collection_name FROM gene_set")
	A <- DBI::dbGetQuery(misgdb,str_glue("
    SELECT DISTINCT gene_set_id as pathway_id, source_id as ensembl_id
    FROM gene_set_gene_symbol, source_member
    WHERE (gene_set_gene_symbol.gene_symbol_id = source_member.gene_symbol_id) AND (namespace_id={namespace_id})
  "))
	list(P=P,A=A)	
}

msigdb_hs <- function() {
	msigdb_generic("data/pathway/msigdb_v2024.1.Hs.db",4)
}

msigdb_mm <- function() {
	msigdb_generic("data/pathway/msigdb_v2024.1.Mm.db",5)
}




# Add pathway lev1/lev2 groups
pathway_level <- function(P,H) {
  library(igraph)
  library(tidygraph)
  G <- tbl_graph(
    nodes = bind_rows(data.frame(pathway_id="root",pathway_name="root"),P),
    edges = bind_rows(data.frame(parent_id="root",child_id=setdiff(H$parent_id,H$child_id)),H),
    node_key = "pathway_id"
  )
  G <- mutate(G,depth = node_distance_from(1,mode="out")) 
  
  V(G)$lev1 <- local({
    D <- t(igraph::distances(G,V(G)[depth==1],mode="out"))
    lev <- max.col(-D)
    lev[!is.finite(D[cbind(seq_along(lev),lev)])] <- NA
    V(G)[depth==1]$pathway_id[lev]
  })
  V(G)$lev2 <- local({
    D <- t(igraph::distances(G,V(G)[depth==2],mode="out"))
    lev <- max.col(-D)
    lev[!is.finite(D[cbind(seq_along(lev),lev)])] <- NA
    V(G)[depth==2]$pathway_id[lev]
  })
  
  P <- left_join(P,select(as_tibble(G),pathway_id,lev1,lev2),by="pathway_id") |>
    left_join(select(P,lev1=pathway_id,lev1_name=pathway_name),by="lev1") |>
    left_join(select(P,lev2=pathway_id,lev2_name=pathway_name),by="lev2")
  P
}



