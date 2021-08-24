# Description:
# run goseq to identify enriched gene ontology biological processes

library(goseq)

# Repeated for ARPE-19 and MIO-M1 analyses
path = '/analysis/arpe19_miom1/results'
dir.create(file.path(path, 'goseq'), showWarnings = F)

# load gene length info
gene.length = read.csv(file.path(path, 'gene_length.csv'), row.names = 1, 
                       stringsAsFactors = F)

# |log2FC| limit for being considered as DE
low_log2fc_cutoff = 0

# FDR cutoff for DE
FDR_cutoff = 0.05

# goseq stat. significance cutoff
goseq_stat_cutoff = 0.05

# biomart was used to extract genes in each GO BP on 2019.06.26,
# which was saved in an object named go.db
load('reference/go/go.db.RData')

GOseq_DE_genes = function(filtered.genes, gene.length.tested, de.genes, log2fc_cutoff, 
                          de.name){
    
    # function that runs goseq with inputs
    # 1. vector of detected genes (filtered.genes)
    # 2. length of detected genes (gene.length.tested)
    # 3. vector of 0/1 indicating whether the gene is to be tested,
    #    typically DE genes (de.genes)
    
    # goseq PWF, wallenius approximation, BH adjust p-values
    pwf=nullp(DEgenes = de.genes, bias.data=gene.length.tested, plot.fit = FALSE)

	# check if at least one DE gene in the db
	if (sum(names(de.genes)[de.genes] %in% names(go.db))>0){
		
		# calculate p-value, adj. p-value, separate GO ID and GO name
		go.wall = goseq(pwf, gene2cat=go.db)
		go.wall$adj.pvalue = p.adjust(go.wall$over_represented_pvalue, method="BH")
		go.wall$proper.category = sapply(go.wall$category, function(x)
			paste0('GO:',x), USE.NAMES = F)
		go.wall$go.name = sapply(go.wall$proper.category, function(x)
			paste0(toupper(substring(x, 12, 12)), 
				   substring(x, 13, nchar(x))), 
			USE.NAMES = F)
		go.wall$go.id = sapply(go.wall$proper.category, function(x)
			substring(x, 1, 10), USE.NAMES = F)
		
		# select columns of interest and rename column names
		go.wall = go.wall[,c('go.name','go.id','adj.pvalue','over_represented_pvalue',
							 'numDEInCat','numInCat')]
		colnames(go.wall) = c('go.name','go.id','adj.pvalue','pvalue',
							  'num.de','num.category')
		
		# select stat. significant GO
		go.wall = go.wall[order(go.wall$adj.pvalue),]
		go.enriched = go.wall[go.wall$adj.pvalue<goseq_stat_cutoff,]
		
		# write to drive
		write.csv(go.enriched, file.path(path, 'goseq', de.name, 
										 paste0('goseq_', log2fc_cutoff, '_',
												go.db.name,'.csv')), 
				  row.names = F)

    }
    
}

Call_goseq = function(de=de, de.name=de.name){
    
    # create comparison folder
    dir.create(file.path(path, 'goseq', de.name), showWarnings = F)
    
    # parse all filtered genes and their length
    filtered.genes = de$id
    gene.length.tested = sapply(filtered.genes, function(x) 
        gene.length$length[rownames(gene.length)==x], 
        USE.NAMES = F)
   
    # run goseq
    de.genes = setNames(de$FDR<FDR_cutoff & abs(de$logFC)>low_log2fc_cutoff, filtered.genes)
    GOseq_DE_genes(filtered.genes = filtered.genes, 
               gene.length.tested = gene.length.tested, 
               de.genes = de.genes, 
               log2fc_cutoff = paste0('log2FC=',low_log2fc_cutoff), 
               de.name = de.name)
}

# loop goseq through each comparison
file.path = list.files(file.path(path, 'edgeR'), pattern = '*.csv$', full.names = T)
file.name = list.files(file.path(path, 'edgeR'), pattern = '*.csv$', full.names = F)
file.name = tools::file_path_sans_ext(file.name)
file.name = sapply(file.name, function(x) substring(x, 10), 
                   USE.NAMES = F)

for (i in seq(file.name)){
    de = read.csv(file.path[i], stringsAsFactors = F)
    de.name = file.name[i]
    Call_goseq(de=de, de.name=de.name)
}


