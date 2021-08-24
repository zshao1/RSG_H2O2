# Description:
# Serves several functions:
    # 1. Take all featureCounts-generated count files and form a count table
    # with sample ID as columns and ensembl ID as rows (minus the version number).
    # 2. generate mapping table
    # 3. generate gene length table
    # 4. generate TPM table

# Repeated for ARPE-19 and MIO-M1 analyses
result_path = '/analysis/arpe19_miom1/results/'

# generate sample mapping csv
count_path = '/data/secondary/arpe19_miom1/star_featureCounts'

# identify .count files generated from featureCounts
files = list.files(count_path, full.names = T, pattern = '*.counts$')
file_name = tools::file_path_sans_ext(list.files(count_path, pattern = '*.counts$'))

condition = c(rep('Control', 6),
              rep('H2O2 + RSG', 6),
              rep('H2O2', 6))

mapping = data.frame(sample_id = file_name, condition = condition,
                     path = files)
write.csv(mapping, paste0(result_path, '/mapping.csv'), row.names = F)


# generate count table

counts = list()
gene_id = list()

for (i in seq(files)) {
    df = read.csv(files[i], stringsAsFactors = F, sep = '\t', skip = 1)
    counts[[file_name[i]]] = df[,7]
    gene_id[[file_name[i]]] = df[,1]
}

# check all count files have the same gene ID in the same order
id.all.true = all(sapply(seq(length(names(gene_id)))[-1], function(x)
                    all.equal(gene_id[[x-1]], gene_id[[x]])))

if (id.all.true==T){
    
    # generate count table, add gene id as row name
    fc_counts = as.data.frame(counts)
    rownames(fc_counts) = gene_id[[1]]
    
    # remove rows with gene ID ending in '_PAR_Y'
    # these genes are duplicate between X and Y chromosome
    keep=sapply(rownames(fc_counts), function(x) 
        substr(x, nchar(x)-5, nchar(x))!='_PAR_Y', USE.NAMES = F)
    fc_counts = fc_counts[keep,]
    
    # remove version number in gene id
    gene.id = sapply(rownames(fc_counts), function(x) 
                        strsplit(x, '\\.')[[1]][1], USE.NAMES = F)
    rownames(fc_counts) = gene.id
    write.csv(fc_counts, paste0(result_path,'/counts.csv'))
} else {
    stop('Not all gene IDs are same or in same order across the .count files')
}

# output gene length csv, limited to genes without 'PAR_Y' ending
genes = gene_id[[1]][keep]
df = read.csv(files[1], stringsAsFactors = F, sep = '\t', skip = 1)
length = sapply(genes, function(x) df$Length[df$Geneid==x], USE.NAMES = F)
genes = sapply(genes, function(x) 
    strsplit(x, '\\.')[[1]][1], USE.NAMES = F)

df_length = data.frame(row.names=genes, length=length)
write.csv(df_length, paste0(result_path, '/gene_length.csv'))

# calculate TPM based on count and gene length tables
# first check that count table and gene length are in same order of genes
if (all(rownames(fc_counts)==rownames(df_length))){
    # TPM_i = 10^6 * (n_i/l_i)/sum(n_j/l_j)
    
    # calculate TPM without appending one to all genes
    TPM_1 = fc_counts/df_length$length
    TPM_2 = colSums(TPM_1)
    TPM_3 = t(t(TPM_1)/TPM_2)
    TPM = 1000000 * TPM_3
    write.csv(TPM, paste0(result_path, '/TPM.csv'))
    
    # recalculate with one appended to all genes
    TPM_1 = (fc_counts+1)/df_length$length
    TPM_2 = colSums(TPM_1)
    TPM_3 = t(t(TPM_1)/TPM_2)
    TPM = 1000000 * TPM_3
    write.csv(TPM, paste0(result_path, '/TPM_nonZero.csv'))
    
} else {
    stop('count and gene length tables are not in the same order of genes')
}


