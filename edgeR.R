# Description:
# run edgeR to identify DE genes

library('edgeR')

# Repeated for ARPE-19 and MIO-M1 analyses
path = '/analysis/arpe19_miom1/results/'
dir.create(file.path(path, 'edgeR'), showWarnings = F)

# change n to number of replicates per condition, remove genes with CPM below 
# threshold in less than n samples
replicate.num = 6

# cpm cutoff used to remove low expressing genes
cpm.cutoff = 4

mapping_all = read.csv(file.path(path, 'mapping.csv'))
counts_all = read.csv(file.path(path, 'counts.csv'), row.names = 1)

keep = rowSums(cpm(counts_all)>cpm.cutoff) >= replicate.num 
counts_all = counts_all[keep,]

DE_et_upgrade = function(df_diff){
    # modify result DE_et table with:
    # 1. add gene name
    # 2. sort by FDR
    # 3. add gene ID as a column
	# 4. change row.name to rank
    
	# df_id_to_name is a biomart generated table mapping ensembl id to gene name
	df_id_to_name = read.csv('/data/reference/human/ensembl_gene_id_to_gene_name/id_to_name.csv', 
                             stringsAsFactors = F)
							 
    df_diff_name = sapply(row.names(df_diff), function(x) {
        name = df_id_to_name$name[df_id_to_name$gene_no_index==x][1];
        if (length(name)==0 | is.na(name)) {
            name = ''
        };
        return(name)}, 
        USE.NAMES = F)
		
    df_diff$name = df_diff_name
    df_diff = df_diff[order(df_diff$FDR),]
    df_diff$id = row.names(df_diff)
    df_diff = df_diff[, c('id','name','logFC','logCPM','PValue','FDR')]
    rownames(df_diff) = seq(1, length(df_diff$id))
	
    return(df_diff)
}

edgeR_test = function(test_ctrl, name, additional.edgeR=NULL){
    # main function that parse out specific samples and run edgeR test
    # additiona.edgeR=T, then the csv output will be saved in a separate
    # folder in the main edgeR folder.
    
    mapping = mapping_all[mapping_all$condition %in% test_ctrl,]
    counts = counts_all[,mapping$sample_id]
    
    if (!(all(colnames(counts)==mapping$sample_id))){
        stop('count.table not in order of condition.table')
    }
	
    group = factor(mapping$condition)
    y = DGEList(counts=counts, group=group)
    y = calcNormFactors(y)
    y = estimateCommonDisp(y)
    y = estimateTagwiseDisp(y)
    
    et = exactTest(y, pair=c(test_ctrl[2], test_ctrl[1]))
    DE_et = topTags(et, n=1000000)$table
    
    df_diff = DE_et_upgrade(DE_et)
    
    if (is.null(additional.edgeR)){
        write.csv(df_diff,file=file.path(path,'edgeR', paste0(name, '.csv')), 
                  row.names = F)
    } else {
        dir.create(file.path(path, 'edgeR', 'additional'), showWarnings = F)
        write.csv(df_diff,file=file.path(path,'edgeR', 'additional',
                                         paste0(name, '.csv')), 
                  row.names = F)
    }
    
}

# parse out two conditions for comparison
# Control vs. Test
edgeR_test(test_ctrl = c('H2O2', 'Control'),
           name = 'H2O2_v_control')

edgeR_test(test_ctrl = c('H2O2 + RSG', 'H2O2'),
           name = 'H2O2_RSG_v_H2O2')

edgeR_test(test_ctrl = c('H2O2 + RSG', 'Control'),
           name = 'H2O2_RSG_v_control')

