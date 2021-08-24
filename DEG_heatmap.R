# Description:
# TPM without zero are used for heatmap
# Takes DEGs from selected comparisons
# Next, TPMs are normalized to mean of each gene
# Next, mean of each condition calculated
# Log2 of normalized expression are then plotted with genes clustered

library(ggplot2)
library(Hmisc)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(scales)

# Repeated for ARPE-19 and MIO-M1 analyses
path = '/analysis/arpe19_miom1/results'

# cap for normalized TPM values
# abs(norm.tpm) > cap will be set to cap
norm.tpm.cap = 0.6

heatmap.width = 5
heatmap.height = 10

# load TPM and mapping
tpm = read.csv(file.path(path, 'TPM_nonZero.csv'), row.names = 1, stringsAsFactors = F)
mapping = read.csv(file.path(path, 'mapping.csv'), stringsAsFactors = F)

# path and name of comparisons
comparisons = list.files(file.path(path, 'edgeR'), pattern = '*.csv$')
comp.name = c('H2O2 vs. Control',
              'H2O2 + RSG vs. H2O2',
              'H2O2 + RSG vs. Control')


# output dir
dir.create(file.path(path, 'heatmap'), showWarnings = F)

deg.hm.cond.only = function(deg, comparison.name) {
    # main heatmap plotting function
    # genes in deg vector will be parsed out of tpm to be plotted
    # comparison.name is the name of the comparison where the DEGs are extracted
    
    # parse out DEGs, divide row-wise by gene average, 
    # take normalized log2 of the values
    tpm.deg = tpm[deg,]
    tpm.deg.avg = rowMeans(tpm.deg)
    tpm.deg.norm.log2 = log2(tpm.deg/tpm.deg.avg)
    
    # calculate mean of each condition for each gene
    tpm.deg.norm.log2.avg = list()
    for (cond in conditions){
        tpm.deg.norm.log2.i = tpm.deg.norm.log2[,mapping$sample_id[mapping$condition==cond]]
        tpm.deg.norm.log2.i.avg = rowMeans(tpm.deg.norm.log2.i)
        tpm.deg.norm.log2.i.avg = as.vector(tpm.deg.norm.log2.i.avg)
        tpm.deg.norm.log2.avg[[cond]] = tpm.deg.norm.log2.i.avg
    }
    df.tpm.deg.norm.log2.avg = as.data.frame(tpm.deg.norm.log2.avg)
        
    # cap norm.tpm to cap on both directions
    for (col in colnames(df.tpm.deg.norm.log2.avg)){
        df.tpm.deg.norm.log2.avg[,col][df.tpm.deg.norm.log2.avg[,col] > norm.tpm.cap] = norm.tpm.cap
        df.tpm.deg.norm.log2.avg[,col][df.tpm.deg.norm.log2.avg[,col] < (-1* norm.tpm.cap)] = (-1 * norm.tpm.cap)
    }
        
    # code that assign colors
    colors = c('#bfbfbf', '#D55E00', '#56B4E9')
    names(colors) = unique(conditions)
    
    # plot heatmap
    ha = HeatmapAnnotation(df = data.frame(Condition=conditions),
                           col = list(Condition=colors), show_annotation_name = F,
                           simple_anno_size = unit(0.65, "cm"))
    f1 = colorRamp2(seq((-1 * norm.tpm.cap), norm.tpm.cap, length = 3), c("blue", "#ffffff", "red"))
    
    hm = Heatmap(df.tpm.deg.norm.log2.avg, col = f1, 
                 cluster_rows = TRUE, cluster_columns = FALSE,
                 name='Norm. \nExp.',
                 show_column_names = F, show_row_names = F,
                 top_annotation = ha, 
                 column_title_gp = gpar(fontsize = 16),
                 row_dend_width = unit(3, "cm"))
    png(file.path(path, 'heatmap', paste0('Heatmap_CondAvg_',comparison.name,'.png')), 
        width=heatmap.width, height=heatmap.height, res=200, units = 'in')
    draw(hm)
    dev.off()
}

# call heatmap function with each comparison's vector of DE genes as input.
    
    
