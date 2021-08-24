# Description:
# PCA analysis with DESeq2 normalized (estimateSizeFactors) and transformed
# (vst) values

library(DESeq2)
library(factoextra)
library(SummarizedExperiment)
library(RColorBrewer)
library(plotly)
library(scales)


# Repeated for ARPE-19 and MIO-M1 analyses
path = '/analysis/arpe19_miom1/results'
dir.create(file.path(path, 'PCA'), showWarnings = F)

# import counts and mapping
counts = read.csv(file.path(path, 'counts.csv'), row.names = 1)
mapping = read.csv(file.path(path, 'mapping.csv'))

rownames(mapping) = mapping$sample_id
dds = DESeqDataSetFromMatrix(countData = counts, colData = mapping, 
                             design = ~ condition)
dds = estimateSizeFactors(dds)
vsd = vst(dds, blind = T)


# num of components to analyze:
num_component = 10

# order of conditions in PCA
cond.order = c('Control', 'H2O2', 'H2O2 + RSG')

# keep top 15000 genes
td = assay(vsd)
td.mean = data.frame(rank=seq(1:length(rownames(td))), mean = rowMeans(td))
td.mean = td.mean[order(td.mean$mean, decreasing = T), ]
td = td[td.mean$rank[1:15000],]
td = t(td)

# calculate pca and %variance
td.pca = prcomp(td)
td.pca.eig = fviz_eig(td.pca)$data

# parse out pca results
pca.cord = as.data.frame(td.pca$x)
pca.cord$cond = sapply(rownames(pca.cord), function(x) 
    mapping$condition[mapping$sample_id==x], USE.NAMES = F)
pca.cord$cond = factor(pca.cord$cond, levels = cond.order)

# 3D plot of dim 1 to 3
colors = c('#bfbfbf', '#D55E00', '#56B4E9')

dim1.label = paste0('Dimension 1 (', format(round(td.pca.eig$eig[1], 2), nsmall = 2),'%)')
dim2.label = paste0('Dimension 2 (', format(round(td.pca.eig$eig[2], 2), nsmall = 2),'%)')
dim3.label = paste0('Dimension 3 (', format(round(td.pca.eig$eig[3], 2), nsmall = 2),'%)')

# plot 3D PCA plot
# some axis reversed due to default plot_ly plotting behavior
fig = plot_ly(pca.cord, x = ~PC1, y = ~PC2, z = ~PC3, color = ~cond, colors = colors)
fig = fig %>% add_markers(symbol = ~cond, symbols = c('circle', 'square', 'diamond'))
fig = fig %>% layout(scene = list(xaxis = list(title = dim1.label),
                                  yaxis = list(title = dim2.label, autorange='reversed', dtick = 2),
                                  zaxis = list(title = dim3.label, dtick = 2)),
                     font = list(size=20))
fig
