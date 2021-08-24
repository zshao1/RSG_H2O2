# Description:
# For two comparisons, plot DEG in first comparison (by default) with 
# FC in comp_1 on x-axis and FC in comp_2 on y-axis.

library(ggplot2)

# Repeated for ARPE-19 and MIO-M1 analyses
path = '/analysis/arpe19_miom1/results'

# a box within axis limit will plotted
axis.limit = 1

# |log2FC| limit for being considered as DE
low_log2fc_cutoff = 0

dir.create(file.path(path, 'DEG_comparison_FC'), showWarnings = F)

plot.fc = function(df1, df2, df1.name, df2.name, selected.DEG=NULL,
                   selected.DEG.name='', use.second.deg=F){
    
	# function will plot DE genes in comp_1, unless selected.DEG is not NULL
	# in that case, DE genes in selected.DEG vector will be plotted instead and
	# selected.DEG.name will need to be specified
	
	# use.second.deg=T to select for the DEG in comp_2 instead of comp_1 (default)
	
	
    if (use.second.deg==F){
        # DE genes in first comparison
        df1.deg = df1$id[df1$FDR<0.05 & abs(df1$logFC)>low_log2fc_cutoff]
        
        # check which df1.deg are not in df2, those not in df2 are removed
        deg.common = df1.deg[df1.deg %in% df2$id]
        print(paste0(length(df1.deg)-length(deg.common), 
                     ' DE genes in comparison 1 are removed as they are not found in comparison 2'))
        
        if (!is.null(selected.DEG)){
            len.deg.common = length(deg.common)
            deg.common = deg.common[deg.common %in% selected.DEG]
            print(paste0(length(deg.common),'/',len.deg.common, ' genes selected'))
        }
    } else if (use.second.deg==T){
        # DE genes in second comparison
        df2.deg = df2$id[df2$FDR<0.05 & abs(df2$logFC)>low_log2fc_cutoff]
        
        # check which df2.deg are not in df1, those not in df1 are removed
        deg.common = df2.deg[df2.deg %in% df1$id]
        print(paste0(length(df2.deg)-length(deg.common), 
                     ' DE genes in comparison 2 are removed as they are not found in comparison 1'))
        
        if (!is.null(selected.DEG)){
            len.deg.common = length(deg.common)
            deg.common = deg.common[deg.common %in% selected.DEG]
            print(paste0(length(deg.common),'/',len.deg.common, ' genes selected'))
        }
    }
    
    
    # FC in the two comparisons
    deg.common.x = sapply(deg.common, function(x) df1$logFC[df1$id==x], 
                              USE.NAMES = F)
    deg.common.y = sapply(deg.common, function(x) df2$logFC[df2$id==x], 
                              USE.NAMES = F)
    
    df.plot = data.frame(gene=deg.common, fc.1=deg.common.x,
                         fc.2=deg.common.y)
    
    # remove genes outside of the plotting window
    len.incl.outlier = length(df.plot$gene)
    df.plot = df.plot[abs(df.plot$fc.1)<axis.limit & abs(df.plot$fc.2)<axis.limit,]
    len.no.outlier = length(df.plot$gene)
    print(paste0((len.incl.outlier-len.no.outlier), '/', len.incl.outlier, 
                 ' genes are removed as they are outside plotting window'))

    # pearson's r correlation
    pearson.r = cor(df.plot$fc.1, df.plot$fc.2, method = "pearson")
    pearson.r = round(pearson.r, 4)
    
    # genes in each quadrant
    # qudarant 1 on top left and numbered by clockwise fashion
    q1 = length(df.plot$gene[df.plot$fc.1<0 & df.plot$fc.2>0])
    q2 = length(df.plot$gene[df.plot$fc.1>0 & df.plot$fc.2>0])
    q3 = length(df.plot$gene[df.plot$fc.1>0 & df.plot$fc.2<0])
    q4 = length(df.plot$gene[df.plot$fc.1<0 & df.plot$fc.2<0])
    
	# specify plotting title
    if (selected.DEG.name==''){
        gg.title = paste0("Pearson's r = ", pearson.r)
    } else {
        gg.title = paste0(selected.DEG.name, "\nPearson's r = ", pearson.r)
    }
    
	# plot
    gg = ggplot(df.plot, aes(x=fc.1, y=fc.2)) + 
        geom_point(alpha=0.3) +
        geom_hline(yintercept=0, linetype="dashed", color = "blue") +
        geom_vline(xintercept=0, linetype="dashed", color = "blue") +
        xlim((-1 * axis.limit), axis.limit) +
        ylim((-1 * axis.limit), axis.limit) +
        xlab(paste0('log2FC - ', df1.name)) +
        ylab(paste0('log2FC - ', df2.name)) +
        ggtitle(gg.title) +
        theme(plot.title = element_text(hjust = 0.5)) +
        annotate("text", label = paste0(q1, ' Genes'), 
			x = (-1*axis.limit+0.25), y = axis.limit, size = 4, colour = "red")+
        annotate("text", label = paste0(q2, ' Genes'), 
			x = (axis.limit-0.25), y = axis.limit, size = 4, colour = "red")+
        annotate("text", label = paste0(q3, ' Genes'), 
			x = (axis.limit-0.25), y = -1*axis.limit, size = 4, colour = "red")+
        annotate("text", label = paste0(q4, ' Genes'), 
			x = (-1*axis.limit+0.25), y = -1*(axis.limit), size = 4, colour = "red") +
        geom_smooth(method = "lm", se = FALSE) +
        theme(axis.text.x = element_text(color='black'),
              axis.text.y = element_text(color='black'))
    
    ggsave(file.path(path, 'DEG_comparison_FC', paste0('FCPlot_df1=',df1.name,
                                                       '_df2=',df2.name,
                                                       selected.DEG.name,'.png')), 
           plot = gg, device = 'png',
           dpi = 300, width = 4, height = 4 )
}

# plot all H2O2 v control DE genes
df1 = read.csv(file.path(path, 'edgeR', '20180731_H2O2_v_control.csv'), 
               stringsAsFactors = F)
df2 = read.csv(file.path(path, 'edgeR', '20180731_H2O2_RSG_v_H2O2.csv'), 
               stringsAsFactors = F)
df1.name = 'H2O2 vs. Control'
df2.name = 'H2O2 + RSG vs. H2O2'

plot.fc(df1, df2, df1.name, df2.name)



# plot genes from selected GO BP with H2O2 vs Control DE genes
GOs = list('GO:0030198'='Extracellular matrix organization',
           'GO:0072359'='Circulatory system development',
           'GO:0019222'='Regulation of metabolic process',
           'GO:0034097'='Response to cytokine',
           'GO:0008219'='Cell death',
           'GO:0008283'='Cell proliferation',
           'GO:0002376'='Immune system process',
           'GO:0007155'='Cell adhesion',
           'GO:0016477'='Cell migration',
           'GO:0050896'='Response to stimulus',
           'GO:0070848'='Response to growth factor',
           'GO:0032502'='Developmental process')

# goseq output file that contains list of DE genes in each GO BP
enriched_GOs = read.csv(file.path(path, 'goseq','H2O2_v_control',
                                  'goseq_log2FC=0_relaxed_both_BP.csv'))

for (i in seq(GOs)) 
{
    GOs_i = GOs[i]
    selected.DEG.name=GOs_i[[1]]
    selected.DEG=enriched_GOs$genes[enriched_GOs$go.id==names(GOs_i)]
    selected.DEG = strsplit(as.character(selected.DEG), split = ';')[[1]]
    plot.fc(df1, df2, df1.name, df2.name, 
            selected.DEG = selected.DEG,
            selected.DEG.name = paste0(selected.DEG.name))
}



# plot 'GO:0042060'='Wound healing' with RSG_H2O2_V_H2O2 DE genes

# goseq output file that contains list of DE genes in each GO BP
enriched_GOs = read.csv(file.path(path, 'goseq','H2O2_RSG_v_H2O2',
                                  'goseq_log2FC=0_relaxed_both_BP.csv'))

selected.DEG.name='Wound healing'
selected.DEG=enriched_GOs$genes[enriched_GOs$go.id=='GO:0042060']
selected.DEG = strsplit(as.character(selected.DEG), split = ';')[[1]]
plot.fc(df1, df2, df1.name, df2.name, 
        selected.DEG = selected.DEG,
        selected.DEG.name = selected.DEG.name,
        use.second.deg=T)

