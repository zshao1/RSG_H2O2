# Description:
# For two comparisons, plot DEG in both comparisons with FC in comp_1 on x-axis
# and FC in comp_2 on y-axis

library(ggplot2)

path = '/analysis/arpe19_miom1/results'
dir.create(file.path(path, 'DEG_comparison_FC'), showWarnings = F)

# a box within axis limit will plotted
axis.limit = 1

# |log2FC| limit for being considered as DE
low_log2fc_cutoff = 0


plot.fc = function(df1, df2, df1.name, df2.name, deg.select){
    
	# deg.select can either be union or intersect of DE genes
    
    print('Plotting FC correlation with DEGs from BOTH COMPARISONS')
    
    # DE genes in first comparison
    df1.deg = df1$id[df1$FDR<0.05 & abs(df1$logFC)>low_log2fc_cutoff]
    df2.deg = df2$id[df2$FDR<0.05 & abs(df2$logFC)>low_log2fc_cutoff]
    
    # combine both DEG lists and select DEGs either union or intersect
    # based on the need of analysis
    if (deg.select=='union'){
        deg.union = union(df1.deg, df2.deg)
        deg.common = deg.union[deg.union %in% df1$id]
        deg.common = deg.common[deg.common %in% df2$id]
        print(paste0(length(deg.union)-length(deg.common), 
     ' DE genes are removed as they are found in only one of two comparisons'))
    } else if (deg.select=='intersect'){
        deg.common = intersect(df1.deg, df2.deg)
        print(paste0(length(deg.common), ' intersect DEGs found in both comparisons'))
    } else {
        warning('Only union or intersection are valid options for selecting DEGs')
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
  
    gg = ggplot(df.plot, aes(x=fc.1, y=fc.2)) + 
        geom_point(alpha=0.3) +
        geom_hline(yintercept=0, linetype="dashed", color = "blue") +
        geom_vline(xintercept=0, linetype="dashed", color = "blue") +
        xlim((-1 * axis.limit), axis.limit) +
        ylim((-1 * axis.limit), axis.limit) +
        xlab(paste0('log2FC - ', df1.name)) +
        ylab(paste0('log2FC - ', df2.name)) +
        ggtitle(paste0(
                       "\nPearson's r = ", pearson.r)) +
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
                                                       '_df2=',df2.name,'.png')), 
           plot = gg, device = 'png',
           dpi = 300, width = 4, height = 4 )
}

# H2O2 vs control between ARPE-19 and MIO-M1
df1 = read.csv('/analysis/arpe19_miom1/results/edgeR/H2O2_v_control_arpe19.csv', 
               stringsAsFactors = F)
df2 = read.csv('/analysis/arpe19_miom1/results/edgeR/H2O2_v_control_miom1.csv', 
               stringsAsFactors = F)
df1.name = 'ARPE-19 - H2O2 vs. Control'
df2.name = 'MIO-M1 - H2O2 vs. Control'

plot.fc(df1, df2, df1.name, df2.name, deg.select='union')


# H2O2 + RSG vs. H2O2 between ARPE-19 and MIO-M1
df1 = read.csv('/analysis/arpe19_miom1/results/edgeR/H2O2_RSG_v_H2O2_arpe19.csv', 
               stringsAsFactors = F)
df2 = read.csv('/analysis/arpe19_miom1/results/edgeR/H2O2_RSG_v_H2O2_miom1.csv', 
               stringsAsFactors = F)
df1.name = 'ARPE-19 - H2O2 + RSG vs. H2O2'
df2.name = 'MIO-M1 - H2O2 + RSG vs. H2O2'

plot.fc(df1, df2, df1.name, df2.name, deg.select='union')


