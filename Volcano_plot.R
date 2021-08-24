# Description:
# Loop through each edgeR comparison and generate Volcano plot consist of
# log2FC on x-axis and FDR on y-axis

library(ggplot2)

# Repeated for ARPE-19 and MIO-M1 analyses
path = '/analysis/arpe19_miom1/results'

# |log2FC| limit for being considered as DE
low_log2fc_cutoff = 0

FC_plot_limit = 1 # set limit of log2FC plotting window, x-axis
FDR_plot_limit = 20 # set limit of -log10FDR plotting window, y-axis

# assign name to each comparison
comparisons = list.files(file.path(path, 'edgeR'), pattern = '*.csv$')
comp.name = c('H2O2\nvs. Control',
              'H2O2 + RSG\nvs. H2O2',
              'H2O2 + RSG\nvs. Control')

for (i in seq(length(comparisons))){
    
    # read in comparison, identify number of DEGs
    df = read.csv(file.path(path, 'edgeR', comparisons[i]), 
                  stringsAsFactors = F)
    
    # statistics of DE
    num.de = sum(df$FDR<0.05 & abs(df$logFC) > low_log2fc_cutoff)
    num.nonde = length(df$id) - num.de
    
    # take -1 * log10 of FDR
    df$plot.FDR = sapply(df$FDR, function(x) (-1 * log10(x)), USE.NAMES = F)
    
    # assign red for DEG, black for others
    df$color = sapply(seq(df$FDR), function(x) {
        if (df$FDR[x] < 0.05 & abs(df$logFC[x]) > low_log2fc_cutoff){
            color='red'
        } else {
            color='black'
        }
        return(color)
    }, USE.NAMES = F)
    
    # set plotting log2FC, where |log2FC|>FC_plot_limit are assigned to limit
    df$plot.log2FC = sapply(df$logFC, function(x){
        if (x > FC_plot_limit){
            FC_i = FC_plot_limit
        } else if (x < (-1 * FC_plot_limit)){
            FC_i = -1 * FC_plot_limit
        } else {
            FC_i = x
        }
        return(FC_i)
    }, USE.NAMES = F)
    
    # set plotting symbol, dot for genes within log2FC limit, bracket of
    # right direction for those outside
    # also set triangle for those with FDR too small
    df$plot.symbol = sapply(seq(dim(df)[1]), function(x){
        if (df$logFC[x] > FC_plot_limit){
            symbol = 'right'
        } else if (df$logFC[x] < (-1 * FC_plot_limit)){
            symbol = 'left'
        } else if (df$plot.FDR[x] > FDR_plot_limit){
            symbol = 'top'
        }
        else {
            symbol = 'in'
        }
        return(symbol)
    }, USE.NAMES = F)
    
    # set df$plot.FDR to FDR limit
    df$plot.FDR = sapply(df$plot.FDR, function(x){
        if (x > FDR_plot_limit){
            FDR_i = FDR_plot_limit
        } else {
            FDR_i = x
        }
        return(FDR_i)
    }, USE.NAMES = F)
    
    # add points outside of view for legend
    # in case there are no DE genes for red legend and genes associated with symbols
    red.dot = data.frame(id='red', name='', logFC=0, logCPM=0, PValue=0,
                         FDR=0, color='red', plot.log2FC=0, plot.symbol='in',
                         plot.FDR=-10)
    left.dot = data.frame(id='red', name='', logFC=0, logCPM=0, PValue=0,
                          FDR=0, color='red', plot.log2FC=0, plot.symbol='left',
                          plot.FDR=-10)
    right.dot = data.frame(id='red', name='', logFC=0, logCPM=0, PValue=0,
                           FDR=0, color='red', plot.log2FC=0, plot.symbol='right',
                           plot.FDR=-10)
    top.dot = data.frame(id='red', name='', logFC=0, logCPM=0, PValue=0,
                         FDR=0, color='red', plot.log2FC=0, plot.symbol='top',
                         plot.FDR=-10)
    df = rbind(df, red.dot, left.dot, right.dot, top.dot)
    
    # set levels
    df$color = factor(df$color, levels = c('red', 'black'))
    df$plot.symbol = factor(df$plot.symbol, levels=c('left', 'in', 'right', 'top'))
    
    # plot results
    gg =  ggplot(df, aes(x=plot.log2FC, y=plot.FDR, color=color, 
                         shape=plot.symbol)) + 
        geom_point(alpha = 0.5, size = 1) +
        scale_shape_manual(values = c(60, 16, 62, 17), guide=FALSE) +
        labs(title=comp.name[i], 
             x='Fold change (Log2FC)', 
             y='Statistical significance (-log10FDR)',
             color=NULL) +
        theme(plot.title = element_text(hjust = 0.5),
              legend.title.align=0.5,
              legend.position='bottom',
              legend.direction = "vertical") +
        geom_hline(yintercept=(-1 * log10(0.05)), linetype='dashed', color='blue') +
        ylim(c(0, FDR_plot_limit)) +
        xlim(c((-1 * FC_plot_limit), FC_plot_limit))
    
    if (low_log2fc_cutoff==0){
        gg = gg + scale_color_manual(values=c("#FF0000", "#000000"),
                                     labels = c(paste0('FDR < 0.05 (', num.de, ' genes)'),
                                                paste0('Non-significant (', num.nonde, ' genes)')))
    } else {
        gg = gg + scale_color_manual(values=c("#FF0000", "#000000"),
                                     labels = c(paste0('FDR < 0.05 & |log2FC| > ',low_log2fc_cutoff, 
                                                       ' (', num.de, ' genes)'),
                                                paste0('Non-significant (', num.nonde, ' genes)'))) 
    }
    
    
    ggsave(filename = file.path(path, 'edgeR', paste0('Volcano_', 
                  tools::file_path_sans_ext(comparisons[i]),'.png')), 
           plot = gg, dpi = 300, device='png', width = 5, height = 4)
    
    warning('Removing 4 rows in ggplot2 is normal')

}