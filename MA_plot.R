# Description:
# Loop through each edgeR comparison and generate MA plot consist of
# logCPM on x-axis and log2FC on y-axis

library(ggplot2)

# Repeated for ARPE-19 and MIO-M1 analyses
path = '/analysis/arpe19_miom1/results'

# |log2FC| limit for being considered as DE
low_log2fc_cutoff = 0

# set MA plotting window
FC_limit = 1 # limit y-axis
log2CPM_min = 0 # limit minimum log2CPM

# assign name to each comparison
comparisons = list.files(file.path(path, 'edgeR'), pattern = '*.csv$')
comp.name = c('H2O2 vs. Control',
              'H2O2 + RSG vs. H2O2',
              'H2O2 + RSG vs. Control')

for (i in seq(length(comparisons))){
    
    # read in comparison, identify number of DEGs
    df = read.csv(file.path(path, 'edgeR', comparisons[i]), 
                  stringsAsFactors = F)
    
    # assign red for DEG, black for others
    df$color = sapply(seq(df$FDR), function(x) {
        if (df$FDR[x] < 0.05 & abs(df$logFC[x]) > low_log2fc_cutoff){
            color='red'
        } else {
            color='black'
        }
        return(color)
    }, USE.NAMES = F)
    
    # set plotting log2FC, where |log2FC|>FC_limit are assigned to limit
    df$plot.log2FC = sapply(df$logFC, function(x){
        if (x > FC_limit){
            FC_i = FC_limit
        } else if (x < (-1 * FC_limit)){
            FC_i = -1 * FC_limit
        } else {
            FC_i = x
        }
        return(FC_i)
    }, USE.NAMES = F)
    
    # set plotting symbol, dot for genes within log2FC limit, triangle of
    # right direction for those outside
    df$plot.symbol = sapply(df$logFC, function(x){
        if (x > FC_limit){
            symbol = 'up'
        } else if (x < (-1 * FC_limit)){
            symbol = 'down'
        } else {
            symbol = 'in'
        }
        return(symbol)
    }, USE.NAMES = F)
    
    # statistics of DE
    num.de = sum(df$FDR<0.05 & abs(df$logFC) > low_log2fc_cutoff)
    num.nonde = length(df$id) - num.de
    
    # add points outside of view for legend
    # in case there are no DE genes for red legend and genes associated with symbols
    red.dot = data.frame(id='red', name='', logFC=0, logCPM=-10, PValue=0,
                         FDR=-10, color='red', plot.log2FC=0, plot.symbol='in')
    up.dot = data.frame(id='red', name='', logFC=0, logCPM=-10, PValue=0,
                        FDR=-10, color='red', plot.log2FC=0, plot.symbol='up')
    down.dot = data.frame(id='red', name='', logFC=0, logCPM=-10, PValue=0,
                        FDR=-10, color='red', plot.log2FC=0, plot.symbol='down')
    df = rbind(df, red.dot, up.dot, down.dot)
    
    # set levels
    df$color = factor(df$color, levels = c('red', 'black'))
    df$plot.symbol = factor(df$plot.symbol, levels=c('in', 'up', 'down'))
    
    # plot results
    gg =  ggplot(df, aes(x=logCPM, y=plot.log2FC, color=color, 
                         shape=plot.symbol, fill = color)) + 
        geom_point(alpha = 0.5, size = 0.5) +
        scale_fill_manual(values=c("#FF0000", "#000000"), guide=FALSE) +
        scale_shape_manual(values = c(21, 24, 25), guide=FALSE) +
        labs(title=comp.name[i], 
             x='Gene Expression (log2CPM)', 
             y='Fold change (Log2FC)',
             color=NULL) +
        xlim(c(log2CPM_min, max(df$logCPM))) +
        ylim(c((-1 * FC_limit), FC_limit)) +
        theme(plot.title = element_text(hjust = 0.5),
              legend.title.align=0.5,
              legend.position = c(0.73, 0.875),
              legend.direction = "vertical",
              legend.background = element_rect(size=0.5, linetype="solid", 
                                               colour ="black")) +
        guides(color = guide_legend(override.aes = list(size=2)))
    
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
    
    ggsave(filename = file.path(path, 'edgeR', paste0('MA_', 
                  tools::file_path_sans_ext(comparisons[i]),'.png')), 
           plot = gg, dpi = 300, device='png', width = 4.5, height = 4)
    
    warning('Removing 3 rows in ggplot2 is normal')
}