# Description: 
# Run ERSSA to test if sample size sufficient for DE discovery

library(ERSSA)

# Repeated for ARPE-19 and MIO-M1 analyses
path = '/analysis/arpe19_miom1/results'
dir.create(file.path(path, 'erssa'), showWarnings = F)

mapping_all = read.csv(file.path(path, 'mapping.csv'))
counts_all = read.csv(file.path(path, 'counts.csv'), row.names = 1)


call_erssa = function(test_ctrl, name){
    
    # function that take edgeR analysis function input and start erssa with edgeR
    
    cond.A = test_ctrl[1]
    cond.B = test_ctrl[2]
    
    mapping = mapping_all[mapping_all$condition %in% 
                              c(cond.A, cond.B),]
    counts = counts_all[,mapping$sample_id]
    condition_table = mapping[,c('sample_id','condition')]
    
    erssa.path = file.path(path, 'manuscript','erssa', name)
    dir.create(erssa.path, showWarnings = F)
    
    # CPM filter set to 4 to match edgeR
    erssa = erssa(count_table = counts, condition_table = condition_table, 
                  DE_ctrl_cond = cond.B, comb_gen_repeat = 50, 
                  DE_cutoff_Abs_logFC = 0.25, num_workers = 22, path = erssa.path, 
                  save_log=T, marginalPlot_stat = 'median', filter_cutoff = 4)
    save(erssa, file = file.path(erssa.path, 'erssa.RData'))
}

set.seed(1)
call_erssa(test_ctrl = c('H2O2', 'Control'),
           name = 'H2O2_v_control')

set.seed(1)
call_erssa(test_ctrl = c('H2O2 + RSG', 'H2O2'),
           name = 'H2O2_RSG_v_H2O2')

set.seed(1)
call_erssa(test_ctrl = c('H2O2 + RSG', 'Control'),
           name = 'H2O2_RSG_v_control')