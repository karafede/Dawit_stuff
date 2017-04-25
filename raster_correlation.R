### correlation 
z <- stack(raste_AVW_chlo_arr,raste_AOD_stack_arr)

r <- calc(z, fun=function(x) cor(x[1:42], x[55:96], method='spearman', use= "pairwise.complete.obs" ))

