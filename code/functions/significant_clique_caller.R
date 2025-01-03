significant_clique_caller=function(output_subdir="/storage10/kwangmoon/MINTsC/results/Ramani2017/L20_lower1_upper40/",
                                   S=3,r=2,cell_type="GM12878",fdr_cutoff=0.05,plot.it=TRUE,
                                   chrlist=c(paste0("chr",c(1:22,"X"))),p.adjust.method="BH"){
  
  library(dplyr)
  num_chr=length(chrlist)
  
  porder_list=list()
  Zorder_list=list()
  for( chr in 1:num_chr){
    setwd(output_subdir)
    tmp_p=tryCatch(qs::qread(paste0('size',S,'/',chrlist[chr],"/","porder_",cell_type,"_DirMult_bandlevel_mle_Aset_filtered.qs")) %>% unlist%>% suppressWarnings
                   ,error=function(e){return(NA)})
    porder_list[[chr]]=tmp_p
    
    
    tmp_z=tryCatch(qs::qread(paste0('size',S,'/',chrlist[chr],"/","Zorder_",cell_type,"_DirMult_bandlevel_mle_Aset_filtered.qs")) %>% unlist%>% suppressWarnings
                   ,error=function(e){return(NA)})
    Zorder_list[[chr]]=tmp_z      
    
  }
  
  porder_allchr=(porder_list) %>% unlist
  Zorder_allchr=(Zorder_list) %>% unlist    
  set.seed(1)
  x=runif(length(porder_allchr))
  param_porder=pbeta(porder_allchr,r,choose(S,2)-r+1)
  fdr=p.adjust(param_porder, method ="BH")#<0.05    
  sig_group=sort(fdr,decreasing=TRUE)<fdr_cutoff
  sig_group_for_plot=sort(fdr,decreasing=TRUE)<fdr_cutoff
  if(plot.it){
    plot(sort(-log(x),na.last = TRUE),sort(-log(param_porder),na.last = TRUE),cex=2,cex.lab=1.5,cex.axis=1.5,cex.main=2,main=paste0(cell_type," p-value (Size=",S,")"),ylab="-log10(Observed)",xlab="-log10(Expected)",col=as.character(ifelse(sig_group_for_plot,"indianred","grey70")))
    legend("topleft", legend=c("FDR<0.05","FDR>0.05"), pch=16, col=c("indianred","grey70"),cex=2)
    abline(0,1,lty=2,col='blue',lwd=2)
    
    
  }
  df=data.frame(clique=names(porder_allchr),pscore=porder_allchr,zscore=Zorder_allchr,pvalue=param_porder,fdr=fdr,multiway_interaction=sig_group)
  rownames(df)=NULL
  return(df)
  
  
}