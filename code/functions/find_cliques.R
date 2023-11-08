


options(scipen = 100, digits = 4)
#library(purrr)
#library(furrr)
#library(parallel)
#library(data.table)
#library(dplyr)

pacman::p_load(purrr, furrr, parallel, data.table, dplyr,stringr,gtools)


find_cliques=function(binsize=500000,
                      data_dir='/afs/cs.wisc.edu/p/keles/Collab_2022/VolumeA/Manybody/data/Ramani2017/tmp/',
                      output_dir='/afs/cs.wisc.edu/p/keles/Collab_2022/VolumeA/Manybody/results/Ramani2017/tmp/',
                      celltypenum=4,
                      corenum=20,
                      chrnum=23,
                      sizefile='hg19.chrom.sizes',
                      ncellsthreshold=3){
  
  
  
  chrlist=c(paste0("chr",c(1:(chrnum-1),"X")))
  
  
  setwd(data_dir)
  cell_type=qs::qread("cell_type.qs")
  ctlists=unique(cell_type[,2])
  
  future::plan(multicore, workers = celltypenum)
  
  
  future_map(ctlists,function(ct){   
    
    
    
    hic_df_ct=lapply(1:chrnum,function(x)qs::qread(paste0(chrlist[x],'/hic_df_',chrlist[x],'_',ct,'.qs')) )
    hic_df_ct=do.call('rbind',hic_df_ct)
    hic_df_ct=split(hic_df_ct,by='cell')    
    
    future::plan(multicore, workers = corenum)
    
    cliques_loop=future_map(hic_df_ct,function(tmp){   
      
      
      nodes<-c(paste0(tmp$chr,"_",tmp$binA),paste0(tmp$chr,"_",tmp$binB)) %>% unique      
      vertices=data.frame(nodes)    
      edges <-tmp[,c('from','to'):=.(paste0(chr,"_",binA),paste0(chr,"_",binB))] #%>% mutate(from=paste0(V1,"_",V2),to=paste0(V3,"_",V4))
      edges<-edges[,c('from','to')] #%>% dplyr::select(from,to,V5)
      g <- igraph::graph_from_data_frame(edges, directed=FALSE, vertices=nodes)
      three_clique<-igraph::cliques(g,min=3)
      
      
      return(three_clique)    
    })
    rm(hic_df_ct)
    gc()
    
    
    qs::qsave(cliques_loop,paste0(output_dir,"/cliques_loop_3tomax_",ct,".qs"))    
    print(paste0(ct," done"))    
    
    
    
    
    
  }           
  
  
  )
  
  
  
  cat('Clique finding for each cell coimplete...')
  cat("\n")
  
  
  cliques_loop_ct=list()
  for(ct in ctlists){
    
    cliques_loop_ct[[ct]]=qs::qread(paste0(output_dir,"/cliques_loop_3tomax_",ct,".qs"))        
    
  }
  qs::qsave(cliques_loop_ct,paste0(output_dir,"/cliques_loop_3tomax.qs"))
  
  
  
  
  cat('Saving cliques')
  cat("\n")
  
  
  
  
#  cliques_loop_ct=qs::qread(paste0(output_dir,"/cliques_loop_3tomax.qs"))
  future::plan(multicore, workers = celltypenum)
#  Q=future_map(1:length(cliques_loop_ct),function(ct){as.character(unlist(lapply(cliques_loop_ct[[ct]],function(y)unlist(lapply(lapply(y,function(x)names(x)),function(x)paste(mixedsort(x),collapse = "-")) )))) }) %>% unlist seems unnecessary.

Q=future_map(1:length(cliques_loop_ct),function(ct){as.character(unlist(lapply(cliques_loop_ct[[ct]],function(y)unlist(lapply(lapply(y,function(x)names(x)),function(x)paste(x,collapse = "-")) )))) }) %>% unlist      
      
      
  qs::qsave(Q,paste0(output_dir,"/Q_notunique.qs"))   
  tableQ=Q %>% table
  
  
  print('Sorting cluques and generating pairwise interactions within each clique...')
  cat("\n")
  
  
  Q_unique=names(tableQ)
  cliquesize<-sapply(strsplit(Q_unique, "-"), length)
  Smax=max(cliquesize)
  
  
  chrlist=c(paste0("chr",c(1:(chrnum-1),"X")))
  options(future.globals.maxSize= Inf)  
  #  for(S in 3:Smax){
  future_map(4:Smax,function(S){    
    assign(paste0("Q",S),Q_unique[cliquesize==S]) 
    
    
    
    qs::qsave(lapply(chrlist,function(x){get(paste0("Q",S))[word(get(paste0("Q",S)),1,sep="_")==x]}),
              paste0(output_dir,'/Q',S,'_filtered_list.qs'))    
  }
  
  )   
  
  Q3=unique(Q[Q%in%Q_unique[cliquesize==3&tableQ>ncellsthreshold]]) 
  qs::qsave(lapply(chrlist,function(x){Q3[word(Q3,1,sep="_")==x]}),
            paste0(output_dir,'/Q3_filtered_list.qs')) 
  
  
  #  for(S in 3:Smax){
  future_map(3:Smax,function(S){     
    tmp_list=list()
    for(chr in 1:chrnum){
      
      tmp=lapply(qs::qread(paste0(output_dir,'/Q',S,'_filtered_list.qs'))[[chr]] ,function(x)strsplit(x, "-"))
      tmp_list[[chr]]=lapply(tmp,function(x)apply(combn(x[[1]],2),2,function(x)paste(x,collapse = "-")))
      
    }
    qs::qsave(tmp_list,paste0(output_dir,'/pairwise_',S,'_filtered_list.qs'))
  }
  #  }
  
  
  
  )
  
  
  
  
  
  
  
  
}
