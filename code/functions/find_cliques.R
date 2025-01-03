options(scipen = 100, digits = 4)


pacman::p_load(purrr, furrr, parallel, data.table, dplyr,stringr,gtools,igraph)


find_cliques=function(binsize=500000,
                      data_dir='/storage10/kwangmoon/MINTsC/data/Ramani2017/',
                      output_dir='/storage10/kwangmoon/MINTsC/results/Ramani2017/',
                      corenum_celltype=1,
                      corenum_cells=20,
                      chrnum=23,
                      sizefile='hg19.chrom.sizes',
                      ncellsthreshold_c0=1,ncellsthreshold_c1=0,Smax=10,
                      chrlist=NULL){
  
  
  if(is.null(chrlist)){chrlist=c(paste0("chr",c(1:(chrnum-1),"X")))}  
  
  
  setwd(data_dir)
  cell_type=qs::qread("cell_type.qs")
  ctlists=unique(cell_type[,2])
  
  future::plan(multicore, workers = corenum_celltype)
  
  
  future_map(ctlists,function(ct){   
    
    
    
    hic_df_ct=lapply(1:chrnum,function(x)qs::qread(paste0(chrlist[x],'/hic_df_',chrlist[x],'_',ct,'.qs'))[,c('cell','chr','binA','binB')] )
    hic_df_ct=do.call('rbind',hic_df_ct)
    hic_df_ct=split(hic_df_ct,by='cell')    
    
    future::plan(multicore, workers = corenum_cells)
    system(paste0('mkdir ',output_dir,"/","'",ct,"'"))       
    future_map(hic_df_ct,function(tmp){   
      
      if(!paste0("cell_clique_",tmp$cell[1],'.qs')%in%list.files(paste0(output_dir,"/",ct))){
        
        
        nodes<-c(paste0(tmp$chr,"_",tmp$binA),paste0(tmp$chr,"_",tmp$binB)) %>% unique      
        vertices=data.frame(nodes)    
        
        edges <-tmp[,.('from'=paste0(chr,"_",binA),'to'=paste0(chr,"_",binB))] 
        
        g <- igraph::graph_from_data_frame(edges, directed=FALSE, vertices=nodes)
        three_clique<-igraph::cliques(g,min=3,max=Smax)
        
        tmpres=unlist(lapply(three_clique,function(x)paste(mixedsort(names(x)),collapse="-")))
        
        qs::qsave(tmpres, paste0(output_dir,"/",ct,"/cell_clique_",tmp$cell[1],'.qs'))        
        
        
        
      }
      
      
    })
    
    #If you do setwd within lapply, it will globally setwd as well. be careful.                     
    cliques_loop=lapply(list.files(paste0(output_dir,"/",ct)),function(x){qs::qread(paste0(output_dir,"/",ct,"/",x))})                     
    
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
  
  
  
  
  
  
  Q=unlist(cliques_loop_ct)
  qs::qsave(Q,paste0(output_dir,"/Q_notunique.qs"))   
  tableQ=Q %>% table
  
  
  print('Sorting cluques and generating pairwise interactions within each clique...')
  cat("\n")
  
  
  Q_unique=names(tableQ)
  cliquesize<-sapply(strsplit(Q_unique, "-"), length)
  chrlabel=word(Q_unique,1,sep="_")                                                                                                                               
  Smax=max(cliquesize)
  
  
  
  options(future.globals.maxSize= Inf)
  future::plan(multicore, workers = length(3:Smax))                                                                                                                                 
  
  
  Qdat=data.table(Q_unique,chrlabel,cliquesize,ncell=as.numeric(tableQ ))
  rm(Q_unique)
  rm(chrlabel)                                                                                                                                 
  rm(tableQ)
  rm(cliquesize)
  gc()
  qs::qsave(Qdat,
            paste0(output_dir,'/Q_summary.qs')) 
  
  
  
  qs::qsave(lapply(chrlist,function(x){Qdat[cliquesize==3&chrlabel==x&ncell>ncellsthreshold_c0]$Q_unique}),
            paste0(output_dir,'/Q3_filtered_list.qs')) 
  
  
  
  future_map(4:Smax,function(S){    
    
    
    
    qs::qsave(lapply(chrlist,function(x){Qdat[cliquesize==S&chrlabel==x&ncell>ncellsthreshold_c1]$Q_unique}),
              paste0(output_dir,'/Q',S,'_filtered_list.qs'))                                                                                                                                  
    
  }
  
  )   
  
  #  Q3=unique(Q[Q%in%Q_unique[cliquesize==3&tableQ>ncellsthreshold]])
  # chrlabel_3=word(Q3,1,sep="_")                                                                                                                                 
  #  qs::qsave(lapply(chrlist,function(x){Q3[chrlabel_3==x]}),
  #            paste0(output_dir,'/Q3_filtered_list.qs')) 
  
  
  
  
  rm(Qdat)
  gc()                                                                                                                                 
  
  
  
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
