options(scipen = 100, digits = 4)
pacman::p_load(splines,dplyr,data.table,stringr,purrr,furrr,COUNT,optimParallel,
               RhpcBLASctl,purrr,doParallel,foreach)

#cell_depth=data.frame(cell_type,depth=summaryfile$depth)[,-2]
#colnames(cell_depth)=c("cell","depth")
#mean_depth=mean(cell_depth[,2])

split_vector_into_chunks <- function(vector, n) {
  total_sum <- sum(vector)
  target_sum <- total_sum / n
  
  chunks <- list()
  current_chunk <- vector[1]
  chunk_index <- 1
  
  for (i in 2:length(vector)) {
    if (sum(current_chunk) + vector[i] <= target_sum || chunk_index == n) {
      current_chunk <- c(current_chunk, vector[i])
    } else {
      chunks[[chunk_index]] <- current_chunk
      current_chunk <- vector[i]
      chunk_index <- chunk_index + 1
    }
  }
  
  # Add the last chunk
  chunks[[chunk_index]] <- current_chunk
  
  # If there are still remaining elements, create additional empty chunks
  while (length(chunks) < n) {
    chunks <- c(chunks, list(vector[0]))
  }
  
  return(chunks)
}



stat_gen_pois_corr_mle=function(r,a_dat,p,lower,upper,binsize=10000,S,epsilon,k){
  
  
  
  
  fstloc=as.numeric(word(word(p,1,sep="-"),2,sep="_"))
  sndloc=as.numeric(word(word(p,2,sep="-"),2,sep="_"))
  
  
  diag=(sndloc-fstloc)/binsize
  
  
  
  #diag=(as.numeric(word(word(p,2,sep="-"),2,sep="_"))-as.numeric(word(word(p,1,sep="-"),2,sep="_")))/binsize
  
  if(!(all(diag<=upper&diag>=lower)&length(unique(r$locipair))==choose(S,2))){return(NA)}        
  if(all(diag<=upper&diag>=lower)&length(unique(r$locipair))==choose(S,2)){     
    
    tmp=data.table(d=diag,p=p)
    
    
    zerotmp=data.table(
      d=tmp$d,
      cc=0,
      cell=rep(cell_depth_dat$cell,each=choose(S,2)),#It was S.. Oct 15.
      depth=rep(cell_depth_dat$depth,each=choose(S,2)),
      locipair=tmp$p
    )
    
    tmpmat<- r[,c("cell","locipair","cc")][zerotmp, on = c("cell","locipair")]    
    
    
    
    
    
    tmpmat$cc[is.na(tmpmat$cc)]=0
    
    
    
    
    tmpmat=tmpdat[,c("cell","d","phat")][tmpmat,on=c("cell","d")]
    
    
    
    
    tmpmat=a_dat[tmpmat,on="locipair"]#This made an error. saving A_dat was a problem . It gave one locuspair
    
    tmpmat$indA[is.na(tmpmat$indA)]=0# Not existing in other cell types but exists in this cell type.
    
    
    tmpmat=alphadat[tmpmat,on="d"] #%>% na.omit #This made an error. saving A_dat was a problem . It gave one locuspair
    
    
    tmpmat=tmpmat[,c('alpha','alpha0') := .(alphaj*indA+epsilon*(1-indA),unique(N_A)*alphaj+unique(N_AC)*epsilon),by=locipair] 
    
    # tmpcovmat=tmpmat[,"halfabscov":=.(sqrt(depth/(1+alpha0)*(alpha0+depth))*phat*(alpha/alpha0))]
    # tmpcovmat=tmpcovmat %>% group_by(cell) %>% reframe(cv=-sum(combn(halfabscov,m=2,FUN=prod),na.rm=T))
    # sumcov=sum(tmpcovmat$cv)
    # 
    tmpmat=tmpmat[,c("E","V") := .(depth*phat*alpha/alpha0,depth*phat*alpha/alpha0*(1-phat*(1+alpha)/(1+alpha0)
                                                                                    +depth*phat*1/(1+alpha0)*(1-alpha/alpha0)))] #%>% group_by(locipair) %>% mutate(alpha=alphaj*indA,alpha0=unique(N_A)*alphaj)    
    
    tmpmat=(tmpmat %>% group_by(locipair) %>% summarise(Zstat=sum(cc-E)/sqrt(sum(V)))) %>% mutate(pvals=pnorm(Zstat,lower.tail=FALSE),ord=rank(pvals))
    
    
    part_loc=unique(c(fstloc[tmpmat$ord<(k+1)],sndloc[tmpmat$ord<(k+1)]))
    n_top_contacting=length(part_loc)
    top_locipair=paste0(chrlist[chr],'_',fstloc[tmpmat$ord<(k+1)],'-',sndloc[tmpmat$ord<(k+1)])
    
    
    
    
    
    
    #return(list(pvals=tmpmat$pvals,n_top_contacting=n_top_contacting))
    return(list(pvals=tmpmat$pvals,Zvals=tmpmat$Zstat,n_top_contacting=n_top_contacting,name_top=paste(paste0(chrlist[chr],"_",part_loc),collapse = "-"),top_locipair=top_locipair))
  }        
  
  
  
  
}




trim_outliers_iqr <- function(x, outlier_factor = 1.5) {
  # Calculate the first quartile (Q1) and third quartile (Q3)
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  
  # Calculate the interquartile range (IQR)
  iqr <- q3 - q1
  
  # Define the lower and upper bounds for trimming outliers
  lower_bound <- q1 - outlier_factor * iqr
  upper_bound <- q3 + outlier_factor * iqr
  
  # Trim the outliers
  x_trimmed <- x[x >= lower_bound & x <= upper_bound]
  
  return(x_trimmed)
}







nlogL=function(alpha,dat,epsilon){
  
  alpha_aug=rep(alpha,dat[,.(Nnz=mean(.N)),by="d"]$Nnz)          
  celln=mean(dat$celln)# %>% unique
  
  
  dat=dat[,c("Nalpha","NCepsilon","alpha"):=.(N_A*alpha_aug,N_AC*epsilon,alpha_aug)]    
  
  
  
  d_dat=dat[,.(Nalpha=mean(Nalpha),indloggammaalpha=sum(indA*lgamma(alpha+cc)),Zloggamma=mean(C_NzA*lgamma(alpha)),Nloggamma=mean(N_A*lgamma(alpha))),by="d"]
  
  res=celln*sum(lgamma(d_dat$Nalpha+d_dat$NCepsilon))+
    sum(d_dat$indloggammaalpha)+ #already summed over k
    sum(d_dat$Zloggamma)-#already summed over k
    celln*sum(d_dat$Nloggamma)-
    sum(lgamma((dat[,.(Ck=mean(depth)),by="cell"] )$Ck+sum(d_dat$Nalpha)+sum(d_dat$NCepsilon)))#already summed over k
  
  
  return(-res)
  
}


#library(data.table)  
gr_nlogL=function(alpha,dat,epsilon){
  
  
  alpha_aug=rep(alpha,dat[,.(Nnz=mean(.N)),by="d"]$Nnz)          
  celln=mean(dat$celln)# %>% unique
  
  
  dat=dat[,c("Nalpha","NCepsilon","alpha"):=.(N_A*alpha_aug,N_AC*epsilon,alpha_aug)]    
  
  
  
  d_dat=dat[,.(N_A=mean(N_A),Nalpha=mean(Nalpha),NCepsilon=mean(NCepsilon),inddigammaalpha=sum(indA*digamma(alpha+cc)),Zdigamma=mean(C_NzA*digamma(alpha)),Ndigamma=mean(N_A*digamma(alpha))),by="d"]
  
  grad=celln*d_dat$N_A*digamma(d_dat$Nalpha+d_dat$NCepsilon)+
    d_dat$inddigammaalpha+ #already summed over k
    d_dat$Zdigamma-#already summed over k
    celln*d_dat$Ndigamma-
    d_dat$N_A*sum(digamma((dat[,.(Ck=mean(depth)),by="cell"] )$Ck+sum(d_dat$Nalpha)+sum(d_dat$NCepsilon)))#already summed over k
  
  
  
  return(-grad)
  
}








################

################


P_value_generate=function(
    L=20,
    lower=1,
    tau=0.005,
    binsize=500000,
    data_dir='/afs/cs.wisc.edu/p/keles/Collab_2022/VolumeA/Manybody/data/Ramani2017/tmp/',
    output_dir='/afs/cs.wisc.edu/p/keles/Collab_2022/VolumeA/Manybody/results/Ramani2017/tmp/',
    Smax=7,Smin=3,
    celltypenum=4,
    corenum=4,
    core_chr=4,
    chrnum=23,
    sizefile='hg19.chrom.sizes',
    chrlevel.exist=TRUE,upper=NULL,epsilon=0.00000001,automate_upper=TRUE,k_default=TRUE,k_lists=NULL){
  
  
  options(scipen = 100, digits = 4)
  setwd(data_dir)
  cell_type=qs::qread("cell_type.qs")
  colnames(cell_type)=c("cell","cluster")
  chrlist=c(paste0("chr",c(1:(chrnum-1),"X")))
  
  
  
  setwd(data_dir)
  if(chrlevel.exist==FALSE){
    for(chr in 1:chrnum){
      
      cell_type=qs::qread("cell_type.qs")
      ctlists=unique(cell_type[,2])
      hic_df=lapply(ctlists,function(x)qs::qread(paste0(chrlist[chr],'/hic_df_',chrlist[chr],'_',x,'.qs')) )
      hic_df=do.call('rbind',hic_df)
      qs::qsave(hic_df,paste0('hic_df_',chrlist[chr],'.qs'))    
    }
    
    
  }
  
  if(automate_upper==TRUE){
    hic_df=qs::qread(paste0('hic_df_',chrlist[chrnum-1],'.qs')) 
    hic_df=hic_df[,c('binA','binB','d') := .(binA/binsize,binB/binsize,binB/binsize-binA/binsize)]#%>%mutate(binA=binA/binsize,binB=binB/binsize,d=binB-binA) 
    upper=(hic_df[,.(zeropoint=min(which(d %>% unique %>% sort  %>% diff!=1))+1),by="cluster"]$zeropoint %>% min)%/%10*10%>%suppressWarnings()
    if(!is.finite(upper)){upper=max(hic_df$d)}
    
  }
  
  
  system(paste0('cd ', output_dir,'; mkdir ',"L",L,"_lower",lower,"_upper",upper))
  S_ind=0
  for(S in Smin:Smax){
    S_ind=S_ind+1  
    if(k_default){k=S-1}
    if(!k_default){k=k_lists[S_ind]}
    
    pairwise_list=qs::qread(paste0(output_dir,'/pairwise_',S,'_filtered_list.qs'))
    Q_list=qs::qread(paste0(output_dir,'/Q',S,'_filtered_list.qs'))
    
    
    size<-fread(sizefile)
    
    
    system(paste0('cd ', output_dir,'; mkdir ',"L",L,"_lower",lower,"_upper",upper,'/size',S))
    
    
    
    
    
    
    
    
    #  for(chr in 1:chrnum){
    #cl_chr <- makeCluster(core_chr)
    cl_chr <-parallelly::makeClusterPSOCK(core_chr)
    registerDoParallel(cl_chr)     
    foreach(chr=c(1:chrnum),.packages=(.packages()), .export=c('nlogL','gr_nlogL','stat_gen_pois_corr_mle'))%dopar%{      
      system(paste0('cd ', output_dir,"/L",L,"_lower",lower,"_upper",upper,'/size',S,'; mkdir ',chrlist[chr]))        
      pairwise_3=pairwise_list[[chr]]
      chrsize=size[V1==chrlist[chr]]$V2%/%binsize
      
      
      
      withinrange=lapply(pairwise_3,function(p){d=(as.numeric(word(word(p,2,sep="-"),2,sep="_"))-as.numeric(word(word(p,1,sep="-"),2,sep="_")))/binsize
      ;return(all(d<=upper&d>=lower))
      
      }
      
      ) %>% unlist    
      
      
      if(mean(withinrange)==0|is.null(withinrange)){
        
        
        for(ct in unique(cell_type[,2])){
          
          saveRDS(NULL,paste0(output_dir,"/L",L,"_lower",lower,"_upper",upper,'/size',S,"/",chrlist[chr],"/","porder_",ct,"_DirMult_bandlevel_mle_Aset_filtered.qs"))                
        }
        
      }  
      
      if(!(mean(withinrange)==0|is.null(withinrange))){
        #pois, Multi
        
        
        options(future.fork.multithreading.enable = FALSE)
        
        
        hic_df=qs::qread(paste0('hic_df_',chrlist[chr],'.qs'))    
        hic_df=hic_df[,c('binA','binB','d') := .(binA/binsize,binB/binsize,binB/binsize-binA/binsize)]
        hic_df=hic_df[d<=upper&d>=lower]
        celln_table=cell_type %>% as.data.frame %>% group_by(cluster) %>% summarise(celln=n()) %>% as.data.table
        hic_df=celln_table[hic_df,on="cluster"] %>% na.omit#merge(hic_df,celln_table,by="cluster")
        
        hic_df=hic_df[,c('prob_nonzero') := .N/celln,by=.(cluster,locipair)]
        
        
        
        hic_df= hic_df[,c('depth') := sum(cc),by=cell]
        hic_df=as.data.table(hic_df)
        
        
        hic_df=hic_df[,c('N') := .( chrsize-d)]
        
        
        
        options(future.globals.maxSize= Inf)    
        future::plan(multicore, workers = celltypenum)    
        future_map(unique(cell_type[,2]),function(x){ 
          celln=as.numeric(table(cell_type[,2])[x])
          A_dat=hic_df[,c('cluster','locipair','d','prob_nonzero','N')]%>% distinct
          A_dat=A_dat[cluster!=x][,.(pA=mean(prob_nonzero>tau),d=unique(d),N=unique(N)),by=locipair] 
          A_dat=A_dat[,'indA' := .(as.numeric(pA!=0))]
          A_dat=A_dat[,-c('pA')]
          A_dat=A_dat[,N_A := .(sum(indA)),by=d] 
          A_dat=A_dat[,N_AC := N-N_A]
          tmp=hic_df[cluster==x]
          C_NnzAtmp=tmp[,c('locipair','d')][,.(C_NnzA=sum(locipair%in%A_dat[indA==1]$locipair)),by=d] %>% arrange(d)
          A_dat=A_dat[C_NnzAtmp,on="d"][,"C_NzA":=celln*N_A-C_NnzA]
          qs::qsave(A_dat,paste0(output_dir,"/L",L,"_lower",lower,"_upper",upper,'/size',S,'/',chrlist[chr],"/A_dat_",x,'_filtered.qs'))    
        })
        
        
        
        
        
        
        
        hic_df=0
        rm(hic_df)    
        gc()    
        
        # clust <- makeCluster(celltypenum)
        # registerDoParallel(clust)
        # foreach(ct = unique(cell_type[,2]),.packages=(.packages()), .export=c('nlogL','gr_nlogL','stat_gen_pois_corr_mle')) %dopar% {        
        future::plan(multicore, workers = celltypenum)    
        options(future.globals.maxSize= Inf)  
        future_map(unique(cell_type[,2]),function(ct){
          
          #          print(nlogL)
          
          
          
          
          hic_df=qs::qread(paste0(chrlist[chr],'/hic_df_',chrlist[chr],'_',ct,'.qs'))    
          
          
          
          hic_df=hic_df[,c('binA','binB','d') := .(binA/binsize,binB/binsize,binB/binsize-binA/binsize)]
          hic_df=hic_df[d<=upper&d>=lower]
          celln_table=cell_type %>% as.data.frame %>% group_by(cluster) %>% summarise(celln=n()) %>% as.data.table
          hic_df=celln_table[hic_df,on="cluster"] %>% na.omit
          
          hic_df=hic_df[,c('prob_nonzero') := .N/celln,by=.(cluster,locipair)]
          
          
          hic_df= hic_df[,c('depth') := sum(cc),by=cell]
          hic_df=as.data.table(hic_df)
          
          
          hic_df=hic_df[,c('N') := .( chrsize-d)]
          
          
          
          
          #summarizing at the band level
          fulldat=hic_df
          
          #above cell does this.
          
          fulldat=fulldat[,c('depth') := sum(cc),by=cell]# %>% group_by(cell) %>% mutate(depth=sum(cc))
          
          
          fulldat=fulldat[,.(bandS=sum(cc)),by=.(d,cell,cluster)] #%>% group_by(d,cell,cluster) %>% summarise(bandS=sum(cc))
          
          fulldat =fulldat[,c('depth') := sum(bandS),by=cell]#fulldat  %>% group_by(cell) %>% mutate(depth=sum(bandS))
          fulldat=as.data.table(fulldat)
          
          
          wholed=(1:chrsize)
          
          
          zerobinlist=lapply(unique(fulldat$cell),function(x)data.table(d=which(!wholed%in%fulldat[cell==x]$d),cell=x,bandS=0,depth=unique(fulldat[cell==x]$depth),cluster=unique(fulldat[cell==x]$cluster)))
          zerobindat=do.call("rbind",zerobinlist)                   
          
          zerobindat=zerobindat[,c("d","cell","cluster","bandS","depth")]
          
          fulldat=rbind(fulldat,zerobindat)
          
          fulldat =fulldat[,c('N') := .( chrsize-d)]%>%arrange(cell,d)  
          
          
          tmpdat=fulldat[(d<=upper&d>=lower&cluster==ct)]
          
          model_fixed_pois <- glm(data=tmpdat,bandS ~ ns(d,L), family="poisson",offset=c(log(depth)))    
          tmpdat=tmpdat[,"muhat" := .(exp(predict(model_fixed_pois,tmpdat)))]
          tmpdat=tmpdat[,"phat" := muhat/depth,by=cell] 
          qs::qsave(tmpdat,paste0(output_dir,"/L",L,"_lower",lower,"_upper",upper,'/size',S,'/',chrlist[chr],"/phat_",ct,'.qs'))
          #tmpdat contains bandlevel prediction even for zero count bands.  
          
          
          
          rawdat=hic_df[cluster%in%c(ct)]
          #rawdat contains raw hic_df for nonzeros
          
          
          
          rawdat=tmpdat[,c("cell","d","bandS")][rawdat,on=c('cell','d')]%>%na.omit
          
          
          
          
          
          
          #rawdat=as.data.table(rawdat)
          cell_depth_dat=distinct(rawdat[cell%in%unique(rawdat$cell)][,c("cell","depth")])
          rawdat=rawdat[,c("cell","locipair","cc")] #%>% dplyr::select(c("cell","locipair","cc"))
          
          setDT(rawdat)
          
          
          
          options(future.globals.maxSize= Inf)  
          future::plan(multicore, workers = corenum)                     
          
          
          tmp <- future_map(pairwise_3[withinrange], function(groups) {                   
            r=rawdat[locipair %in% groups]
            if(length(unique(r$locipair))!=choose(S,2)){return(NA)}
            if(length(unique(r$locipair))==choose(S,2)){return(r)}
          })
          
          rawdat=tmp
          
          #          rm(tmp)                   
          #          gc()
          
          setDT(tmpdat)
          
          nd=length(lower:upper)#length(unique(rawdat$d))
          
          
          
          
          ###A_dat result come in: 
          #A_dat contains summaries for each locipair. Whether or not it is contained in Aset.             
          #N_A per d summarizes how many of the loci pairs have indA==1
          
          
          A_dat=qs::qread(paste0(output_dir,"/L",L,"_lower",lower,"_upper",upper,'/size',S,'/',chrlist[chr],"/A_dat_",ct,'_filtered.qs'))
          
          tmp=A_dat[hic_df,on=c("locipair","d","N")]%>% na.omit
          
          #tmp =tmp[,NnzAc := .(N_A - NnzA) ,by=d] 
          
          tmp=tmp %>% arrange(d)
          celln=unique(tmp$celln)
          
          A_dat=A_dat[,-c("N","d")] 
          
          
          setDT(A_dat)
          
          
          
          
          
          
          
          future::plan(multicore, workers = corenum)                   
          tmp2 <- future_map(pairwise_3[withinrange], function(groups) {                   
            A_dat[locipair %in% groups] 
          })
          
          
          A_dat=tmp2                   
          
          #          rm(tmp2)                   
          #          gc()                   
          
          #nzeroj_ind=which(arrange(tmp[,.(s=sum(N_A)),by="d"],d)$s!=0)
          tt=arrange(tmp[,.(s=sum(N_A)),by="d"],d)
          nzeroj_ind=(lower:upper)[(lower:upper)%in%tt[(tt$s!=0),]$d]
          
          # cl <- makeCluster(corenum)     
          
          cl <-parallelly::makeClusterPSOCK(corenum)
          
          
          clusterEvalQ(cl, library("data.table")) %>% invisible
          # upper = 1,
          setDefaultCluster(cl=cl)
          set.seed(1)  
          init=sort(runif(length(lower:upper),min=0,max=0.01),decreasing=TRUE)#rep(0.5,length(lower:upper))
          opt_res=optimParallel(par = init[nzeroj_ind], fn = nlogL,gr = gr_nlogL, dat = tmp[,c("C_NzA","N_A","N_AC","d","cell","depth","celln","cc","indA")][d%in%nzeroj_ind],lambda=0, method = "L-BFGS-B",
                                lower = rep(1e-15,length(nzeroj_ind)),epsilon=epsilon,
                                control = list(trace = 3,maxit = 500,
                                               #                                               maxit = 10, factr = 1e12,lmm=2),
                                               lmm = 2),
                                parallel=list(loginfo=TRUE))
          
          
          stopCluster(cl)  
          
          
          
          
          
          
          param=rep(1e-15,length(lower:upper))
          param[which(tmp[,.(s=sum(N_A)),by="d"]$s!=0)]=opt_res$par
          alphadat=data.table(d=lower:upper,alphaj=param)    
          qs::qsave(alphadat,paste0(output_dir,"/L",L,"_lower",lower,"_upper",upper,'/size',S,'/',chrlist[chr],"/alphahat_",ct,'.qs'))
          #alphadat=data.table(d=lower:upper,alphaj=opt_res$par)    
          rm(opt_res)
          gc()                            
          
          future::plan(multicore, workers = corenum)    
          
          
          options(future.globals.maxSize= Inf)  
          start_time <- Sys.time()
          res=future_pmap(list(rawdat[!is.na(rawdat)],A_dat[!is.na(rawdat)],pairwise_3[withinrange][!is.na(rawdat)]),stat_gen_pois_corr_mle,lower=lower,upper=upper,binsize=binsize,S=S,epsilon=epsilon,k=k)
          end_time <- Sys.time()
          end_time - start_time
          
          
          
          
          
          
          
          #zscores=res %>% unlist
          porder=lapply(lapply(res,function(x)x$pvals),function(x)sort(x)[k]) %>% unlist
          Zvals=lapply(lapply(res,function(x)x$Zvals),function(x)sort(x,decreasing = TRUE)[k]) %>% unlist
          names(porder)=Q_list[[chr]][withinrange][!is.na(rawdat)]#pairwise_3[withinrange]       
          
          porder=porder[which(lapply(res,function(x)x$n_top_contacting>=S) %>% unlist)]
          Zvals=Zvals[which(lapply(res,function(x)x$n_top_contacting>=S) %>% unlist)]
          top_locipairs=lapply(res,function(x)x$top_locipair)[which(lapply(res,function(x)x$n_top_contacting>=S) %>% unlist)] 
          qs::qsave(Zvals,paste0(output_dir,"/L",L,"_lower",lower,"_upper",upper,'/size',S,"/",chrlist[chr],"/","Zorder_",ct,"_DirMult_bandlevel_mle_Aset_filtered.qs"))
          qs::qsave(porder,paste0(output_dir,"/L",L,"_lower",lower,"_upper",upper,'/size',S,"/",chrlist[chr],"/","porder_",ct,"_DirMult_bandlevel_mle_Aset_filtered.qs"))
          qs::qsave(top_locipairs,paste0(output_dir,"/L",L,"_lower",lower,"_upper",upper,'/size',S,"/",chrlist[chr],"/","top_locipairs_",ct,"_DirMult_bandlevel_mle_Aset_filtered.qs"))
          rm(porder)
          gc()                     
          print(ct)                   
          
        }
        #stopCluster(cl_chr)  
        )
        
      }
      
      
      
      
      
      
    }
    
    stopCluster(cl_chr)  
    
    
    
    
    
    
  }
  
  
  
  
  
  
  
  
  
}
