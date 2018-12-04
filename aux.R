load_libs<-function(pg_choice,Rlib){
  if(pg_choice=="RaceID3"){library(RaceID,lib.loc=Rlib)
    library(ggplot2,lib.loc=Rlib) } else if (pg_choice == "Monocle2"){
      library(monocle,lib.loc=Rlib)
    } else if (pg_choice == "Seurat"){
      library(Seurat,lib.loc=Rlib)}
}

get_cluinit<-function(pg_choice,sc){
   if(pg_choice=="RaceID3"){
    sc@cpart<-sc@cluster$kpart
    cluinit<-max(sc@cluster$kpart)}
   else if (pg_choice == "Monocle2"){
    cluinit<-max(pData(sc)$Cluster)}
  else if (pg_choice == "Seurat"){
    cluinit<-max(sc@ident)}
  return(cluinit)
}

recluster_plot_tsne<-function(pg_choice,sc,numclu){
  sc@cpart<-sc@cluster$kpart
  scnew<-sc
  if(pg_choice=="RaceID3"){
    if(numclu!=max(sc@cluster$kpart)){   
      scnew<-clustexp(sc,rseed=314,FUNcluster="kmedoids",sat=FALSE,cln=numclu)
      scnew<-findoutliers(scnew)
      scnew@cpart<-scnew@cluster$kpart
      }else if (pg_choice == "Monocle2"){
        if(numclu!=max(pData(sc)$Cluster)){
        scnew <- clusterCells(sc, num_clusters=numclu)}
      }else if (pg_choice == "Seurat"){
        if(numclu!=max(sc@ident)){
        scnew<-FindClusters(object =sc)}
      }
  }
  
  return(scnew)
} 

get_clu_plot<-function(pg_choice,sc){
  if(pg_choice=="RaceID3"){
    plotmap(sc,final=FALSE)
  }else if (pg_choice == "Monocle2"){
    plot_cell_clusters(sc,1, 2, color="Cluster")
  }else if (pg_choice == "Seurat"){
    TSNEPlot(object = sc)
  }
}

get_top10<-function(pg_choice,sc){
  if(pg_choice=="RaceID3"){
    res10L<-lapply(unique(sc@cpart),function(X){
        dg<-clustdiffgenes(sc,X,pvalue=.01)
        dgsub<-dg[dg$fc>=2&dg$padj<0.05,]
        if(nrow(dgsub)>0){
          dg<-dgsub
          dg<-head(dg,n=10)
          dg$Cluster<-X
          dg$Gene<-rownames(dg)}else{dg<-NULL}
        return(dg)})
      top10<-as.data.frame(do.call(rbind,res10L))
      #top10<-top10[top10$padj<0.05,]
      top10<-top10[with(top10, order(Cluster, padj)),]
   }else if (pg_choice == "Monocle2"){
     pData(sc)$Cluster<-as.character(pData(sc)$Cluster)
     res10L<-lapply(unique(pData(sc)$Cluster),function(X){
       pData(sc)$Cluster[!pData(sc)$Cluster %in% X]<-"0"
       pData(sc)$Cluster<-factor(pData(sc)$Cluster)
       pData(sc)$Cluster<-relevel(pData(sc)$Cluster,ref="0")#probably not really necessary since 0<any cluster number
       res<-differentialGeneTest(sc,fullModelFormulaStr = "~Cluster",cores=4)
       res.filt<-res[res$qval<0.05,]
       if(nrow(res.filt)>0){
         res<-res.filt
         res<-res[order(res$qval),]
         res$Cluster<-X
       }else{res<-NULL}
       return(res)
     })
     top10<-as.data.frame(do.call(rbind,res10L))
     top10<-top10[with(top10, order(Cluster, qval)),]
   }else if (pg_choice == "Seurat"){
     top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
     colnames(top10)[colnames(top10) %in% "cluster"]<-"Cluster"
   }
  return(top10)
}

get_marker_plot<-function(pg_choice,sc,genes){
  if(pg_choice=="RaceID3"){
    plotmarkergenes(sc,genes)
  }else if (pg_choice == "Monocle2"){
    
  }else if (pg_choice == "Seurat"){
    DoHeatmap(object = sc,genes.use = genes,slim.col.label = TRUE,remove.key = TRUE)
  }
}

render_data_head<-function(pg_choice,sc){
  if(pg_choice=="RaceID3"){
    ntemp<-as.data.frame(as.matrix(sc@ndata)*5000,stringsAsFactors=FALSE)
  }else if (pg_choice == "Monocle2"){
    ntemp<-as.data.frame(t(t(exprs(sc)) /  pData(sc)[, 'Size_Factor']),stringsAsFactors=FALSE)
  }else if (pg_choice == "Seurat"){
    
  }
  return(ntemp)
}

get_feature_plot<-function(pg_choice,sc,nv,nt,tsnelog){
  if(pg_choice=="RaceID3"){
    plotexpmap(sc,nv,n=nt,logsc=tsnelog)
  }else if (pg_choice == "Monocle2"){
    
  }else if (pg_choice == "Seurat"){
    
  }
  
}