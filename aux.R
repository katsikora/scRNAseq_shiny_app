load_libs<-function(pg_choice,Rlib){
  if(pg_choice=="RaceID3"){
    library(Matrix,lib.loc=Rlib)
    library(RaceID,lib.loc=Rlib)
    library(ggplot2,lib.loc=Rlib) } else if (pg_choice == "Monocle2"){
      library(VGAM, lib.loc=Rlib)
      library(irlba, lib.loc=Rlib)
      library(DDRTree, lib.loc=Rlib)
      library(monocle,lib.loc=Rlib)
      library(scales,lib.loc=Rlib)
    } else if (pg_choice == "Seurat"){
      library(Seurat,lib.loc=Rlib)
      } #,lib.loc=Rlib
}

get_cluinit<-function(pg_choice,sc){
   if(pg_choice=="RaceID3"){
    sc@cpart<-sc@cluster$kpart
    cluinit<-max(sc@cluster$kpart)}
   else if (pg_choice == "Monocle2"){
    cluinit<-max(as.numeric(pData(sc)$Cluster))}
  else if (pg_choice == "Seurat"){
    cluinit<-max(as.numeric(as.character(sc@ident)))}
  return(cluinit)
}

recluster_plot_tsne<-function(pg_choice,sc,numclu){
  if(pg_choice=="RaceID3"){
    sc@cpart<-sc@cluster$kpart
    scnew<-sc
    if(numclu!=max(sc@cluster$kpart)){   
      scnew<-clustexp(sc,rseed=314,FUNcluster="kmedoids",sat=FALSE,cln=numclu)
      scnew<-findoutliers(scnew)
      scnew@cpart<-scnew@cluster$kpart
    }}else if (pg_choice == "Monocle2"){
       scnew<-sc
        if(numclu!=max(as.numeric(pData(sc)$Cluster))){
        scnew <- clusterCells(sc, num_clusters=numclu)}
    }else if (pg_choice == "Seurat"){
       scnew<-sc
        if(numclu!=max(as.numeric(as.character(sc@ident)))){
          scnew<-FindNeighbors(object = UpdateSeuratObject(sc))
          scnew<-FindClusters(object =scnew)}
      }
  return(scnew)
} 

get_clu_plot<-function(pg_choice,sc){
  if(pg_choice=="RaceID3"){
    plotmap(sc,final=FALSE)
  }else if (pg_choice == "Monocle2"){
    plot_cell_clusters(sc,1, 2, color="Cluster")
  }else if (pg_choice == "Seurat"){
    if(!class(sc) %in% "Seurat"){
    DimPlot(object = UpdateSeuratObject(sc),reduction="tsne")
  }else{DimPlot(object = sc,reduction="tsne")}}
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
         res$Gene<-rownames(res)
         res<-subset(res,select=c("Gene","qval","num_cells_expressed","Cluster"))
       }else{res<-NULL}
       return(res)
     })
     top10<-as.data.frame(do.call(rbind,res10L))
     top10<-top10[with(top10, order(Cluster, qval)),]
   }else if (pg_choice == "Seurat"){
     markers<-FindAllMarkers(object = UpdateSeuratObject(sc),only.pos = TRUE,min.pct = 0.25,thresh.use = 0.25)
     top10 <- markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(10, avg_logFC)
     colnames(top10)[colnames(top10) %in% "cluster"]<-"Cluster"
     colnames(top10)[colnames(top10) %in% "gene"]<-"Gene"
     top10<-as.data.frame(top10,stringsAsFactors=FALSE)
   }
  return(top10)
}

get_marker_plot<-function(pg_choice,sc,topn){
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  set.seed(123)
  genes<-unique(topn$Gene)
  if(pg_choice=="RaceID3"){
    #plotmarkergenes(sc,genes)
    plotdat<-as.data.frame(as.matrix(sc@ndata[rowSums(as.matrix(sc@ndata))>0,]),stringsAsFactors=FALSE)
    plotdat<-subset(plotdat,subset=rownames(plotdat) %in% genes)
    plotdat2<-as.matrix(log2(plotdat+0.01))
    #rownames(plotdat2)<-rownames(plotdat)
    plotdat2<-plotdat2[match(genes,rownames(plotdat2)),order(as.numeric(sc@cpart))]
    colv<-sample(col_vector,max(as.numeric(sc@cpart)))[sort(as.numeric(sc@cpart))]
    heatmap.2(plotdat2, scale="column", trace="none", dendrogram="none",
              col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255),labCol="",ColSideColors=colv,Colv=FALSE,Rowv=FALSE,
              main="Gene Selection",margins=c(10,12))
  }else if (pg_choice == "Monocle2"){
    plotdat<-as.data.frame(t(t(exprs(sc)) /  pData(sc)[, 'Size_Factor']),stringsAsFactors=FALSE)
    plotdat<-subset(plotdat,subset=rownames(plotdat) %in% genes)
    plotdat2<-as.matrix(log2(plotdat+0.01))
    #rownames(plotdat2)<-rownames(plotdat)
    plotdat2<-plotdat2[match(genes,rownames(plotdat2)),order(as.numeric(pData(sc)$Cluster))]
    colv<-sample(col_vector,max(as.numeric(pData(sc)$Cluster)))[sort(as.numeric(pData(sc)$Cluster))]
    heatmap.2(plotdat2, scale="column", trace="none", dendrogram="none",
              col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255),labCol="",ColSideColors=colv,Colv=FALSE,Rowv=FALSE,
              main="Gene Selection",margins=c(10,12))  
  }else if (pg_choice == "Seurat"){
    DoHeatmap(object = UpdateSeuratObject(sc),features=genes)
  }
}

render_data_head<-function(pg_choice,sc){
  if(pg_choice=="RaceID3"){
    ntemp<-as.data.frame(as.matrix(sc@ndata)*5000,stringsAsFactors=FALSE)
  }else if (pg_choice == "Monocle2"){
    ntemp<-as.data.frame(t(t(exprs(sc)) /  pData(sc)[, 'Size_Factor']),stringsAsFactors=FALSE)
  }else if (pg_choice == "Seurat"){
   ntemp<-as.data.frame(expm1(sc@data)) 
  }
  return(ntemp)
}

get_feature_plot<-function(pg_choice,sc,nv,nt,tsnelog){
  if(pg_choice=="RaceID3"){
    plotexpmap(sc,nv,n=nt,logsc=tsnelog)
  }else if (pg_choice == "Monocle2"){
    plotdat<-as.data.frame(t(sc@reducedDimA),stringsAsFactors=FALSE)
    ndata<-as.data.frame(t(t(exprs(sc)) /  pData(sc)[, 'Size_Factor']),stringsAsFactors=FALSE)
    l<-apply(ndata[nv,]-.1,2,sum)+.1
    if (tsnelog) {
      f <- l == 0
      l <- log2(l)
      l[f] <- NA
    }
    plotdat$label<-l
    ggplot(plotdat %>% dplyr::arrange(label), aes(x = V1, y = V2, color = label))+geom_point(size = 2)+
      scale_colour_continuous(low = "steelblue3", high ="darkorange", space = "Lab", na.value = "grey50",                                                               guide = "colourbar",name=ifelse(tsnelog==FALSE,"Counts","Log2Counts"))+xlab("Dim1")+ylab("Dim2")+theme(axis.text=element_text(size=14),axis.title=element_text(size=16),strip.text=element_text(size=14))+ggtitle(nt)
  }else if (pg_choice == "Seurat"){
    if(length(nv)==1){
      FeaturePlot(UpdateSeuratObject(sc),nv,cols = c("lightgrey", "blue"),ncol = 1) 
    }
  }
  
}
