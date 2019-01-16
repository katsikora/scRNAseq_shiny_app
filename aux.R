load_libs<-function(pg_choice,Rlib){
  if(pg_choice=="RaceID3"){
    library(Matrix,lib.loc=Rlib)
    library(RaceID,lib.loc=Rlib)
    library(ggplot2,lib.loc=Rlib) } else if (pg_choice == "Monocle2/Seurat"){
      library(VGAM, lib.loc=Rlib)
      library(irlba, lib.loc=Rlib)
      library(DDRTree, lib.loc=Rlib)
      library(monocle,lib.loc=Rlib)
      library(scales,lib.loc=Rlib)
      library(Seurat,lib.loc=Rlib)
    } 
}

get_cluinit<-function(pg_choice,sc){
   if(pg_choice=="RaceID3"){
    sc@cpart<-sc@cluster$kpart
    cluinit<-max(sc@cluster$kpart)}
   else if (pg_choice == "Monocle2/Seurat"){
    cluinit<-max(as.numeric(pData(sc)$Cluster))}
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
    }}else if (pg_choice == "Monocle2/Seurat"){
       scnew<-sc
        if(numclu!=max(as.numeric(pData(sc)$Cluster))){
        scnew <- clusterCells(sc, num_clusters=numclu)}
    }
  return(scnew)
} 

get_clu_plot<-function(pg_choice,sc){
  if(pg_choice=="RaceID3"){
    plotmap(sc,final=FALSE)
  }else if (pg_choice == "Monocle2/Seurat"){
    plot_cell_clusters(sc,1, 2, color="Cluster")
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
   }else if (pg_choice == "Monocle2/Seurat"){
     #convert seurat object from monocle... 
     seuset <- Seurat::CreateSeuratObject(counts=exprs(sc), 
                                              min.cells=4,
                                              min.features=0,
                                              project = "exportCDS",
                                              meta.data = pData(sc),
                                              assay="RNA")
     
     seuset@meta.data <- pData(sc)
     seuset <- NormalizeData(object = seuset)
     seuset <- FindVariableFeatures(object = seuset)
     seuset <- ScaleData(object = seuset)
     seuset <- RunPCA(object = seuset)
     seuset <- FindNeighbors(object = seuset)
     seuset <- FindClusters(object = seuset)
     seuset<-SetIdent(seuset, value = pData(sc)$Cluster)
     markers<-FindAllMarkers(object = seuset,only.pos = TRUE,min.pct = 0.25,thresh.use = 0.25)
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
  }else if (pg_choice == "Monocle2/Seurat"){
    seuset <- Seurat::CreateSeuratObject(counts=exprs(sc), 
                                         min.cells=4,
                                         min.features=0,
                                         project = "exportCDS",
                                         meta.data = pData(sc),
                                         assay="RNA")
    
    seuset@meta.data <- pData(sc)
    seuset <- NormalizeData(object = seuset)
    seuset <- FindVariableFeatures(object = seuset)
    seuset <- ScaleData(object = seuset)
    seuset <- RunPCA(object = seuset)
    seuset <- FindNeighbors(object = seuset)
    seuset <- FindClusters(object = seuset)
    #seuset <- RunTSNE(object = seuset)
    seuset<-SetIdent(seuset, value = pData(sc)$Cluster)
    DoHeatmap(object = seuset,features=genes)
  }
}

render_data_head<-function(pg_choice,sc){
  if(pg_choice=="RaceID3"){
    ntemp<-as.data.frame(as.matrix(sc@ndata)*5000,stringsAsFactors=FALSE)
  }else if (pg_choice == "Monocle2/Seurat"){
    ntemp<-as.data.frame(t(t(exprs(sc)) /  pData(sc)[, 'Size_Factor']),stringsAsFactors=FALSE)
  }
  return(ntemp)
}

get_feature_plot<-function(pg_choice,sc,nv,nt,tsnelog){
  if(pg_choice=="RaceID3"){
    plotexpmap(sc,nv,n=nt,logsc=tsnelog)
  }else if (pg_choice == "Monocle2/Seurat"){
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
  }
}
