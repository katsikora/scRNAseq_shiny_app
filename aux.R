#debug_path="/var/log/shiny-server"
#debug_path="/data/manke/sikora/shiny_apps/debug"
set.seed(314)
load_libs<-function(pg_choice,Rlib){
  .libPaths(Rlib)
  if(pg_choice=="RaceID3"){
    library(RaceID,lib.loc=Rlib)
     } else if (pg_choice == "Monocle2"){
      library(monocle,lib.loc=Rlib)
      library(Seurat,lib.loc=Rlib)
           } else if (pg_choice == "Seurat3"){
             library(Seurat,lib.loc=Rlib)
           }
  # 
  library(ggplot2)
  library(gplots)
  library(RColorBrewer)
  library(cluster,lib.loc=Rlib)
}


check_class<-function(pg_choice,sc){
  if(pg_choice=="RaceID3"){
    #check class
    if (class(sc)[1]=="SCseq"){
      res<-TRUE
    }else{
      res<-FALSE
    }
  } else if (pg_choice == "Monocle2"){
    if(class(sc)[1]=="CellDataSet"){
      res<-TRUE
    }else{
      res<-FALSE
      }
    }else if (pg_choice == "Seurat3"){
      if(class(sc)[1]=="Seurat"){
        res<-TRUE
      }else{
        res<-FALSE
      }
  }
  return(res)
}

check_slots<-function(pg_choice,sc){
  if(pg_choice=="RaceID3"){
    #check class
    if (all(isTruthy(sc@ndata),isTruthy(sc@tsne),isTruthy(sc@cluster$kpart),isTruthy(sc@distances),isTruthy(sc@cpart))){
      res<-TRUE
    }else{
      res<-FALSE
    }
  } else if (pg_choice == "Monocle2"){
    if(all(isTruthy(sc@reducedDimA),isTruthy(pData(sc)$Cluster),isTruthy(pData(sc)$rho))){
      res<-TRUE
    }else{
      res<-FALSE
    }
  }else if (pg_choice == "Seurat3"){
    if(all(isTruthy(length(sc@commands)>=7),isTruthy(sc@active.assay %in% "RNA"),isTruthy(startsWith(as.character(sc@version),"3")),isTruthy("FindClusters" %in% names(sc@commands)),isTruthy("RunTSNE" %in% names(sc@commands)))){
      res<-TRUE
    }else{
      res<-FALSE
    }
  }
  return(res)
}


get_cluinit<-function(pg_choice,sc){
   if(pg_choice=="RaceID3"){
    sc@cpart<-sc@cluster$kpart
    cluinit<-max(sc@cluster$kpart)}
   else if (pg_choice == "Monocle2"){
    cluinit<-max(as.numeric(pData(sc)$Cluster))}
   else if (pg_choice == "Seurat3"){
     #return resolution rather than the number of clusters
    #cluinit<-max(as.numeric(sc[[]]$seurat_clusters))}
     cluinit<-sc@commands$FindClusters@params[["resolution"]]
   }
   return(cluinit)
}

#seuset@commands$FindClusters@params[["resolution"]]

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
        if(numclu+1!=max(as.numeric(pData(sc)$Cluster))){
        scnew <- clusterCells(sc, num_clusters=(numclu+1))}
    }else if (pg_choice == "Seurat3"){
      #the relevant parameter is resolution rather than number of clusters
      res<-numclu
      scnew<-FindClusters(sc,resolution=as.numeric(res),random.seed =314)
    }
  return(scnew)
} 

get_clu_plot<-function(pg_choice,sc){
  if(pg_choice=="RaceID3"){
    plotmap(sc,final=FALSE)
  }else if (pg_choice == "Monocle2"){
    plot_cell_clusters(sc,1, 2, color="Cluster")
  }else if (pg_choice == "Seurat3"){
    DimPlot(object = sc, reduction = "tsne")
  }
}

plot_silhouette<-function(pg_choice,sc){
  if(pg_choice=="RaceID3"){
    m <- sc@cluster$kpart
    dp  <- as.dist(sc@distances)
      }else if (pg_choice=="Monocle2"){
    tsne_data<-reducedDimA(sc)
    dp<-as.matrix(dist(t(tsne_data)))
    m<-as.integer(pData(sc)$Cluster)
    names(m)<-colnames(dp)
      }else if (pg_choice=="Seurat3"){
        tsne_data<-Embeddings(sc,reduction = "tsne")
        dp<-as.matrix(dist(tsne_data))
        m<-as.integer(sc[[]]$seurat_clusters)
        names(m)<-colnames(dp)
      }
  si<-silhouette(m,dp)
  plot(si,col=1:max(m), border=NA,main=sprintf("Silhouette plot for %s clusters",max(m)))
}

plot_clu_separation<-function(pg_choice,sc){
  if(pg_choice=="RaceID3"){
    plotsaturation(sc,disp=TRUE)
  }else if (pg_choice=="Monocle2"){
    plot_rho_delta(sc)
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
      top10<-top10[with(top10, order(Cluster, padj)),]
      seuset<-NULL
   }else if (pg_choice == "Monocle2"){
     #convert seurat object from monocle... 
     seuset <- Seurat::CreateSeuratObject(counts=Biobase::exprs(sc), 
                                              min.cells=4,
                                              min.features=0,
                                              project = "exportCDS",
                                              meta.data = pData(sc),
                                              assay="RNA")
     
     seuset@meta.data <- pData(sc)
     seuset <- NormalizeData(object = seuset)
     seuset <- FindVariableFeatures(object = seuset)#,set.var.genes = FALSE
     seuset <- ScaleData(object = seuset,features = NULL)
     seuset <- RunPCA(object = seuset)
     seuset <- FindNeighbors(object = seuset)
     seuset <- FindClusters(object = seuset)
     seuset<-SetIdent(seuset, value = pData(sc)$Cluster)
     markers<-FindAllMarkers(object = seuset,only.pos = TRUE,min.pct = 0.25,thresh.use = 0.25)
     top10 <- markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(10, avg_logFC)
     colnames(top10)[colnames(top10) %in% "cluster"]<-"Cluster"
     colnames(top10)[colnames(top10) %in% "gene"]<-"Gene"
     top10<-as.data.frame(top10,stringsAsFactors=FALSE)
   }else if (pg_choice == "Seurat3"){
     markers<-FindAllMarkers(object = sc,only.pos = TRUE,min.pct = 0.25,thresh.use = 0.25)
     top10 <- markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(10, avg_logFC)
     colnames(top10)[colnames(top10) %in% "cluster"]<-"Cluster"
     colnames(top10)[colnames(top10) %in% "gene"]<-"Gene"
     top10<-as.data.frame(top10,stringsAsFactors=FALSE)
     seuset<-NULL
   }
  return(list(top10,seuset))
}

get_marker_plot<-function(pg_choice,sc,topn,seuset){
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  set.seed(314)
  genes<-unique(topn$Gene)
  #genes<-gsub("--","__",genes)
  if(pg_choice=="RaceID3"){
    plotdat<-as.data.frame(as.matrix(sc@ndata[rowSums(as.matrix(sc@ndata))>0,]),stringsAsFactors=FALSE)
    plotdat<-subset(plotdat,subset=rownames(plotdat) %in% genes)
    plotdat2<-as.matrix(log2(plotdat+0.01))
    plotdat2<-plotdat2[match(genes,rownames(plotdat2)),order(as.numeric(sc@cpart))]
    colv<-sample(col_vector,max(as.numeric(sc@cpart)))[sort(as.numeric(sc@cpart))]
    heatmap.2(plotdat2, scale="column", trace="none", dendrogram="none",
              col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255),labCol="",ColSideColors=colv,Colv=FALSE,Rowv=FALSE,
              main="Gene Selection",margins=c(10,12))
  }else if (pg_choice == "Monocle2"){
    if(!is.null(VariableFeatures(seuset))){
    VariableFeatures(seuset)<-unique(c(VariableFeatures(seuset),genes))
    seuset <- ScaleData(object = seuset)
    DoHeatmap(object = seuset,features=genes)}
  }else if (pg_choice == "Seurat3"){
    VariableFeatures(sc)<-unique(c(VariableFeatures(sc),genes))
    seuset <- ScaleData(object = sc)
    DoHeatmap(object = seuset,features=genes)
  }
}

render_data_head<-function(pg_choice,sc){
  if(pg_choice=="RaceID3"){
    ntemp<-as.data.frame(as.matrix(sc@ndata)*5000,stringsAsFactors=FALSE)
  }else if (pg_choice == "Monocle2"){
    ntemp<-as.data.frame(t(t(Biobase::exprs(sc)) /  pData(sc)[, 'Size_Factor']),stringsAsFactors=FALSE)
  }else if (pg_choice == "Seurat3"){
    ntemp<-as.data.frame(expm1(as.matrix(seuset[["RNA"]]@data)),stringsAsFactors=FALSE)
  }
  return(ntemp)
}

get_feature_plot<-function(pg_choice,sc,nv,nt,tsnelog){
  if(pg_choice=="RaceID3"){
    plotexpmap(sc,nv,n=nt,logsc=tsnelog)
  }else if (pg_choice == "Monocle2"){
    plotdat<-as.data.frame(t(sc@reducedDimA),stringsAsFactors=FALSE)
    ndata<-as.data.frame(t(t(Biobase::exprs(sc)) /  pData(sc)[, 'Size_Factor']),stringsAsFactors=FALSE)
    l<-apply(ndata[nv,]-.1,2,sum)+.1
    if (tsnelog) {
      f <- l == 0
      l <- log2(l)
      l[f] <- NA
    }
    plotdat$label<-l
    ggplot(plotdat %>% dplyr::arrange(label), aes(x = V1, y = V2, color = label))+geom_point(size = 2)+
      scale_colour_continuous(low = "steelblue3", high ="darkorange", space = "Lab", na.value = "grey50",                                                               guide = "colourbar",name=ifelse(tsnelog==FALSE,"Counts","Log2Counts"))+xlab("Dim1")+ylab("Dim2")+theme(axis.text=element_text(size=14),axis.title=element_text(size=16),strip.text=element_text(size=14))+ggtitle(nt)
  } else if (pg_choice == "Seurat3"){
    plotdat<-as.data.frame(Embeddings(sc,reduction = "tsne"),stringsAsFactors=FALSE)
    ndata<-as.data.frame(expm1(as.matrix(seuset[["RNA"]]@data)),stringsAsFactors=FALSE)
    l<-apply(ndata[nv,]-.1,2,sum)+.1
    if (tsnelog) {
      f <- l == 0
      l <- log2(l)
      l[f] <- NA
    }
    plotdat$label<-l
    ggplot(plotdat %>% dplyr::arrange(label), aes(x = tSNE_1, y = tSNE_2, color = label))+geom_point(size = 2)+
      scale_colour_continuous(low = "steelblue3", high ="darkorange", space = "Lab", na.value = "grey50",                                                               guide = "colourbar",name=ifelse(tsnelog==FALSE,"Counts","Log2Counts"))+xlab("Dim1")+ylab("Dim2")+theme(axis.text=element_text(size=14),axis.title=element_text(size=16),strip.text=element_text(size=14))+ggtitle(nt)
  }
}
