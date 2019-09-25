## app.R ##
ver_sion<-"1.0.0"
Rlib="/data/manke/sikora/shiny_apps/Rlibs3.5.0_bioc3.7"
#debug_path="/var/log/shiny-server"
#debug_path="/data/manke/sikora/shiny_apps/debug"
#debug_path="/root/container-logs"
.libPaths(Rlib)
set.seed(314)

options(shiny.maxRequestSize = 1000*1024^2)

library(Rcpp,lib.loc=Rlib)
library(shinydashboard)#,lib.loc=Rlib,verbose=TRUE
library(rhandsontable) #,lib.loc=Rlib,verbose=TRUE
library(DT) #,lib.loc=Rlib,verbose=TRUE


ui <- function(request) {dashboardPage(
    dashboardHeader(title = "Dataset selection"),
    ## Sidebar content
    dashboardSidebar(tagList(

      #selectInput(inputId="genome", label="Select organism", choices=c("PLEASE SELECT A GENOME","Zebrafish [zv10]","Fission yeast","Fruitfly [dm6]","Fruitfly [dm3]","Human [hg37]","Human [hg38]","Mouse [mm9]","Mouse [mm10]"), selected = NULL),#"PLEASE SELECT A GENOME",, selected = NULL
      selectInput(inputId="selectformat",label="Select R package",choices=c("Please select a package","RaceID3","Monocle2","Seurat3"), selected = NULL),
      fileInput('file1', 'Choose file to upload',accept = c('.RData','.RDS')),
      uiOutput("adddataset"),
      #actionButton(inputId="adddataset", label="Submit dataset"),
      tags$footer(list(textOutput("version"),"Copyright 2018 MPI-IE Freiburg Bioinfo Core Unit",imageOutput("logo")),style = "position:absolute;bottom:0")
            )),
        
    dashboardBody(
        h2("Single cell RNAseq analysis"),
        uiOutput("resultPanels")
    )
               
    

 )}


server <- function(input, output, session) {
  
    
    output$walkThrough<-renderUI(HTML("<ul><li>1.Choose  R package used to produce your single cell object. Upload your single cell object in \".RData\" or \".RDS\" format. Click on retrieve dataset. Your data will appear in the InputData tab.</li><li>2.You can visualize the clusters in your dataset as well as change their number in the \"Cell map and clustering\" tab. You can get up to 10 marker genes per cluster and visualize them on a heatmap. </li><li>3.Provide semicolon-separated Gene IDs to calculate aggregate expression for. If your genes are expressed under the filtering criteria, you can visualize their expression on a tsne plot in tab \"Marker Gene Visualization\". </li><li>4.To  list top10 correlated genes as well as to plot pairwise gene expression of genes of interest, go to the tab \"Correlation Analyses\"</li></ul>"))
    output$FAQ<-renderText("Currently, no uniform gene naming system is prerequisite. You have to provide Gene IDs consistent with the naming used to produce your dataset.\n For questions, bug reports or feature requests, contact sikora@ie-freiburg.mpg.de.\n For reporting issues or pull requests on GitHub, go to https://github.com/maxplanck-ie/scRNAseq_shiny_app .")
    
    output$fileDescription<-renderText("Please provide a semicolon-separated list of Gene IDs you would like to obtain results for.")
    output$fileDescription2<-renderText("Please provide a semicolon-separated list of Gene IDs you would like to obtain results for.")
    
    output$logo<-renderImage({list(src="/data/manke/sikora/shiny_apps/userIN_to_yaml/MPIIE_logo_sRGB.jpg",width=100,height=100)},deleteFile =FALSE)
    
    output$version<-renderText(sprintf("App version %s",ver_sion))


################################
    #check if genome has been selected and alert if not the case

    source("/data/manke/sikora/shiny_apps/scRNAseq_shiny_app/aux.R")
    
    values<-reactiveValues()
    values$rowsSel<-""
    values$cList<-"All"
    values$inGenes=""
    
    values$init_checks_passed<-reactive({all(c(input$selectformat!="Please select a package",!is.null(input$selectformat)))})
    output$init_checks_passed<-reactive({values$init_checks_passed()})
    outputOptions(output, "init_checks_passed", suspendWhenHidden = FALSE)
    
    observeEvent(input$selectformat,{
    #imports depend on selected format!
    #load packages in function of the input format (or use namespace loading...)
    load_libs(input$selectformat,Rlib)},ignoreInit=TRUE)#
    output$sessionInfo <- renderPrint({capture.output(sessionInfo())})
    
    output$downloadSessionInfo <- downloadHandler(
      filename = "sessionInfo.txt",
      content = function(con) {
        sink(con)
        print(sessionInfo())
        sink()
      }
    )

    
    

################################
    observe({
        req(isTruthy(input$file1))
        output$adddataset<-renderUI({actionButton(inputId="adddataset", label="Submit dataset")})
        #}
    })
    
    
    observeEvent(input$adddataset,{
      ######################################################################################################
      if(values$init_checks_passed()){

      inFormat<-isolate(input$selectformat)
      
       if (!is.null(input$file1)){values$datpath<-isolate(input$file1)$datapath}
      datPath<-isolate(values$datpath)
      
      if(grepl("rds$",datPath,ignore.case=TRUE)){
            values$sc<-readRDS(datPath)}
      else if (grepl("rdata$",datPath,ignore.case=TRUE)){
           myEnv<-environment()
           sctmp<-load(datPath, envir = myEnv)
           
           values$sc <- myEnv[[sctmp]]
           
      }
     
       sc<-values$sc
       
       class_ok<-check_class(inFormat,sc)
       if(!isTruthy(class_ok)){showModal(modalDialog(title = "DATASET CLASS INCORRECT",
         "Please provide a dataset of the class matching your R package selection!",
         easyClose = TRUE
       ))}
       req(isTruthy(class_ok))
       slots_ok<-check_slots(inFormat,sc)
       if(!isTruthy(slots_ok)){showModal(modalDialog(title = "DATASET DOESN'T HAVE ALL REQUIRED SLOTS POPULATED",
                                                     "Please provide a fully processed dataset containing clustering information and tsne coordinates!",
                                                     easyClose = TRUE
       ))}
       req(isTruthy(slots_ok))
       
       
    ###########################################################################################################   
       cluinit<-get_cluinit(input$selectformat,sc)
       output$CluCtrl<-renderUI({tagList(sliderInput("numclu", ifelse(input$selectformat=="Seurat3","Resolution","Number of clusters"),min=ifelse(input$selectformat=="Seurat3",0,1),max=ifelse(input$selectformat=="Seurat3",1,2*cluinit),value=cluinit,round=TRUE))})
       output$selectdimred<-renderUI({
          if(input$selectformat=="Seurat3"){
            aa<-sc@active.assay
            if(isTruthy(grepl(paste0("RunUMAP.",aa),names(sc@commands)))) {
              tagList(selectInput("selectdimred","Select dimensionality reduction method.",choices=c("tSNE","UMAP")))}else{tagList(selectInput("selectdimred","Select dimensionality reduction method.",choices=c("tSNE")))}
            }else{tagList(selectInput("selectdimred","Select dimensionality reduction method.",choices=c("tSNE")))}})
       output$selectdimred2<-renderUI({
         if(input$selectformat=="Seurat3"){
           aa<-sc@active.assay
           if(isTruthy(grepl(paste0("RunUMAP.",aa),names(sc@commands)))){tagList(selectInput("selectdimred2","Select dimensionality reduction method.",choices=c("tSNE","UMAP")))}else{tagList(selectInput("selectdimred","Select dimensionality reduction method.",choices=c("tSNE")))}
           }else{tagList(selectInput("selectdimred","Select dimensionality reduction method.",choices=c("tSNE")))}}) 
       output$cluSep<-renderPlot({
         sc<-values$sc
         plot_clu_separation(input$selectformat,sc)},height=700)
       output$tsneClu<-renderPlot({
         sc<-values$sc
         get_clu_plot(input$selectformat,sc,input$selectdimred)})
       output$silhPlot<-renderPlot({
         sc<-values$sc
         plot_silhouette(input$selectformat,sc)})
    ###########################################################################################################      
       observeEvent(input$plotclu, {
         showModal(modalDialog(title = "YOUR REQUEST IS BEING PROCESSED",
                               "The cells are being reclustered. Please allow (up to) some minutes.",
                               easyClose = TRUE)) 
           sc<-values$sc
           values$sc<-recluster_plot_tsne(input$selectformat,sc,input$numclu)
           output$tsneClu<-renderPlot({
             sc<-values$sc
             get_clu_plot(input$selectformat,sc,input$selectdimred)})
           output$cluSep<-renderPlot({
             sc<-values$sc
             plot_clu_separation(input$selectformat,sc)},height=700)#,width=600
           output$silhPlot<-renderPlot({
             sc<-values$sc
             plot_silhouette(input$selectformat,sc)})
           #})
       },ignoreInit=TRUE)#end observe plotclu
       
    ###########################################################################################################   
       
        observeEvent(input$getmkrs, {
          showModal(modalDialog(title = "YOUR REQUEST IS BEING PROCESSED",
                                "The markers for the new clusters are being extracted. Please allow (up to) some minutes.",
                                easyClose = TRUE))  
         sc<-values$sc
         top10_seuset<-get_top10(input$selectformat,sc)
         top10<-top10_seuset[[1]]
         
         seuset<-top10_seuset[[2]]#
    ######################################################################
         observeEvent(input$numDEGs, { 
             req(input$getmkrs)            
             resnL<-lapply(unique(top10$Cluster),function(X){
               head(top10[top10$Cluster %in% X,],n=input$numDEGs)})
             topn<-as.data.frame(do.call(rbind,resnL))
             mdict<-c("RaceID3"="padj","Monocle2"="p_val_adj","Seurat3"="p_val_adj")
             topn<-topn[with(topn, order(Cluster, eval(as.name(mdict[input$selectformat])))),]
             output$topn<-renderTable({topn})
             values$topn<-topn
             sc<-values$sc
             output$geneheatmap<-renderPlot({get_marker_plot(input$selectformat,sc,topn,seuset)})
             
           })#end of observe numDEGs
          
           
       },ignoreInit=TRUE)#end observe getmkrs
       
       
       
    ####################################################################   
        
        output$downloadTable <- downloadHandler(
          filename = function() {
            paste("sc.clu",input$numclu,".top",input$numDEGs, ".tsv", sep = "")
          },
          content = function(file) {
            line<-paste0("#",datPath)
            write(line, file=file, append=FALSE )
            write.table(values$topn, file, row.names = FALSE,sep="\t",quote=FALSE,append=TRUE)
          }
        )
    #####################################################################    
        
        #render the head
        sc<-values$sc
        ntemp<-render_data_head(input$selectformat,sc)
        values$ndata<-ntemp[rowSums(ntemp)>0,]
        ndata<-values$ndata
        output$datHead<-renderTable({ndata[1:10,1:min(8,ncol(ndata))]},caption="Normalized data",caption.placement = getOption("xtable.caption.placement", "top"),include.rownames=TRUE)
        output$dataDims<-renderText({sprintf("Your input has %s nonzero rows and %s columns.",nrow(ndata),ncol(ndata))})
        values$summaryTPC<-summary(colSums(as.matrix(ndata)))
        output$summaryTPC<-renderPrint({values$summaryTPC})
        
        
        }else{ #end of if init_checks_passed
        #waiting_for_click(1)
        
        showModal(modalDialog(
          title = "NO R PACKAGE SELECTED",
          "Please provide the missing selection and re-submit dataset before continuing!",
          easyClose = TRUE
        ))}
        
       },ignoreInit=TRUE)#end of observe input$submitinput  
        #}) #end of observe
      
      
          observeEvent(input$selectgenes,{
          inGenesL<-isolate(input$geneid)
          if(inGenesL!=""){
             inGenes<-unique(unlist(strsplit(inGenesL,split=";")))}
          if(!input$selectformat %in% "Seurat3"){inGenes<-gsub("--","__",inGenes)}
          inGenes<-trimws(inGenes)
          values$inGenes<-inGenes
          output$genesSel<-renderText({paste0("Selected genes: ",paste0(inGenes,collapse=" "))})
          ndata<-isolate(values$ndata)
          
          nv<-inGenes[inGenes %in% rownames(ndata)]
          if(!isTruthy(nv)){showModal(modalDialog(title = "NO EXPRESSED GENES IN SELECTION!","The GeneIDs you provided either don't match your countdata gene identifiers or are not expressed in at least 1 cell.",easyClose = TRUE))}
          req(nv)
          output$genesExpr<-renderText({paste0("Expressed genes: ",paste0(nv,collapse=" "))})

          
        },ignoreInit=TRUE) #end observe selectgenes
        
          
          
       observeEvent(input$plottsne,{
         
            inGenes<-isolate(values$inGenes)
            
            sc<-isolate(values$sc)
            ndata<-isolate(values$ndata)

            nv<-inGenes[inGenes %in% rownames(ndata)]
            
            if(length(nv)>0){
                nt<-isolate(input$tsnetit)
                output$tsneAgg<-renderPlot({get_feature_plot(input$selectformat,sc,nv,nt,as.logical(input$tsnelog),input$selectdimred2)})  
            
            }#fi
       },ignoreInit=TRUE)#end of observe plottsne
       
       
       observeEvent(input$selectgenes2,{
         inGenesL<-isolate(input$geneid2)
         if(inGenesL!=""){
           inGenes<-unique(unlist(strsplit(inGenesL,split=";")))}
         if(!input$selectformat %in% "Seurat3"){inGenes<-gsub("--","__",inGenes)}
         inGenes<-trimws(inGenes)
         values$inGenes2<-inGenes
         output$genesSel2<-renderText({paste0("Selected genes: ",paste0(inGenes,collapse=" "))})
         ndata<-isolate(values$ndata)
         
         nv<-inGenes[inGenes %in% rownames(ndata)]
         if(!isTruthy(nv)){showModal(modalDialog(title = "NO EXPRESSED GENES IN SELECTION!","The GeneIDs you provided either don't match your countdata gene identifiers or are not expressed in at least 1 cell.",easyClose = TRUE))}
         req(nv)
         output$genesExpr2<-renderText({paste0("Expressed genes: ",paste0(nv,collapse=" "))})
         
         
         inGenes<-isolate(values$inGenes2)
         ndata<-isolate(values$ndata)
         
         nv<-inGenes[inGenes %in% rownames(ndata)]
            ###produce top correlated genes for aggregated selected gene(s)
            cor.log2<-cor(x=log2(colSums(ndata[rownames(ndata) %in% nv,])+0.1),y=t(log2(ndata+0.1))) 
            cor.log2T<-as.data.frame(t(cor.log2),stringsAsFactors=FALSE)
            colnames(cor.log2T)<-"cor"
            cor.log2T$abscor<-abs(cor.log2T$cor)
            cor.log2T<-cor.log2T[order(cor.log2T$abscor,decreasing=TRUE),]
            output$corlog2<-renderPlot({ggplot(cor.log2T)+geom_violin(aes(x="all",y=cor))})
            
            output$top10cor<-renderTable({head(cor.log2T[,"cor",drop=FALSE],n=10)},caption="Top 10 correlated genes",caption.placement = getOption("xtable.caption.placement", "top"),include.rownames=TRUE)
            
       },ignoreInit=TRUE)#end of observe selectgenes2
 


        observeEvent(input$plotpwcor,{
            sc<-isolate(values$sc)
            ndata<-isolate(values$ndata)
            if((input$pwselX!="")&(input$pwselY!="")){
            inpwselXL<-isolate(input$pwselX)
            inpwselYL<-isolate(input$pwselY)}
            if(!input$selectformat %in% "Seurat3"){inpwselX<-trimws(gsub("--","__",unlist(strsplit(inpwselXL,split=";"))))
            inpwselY<-trimws(gsub("--","__",unlist(strsplit(inpwselYL,split=";"))))}else{
              inpwselX<-trimws(unlist(strsplit(inpwselXL,split=";")))
              inpwselY<-trimws(unlist(strsplit(inpwselYL,split=";")))
            }
            xtest<-inpwselX[inpwselX %in% rownames(ndata)]
            ytest<-inpwselY[inpwselY %in% rownames(ndata)]
            if(any(!isTruthy(xtest),!isTruthy(ytest))){
              showModal(modalDialog(title = "NO EXPRESSED GENES IN SELECTION!","The GeneIDs you provided either don't match your countdata gene identifiers or are not expressed in at least 1 cell.",easyClose = TRUE))
            }
            req(xtest,ytest)
            plotdat<-as.data.frame(cbind(colSums(ndata[rownames(ndata) %in% inpwselX,]),colSums(ndata[rownames(ndata) %in% inpwselY,])),stringsAsFactors=FALSE)
            colnames(plotdat)<-c("X","Y")
            
            corv<-round(cor(x=log2(plotdat$X+0.1),y=log2(plotdat$Y+0.1)),digits=2)
            pt<-isolate(input$pwcortit)
            output$pwplot<-renderPlot({ggplot(data=plotdat)+geom_point(aes(x=X,y=Y))+ggtitle(paste(pt,"cor=",corv,sep=" "))})

       },ignoreInit=TRUE)#end of observe input$plotpwcor 
        
        cludesc<-c("RaceID3"="Kmedoids clustering was run on logpearson distances between cells.","Monocle2"="Density peak clustering was run on distances between cells.","Seurat3"="Louvain partitioning of the shared nearest neighbourgh graph was applied.")
        
         output$get_vignette <- downloadHandler(
           filename = "scRNAseq_app_vignette.html",
           content = function(con) {
             file.copy(from="/data/manke/sikora/shiny_apps/scRNAseq_docs/scRNAseq_vignette.html", to=con, overwrite =TRUE)
               }
         )
         
         output$get_raceid <- downloadHandler(
           filename = "raceid3_dataset.RData",
           content = function(con) {
             file.copy(from="/data/processing/scRNAseq_shiny_app_example_data/GSE81076_raceid.workspaceR/sc.minT1000.RData", to=con, overwrite =TRUE)
           }
         )
         
         output$get_monocle <- downloadHandler(
           filename = "monocle3_dataset.RData",
           content = function(con) {
             file.copy(from="/data/processing/scRNAseq_shiny_app_example_data/GSE81076_monocle.workspaceR/minT5000.mono.set.RData", to=con, overwrite =TRUE)
           }
         )
         
         output$get_seurat <- downloadHandler(
           filename = "seurat3_dataset.RDS",
           content = function(con) {
             file.copy(from="/data/processing/scRNAseq_shiny_app_example_data/GSE75478_seuset.umap.RDS", to=con, overwrite =TRUE)
           }
         )
         
         output$downloadData <- downloadHandler(
           filename = "MyData.RDS",
           content = function(con) {
             base::saveRDS(isolate(values$sc),file=con)
           }
         )
        
         output$geneid<-renderUI({tagList(textInput(inputId="geneid", label="GeneID", value="",placeholder="TYPE IN GENE ID"),
                                  actionButton("selectgenes", "Select genes",style = "color: black;background-color:#6495ED"))})
         
         output$geneid2<-renderUI({tagList(textInput(inputId="geneid2", label="GeneID", value="",placeholder="TYPE IN GENE ID"),
                                          actionButton("selectgenes2", "Select genes",style = "color: black;background-color:#6495ED"))})
         
         output$summary_produced<-reactive({isTruthy(values$summaryTPC)})
         output$waitmssg<-renderText({"Please wait until data loading is completed. This may take some minutes."})
         output$initmssg<-renderText({"No data has been loaded. Please upload a file and click the 'Submit dataset' button."})

         outputOptions(output, "summary_produced", suspendWhenHidden = FALSE)
         
         # 
############################
    output$resultPanels<-renderUI({myTabs<-list(tabPanel(title="Walkthrough",
                                                      fluidPage(
                                                        fluidRow(
                                                          box(title="Walkthrough",uiOutput("walkThrough")),
                                                          box(title="Miscellaneous information",textOutput("FAQ")),
                                                          downloadButton("get_vignette", label = "Download vignette html")),
                                                        fluidRow(
                                                          box(title="Example datasets",width=12,
                                                              downloadButton("get_raceid",label="Download RaceID3 example dataset"),
                                                              downloadButton("get_monocle",label="Download Monocle2 example dataset"),
                                                              downloadButton("get_seurat",label="Download Seurat3 example dataset"))
                                                               )
                                                      )
                                                          ),
                                                  tabPanel(title="Input Data",
                                                           conditionalPanel(condition=("output.summary_produced"),
                                                      fluidPage(
                                                          fluidRow(
                                                            tableOutput("datHead"),
                                                            box(textOutput("dataDims"),width=4,title="Dimensions of your data."),
                                                            box(verbatimTextOutput("summaryTPC"),width=4,title="Summary of transcript per cell (TPC) in your data.")
                                                                   )
                                                             ) 
                                                            ),
                                                     # conditionalPanel(condition="(input.adddataset>0 ||$('html').hasClass('shiny-busy') || $('summaryTPC').hasClass('recalculating'))&&((!output.summary_produced)&&(output.init_checks_passed)", # &&(!output.waiting_for_click))
                                                      #                 fluidPage(
                                                       #                  fluidRow(
                                                        #                   box(textOutput("waitmssg")),
                                                         #                  renderImage({list(src="/data/manke/sikora/shiny_apps/Agmh.gif",width=200,height=200)},deleteFile =FALSE)
                                                          #              )
                                                           #           )
                                                            #           ),
                                                      conditionalPanel(condition="(input.adddataset==0) || !output.init_checks_passed ", #|| output.waiting_for_click
                                                                       fluidPage(
                                                                         fluidRow(
                                                                           box(textOutput("initmssg"))
                                                                         )
                                                                       )
                                                      )
                                                              
                                                          ),
                                                ##
                                                tabPanel(title="Cell map and clustering",
                                                         fluidPage(
                                                           fluidRow(
                                                             box(title="Cell map with cluster assignment",plotOutput("tsneClu"),width=5,height=600),
                                                             box(title="Method Description",renderText(cludesc[input$selectformat])),
                                                             box(title = "Plot controls",uiOutput("CluCtrl"),uiOutput("selectdimred"),actionButton(inputId="plotclu",label="Update cluster plots",style = "color: black;background-color:#6495ED"))),
                                                           fluidRow(
                                                             box(title="Metrics for cluster number selection",plotOutput("cluSep"),height=800),
                                                             box(title="Silhoutte Plot",plotOutput("silhPlot"))
                                                           )
                                                           
                                                         )

                                                        
                                                ),
                                                
                                                #tabPanel(title="Annotation.Table",
                                                         #fluidRow(
                                                           #column(4,uiOutput("configurator"))
                                                           #),
                                                         
                                                         #DTOutput("gtf"),
                                                         #actionButton(inputId="selGenesFromTab",label="Select Gene IDs from table"),
                                                         #actionButton(inputId="clearRowSel",label="Clear row selection")
                                                         
                                                #),
                                                tabPanel(title="Marker Gene Calculation",
                                                         fluidPage(
                                                           fluidRow(
                                                             box(title="Heatmap of top marker genes",plotOutput("geneheatmap"),width=6),
                                                             box(title="Marker genes",actionButton(inputId="getmkrs",label="Get marker genes",style = "color: black;background-color:#F5DEB3"),sliderInput("numDEGs", "Number of top markers to show",min=1,max=10,value=2,round=TRUE),tableOutput("topn"),width=6)),
                                                           fluidRow(box(downloadButton(outputId="downloadTable", label="Download table"),width=5))
                                                         )
                                                
                                                ),
                                                   tabPanel(title="Marker Gene Visualization",
                                                      fluidPage(
                                                          box(plotOutput("tsneAgg"),width=5),
                                                          box(title = "Plot controls",selectInput("tsnelog", "Log scale",choices=c("TRUE","FALSE"),selected="TRUE"),textInput("tsnetit","Plot title",value="Selected genes",placeholder="TYPE IN PLOT TITLE"),uiOutput("selectdimred2"),actionButton(inputId="plottsne",label="Plot cell map",width="200px",style = "color: black;background-color:#6495ED")),
                                                          box(title="Method Description",renderText("Counts were aggregated over selected genes as sum and the resulting (log2) expression levels were colour-coded on the cell map.")),
                                                          box(uiOutput("geneid",width=4),textOutput("fileDescription")),
                                                          box(title="Genes used",textOutput("genesSel"),textOutput("genesExpr"),width=5)
                                                          
                                                                )
                                                              ),
                                                   tabPanel(title="Correlation Analyses",
                                                      fluidPage(
                                                        fluidRow(
                                                          box(plotOutput("corlog2"),width=4),
                                                          box(tableOutput("top10cor"),width=4,height=420),
                                                          box(uiOutput("geneid2"),textOutput("fileDescription2"),width=4)
                                                          ),
                                                        fluidRow(
                                                          box(title="Method Description",renderText("Pearson correlation was calculated between log2-transformed aggregated counts (sum) for gene selection and all log2-transformed genes in the ndata slot of the sc object. Top 10 genes are listed.")),
                                                          box(title="Genes used",textOutput("genesSel2"),textOutput("genesExpr2"))
                                                          ),
                                                       fluidRow(
                                                          box(plotOutput("pwplot"),width=4),
                                                          box(title = "Select gene(s) X",textInput(inputId="pwselX", label="Gene symbol(s) for X axis (semicolon-separated)",value="")),
                                                          box(title = "Select gene(s) Y",textInput(inputId="pwselY", label="Gene symbol(s) for Y axis (semicolon-separated)",value="")),
                                                          box(textInput("pwcortit","Plot title",value="Selected genes",placeholder="TYPE IN PLOT TITLE"))
                                                           ),
                                                       fluidRow(actionButton(inputId="plotpwcor",label="Plot expression",width="200px",style = "color: black;background-color:  	 	#6495ED"),
                                                                box(title="Method Description",renderText("Pairwise plot of normalized counts, aggregated by summing over provided Gene IDs for each axis.")))
                                                          )
                                                   ),
                                                  
                                                  tabPanel(title="Session Info",
                                                      fluidPage(
                                                          verbatimTextOutput("sessionInfo"),
                                                          downloadButton(outputId="downloadSessionInfo", label="Download session info"),
                                                          downloadButton(outputId="downloadData", label="Download your data")
                                                               )
                                                          )

                                                 )

            do.call(tabsetPanel, myTabs)})


}

shinyApp(ui, server)
