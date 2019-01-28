## app.R ##
#Rlib="/data/boehm/group/shiny_apps/Rlibs3.5.0"
Rlib="/data/manke/sikora/shiny_apps/Rlibs3.5.0_bioc3.7"
#Rlib="/data/manke/sikora/shiny_apps/Rlibs3.5.0_bioc3.6"

#library(Rcpp,lib.loc=Rlib,verbose=TRUE)#,lib.loc=Rlib
#library(scales,lib.loc=Rlib,verbose=TRUE)
#library(shiny,lib.loc=Rlib,verbose=TRUE)
sink("/var/log/shiny-server/sessionInfo.txt")
print(sessionInfo())
sink()
library(shinydashboard,lib.loc=Rlib,verbose=TRUE)
library(rhandsontable,lib.loc=Rlib,verbose=TRUE)
library(DT,lib.loc=Rlib,verbose=TRUE)

#options(shiny.error = browser)

ui <- function(request) {dashboardPage(
    dashboardHeader(title = "Dataset selection"),
    ## Sidebar content
    dashboardSidebar(

      selectInput(inputId="genome", label="Select organism", choices=c("PLEASE SELECT A GENOME","Zebrafish [zv10]","Fission yeast","Fruitfly [dm6]","Fruitfly [dm3]","Human [hg37]","Human [hg38]","Mouse [mm9]","Mouse [mm10]")),#"PLEASE SELECT A GENOME",, selected = NULL
      selectInput(inputId="selectformat",label="Select input file format",choices=c("Please select format","RaceID3","Monocle")),
      textInput(inputId="group", label="Group", value = "", width = NULL, placeholder = NULL),
      textInput(inputId="owner", label="Project Owner", value = "", width = NULL, placeholder = NULL),
      textInput(inputId="projectid", label="Project ID", value = "", width = NULL, placeholder = NULL),
      textInput(inputId="pathtodata", label="Data path", value = "", width = NULL, placeholder = NULL),
      actionButton(inputId="adddataset", label="Select dataset"),
      textInput(inputId="geneid", label="GeneID", value="",placeholder="TYPE IN GENE ID"),
      actionButton("selectgenes", "Select genes"),
      textOutput("fileDescription"),
      imageOutput("logo"),
      tags$footer("Copyright 2018 MPI-IE Freiburg Bioinfo Core Unit"),
      bookmarkButton()
        ),
        
    dashboardBody(
        h2("Single cell RNAseq analysis"),
        uiOutput("resultPanels")
               
    )

 )}


server <- function(input, output, session) {
    
    output$walkThrough<-renderUI(HTML("<ul><li>1.Provide group and project ID information to retrieve a serialized R object containing your dataset. Click on retrieve dataset. Your data will appear in the InputData tab.</li><li>2.Provide semicolon-separated Gene IDs to plot aggregate expression for.</li><li>3.#This is a currently not implemented## Provide rules for cell assignment to a known class that will be used to facet the plots.Click on Run analysis.#End of not implemented# Your results will appear in the corresponding tabs.</li><li>The order of providing the information matters!</li></ul>"))
    output$FAQ<-renderText("Currently, no uniform gene naming system is prerequisite. You have to provide Gene IDs consistent with the naming used to produce your dataset.\n Merging data from multiple datasets or batch effect removal are currenlty not supported.\n For questions, bug reports or feature requests, contact sikora@ie-freiburg.mpg.de.\n For reporting issues or pull requests on GitHub, go to https://github.com/katsikora/scRNAseq_shiny_app .")
    
    output$fileDescription<-renderText("GeneID: Please provide a semicolon-separated list of Gene IDs you would like to obtain results for.")
    
    output$logo<-renderImage({list(src="/data/manke/sikora/shiny_apps/userIN_to_yaml/MPIIE_logo_sRGB.jpg",width=100,height=100)},deleteFile =FALSE)


################################
    source("/data/manke/sikora/shiny_apps/scRNAseq_shiny_app/aux.R")
    #imports depend on selected format!
    #load packages in function of the input format (or use namespace loading...)
    observeEvent(input$selectformat,{try(load_libs(input$selectformat,Rlib),outFile="/var/log/shiny-server/library.err")},ignoreInit=TRUE)#
    

    output$sessionInfo <- renderPrint({capture.output(sessionInfo())})

    values<-reactiveValues()
    values$rowsSel<-""
    values$cList<-"All"
    values$inGenes=""
    

################################
    observeEvent(input$adddataset, {
      ######################################################################################################      
      psel<-c("Monocle"="*.mono.set.RData$","RaceID3"="^sc.minT*.RData$") 
      inFormat<-isolate(input$selectformat)
      
      if((input$group!="")&(input$owner!="")&(input$projectid!="")&(input$pathtodata=="")){
        inGroup<-isolate(input$group)
        inOwner<-isolate(input$owner)
        inProjectID<-isolate(input$projectid)
  
        values$datdir<-system(sprintf("find /data/%s/sequencing_data -name Analysis_%s_%s_%s -type d | sort",tolower(gsub("-.+","",inGroup)),inProjectID,inOwner,inGroup),intern=TRUE) 
        details<-file.info(dir(datdir,pattern="*.mono.set.RData$",full.names=TRUE,recursive=TRUE))
        details = details[with(details, order(as.POSIXct(mtime),decreasing=TRUE)), ]
        
        values$datpath<-rownames(details)[1]
              }
      
      else if ((input$group=="")&(input$owner=="")&(input$projectid=="")&(input$pathtodata!="")){
        values$datpath<-isolate(input$pathtodata)
      }  
      datPath<-isolate(values$datpath)
      
    
      if(grepl("rds$",datPath,ignore.case=TRUE)){
            values$sc<-readRDS(datPath)}
      else if (grepl("rdata$",datPath,ignore.case=TRUE)){
           myEnv<-environment()
           try(sctmp<-load(datPath, envir = myEnv),outFile="/var/log/shiny-server/RData.err")
           values$sc <- myEnv[[sctmp]]
             }
       sc<-values$sc
       
    ###########################################################################################################   
       cluinit<-try(get_cluinit(input$selectformat,sc),outFile="/var/log/shiny-server/cluinit.err")
       output$CluCtrl<-renderUI({tagList(sliderInput("numclu", "Number of clusters",min=1,max=2*cluinit,value=cluinit,round=TRUE))})
    ###########################################################################################################      
       observeEvent(input$plotclu, {
           values$sc<-recluster_plot_tsne(input$selectformat,sc,input$numclu)
           sc<-values$sc
       output$tsneClu<-try(renderPlot({get_clu_plot(input$selectformat,sc)}),outFile="/var/log/shiny-server/get_clu_plot.err")
       },ignoreInit=TRUE)#end observe plotclu
       
    ###########################################################################################################   
       
        observeEvent(input$getmkrs, {
         sc<-values$sc
         top10<-get_top10(input$selectformat,sc) #
    ######################################################################
         observeEvent(input$numDEGs, { 
             req(input$getmkrs)            
             resnL<-lapply(unique(top10$Cluster),function(X){
               head(top10[top10$Cluster %in% X,],n=input$numDEGs)})
             topn<-as.data.frame(do.call(rbind,resnL))
             mdict<-c("RaceID3"="padj","Monocle"="p_val_adj")
             topn<-topn[with(topn, order(Cluster, eval(as.name(mdict[input$selectformat])))),]
             output$topn<-renderTable({topn})
             values$topn<-topn
             output$geneheatmap<-renderPlot({get_marker_plot(input$selectformat,sc,topn)})
             
           })#end of observe numDEGs
          
           
       },ignoreInit=TRUE)#end observe getmkrs
       
       
       
    ####################################################################   
        
        output$downloadTable <- downloadHandler(
          filename = function() {
            paste("sc.clu",input$numclu,".top",input$numDEGs, ".tsv", sep = "")
          },
          content = function(file) {
            write.table(values$topn, file, row.names = FALSE,sep="\t",quote=FALSE)
          }
        )
    #####################################################################    
        
        #render the head
        sc<-values$sc
        sink("/var/log/shiny-server/pData_head.txt")
        ntemp<-try(render_data_head(input$selectformat,sc),outFile="/var/log/shiny-server/ntemp.err")
        values$ndata<-ntemp[rowSums(ntemp)>0,]
        ndata<-values$ndata
        output$datHead<-renderTable({ndata[1:10,1:min(8,ncol(ndata))]},caption="Normalized data",caption.placement = getOption("xtable.caption.placement", "top"),include.rownames=TRUE)
         orgv<-c("Zebrafish [zv10]"="GRCz10","Fission yeast"="SchizoSPombe_ASM294v2","Fruitfly [dm6]"="dm6","Fruitfly [dm3]"="dm3","Human [hg37]"="hs37d5","Human [hg38]"="GRCh38","Mouse [mm9]"="GRCm37","Mouse [mm10]"="GRCm38")
        #observeEvent(input$genome,{
            ens_dir<-dir(path=sprintf("/data/repository/organisms/%s_ensembl/ensembl",orgv[input$genome]),pattern="genes.gtf",full.names=TRUE,recursive=TRUE)
            gtf_path<-ens_dir[length(ens_dir)]
            gtf<-as.data.frame(rtracklayer::import(gtf_path))
            gtf<-unique(gtf[,c(1,5,10:16)])
            gtf$GeneSym<-paste0(gtf$gene_name,"__chr",gtf$seqnames)
        
            values$dat <- gtf
        
     
            output$configurator<-renderUI({tagList(selectInput(inputId="gene_biotype",label="Gene Biotype:",choices=c("All",unique(as.character(gtf$gene_biotype))),selected="All"),
                                                 selectInput(inputId="seqnames",label="Chromosome:",choices=c("All",unique(as.character(gtf$seqnames))),selected="All")) })  
        output$gtf<-renderDT({
            
          dat<-values$dat

        if (!input$gene_biotype %in% "All" & !is.na(input$gene_biotype)) {
          dat <- dat[dat$gene_biotype %in% input$gene_biotype,]
        }
        if (!input$seqnames %in% "All"& !is.na(input$seqnames)) {
          dat <- dat[dat$seqnames %in% input$seqnames,]
        }
         
          values$dat2<-dat
          dat},server=TRUE,options = list(autoWidth = TRUE,scrollX=TRUE), filter = "bottom")#end of renderDT
       # },ignoreInit=TRUE)#end of observe input$genome
        
        
       },ignoreInit=TRUE)#end of observe input$submitinput   
    
    
        misc<-observe({req(input$gtf_rows_selected)
                      values$rowsSel<-input$gtf_rows_selected})

      observeEvent(input$selectgenes,{
          inGenesL<-isolate(input$geneid)
          if(inGenesL!=""){
             inGenes<-unique(unlist(strsplit(inGenesL,split=";")))}
          inGenes<-gsub("--","__",inGenes)
          values$inGenes<-inGenes
          output$genesSel<-renderText({paste0("Selected genes: ",paste0(inGenes,collapse=" "))})
          output$genesSel2<-renderText({paste0("Selected genes: ",paste0(inGenes,collapse=" "))})
          ndata<-isolate(values$ndata)
          
          nv<-inGenes[inGenes %in% rownames(ndata)]
          output$genesExpr<-renderText({paste0("Expressed genes: ",paste0(nv,collapse=" "))})
          output$genesExpr2<-renderText({paste0("Expressed genes: ",paste0(nv,collapse=" "))})
          
          
        },ignoreInit=TRUE)
        
        
        observeEvent(input$selGenesFromTab,{
              dat2<-values$dat2
          
                  if(!values$rowsSel %in% ""){
                    dat2<-dat2[values$rowsSel,]}
                ##test which column to use
                   z<-apply(dat2[,c("gene_id","gene_name","GeneSym")],2,function(X) sum(rownames(values$ndata) %in% X ))
                   csel<-names(which.max(z))
              
                inGenes<-unique(dat2[,csel])
                values$inGenes<-inGenes
                output$genesSel<-renderText({paste0("Selected genes: ",paste0(inGenes,collapse=" "))})
                output$genesSel2<-renderText({paste0("Selected genes: ",paste0(inGenes,collapse=" "))})
                ndata<-isolate(values$ndata)
                
                nv<-inGenes[inGenes %in% rownames(ndata)]
                output$genesExpr<-renderText({paste0("Expressed genes: ",paste0(nv,collapse=" "))})
                output$genesExpr2<-renderText({paste0("Expressed genes: ",paste0(nv,collapse=" "))})
                
        },ignoreInit=TRUE)#end of observe input$selectgenesfromtab
        
        
       observeEvent(input$plottsne,{
         
            inGenes<-isolate(values$inGenes)
            
            sc<-isolate(values$sc)
            ndata<-isolate(values$ndata)

            nv<-inGenes[inGenes %in% rownames(ndata)]
            
            if(length(nv)>0){
                nt<-isolate(input$tsnetit)
                output$tsneAgg<-renderPlot({get_feature_plot(input$selectformat,sc,nv,nt,as.logical(input$tsnelog))})  
            
            }#fi
            ###produce top correlated genes for aggregated selected gene(s)
            cor.log2<-cor(x=log2(colSums(ndata[rownames(ndata) %in% nv,])+0.1),y=t(log2(ndata+0.1))) 
            output$corlog2<-renderPlot({boxplot(t(cor.log2))})
            cor.log2T<-as.data.frame(t(cor.log2),stringsAsFactors=FALSE)
            colnames(cor.log2T)<-"cor"
            cor.log2T$abscor<-abs(cor.log2T$cor)
            cor.log2T<-cor.log2T[order(cor.log2T$abscor,decreasing=TRUE),]
            output$top10cor<-renderTable({head(cor.log2T[,"cor",drop=FALSE],n=10)},caption="Top 10 correlated genes",caption.placement = getOption("xtable.caption.placement", "top"),include.rownames=TRUE)
 

       },ignoreInit=TRUE)#end of observe plottsne
       
       

        observeEvent(input$plotpwcor,{
            sc<-isolate(values$sc)
            ndata<-isolate(values$ndata)
            if((input$pwselX!="")&(input$pwselY!="")){
            inpwselXL<-isolate(input$pwselX)
            inpwselYL<-isolate(input$pwselY)}
            inpwselX<-unlist(strsplit(inpwselXL,split=";"))
            inpwselY<-unlist(strsplit(inpwselYL,split=";"))
            plotdat<-as.data.frame(cbind(colSums(ndata[rownames(ndata) %in% inpwselX,]),colSums(ndata[rownames(ndata) %in% inpwselY,])),stringsAsFactors=FALSE)
            colnames(plotdat)<-c("X","Y")
            
            corv<-round(cor(x=log2(plotdat$X+0.1),y=log2(plotdat$Y+0.1)),digits=2)
            pt<-isolate(input$pwcortit)
            output$pwplot<-renderPlot({ggplot(data=plotdat)+geom_point(aes(x=X,y=Y))+ggtitle(paste(pt,"cor=",corv,sep=" "))})

       },ignoreInit=TRUE)#end of observe input$plotpwcor 

   


############################
    output$resultPanels<-renderUI({myTabs<-list(tabPanel(title="WalkThrough",
                                                      fluidPage(
                                                          box(title="Walkthrough",uiOutput("walkThrough")),
                                                          box(title="Miscellaneous information",textOutput("FAQ"))                                                          
                                                               )
                                                          ),
                                                  tabPanel(title="InputData",
                                                      fluidPage(
                                                          fluidRow(
                                                            tableOutput("datHead"),
                                                              tableOutput("countDatHead")
                                                                   ),
                                                          #fluidRow(
                                                              rHandsontableOutput("hot"),
                                                              
                                                            #      ),
                                                          fluidRow(
                                                              box(renderText("Please complete missing sample information. Group variable will be used for plot faceting.")),
                                                              tableOutput("inTabHead"),
                                                              box(title="Debug",
                                                                  textOutput("debug"),
                                                                  textOutput("ingenes"))
                                                                  ) 
                                                              )
                                                          ),
                                                ##
                                                tabPanel(title="Cluster.Number",
                                                         fluidPage(
                                                           box(plotOutput("tsneClu"),width=5,height=600),
                                                           box(title = "Plot controls",uiOutput("CluCtrl")),
                                                           box(title="Method Description",renderText("Kmedoids clustering was run on logpearson distances between cells.")),
                                                           box(actionButton(inputId="plotclu",label="Plot clusters on tsne"),
                                                           actionButton(inputId="getmkrs",label="Get marker genes"),width=6),
                                                           box(plotOutput("geneheatmap"),width=6),
                                                           box(downloadButton(outputId="downloadTable", label="Download table"),width=5),
                                                           box(title="Marker genes",sliderInput("numDEGs", "Number of top markers",min=1,max=10,value=2,round=TRUE),tableOutput("topn"),width=5)
                                                         
                                                         )

                                                ),##
                                                
                                                tabPanel(title="Annotation.Table",
                                                         fluidRow(
                                                           column(4,uiOutput("configurator"))
                                                           ),
                                                         #fluidPage(
                                                         DTOutput("gtf"),
                                                         actionButton(inputId="selGenesFromTab",label="Select Gene IDs from table")
                                                         #)
                                                ),
                                                   tabPanel(title="Tsne.Map",
                                                      fluidPage(
                                                          box(plotOutput("tsneAgg"),width=5),
                                                          box(title = "Plot controls",selectInput("tsnelog", "Log scale",choices=c("TRUE","FALSE"),selected="TRUE"),textInput("tsnetit","Plot title",value="Selected genes",placeholder="TYPE IN PLOT TITLE")),
                                                          box(title="Method Description",renderText("(Log) counts were aggregated over selected genes and the expression levels were colour-coded on the tsne map.")),
                                                          box(title="Genes used",textOutput("genesSel"),textOutput("genesExpr")),
                                                          actionButton(inputId="plottsne",label="Plot tsne map")
                                                                )
                                                              ),
                                                   tabPanel(title="Top.Correl.Genes",
                                                      fluidPage(
                                                          box(plotOutput("corlog2"),width=4),
                                                          box(tableOutput("top10cor")),
                                                          box(title="Method Description",renderText("Pearson correlation was calculated between log2-transformed aggregated counts for gene selection and all log2-transformed genes in the ndata slot of the sc object. Top 10 genes are listed.")),
                                                          box(title="Genes used",textOutput("genesSel2"),textOutput("genesExpr2"))
                                                               )
                                                          ),
                                                  tabPanel(title="Pairwise.Expression",
                                                      fluidPage(
                                                          box(plotOutput("pwplot"),width=4),
                                                          box(title = "Select gene(s) X",textInput(inputId="pwselX", label="Gene symbol(s) for X axis",value="")),
                                                          box(title = "Select gene(s) Y",textInput(inputId="pwselY", label="Gene symbol(s) for Y axis",value="")),
                                                          box(textInput("pwcortit","Plot title",value="Selected genes",placeholder="TYPE IN PLOT TITLE")),
                                                          actionButton(inputId="plotpwcor",label="Plot expression"),
                                                          box(title="Method Description",renderText("Pairwise plot of normalized counts."))
                                                          
                                                               )
                                                          ),
                                                  tabPanel(title="sessionInfo",
                                                      fluidPage(
                                                          verbatimTextOutput("sessionInfo")                                                          
                                                               )
                                                          )

                                                 )

            do.call(tabsetPanel, myTabs)})


}

shinyApp(ui, server)
