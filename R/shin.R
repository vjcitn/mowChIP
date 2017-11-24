#' visualize coverage profiles for selected genes and MOWChIP protocols
#' @import shiny
#' @note This function sets up references to BigWig files in a public
#' AWS S3 bucket and will use user selections to extract relevant coverage
#' scores.
#' @export
mowchiBrowse = function() {
   requireNamespace("GenomeInfoDb")
   requireNamespace("GenomicFiles")
   requireNamespace("EnsDb.Hsapiens.v75")
   requireNamespace("ensembldb")

# data source
   data(caoChIP)
   S3paths = mowChInS3()
   remoteMow= GenomicFiles(files=S3paths, 
       colData=colData(caoChIP[,seq(1,15,2)]))

# gene symbols and ranges
   allg = genes(EnsDb.Hsapiens.v75)
   accty = c("protein_coding", "IG_C_gene", "IG_V_gene", "IG_D_gene",
        "IG_J_gene")
   allpc = allg[which(allg$gene_biotype %in% accty),] #"protein_coding"),] 
   seqlevelsStyle(allpc) = "UCSC"
   allsym = sort(unique(mcols(allpc)$symbol))
   
   ui = fluidPage(
    sidebarLayout(
     sidebarPanel(
      helpText("MOWChIP browser: Exploring ChIP-seq with as few as 100 cells ",
            a(href='https://www.ncbi.nlm.nih.gov/pubmed/?term=26214128',"(Cao et al. Nat Meth 2015)")),
      helpText("Select gene symbol (can delete and enter choice) to visualize MOWChIP coverage"),
      selectInput("sym", "gene", choices=allsym, selected="SPI1"),
      numericInput("rad", "radius(bp)", value=50000, min=5000, max=2e6, step=5000),
      checkboxGroupInput("mark", "mark", choices=c("H3K4me3", "H3K27ac"),
        selected="H3K4me3"),
      helpText("Figure legend: 100c: K4me3 denotes 100 cell protocol, H3K4me3 mark"),
      width=2),
     mainPanel(plotOutput("viewByGene", height="600px"))
     )
    )
   server = function(input, output) {
    output$viewByGene = renderPlot({
     marks = input$mark
     r = NULL
     if ("H3K4me3" %in% marks) r = which(colData(remoteMow)$mark == "K4me3")
     if ("H3K27ac" %in% marks) r = c(r, which(colData(remoteMow)$mark == "K27ac"))
     validate(need(!is.null(r), "please choose mark"))
     showNotification("retrieving coverage from AWS S3")
     viewByGene( remoteMow[, r], sym = input$sym, 
         radius = input$rad, gstr=allpc, 
         namegen = function(x) paste(x$numCells[1], "c: ", x$mark[1], collapse="", sep="") )
     })
   }
shinyApp(ui=ui, server=server)
}
     
