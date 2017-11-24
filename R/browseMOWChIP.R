#' create HTTP paths to S3 bucket holding selected MOWChIP bigwig files
#' @examples
#' mowChInS3()
#' @export
mowChInS3 = function() {
#c("http://s3.amazonaws.com/bcfound-mow/GSM1599147_GM12878_K4me3_10K_Rep1.bw",
#"http://s3.amazonaws.com/bcfound-mow/GSM1599149_GM12878_K4me3_1K_Rep1.bw",
#"http://s3.amazonaws.com/bcfound-mow/GSM1599151_GM12878_K4me3_600_Rep1.bw",
#"http://s3.amazonaws.com/bcfound-mow/GSM1599153_GM12878_K4me3_100_Rep1.bw")
k4n = c("GSM1599147_GM12878_K4me3_10K_Rep1.bw",
"GSM1599149_GM12878_K4me3_1K_Rep1.bw",
"GSM1599151_GM12878_K4me3_600_Rep1.bw",
"GSM1599153_GM12878_K4me3_100_Rep1.bw")
acn = c("GSM1599155_GM12878_K27ac_10K_Rep1.bw",
"GSM1599157_GM12878_K27ac_1K_Rep1.bw",
"GSM1599159_GM12878_K27ac_600_Rep1.bw",
"GSM1599161_GM12878_K27ac_100_Rep1.bw")
sprintf("http://s3.amazonaws.com/bcfound-mow/%s", c(k4n,acn))
}

#' @importFrom GenomeInfoDb "seqlevelsStyle<-"
filterToGene = function( gf, geneSymbol, 
     radius, anno=EnsDb.Hsapiens.v75, 
     slStyle = "UCSC", stdSeqs=TRUE ) {
  mod = try(genemodel(geneSymbol, annoResource=anno))
  if (inherits(mod, "try-error")) stop(sprintf("genemodel fails for %s with given annoResource", geneSymbol))
  gm = range(genemodel(geneSymbol,
        annoResource=anno))+radius
  gm = keepStandardChromosomes(gm, pruning.mode="coarse")
  seqlevelsStyle(gm) = slStyle
  rowRanges(gf) = gm
  reduceByRange(gf, MAP=function(r,f)
            import.bw(f, which=r))
}


#' multimodal Gviz-based visualization of MOWChIP tracks
#' @import Gviz
#' @import EnsDb.Hsapiens.v75 
#' @import GenomicFiles 
#' @import ensembldb 
#' @import erma
#' @import GenomicRanges
#' @import IRanges
#' @param gf GenomicFiles instance giving paths to .bw files
#' @param sym gene symbol
#' @param radius number of bp around coding region to display
#' @param gstr a GRanges with gene annotation information, e.g., genes(EnsDb.Hsapiens.v75)
#' @param namegen a function of one argument assumed to be a row of \code{colData(gf)} that will be used as the \code{name} parameter for the
#' \code{\link[Gviz]{DataTrack}} derived from each sample
#' @param \dots passed to \code{\link[Gviz]{plotTracks}} exclusive of type which is set to "polygon"
#' @examples
#' data(caoChIP)
#' S3paths = mowChInS3()
#' if (interactive()) {
#'  remoteMowK4me3 = GenomicFiles(files=S3paths[1:4], 
#'     colData=colData(caoChIP[,c(1,3,5,7)]))
#'  remoteMowK27ac = GenomicFiles(files=S3paths[5:8], 
#'     colData=colData(caoChIP[,c(9,11,13,15)]))
#'  require(EnsDb.Hsapiens.v75)
#'  gg = genes(EnsDb.Hsapiens.v75)
#'  require(GenomeInfoDb)
#'  seqlevelsStyle(gg) = "UCSC"
#'  viewByGene(gf=remoteMowK4me3, sym="IGHA2", radius=80000, gstr=gg)
#'  viewByGene(gf=remoteMowK27ac, sym="SPI1", radius=80000, gstr=gg)
#'  }
#' @export
viewByGene = function(gf=caoChIP[,1:6], sym="IGHA2",
    radius=100000, gstr, namegen = function(x)x$numCells[1], ...) {
  jj = try(filterToGene(gf, sym, radius)[[1]])
  if (inherits(jj, "try-error")) stop("filterToGene failed")
  nms = vapply(1:length(jj), function(i)
                  namegen(colData(gf)[i,]), character(1))
  zz = lapply(1:length(jj), function(i) DataTrack(jj[[i]], name=nms[i]))
  sn = as.character(seqnames(jj[[1]])[1])
  grt = genesInRegion(sn, start=min(start(jj[[1]])), # jj[[1]] already has radius
          end=max(end(jj[[1]])), gstr=gstr, radius=0)
  plotTracks(c(zz, grt, GenomeAxisTrack(littleTicks=TRUE)), type="polygon", ...)
}

biotypes_ig = function() sprintf("IG_%s_gene", c("C", "D", "V", "J"))

genesInRegion = function(chr, start, end, gstr, radius=0,
   biotypes = c("protein_coding", biotypes_ig()),
   asGRT = TRUE) {
  ans = subsetByOverlaps(gstr, GRanges(chr, IRanges(start,end))+radius)
  ans = ans[which(ans$gene_biotype %in% biotypes)]
  seqlevels(ans) = chr
  if (asGRT) return(GeneRegionTrack(ans, showId=TRUE, name=chr[1]))
  ans
}


getStandardModel = function(sym) {
  mymod = genemodel(sym, anno=EnsDb.Hsapiens.v75)
  seqlevelsStyle(mymod) = "UCSC"
  keepStandardChromosomes(mymod, pruning.mode="coarse")
}
getGRT = function(syms) {
  m = lapply(syms, getStandardModel)
  do.call(GeneRegionTrack, m)
}
