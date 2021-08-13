#Upload Ensembl data for exon breakpoints

makeGeneModel <- function(snvs, cnvs) {
  
  ensembl_annotation_gene <- fread("figures/Figure3/MECP2.gene_model.bed")
  
  ## First step is to setup the appropriate gene model with "tiny" introns so that the reader can focus on coding sequence
  # Create initial exon breakpoints
  coding.start <- ensembl_annotation_gene[,cdsStart]
  coding.end <- ensembl_annotation_gene[,cdsEnd]
  
  exon.starts <- as.integer(unlist(str_split(ensembl_annotation_gene[,exonStarts],",")))
  exon.starts <- exon.starts[1:length(exon.starts)-1]
  exon.lengths <-  as.integer(unlist(str_split(ensembl_annotation_gene[,exonLengths],",")))
  exon.lengths <- exon.lengths[1:length(exon.lengths)-1]
  for (i in 1:length(exon.starts)) {
    exon.starts[i] <- exon.starts[i] + ensembl_annotation_gene[,start]
  }
  exon.ends <- c()
  for (i in 1:length(exon.lengths)) {
    exon.ends <- c(exon.ends,exon.starts[i] + exon.lengths[i])
  }
  
  exon.widths <- abs(exon.ends - exon.starts)
  exon.count <- as.numeric(ensembl_annotation_gene[,numExons])
  strand <- as.character(ensembl_annotation_gene[,strand])
  
  ## Adjust to small introns:
  exon.starts.rel <- exon.starts - exon.starts[1]
  exon.ends.rel <- exon.ends - exon.starts[1]
  
  coding.start.rel <- coding.start - exon.starts[1]
  coding.end.rel <- coding.end - exon.starts[1]
  
  ## Split exons based on coding start into coding and noncoding:
  ## This just sets the relative hight of coding exons versus UTR in the final plot:
  height.cds <- 0.1
  height.utr <- 0.04
  
  exon.sizes <- c()
  exon.starts.cds <- c()
  exon.ends.cds <- c()
  exon.widths <- c()
  
  for (i in 1:length(exon.starts.rel)) {
    curr.start <- exon.starts.rel[i]
    curr.end <- exon.ends.rel[i]
    
    ## Exon is completely prior to CDS start
    if (curr.start < coding.start.rel & curr.end < coding.start.rel) {
      exon.sizes <- c(exon.sizes,height.utr)
      exon.starts.cds <- c(exon.starts.cds, curr.start)
      exon.ends.cds <- c(exon.ends.cds, curr.end)
      exon.widths <- c(exon.widths,curr.end - curr.start)
      ## Exon is completely after CDS end
    } else if (curr.start >= coding.end.rel & curr.end >= coding.end.rel) {
      exon.sizes <- c(exon.sizes,height.utr)
      exon.starts.cds <- c(exon.starts.cds, curr.start)
      exon.ends.cds <- c(exon.ends.cds, curr.end)
      exon.widths <- c(exon.widths,curr.end - curr.start)
    } else if (curr.start < coding.start.rel & curr.end > coding.start.rel) {
      # Have to split the exon in twain
      ex.1.start <- curr.start
      ex.1.end <- coding.start.rel - 1
      ex.2.start <- coding.start.rel
      ex.2.end <- curr.end
      
      exon.sizes <- c(exon.sizes,height.utr)
      exon.sizes <- c(exon.sizes,height.cds)
      
      exon.starts.cds <- c(exon.starts.cds, ex.1.start, ex.2.start)
      exon.ends.cds <- c(exon.ends.cds, ex.1.end, ex.2.end)
      
      exon.widths <- c(exon.widths,ex.1.end - ex.1.start)
      exon.widths <- c(exon.widths,ex.2.end - ex.2.start)
      
    } else if (curr.start < coding.end.rel & curr.end > coding.end.rel) {
      # Have to split the exon in twain
      ex.1.start <- curr.start
      ex.1.end <- coding.end.rel
      ex.2.start <- coding.end.rel + 1
      ex.2.end <- curr.end
      
      exon.sizes <- c(exon.sizes,height.cds)
      exon.sizes <- c(exon.sizes,height.utr)
      
      exon.starts.cds <- c(exon.starts.cds, ex.1.start, ex.2.start)
      exon.ends.cds <- c(exon.ends.cds, ex.1.end, ex.2.end)
      
      exon.widths <- c(exon.widths,ex.1.end - ex.1.start)
      exon.widths <- c(exon.widths,ex.2.end - ex.2.start)
    } else {
      exon.sizes <- c(exon.sizes,height.cds)
      
      exon.starts.cds <- c(exon.starts.cds, curr.start)
      exon.ends.cds <- c(exon.ends.cds, curr.end)
      
      exon.widths <- c(exon.widths,curr.end - curr.start)
    }
    
  }
  
  ## Next is to actually make the small introns themselves:
  exon.count <- length(exon.starts.cds)
  intron.dist <- c(2:exon.count)
  
  for (i in c(2:exon.count)) {
    
    dist <- exon.starts.cds[i] - exon.ends.cds[i-1]
    if (dist <= 200) {
      intron.dist[i-1] <- dist
    } else {
      intron.dist[i-1] <- 200
    }
    
  }
  
  exon.starts.rel.2 <- exon.starts.cds
  exon.ends.rel.2 <- exon.ends.cds
  
  for (i in c(2:exon.count)) {
    
    exon.starts.rel.2[i] <- exon.ends.rel.2[i-1] + intron.dist[i-1]
    exon.ends.rel.2[i] <- exon.starts.rel.2[i] + exon.widths[i]
    
  }
  
  ## This bit builds the final "mapping" key to translate the original variant positions into the new coordinate system I've designed above
  key <- data.table(original=c(0:exon.ends.rel[length(exon.ends.rel)]))
  
  for (i in 1:length(exon.starts.rel.2)) {
    
    x <- c(exon.starts.cds[i]:exon.ends.cds[i])
    y <- c(exon.starts.rel.2[i]:exon.ends.rel.2[i])
    
    for (z in 1:length(x)) {
      
      key[original == x[z],new:=y[z]]
      
      
    }
    
    # Now do intron mapping:
    if (i != length(exon.starts.rel.2)) {
      intron.left <- exon.ends.rel.2[i] + 20
      intron.right <- exon.starts.rel.2[i+1] - 20
      mid <- (exon.ends.cds[i] + exon.starts.cds[i+1]) / 2
      for (z in c((exon.ends.cds[i]+1):(exon.starts.cds[i+1]-1))) {
        if (z <= mid) {
          key[original == z,new:=intron.left]
        } else {
          key[original == z,new:=intron.right]
        }
      }
    }
  }
  
  ## This function does the actual translating:
  lookup <- function(og) {
    
    if (og < 0) {
      return(as.integer(og))
    } else {
      val <- key[original==og,new]
      return(val)
    }
  }
  
  # This is the final GRanges object of the actual gene model
  chromosome <- ensembl_annotation_gene[,chrom]
  features.exons <- GRanges(seqnames = chromosome, strand = strand, IRanges(start = exon.starts.rel.2, width = exon.widths))
  features.exons$fill <- c(col.palette[1])

  ## Set the heights according to what we did above with coding/noncoding 
  features.exons$height <- exon.sizes
   
  snvs[,Position.rel:=start - exon.starts[1]]
  snvs[,Position.rel.f:=lookup(Position.rel),by=1:nrow(snvs)]
  
  cnvs[,Position.rel:=start - exon.starts[1]]
  cnvs[,Position.rel.f:=lookup(Position.rel),by=1:nrow(cnvs)]
  
  cnvs[,End.rel:=end - exon.starts[1]]
  cnvs[,End.rel.f:=lookup(End.rel),by=1:nrow(cnvs)]
  
  ## Set scores:
  snvs[,dummy:=1]
  counts <- snvs[,sum(dummy),by="start"]
  
  mutations.gr <- data.table()
  
  for (i in 1:nrow(counts)) {
    
    collected.sites <- snvs[start == counts[i,start]]
    num <- counts[i,V1]
    collected.sites[,score:=.I]
    
    mutations.gr <- bind_rows(mutations.gr,collected.sites[score==num])
      
  }
  
  positions <- mutations.gr[,Position.rel.f]
  mutations.gr <- GRanges(seqnames = chromosome, IRanges(positions, width=1), score = mutations.gr[,score], color = as.character(mutations.gr[,colour]), shape = as.character("circle"))

  return(list(features = features.exons, snvs = mutations.gr, cnvs = cnvs))
  
}
