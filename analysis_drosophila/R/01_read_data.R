#!/usr/bin/Rscript --no-site-file
library(coala)

# Analysis Options
locus_length <- 150000

# Global Config
options(stringsAsFactors = FALSE)

# Read individuals & populations
individuals <- read.delim("data/indivs.txt")

ind_africa <- c("CK1", "CO10N", "EA119", "EB132", "ED10N", 
                "EF120", "EM239", "ER11",  "GA125", "GU10", 
                "KN133N", "NG10N", "RC1",   
                "SB10",  "SD145", "SE16",  "SF110", "SP173",   
                "TZ10",  "UG19",  "ZI103", "ZS11")

ind_europe <- c('FR14', 'FR151', 'FR207', 'FR310',
                'NL1', 'NL2', 'NL12', 'NL14', 'NL17', 'NL18', 'NL19')


for (i in 1:5) stopifnot(all(c(ind_africa, ind_europe) %in% individuals[, i]))


read_file <- function(file) {
  # Read the sequence file
  seq_file <- gzfile(file)
  message("Loading: ", file)
  snp_table <- read.delim(seq_file, header = FALSE, stringsAsFactors = FALSE)

  chr <- snp_table[1, 1]
  message("Processing Chromosome: ", chr)

  # Determine which individuals we need
  ind_names <- as.character(individuals[[paste0("Chr", chr)]])
  ind_names <- ind_names[ind_names != ""]
  pop1_mask <- ind_names %in% ind_africa
  pop2_mask <- ind_names %in% ind_europe 
  populations <- list(africa = ind_names[pop1_mask], 
                      europe = ind_names[pop2_mask])
  req_inds_idx <- c(which(pop1_mask), which(pop2_mask))
  
  # Create a data.frame with a row for each individual
  variants_raw <- do.call(rbind.data.frame, lapply(snp_table$V3, function(x) {
    strsplit(x, '')[[1]][req_inds_idx]
  }))
  colnames(variants_raw) <- ind_names[req_inds_idx]
  
  variants <- data.frame(chr = snp_table$V1,
                         pos = snp_table$V2,
                         variants_raw,
                         ref = do.call(rbind, strsplit(snp_table$V4, '')))
  rm(snp_table, variants_raw)
  
  # Remove sites that
  # - have missing data ("N") or
  # - are none-segregating, 
  # - have a disagreeing outgroup sequences
  is_valid <- apply(variants, 1, function(x) {
    if (any(is.na(x) | x == "N")) return(1)
    if (length(unique(x[1:length(req_inds_idx) + 2])) == 1) return(2)
    if (x['ref.1'] != x['ref.2']) return(3)
    0
  })
  message("* filtering: Removing ", sum(is_valid > 0), " sites, ", 
          sum(is_valid == 0), " left")
  variants <- variants[is_valid == 0, -ncol(variants)]
  
  for (i in 3:(ncol(variants) - 1)) {
    variants[, i] <- as.numeric(variants[, i] != variants$ref.1)
  }
  
  variants$locus <- floor(variants$pos / locus_length)
  
  segsites <- by(variants, variants$locus, function(x) {
    ss <- create_segsites(snps = t(as.matrix(x[, 3:(ncol(variants) - 2)])),
                          positions = (x$pos %% locus_length) / locus_length)
    attr(ss, "chr") <- unique(x$chr)
    locus <- unique(x$locus)
    attr(ss, "locus") <- locus
    attr(ss, "pos") <- c(start = locus * locus_length, 
                         end = (locus + 1) * locus_length - 1,
                         min = min(x$pos),
                         max = max(x$pos))
    ss
  })
  
  # Remove first and last 3 of loci which have data
  # because the quality might still be bad there
  segsites <- tail(segsites, -3)
  segsites <- head(segsites, -3)
  
  # Remove loci with less than 1000 SNPs
  # Even sweeps in Pop2 should have polymorphic sites in Pop1
  segsites <- segsites[sapply(segsites, ncol) > 1000]
  
  message("* imported ", length(segsites), " loci with an total of ", 
          sum(sapply(segsites, ncol)), " SNPs ")
  message("")
  
  list(chr = chr, populations = populations, segsites = segsites)
}

snp_data <- lapply(list.files("data", "*.txt.gz", full.names = TRUE), read_file)
dir.create("cache", FALSE)
saveRDS(snp_data, file = "cache/snp_data.Rds")
saveRDS(locus_length, file = "cache/locus_length.Rds")
