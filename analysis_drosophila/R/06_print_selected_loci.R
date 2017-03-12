# This script creates a table of the selected loci

locus_list <- readRDS("cache/locus_list.Rds")
select_loci <- readRDS("cache/selected_loci.Rds")
class_probs <- readRDS("cache/class_probs.Rds")

# Create a data.frame with all important information
locus_table <- do.call(rbind, lapply(locus_list, function(ss) {
  data.frame(chr = attr(ss, "chr"),
             locus = attr(ss, "locus"),
             start_pos = attr(ss, "pos")["start"],
             end_pos = attr(ss, "pos")["end"],
             first_snp = attr(ss, "pos")["min"],
             last_snp = attr(ss, "pos")["max"])
}))

locus_table$selected <- seq_len(nrow(locus_table)) %in% select_loci
locus_table$confidence <- attr(class_probs,  "probabilities")[ , "TRUE"]

table(locus_table$chr, locus_table$selected)

dir.create("results", showWarnings = FALSE)
write.csv(locus_table, file = "results/selected_loci.csv", row.names = FALSE)

# Create latex tables
library(stargazer)
rows <- list(1:50, 51:100, 101:150, 151:200, 201:250, 250:266)
for (i in seq_along(rows)) {
  tex_table <- stargazer(locus_table[rows[[i]], ],
                         summary=FALSE, rownames = FALSE, float = FALSE)
  write(tex_table, file = paste0("results/selected_loci_", i, ".tex"))
}

