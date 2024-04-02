# packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))


# main
main = function(merge_file, output_dir = "data/processed/") {
  
  # process sequences
  nanobodies = load_merged_reads(merge_file) %>% 
    filter_reads_by_primers(F_seq = "GAGGTGCAGCTG", R_seq = "ACCGTGAGCAGC") %>%
    translate_reads() %>%
    filter_translations_by_library_design(pattern = "EVQLV[^X\\*]{78,}VTVSS") %>%
    convert_StringSet_to_df() %>%
    mutate(length = nchar(sequence)) %>%
    filter(length == 118)
  
  # calculate AA frequencies
  residues = c(27:33, 55, 97:103, 105:106)
  all_counts = as.data.frame(t(consensusMatrix(AAStringSet(nanobodies$sequence))), row.names = 1:118)
  variable_position_counts = all_counts[residues, ]
  variable_position_proportions = round(variable_position_counts / nrow(nanobodies), 5)
  
  # write data
  prefix = gsub("(.*)NGmerge.fastq.gz", "\\1", basename(merge_file))
  write.csv(variable_position_counts, paste0(output_dir, prefix, "counts.csv"))
  write.csv(variable_position_proportions, paste0(output_dir, prefix, "proportions.csv"))
  write.csv(all_counts, paste0(output_dir, prefix, "counts_all_positions.csv"))
  return(invisible())
} 


# load merged fastq file
load_merged_reads = function(file, format = "fastq") {
  reads = readDNAStringSet(filepath = file, use.names = FALSE, format = format)
  return(reads)
}


# filter DNAStringSet for F and R sequence matching
filter_reads_by_primers = function(reads, F_seq = "", R_seq = "") {
  
  # find reads that match primers (either forward or reverse complement)
  matches = c(reads[vcountPattern(F_seq, reads) > 0 & vcountPattern(R_seq, reads) > 0],
              reverseComplement(reads)[vcountPattern(F_seq, reverseComplement(reads)) > 0 & vcountPattern(R_seq, reverseComplement(reads)) > 0])
  
  # remove any flanking characters outside primers
  matches = DNAStringSet(sub(paste0(".*(", F_seq, ".*", R_seq, ").*"),  "\\1", matches))
  return(matches)
}


# translate DNAstringSet
translate_reads = function(reads) {
  translations = translate(reads, no.init.codon = TRUE, if.fuzzy.codon = "solve")
  return(translations)
}


# filter AAStringSet for library-specific regex
filter_translations_by_library_design = function(translations, pattern = "") {
  translations = translations[grepl(pattern, translations)]
  return(translations)
}


# convert StringSet to df
convert_StringSet_to_df = function(stringset) {
  sequence_df = as_tibble(as.character(stringset))
  names(sequence_df) = "sequence"
  return(sequence_df)
}