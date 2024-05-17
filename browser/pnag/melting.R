library(rmelting)
print("rmalting loaded")
library(dplyr)
print("dplyr loaded")
library(ggplot2)
print("ggplot2 loaded")

args <- commandArgs(trailingOnly = T)

print("Run RSCRIPT MELTING")

df_asos <- read.delim(args)  #read.delim("./Documents/mason/browser/pnag/static/data/188_2/outputs/result_table.tsv")
df_asos$Tm <- sapply(df_asos$ASO_seq, function(x) {
  melting(sequence = x, nucleic.acid.conc = 0.000008, hybridisation.type = "rnarna", Na.conc = 0.1)$Results$`Melting temperature (C)`
})
rownames(df_asos) <- df_asos$ASO
df_asos$ASO <- gsub(".*_(ASO.*)", "\\1", rownames(df_asos))

  
g <- df_asos %>% ggplot(aes(x=ASO, y=Tm)) + geom_col(fill="steelblue") + theme_minimal() +
  coord_cartesian(ylim=c(min(df_asos$Tm)-5, max(df_asos$Tm)+5)) +  ylab("Tm (°C)") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size =12),
        panel.grid.major.x = element_blank()) 

# also add molecular weight calculator:
# Define a function to calculate the molecular weight of an n-mer of nucleobases connected with a peptide bond. The
# peptid bond contains 9 H, 2 N, 2 O and 4 C atoms.
calculate_pna_mw <- function(seq) {
  # Define the molecular weight of each nucleotide (nucleobase only)
  mw_nucleotide <- c(A=135.13, T=125.06, G=152.12, C=111.07)
  # Define the molecular weight of the peptide bond
  mw_peptidebond <- 9*1.01 + 2*14.01 + 2*16.00 + 4*12.01
  # Split the DNA sequence into individual nucleotides
  nucleotides <- strsplit(seq, "")[[1]]
  # Calculate the molecular weight of the sequence, accounting for the peptide bond
  mw <- sum(mw_nucleotide[nucleotides]) + mw_peptidebond * (length(nucleotides) - 1) + 1.01*2
  return(mw)
}

df_asos$Mw <- sapply(df_asos$ASO_seq, calculate_pna_mw)


df_asos$Tm <- format(round(df_asos$Tm, 2), nsmall = 2)
df_asos$Mw <- format(round(df_asos$Mw, 2), nsmall = 0)
write.table(df_asos, args, sep = "\t", quote = F)

sname <- paste0(gsub("result_table.tsv","", args), "tm.png")
ggsave(sname,g)




