library(rmelting)
library(Biostrings)
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = T)
print(args)
print("abcde")
fasta <- as.data.frame(readDNAStringSet(args))
colnames(fasta) <- "PNA_sequence"
fasta$Tm_ASOs <- sapply(fasta$PNA_sequence, function(x) melting(sequence = x, nucleic.acid.conc = 0.000008, hybridisation.type = "rnarna", Na.conc = 0.1)$Results$`Melting temperature (C)`)
fasta$ASO <- gsub(".*;[^_]*_", "", rownames(fasta))

g <- fasta %>% ggplot(aes(x=ASO, y=Tm_ASOs)) + geom_col(fill="steelblue") + theme_minimal() + 
  coord_cartesian(ylim=c(min(fasta$Tm_ASOs)-10, max(fasta$Tm_ASOs)+10)) +  ylab("Tm (°C)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 110/length(fasta$Tm_ASOs)),
        axis.text.y = element_text(size =14), axis.title.x = element_text(size =16),
        axis.title.y = element_text(size =16),
        panel.grid.major.x = element_blank()) 

sname <- paste0(gsub("reference_sequences/aso_sequences.fasta","", args), "outputs/tm.png")
                
print(sname)
ggsave(sname,g)
print(fasta)