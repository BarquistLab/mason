library(rmelting)
library(Biostrings)
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = T)

df_asos <- read.delim(args)  #read.delim("./Documents/mason/browser/pnag/static/data/188_2/outputs/result_table.tsv")
df_asos$melting_temperature <- sapply(df_asos$ASO_sequence, function(x) {
  melting(sequence = x, nucleic.acid.conc = 0.000008, hybridisation.type = "rnarna", Na.conc = 0.1)$Results$`Melting temperature (C)`
})
rownames(df_asos) <- df_asos$ASO_name
df_asos$ASO_name <- gsub(".*_(ASO.*)", "\\1", rownames(df_asos))


  
g <- df_asos %>% ggplot(aes(x=ASO_name, y=melting_temperature)) + geom_col(fill="steelblue") + theme_minimal() + 
  coord_cartesian(ylim=c(min(df_asos$melting_temperature)-5, max(df_asos$melting_temperature)+5)) +  ylab("Tm (°C)") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 80/length(df_asos$melting_temperature)),
        axis.text.y = element_text(size = 80/length(df_asos$melting_temperature)), 
        axis.title.y = element_text(size =12),
        panel.grid.major.x = element_blank()) 

df_asos$melting_temperature <- format(round(df_asos$melting_temperature, 2), nsmall = 2)
write.table(df_asos, args, sep = "\t", quote = F)

sname <- paste0(gsub("result_table.tsv","", args), "tm.png")
ggsave(sname,g)




