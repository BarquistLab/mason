



library(rmelting)
print("rmalting loaded")

library(dplyr)
print("dplyr loaded")
library(ggplot2)
print("ggplot2 loaded")

# Source shared utilities
source("./pnag/r_utils.R")

args <- commandArgs(trailingOnly = T)

print("Run RSCRIPT MELTING")

df_asos <- read.delim(args)  #read.delim("./Documents/mason/browser/pnag/static/data/188_2/outputs/result_table.tsv")
df_asos$Tm <- sapply(df_asos$ASO_seq, function(x) {
  melting(sequence = x, nucleic.acid.conc = 0.000008, hybridisation.type = "rnarna", Na.conc = 0.1)$Results$`Melting temperature (C)`
})
rownames(df_asos) <- df_asos$ASO
df_asos$ASO <- gsub(".*_(ASO.*)", "\\1", rownames(df_asos))


g <- df_asos %>% ggplot(aes(x = ASO, y = Tm)) +
  geom_bar(stat = "identity", fill = "#31688e") +
  ggtitle("Predicted melting temperature (Tm)") +
  labs(x = "ASO sequence", y = "Tm (\u00B0C)") +
  theme_classic() +
  coord_cartesian(ylim = c(min(df_asos$Tm) - 5, max(df_asos$Tm) + 5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, face = "bold")) +
  geom_text(aes(label = round(Tm, 1)), vjust = -0.25, size = 7)

df_asos$Mw <- sapply(df_asos$ASO_seq, calculate_pna_mw)


df_asos$Tm <- format(round(df_asos$Tm, 2), nsmall = 2)
df_asos$Mw <- format(round(df_asos$Mw, 2), nsmall = 0)
write.table(df_asos, args, sep = "\t", quote = F)

wplot <- nrow(df_asos) + 5
sname <- paste0(gsub("result_table.tsv", "", args), "tm.png")
ggsave(sname, g, width = wplot, limitsize = FALSE)
sname_svg <- paste0(gsub("result_table.tsv", "", args), "tm.svg")
ggsave(sname_svg, g, width = wplot, limitsize = FALSE)
