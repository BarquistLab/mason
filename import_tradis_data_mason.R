
# get table from downloads (excel file) and modify, ending with table of ess genes of several strains.

library(readxl)
library(writexl)
library(dplyr)

# Read the input Excel file
tradis_data <- "~/Downloads/mbio.01798-24-s0002.xlsx"
tradis_da <- read_excel(tradis_data)

#retain only columns containing "TraDIS Essentiality:" or "Locus:" or "Keio gene name" in their name
tradis_d <- tradis_da %>%
  select(contains("Keio gene name"), contains("TraDIS Essentiality:"), contains("Locus:")) %>%
  # make all "NA" NA
    mutate(across(everything(), ~ na_if(., "NA")))


  # for all columns containing "TraDIS Essentiality:", replace values <=0 with "essential" and values >0 with "non-essential". keep NAs as they are
tradis_d_ess <- tradis_d %>%
    mutate(across(contains("TraDIS Essentiality:"), ~ case_when(. <= 0 ~ "essential",
                                      . > 0 ~ "non-essential",
                                      TRUE ~ as.character(.))))


# get each organism from headers
orgnames <- unique(gsub("TraDIS Essentiality: (.*)", "\\1", colnames(tradis_d_ess)[grepl("TraDIS Essentiality:", colnames(tradis_d_ess))]))


# ok make a table for each org with all their essenial genes, and keio names, and locus tags. then save as excel file with one sheet per org
output_list <- list()
for (org in orgnames) {
    org_table <- tradis_d_ess %>%
        select(contains(paste0("TraDIS Essentiality: ", org)), contains("Keio gene name"), paste0("Locus: ", org)) %>%
        rename(Essentiality = contains(paste0("TraDIS Essentiality: ", org))) %>%
        filter(Essentiality == "essential") %>%
        select(-Essentiality) %>%
        distinct()
    # change column names to keio_gene_name and locus_tag
    colnames(org_table) <- c("Keio_gene_name", "Locus_tag")
    output_list[[org]] <- org_table
  # if no keio gene name, fill with locus tag
    org_table <- org_table %>%
        mutate(Keio_gene_name = ifelse(is.na(Keio_gene_name), Locus_tag, Keio_gene_name))
  # save the table as excel file with sheet name as org
    write_xlsx(org_table, paste0("~/Downloads/", org, "_essential_genes.xlsx"))
}

# add one more for the keio list.
ecogene_list <- tradis_da %>%
    select(contains("Keio gene name"), contains("EcoGene Essentiality: "), contains("Locus: Escherichia coli BW25113 (Keio)")) %>%
   # filter for 1 in "EcoGene Essentiality: Escherichia coli BW25113" column
    filter(`EcoGene Essentiality: Escherichia coli BW25113` == 1) %>%
  # just get the keio gene name and locus tag, and rename them to keio_gene_name and locus_tag
    select(contains("Keio gene name"), contains("Locus: Escherichia coli BW25113 (Keio)")) %>%
    rename(Keio_gene_name = contains("Keio gene name"), Locus_tag = contains("Locus: Escherichia coli BW25113 (Keio)")) %>%
    distinct()

# ok now add ecogene list excel:
write_xlsx(ecogene_list, "~/Downloads/EcoGene_BW25113_essential_genes.xlsx")


