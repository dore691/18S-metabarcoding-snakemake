# Assignation taxonomique PR2 avec Snakemake

# 1. Gestion des logs
dir.create(dirname(snakemake@log[[1]]),
           recursive = TRUE, showWarnings = FALSE)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# 2. Charger DADA2
library(dada2)

# 3. Recuperer les fichiers depuis Snakemake
seqtab.nochim <- readRDS(snakemake@input[["asv"]])
database <- snakemake@input[["database"]]

print(dim(seqtab.nochim))

# 4. Assignation taxonomique
taxa_pr2 <- assignTaxonomy(
    seqtab.nochim,
    database,
    taxLevels = as.character(
        snakemake@params[["taxLevels"]]
    ),
    minBoot = as.numeric(
        snakemake@params[["minBoot"]]
    ),
    multithread = snakemake@threads
)

# 5. Sauvegarder les resultats
saveRDS(
    taxa_pr2,
    snakemake@output[["taxonomy"]]
)

print("Assignation taxonomique terminee")

# 6. Fermer le fichier de log
sink(type = "message")
sink(type = "output")
close(log_file)
