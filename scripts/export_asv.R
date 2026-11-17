# Charger la table des ASV sans chimeres
x <- readRDS(snakemake@input[["asv"]])

# Creer les identifiants des ASV
ids <- paste0("ASV_", seq_len(ncol(x)))

# Recuperer les sequences nucleotidiques
seqs <- colnames(x)

# Transposer la table d'abondance
asv_table <- t(x)
rownames(asv_table) <- ids

# Exporter la table d'abondance au format CSV
write.csv(
  asv_table,
  snakemake@output[["table"]]
)

# Exporter les sequences au format FASTA
writeLines(
  as.vector(rbind(paste0(">", ids), seqs)),
  snakemake@output[["fasta"]]
)

# Afficher un resume
cat("Nombre d'echantillons :", nrow(x), "\n")
cat("Nombre d'ASV :", ncol(x), "\n")
cat("Exportation terminee avec succes.\n")
