# Gestion des logs Snakemake
dir.create(dirname(snakemake@log[[1]]),
           recursive = TRUE, showWarnings = FALSE)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

#pour savoir le repertoire a installer la version rlang nécessaire


library(dada2);packageVersion("dada2")
library('dplyr')
library('ggplot2')
library('reshape2')
library('ShortRead')
library("Biostrings")

#importation des séquences trimmés de mon dossier

# Recuperation des fichiers depuis Snakemake

path <- file.path("results", "trimmed")

fnFs <- as.character(snakemake@input[["r1"]])
fnRs <- as.character(snakemake@input[["r2"]])

stopifnot(length(fnFs) > 0)
stopifnot(length(fnRs) > 0)
stopifnot(length(fnFs) == length(fnRs))

# Extract sample names
sample.namesF <- basename(fnFs)
sample.namesF <- sub("_R1.*", "", sample.namesF)
sample.namesF <- sub("_S[0-9]+$", "", sample.namesF)

sample.namesR <- basename(fnRs)
sample.namesR <- sub("_R2.*", "", sample.namesR)
sample.namesR <- sub("_S[0-9]+$", "", sample.namesR)

stopifnot(identical(sample.namesF, sample.namesR))
stopifnot(length(unique(sample.namesF)) == length(sample.namesF))
stopifnot(length(unique(sample.namesR)) == length(sample.namesR))

#filtrage qualité
############################################################

filt_path <- file.path(path, "filtered")
dir.create(filt_path, showWarnings = FALSE)

outdir <- file.path(path, "results_dada2")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

filtFs <- file.path(filt_path, paste0(sample.namesF, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.namesR, "_R_filt.fastq.gz"))
names(filtFs) <- sample.namesF
names(filtRs) <- sample.namesR

print(sample.namesF)
# Vérification finale unicité
stopifnot(length(unique(filtFs)) == length(filtFs))
stopifnot(length(unique(filtRs)) == length(filtRs))

  ############################################################
  # 4) Filtrage
  ############################################################
  
out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = as.integer(snakemake@params[["truncLen"]]),
  maxN = as.integer(snakemake@params[["maxN"]]),
  maxEE = as.numeric(snakemake@params[["maxEE"]]),
  truncQ = as.integer(snakemake@params[["truncQ"]]),
  rm.phix = as.logical(snakemake@params[["rm_phix"]]),
  compress = as.logical(snakemake@params[["compress"]]),
  multithread = snakemake@threads
)
  
saveRDS(out, file.path(outdir, "out.rds"))
  
  ############################################################
  # 5) Modélisation des erreurs
  ############################################################
  
errF <- learnErrors(filtFs, multithread = snakemake@threads)
errR <- learnErrors(filtRs, multithread = snakemake@threads)
  
saveRDS(errF, file.path(outdir, "errF.rds"))
saveRDS(errR, file.path(outdir, "errR.rds"))
  
pdf(file.path(outdir, "ErrorLearning.pdf"), width = 8, height = 10)
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()
  
  ############################################################
  # 6) Déréplication
  ############################################################
  
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
  
names(derepFs) <- sample.namesF
names(derepRs) <- sample.namesR
  
  ############################################################
  # 7) Inférence des ASVs
  ############################################################
  
dadaFs <- dada(derepFs, err = errF, multithread = snakemake@threads)
dadaRs <- dada(derepRs, err = errR, multithread = snakemake@threads)
  
saveRDS(dadaFs, file.path(outdir, "dadaFs.rds"))
saveRDS(dadaRs, file.path(outdir, "dadaRs.rds"))
  
  ############################################################
  # 8) Fusion paired-end
  ############################################################
  
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
saveRDS(mergers, file.path(outdir, "mergers.rds"))
  
merged_reads <- sapply(mergers, function(x) {
  if (is.null(x) || nrow(x) == 0) return(0)
  nrow(x)
})
  
cat("Nombre de reads mergés par échantillon:\n")
print(summary(merged_reads))
  
  ############################################################
  # 9) Table de séquences
  ############################################################
  
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, file.path(outdir, "seqtab.rds"))
  
  ############################################################
  # 10) Suppression des chimères
  ############################################################
  
seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = snakemake@threads,
  verbose = TRUE
)
  
saveRDS(seqtab.nochim, file.path(outdir, "seqtab.nochim.rds"))
  
seq_lengths <- nchar(getSequences(seqtab.nochim))
  
pdf(file.path(outdir, "ASV_length_distribution.pdf"))
hist(seq_lengths,
      breaks = 50,
      main = "Distribution des longueurs des ASVs",
      xlab = "Longueur des séquences (bp)")
dev.off()
  
cat("Résumé longueurs ASVs :\n")
print(summary(seq_lengths))
  
cat("Proportion non chimérique:\n")
print(sum(seqtab.nochim) / sum(seqtab))
  
  ############################################################
  # 11) Track des reads
  ############################################################
  
getN <- function(x) {
  if (is.null(x)) return(0)
  sum(getUniques(x))
}
  
getMergedN <- function(x) {
  if (is.null(x) || nrow(x) == 0) return(0)
  sum(x$abundance)
}
  
track <- cbind(
  input = out[,1],
  filtered = out[,2],
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged = sapply(mergers, getMergedN),
  nonchim = rowSums(seqtab.nochim)
)
  
rownames(track) <- sample.namesF
  
write.table(track, file.path(outdir, "track.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
retention_percent <- sweep(track, 1, track[, "input"], FUN = "/") * 100
write.csv(retention_percent, file.path(outdir, "retention_percent.csv"))
  
retention_df <- as.data.frame(retention_percent)
retention_df$Sample <- rownames(retention_df)
retention_melted <- melt(retention_df, id.vars = "Sample")
mean_retentions <- aggregate(value ~ variable, data = retention_melted, mean)
  
pdf(file.path(outdir, "retention_track.pdf"), width = 12, height = 10)
ggplot(mean_retentions, aes(x = variable, y = value)) +
  geom_bar(stat = "identity", fill = "#0072B2") +
  geom_text(aes(label = sprintf("%.1f%%", value)), vjust = -0.5) +
  ylim(0, 100) +
  theme_minimal() +
  labs(x = "Étape du pipeline",
       y = "Taux moyen de rétention (%)",
       title = "Taux moyen de rétention par étape") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"))
dev.off()
  


# Fermeture du fichier de log
sink(type = "message")
sink(type = "output")
close(log_file)
