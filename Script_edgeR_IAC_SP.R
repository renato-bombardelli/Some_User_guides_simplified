### Preparo dos Dados ###
count_table <- read.table("NEW_Quantified_all_Ssci_reads_greater_than_50bp.txt", header = TRUE, row.names = "Geneid")
gene_length <- count_table[,5]
names(gene_length) <- rownames(count_table)

count_table <- count_table[, -(1:5)]


colnames(count_table) <- c("SPinoc_rep1", "SPinoc_rep2", "SPinoc_rep3",
                           "IACinoc_rep1", "IACinoc_rep2", "IACinoc_rep3")

condition <- c("Resistent", "Resistent", "Resistent", "Susceptible",
               "Susceptible", "Susceptible")

sample_info <- data.frame(condition = factor(condition))
rownames(sample_info) <- colnames(count_table)
group <- factor(paste(condition))

### Filtragem e seleção das amostras (todas) ###
library("edgeR")

y <- DGEList(counts = count_table, group = group)
cpms <- cpm(y)
keep <- rowSums(cpms > 1) >= 3
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

plotMDS(y, col = c("blue", "blue", "blue", "red", "red", "red"))

### Delineamento Experimental ###
design <- model.matrix(~ 0 + group, data = y$samples)
colnames(design) <- levels(y$samples$group)
y <- estimateDisp(y, design, robust = TRUE)
plotBCV(y)
y$design

### Teste da Expressão Diferencial - LRT(Likelihood ratio test) ###
fit <- glmFit(y, design, robust = TRUE)
# Podemos alterar o contraste a ser testado. Neste caso utiliza-se a interação
my.contrast <- makeContrasts("Resistent - Susceptible", levels = design)
lrt <- glmLRT(fit, contrast = my.contrast)
topTags(lrt)
DE <- decideTestsDGE(lrt, adjust.method = "none", p.value = 0.05)
summary(DE)
#plotSmear(lrt)

DE_tags <- rownames(y)[as.logical(DE)]
plotSmear(lrt, de.tags=DE_tags, cex = 0.3)
abline(h=c(-1, 1), col="blue")

write.table(file="DE_tags_cuadapt_greater_50bp.csv", DE_tags)
write.table(file="DEs_cuadapt_greater_50bp.csv", lrt$table)
write.table(file="DEs_cuadapt_greater_50bp_FDR.csv", topTags(lrt, n = 6000), sep='\t')
