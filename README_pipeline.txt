#pipeline utilisado na análise dos dados do fungo SP_IAC 48 hai (Ex comandos com os parâmetros utilizados)

#FastqC - verificar qualidade das leituras pré-processadas

fastqc -t 2 "arquivo".fastq -o "output_directory"

#filtragem com cutadapt

-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 20 --max-n 0 -q 20 -j 2

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 20 --max-n 0 -q 20 -j 2 -o ./cutadapt/"$output1"_cut_PE1.fastq -p ./cutadapt/"$output1"_cut_PE2.fastq ./"$i" ./"$file2" > ./cutadapt/"$output1"_summary.txt

#FastqC - verificar a filtragem com o cutadapt

fastqc -t 2 "arquivo_processado".fastq -o "output_directory"

#Alinhamento com hisat2

##criação do index da referência

hisat2-build "Reference.fasta" "prefix"

hisat2-build ./References/Sscitamineum/NCBI_Ssci_genome.fa ./index_Sscitamineum_full_genome/Sscitamineum_full_genome.index

##Alinhamento contra a referência (index)

-p 2 --rg-id "Treatment" --rg SM:"Sample_identification" -x index -1 "pair-end1" -2 "pair-end2" -S "output.sam"

hisat2 -p 2 --rg-id SPinoc --rg SM:SPinoc_rep"$rep" --summary-file ./SPinoc/cutadapt/alignment/summary_align_SPinoc_rep"$rep"_cutadapt.txt -x "$index" -1 "$i" -2 "$file2" -S ./SPinoc/cutadapt/alignment/SPinoc_rep"$rep"_cutadapt_Ssci.sam

#ordenação (sort) e compactação (BAM) com samtools

samtools sort "output_hisat2".sam > "arquivo".bam

samtools index "arquivo".bam

#Preparo da anotação GTF com "gffred" do cufflinks
## Anotação -> NCBI_Ssci_genome.gff3

gffread NCBI_Ssci_genome.gff3 -T -o anotation.gtf

#Contagem das leituras com featureCounts (do pacote subread)

-s 0 -p -T 2 -t CDS -g gene_id -a "ref_anotation" -o "output" "paired-end1" "paired-end2" ### muito importante adicionar a opção '-p' para indicar que são reads paired-end

featureCounts -s 0 -p -T 2 -t CDS -g gene_id -a ./References/Sscitamineum/anotation_Ssci.gtf -o ./counts_cutadapt/Quantified_all_samples.txt $(ls ./SPinoc/cutadapt/alignment/*.bam) $(ls ./IACinoc/cutadapt/alignment/*.bam)

##Contagem das leituras alinhadas ao transcriptoma (comps, GGs) da cana

#######genereting a GTF file from the comps_plus_GGs.fas sugarcane reference
#/media/renato/Backup Files/References/sugarcane/Generating_GTF/annotation_sugarcane.py ##python script for that
#annotation.gtf -> containing the transcripts as if it were chromossomes, with start position 1 and end position the length of the transcript, no strand (.) and transcript_id and gene_id the name of the transcripts _ exon as meta-feature

-s 0 -p -T 2 -t exon -g gene_id -a "ref_anotation" -o "output" "paired-end1" "paired-end2"

featureCounts -s 0 -p -T 2 -t exon -g gene_id -a ./References/sugarcane/Generating_GTF/annotation_sugarcane.gtf -o ./counts_sugarcane/Quantified_all_samples.txt $(ls ./SPinoc/cutadapt/alignment_sugarcane/*.bam) $(ls ./SPcontrol/cutadapt/alignment_sugarcane/*.bam) $(ls ./IACinoc/cutadapt/alignment_sugarcane/*.bam) $(ls ./IACcontrol/cutadapt/alignment_sugarcane/*.bam)

######################################################################################################################################################################################################################

#detecção de genes diferencilamente expressos (DEs) com edgeR (script.R)

Ex: Comparação SP vc IAC (inoculadas - RNA-Seq fungo)

### Preparo dos Dados ###
count_table <- read.table("Quantified_all_samples.txt", header = TRUE, row.names = "Geneid")
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

#escreve arquivos, necessita processamento no Excel
write.table(file="DE_tags.csv", DE_tags)
write.table(file="DEs.csv", lrt$table)

##########################################################################################################################################################################

#After detecting the diferencial expres genes (DE) we need to get the fasta file containing only the DEs for the GO terms enrichment (Blast2GO)

##/media/renato/Backup Files/Getting_DE_sequence_fasta/DEs_sequence_only_cutadapt.py



