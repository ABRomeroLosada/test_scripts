## Encontrar secuencias de transcritos y construir el indice
## kallisto index -i ostreococcus_tauri.idx ostreococcus_tauri_transcripts.fasta

## allisto quant -i ostreococcus_tauri.idx -b 100 -o kallisto_out --genomebam --gtf ostreococcus_tauri.gtf --chromosomes chrom_ostreococcus_tauri.txt --single -l 200 -s 10 sd_zt20_7.fq.gz

library(seqinr)

microalgae <- c("bprasinos", "creinhardtii", "csubellipsoidea",
                "czofingiensis", "dsalina", "hlacustris", "knitens",
                "mendlicherianum", "mpusilla", "ngaditana",
                "ptricornutum","smuscicola", "vcarteri")
#i<-6
for (i in 1:length(microalgae))
{
  genome.data <- read.fasta(file = paste("./genomes/",microalgae[i],".fa.gz", sep=""),seqtype = "DNA")
  genome.seqs <- getSequence(genome.data)
  chrom.data <- data.frame(getName(genome.data),sapply(X = genome.seqs,FUN = length))
  write.table(x = chrom.data,file =paste("./genomes/chrom_",microalgae[i],".txt", sep=""),
              quote = F,row.names = F,col.names = F, sep = "\t")
  
}

ostta.genome.data <- read.fasta(file = "ostreococcus_tauri.fa.gz",seqtype = "DNA")
ostta.genome.seqs <- getSequence(ostta.genome.data)

chrom.data <- data.frame(getName(ostta.genome.data),sapply(X = ostta.genome.seqs,FUN = length))
write.table(x = chrom.data,file = "chrom_ostreococcus_tauri.txt",quote = F,row.names = F,col.names = F, sep = "\t")


chrom.data <- data.frame(getName(ostta.genome.data),sapply(X = ostta.genome.seqs,FUN = length))
write.table(x = chrom.data,file = "chrom_ostreococcus_tauri.txt",quote = F,row.names = F,col.names = F, sep = "\t")
