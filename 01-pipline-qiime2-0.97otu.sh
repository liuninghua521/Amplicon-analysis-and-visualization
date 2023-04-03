
source activate qiime2-2021.2  

#### 1. Data Import####
time qiime tools import  \
  --type 'SampleData[PairedEndSequencesWithQuality]'  \
  --input-path sample-metadata.tsv \
  --output-path demux.qza  \
  --input-format PairedEndFastqManifestPhred33V2

time qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

#### 2. Generate feature tables and representative sequences####
##2.1
time qiime dada2 denoise-paired \
--i-demultiplexed-seqs demux.qza \
--p-n-threads 0 \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 0 \
--p-trunc-len-r 0 \
--o-table table.qza \
--o-representative-sequences rep-seqs.qza \
--o-denoising-stats denoising-stats.qza

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file sample-metadata.tsv

##2.2 
time qiime vsearch cluster-features-de-novo \
  --i-table table.qza \
  --i-sequences rep-seqs.qza \
  --p-perc-identity 0.97 \
  --o-clustered-table table-dn-97.qza \
  --o-clustered-sequences rep-seqs-dn-97.qza

qiime feature-table summarize \
  --i-table table-dn-97.qza \
  --o-visualization table97.qzv \
  --m-sample-metadata-file sample-metadata.tsv

####3 Resampling of feature tables and filtering of low abundance OTUs####
##3.1
qiime tools export \
  --input-path table-dn-97.qza \
  --output-path exported-feature-table

biom convert -i exported-feature-table/feature-table.biom \
    -o exported-feature-table/feature-table97.txt \
    --to-tsv
    
sed -i '1d' exported-feature-table/feature-table97.txt
sed -i 's/#//' exported-feature-table/feature-table97.txt

##3.2
# ======================== R =========================
library(vegan)
library(tidyverse)
path1 = '/home/LNH/lnh-data/seamount_data/K_bacteria'

otu = read.table(paste0(path1, "/exported-feature-table/feature-table97.txt"),
                 header=T, sep="\t", quote = "", row.names=1, comment.char="#",stringsAsFactors = FALSE)

otu_rare1 = otu
sum(otu)
abundance=0.0001
otu_rare1 <- otu[rowSums(otu)/sum(otu)>=(abundance/100),]
colSums(otu_rare1)
otu_rare = as.data.frame(t(rrarefy(t(otu_rare1), min(colSums(otu_rare1)))))
colSums(otu_rare)

write.table(otu_rare, paste0(path1, "/exported-feature-table/rare-feature-table97.txt"),
            sep = "\t", row.names = T, col.names = T, quote = FALSE)
# ======================== R end=========================

##3.3
# Add OTU_ID and a tab \t at the beginning of the first line。
vim exported-feature-table/rare-feature-table97.txt

biom convert -i exported-feature-table/rare-feature-table97.txt \
  -o exported-feature-table/rare-table97.biom \
  --table-type="OTU table" --to-hdf5

qiime tools import \
  --input-path exported-feature-table/rare-table97.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path rare-feature-table97.qza

qiime feature-table summarize \
  --i-table rare-feature-table97.qza \
  --o-visualization rare-table97.qzv \
  --m-sample-metadata-file sample-metadata.tsv

qiime feature-table filter-seqs \
  --i-data rep-seqs-dn-97.qza \
  --i-table rare-feature-table97.qza \
  --o-filtered-data rare-rep-seqs-dn-97.qza

##3.4
# change ID
mkdir change_id

qiime tools export \
  --input-path rare-feature-table97.qza \
  --output-path change_id

biom convert -i change_id/feature-table.biom \
    -o change_id/rare-feature-table97.txt \
    --to-tsv
    
sed -i '1d' change_id/rare-feature-table97.txt
sed -i 's/#OTU ID/otuID/' change_id/rare-feature-table97.txt

qiime tools export \
  --input-path rare-rep-seqs-dn-97.qza \
  --output-path change_id

seqkit fx2tab change_id/dna-sequences.fasta > change_id/ref.seq

# ========================R=========================

otu2 <- paste0(path1, "/change_id/rare-feature-table97.txt") %>%
  read.delim(check.names = FALSE, header = T,sep="\t")

otu2$OTUID <- paste0("OTU",seq_len(nrow(otu2)))

otu2 <- select(otu2, OTUID, everything())

rep <- paste0(path1, "/change_id/ref.seq") %>%
  read.delim(check.names = FALSE, sep="\t", header = F) %>%
  select(-V3)
colnames(rep) <- c("ID", "SEQ")

bind_data <- otu2 %>% left_join(rep , by=c("otuID" = 'ID'))

otu.final <- select(bind_data, -otuID,  -SEQ)
seq.final <- select(bind_data, OTUID, SEQ)

write.table(otu.final, paste0(path1, "/change_id/otu.final.table"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(seq.final, paste0(path1, "/change_id/final.seqs"), sep = "\t", quote = F, col.names = F, row.names = F)
# ========================R end =========================

seqkit tab2fx change_id/final.seqs > change_id/seq.final.seq
rm -rf rare-feature-table97.qza rare-rep-seqs-dn-97.qza

biom convert -i change_id/otu.final.table \
  -o change_id/rare-table97.biom \
  --table-type="OTU table" --to-hdf5

qiime tools import \
  --input-path change_id/rare-table97.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path rare-feature-table97.qza

time qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path change_id/seq.final.seq \
  --output-path rare-rep-seqs-dn-97.qza

qiime feature-table summarize \
  --i-table rare-feature-table97.qza \
  --o-visualization rare-table97-2.qzv \
  --m-sample-metadata-file sample-metadata.tsv
  
time qiime feature-table tabulate-seqs \
  --i-data rare-rep-seqs-dn-97.qza \
  --o-visualization rare-rep-seqs.qzv


#### 4. Alpha diversity analysis and species annotation####
##4.1
time qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rare-rep-seqs-dn-97.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
  
##4.2
time qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table rare-feature-table97.qza \
  --p-sampling-depth 60734 \
  --m-metadata-file sample-metadata.tsv \
  --output-dir core-metrics-results

time qiime diversity alpha-rarefaction \
  --i-table rare-feature-table97.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 60734 \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization alpha-rarefaction.qzv

##4.4
time qiime feature-classifier classify-sklearn \
	  --i-classifier /home/LNH/lnh-data/amplion_classfier/classfier/bacteria_classifier_gg_13_8_97_V4-V5.qza \
	  --i-reads rare-rep-seqs-dn-97.qza \
	  --o-classification taxonomy.qza

	qiime metadata tabulate \
	  --m-input-file taxonomy.qza \
	  --o-visualization taxonomy.qzv

	qiime taxa barplot \
	  --i-table rare-feature-table97.qza \
	  --i-taxonomy taxonomy.qza \
	  --m-metadata-file sample-metadata.tsv \
	  --o-visualization taxa-bar-plots.qzv

