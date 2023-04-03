suppressPackageStartupMessages({
  library(MicrobiotaProcess) 
  library(phyloseq) 
  library(ggplot2) 
  library(plyr)
  library(ggpubr)
  library(tidyverse) 
  library(vegan) 
  library(reshape2) 
  library(ggnewscale) 
  library(paletteer) 
  library(ggtree)
  library(aplot)
  library(aplot)
  library(ggplotify)
  library(microeco)
  library(file2meco)
  library(magrittr)
  library(mecodev)
})

####0. Read data####
ps <- import_qiime2(otuqza = "rare-feature-table97.qza", 
                    taxaqza = "taxonomy.qza", 
                    mapfilename = "sample-metadata.tsv",
                    #refseqqza = "rare-rep-seqs-dn-97.qza", 
                    treeqza = "rooted-tree.qza"
                    )
#sample_order
sample_order <- c("C1.0",	"C1.30",	"C1.50",	"C1.DCM",	"C1.200",	
                  "C1.300",	"C1.500",	"C1.1000",	"C1.2300",	"C2.0",
                  "C2.30",	"C2.50",	"C2.DCM",	"C2.200",	"C2.300",	
                  "C2.500",	"C2.1000",	"C2.1900",	"C3.0",	"C3.30",
                  "C3.50",	"C3.DCM",	"C3.200",	"C3.300",	"C3.500",	
                  "C4.0",	"C4.30",	"C4.50",	"C4.90",	"C5.0",	"C5.30",
                  "C5.50",	"C5.DCM",	"C5.200",	"C5.300",	"C5.600",	
                  "C6.0",	"C6.30",	"C6.50",	"C6.DCM",	"C6.200",	
                  "C6.300",	"C6.500",	"C6.1000",	"C6.1450",	"C7.0",
                  "C7.30",	"C7.50",	"C7.DCM",	"C7.200",	"C7.300",	
                  "C7.500",	"C7.1000",	"C7.2000",	"C7.2400")
# depth_level
depth_level <- c("0",	"30",	"50",	"DCM",	"200",	"300",	"500",	"1000+",	"2000+")
# site_level
site_level <- c("C1",	"C2",	"C3",	"C4",	"C5",	"C6",	"C7")


#====================================================================#
####1. rarefaction curves####
set.seed(1024)
p_rare1 <- ggrarecurve(obj = ps,
                      indexNames = c("Observe","Chao1","ACE"),
                      chunks = 100, 
                      linesize = 1,
                      factorLevels = list()
                      ) +  
  theme_bw() +
  theme(axis.text = element_text(size = 8), 
        panel.grid = element_blank(),
        strip.background = element_rect(colour = NA,fill="grey"),
        strip.text.x = element_text(face = "bold")) 
p_rare1

p_rare1 + scale_color_paletteer_d("ggsci::default_igv")


#====================================================================#
####2. alpha diversity####
set.seed(1024)
alphaobj <- get_alphaindex(ps)

## line chart##
colnames(as.data.frame(alphaobj))[1:6]

ala <- alphaobj@alpha
ala$sample <- rownames(ala)
ala$sample <- factor(ala$sample, levels = sample_order)
ala$group <- rep("bacteria", nrow(ala))

# Diversity
ala_zhe1 <- ggplot(ala, aes(x=sample, y=Shannon, group=group)) + 
  geom_point(size=4,shape=20) +
  geom_line(size = 0.8) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14, color = "black")
  ) +
  ylab("Shannon-wiener index")
ala_zhe1

# Richness
ala_zhe2 <- ggplot(ala, aes(x=sample, y=Observe, group=group)) + 
  geom_point(size=4,shape=20) +
  geom_line(size = 0.8) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14, color = "black")
  ) +
  ylab("Observe OTU")
ala_zhe2

# Evenness
ala_zhe3 <- ggplot(ala, aes(x=sample, y=Pielou, group=group)) + 
  geom_point(size=4,shape=20) +
  geom_line(size = 0.8) +
  theme_bw() +
  xlab("samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, color = "black"),
        axis.title = element_text(size = 14, color = "black")) +
  ylab("Pielou’s evenness index")
ala_zhe3

ala_zhe <- ala_zhe1 %>%
  insert_bottom(ala_zhe2) %>%
  insert_bottom(ala_zhe3)
ala_zhe

## Violin chart
alphaobj@sampleda$site <- factor(alphaobj@sampleda$site, levels = site_level)

p_alpha1 <- ggbox(alphaobj, 
                  geom = "violin", 
                  factorNames = "site",
                  indexNames = c("Observe","Shannon", "Pielou"),
                  compare = TRUE, 
                  p_textsize = 4, 
                  boxwidth = 0.2, 
                  factorLevels = site_level,
                  signifmap = TRUE, 
                  comparelist=list(c("C1", "C7"), c("C2", "C6"), c("C3", "C5")
                  )
                  ) + 
  theme(strip.background = element_rect(colour = NA, fill = "grey")) 
  
p_alpha1 + scale_fill_paletteer_d("ggsci::nrc_npg")

#====================================================================#
####3. Species composition analysis####

phytax <- get_taxadf(obj = ps,
                     taxlevel = 2, 
                     type = "species" 
                     )

phytax@sam_data$depth <- factor(phytax@sam_data$depth, levels = depth_level)
phytax@sam_data$site <- factor(phytax@sam_data$site, levels = site_level)

phybar <- ggbartax(obj = phytax,
                   position = "stack",
                   topn = 9, 
                   width = 0.9, 
                   count = FALSE, 
                   facetNames = "depth", 
                   sampleLevels = sample_order, 
                   plotgroup = FALSE, 
                   ) +
  xlab(NULL) + 
  ylab("Relative abundance (%)") +
  guides(fill=guide_legend(override.aes = list(size=10))) +
  scale_fill_manual(values = paletteer_d("ggsci::category10_d3", 10),
                    breaks=c("p__Proteobacteria","p__Cyanobacteria", "p__Bacteroidetes",
                             "p__Crenarchaeota", "p__Actinobacteria","p__SAR406",
                             "p__Planctomycetes","p__Euryarchaeota","p__Chloroflexi",
                              "Others"),
                    labels=c("Proteobacteria","Cyanobacteria", "Bacteroidetes",
                             "Crenarchaeota", "Actinobacteria", "SAR406",
                             "Planctomycetes", "Euryarchaeota", "Chloroflexi", "Others")
  ) +
  scale_fill_paletteer_d("ggsci::category10_d3") +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "#FFFFF0", colour = "black"),
    strip.text = element_text(size = 14,colour = "black"),
    axis.text.x = element_text(angle = -90, hjust = 1, size = 12,colour = "black"),
    axis.text.y = element_text(size = 12,colour = "black"),
    axis.title = element_text(size = 16,colour = "black"),
    legend.text = element_text(size = 14,colour = "black")
    
  ) 

phybar


#====================================================================#
####4. veen chart####

D_DCM <- subset_samples(ps, depth %in% c("DCM"))

vennlist <- get_vennlist(obj = D_DCM, factorNames = "site")

library(VennDiagram)
vennp <- venn.diagram(vennlist,
                      height = 5,
                      width = 5, 
                      filename = NULL, 
                      alpha = 0.85,
                      col = "white", 
                      #fill = c("#7FC97F", "#386CB0", "#BEAED4", "#FDC086", "#FFFF99"),
                      #          A3          A5         A7         B3          B6
                      fill = brewer.pal(7, "Accent"),
                      fontfamily = "serif",
                      fontface = "bold",
                      cex = 1.3,
                      cat.col= "black",
                      cat.fontface = "bold",
                      cat.fontfamily = "serif",
                      cat.cex = 1.3,
                      cat.default.pos = "outer",
                      #cat.dist = 0.1,
                      margin = 0.05, 
                      lwd = 2, 
                      lty ='dashed', 
                      imagetype = "svg")
grid::grid.draw(vennp)


####5 share OTU####

new1 <- data.frame(C1=0, C2=0, C3=0, C4=0, C5=0, C6=0, C7=0, unique=0)
new1 <- new1[-1,]


for (i in 1:length(depth_level)){
  otu.sub <- subset_samples(ps, depth %in% depth_level[i])
  vennlist.sub <- get_vennlist(obj = otu.sub, factorNames = "site")
  sub.d <- as.data.frame(t((lengths(vennlist.sub))))
  rownames(sub.d) <- depth_level[i]
  sub.d$unique <- length(Reduce(intersect, vennlist.sub))
  new1 <- rbind.fill(new1, sub.d)
  
}

new1$site <- c("0",	"30",	"50",	"DCM",	"200",	"300",	"500",	"1000",	"2000")
new1 <- column_to_rownames(.data = new1, var = "site")

bili.d <- new1
bili.d <- apply(bili.d, 2, function(x) bili.d[, ncol(bili.d)] /x)
bili.d <- as.data.frame(bili.d)

row_mean = as.data.frame(apply(bili.d[, 1:(ncol(bili.d)-1)] ,1, mean))
colnames(row_mean) <- "bili"
bili.d <- cbind(bili.d,  row_mean)

bili.d2 <- rownames_to_column(bili.d, var = "depth")
bili.d2[bili.d2 == 'DCM'] <- 100
bili.d2$depth <- as.numeric(bili.d2$depth)
bili.d2 <- melt(bili.d2,
                id.vars = c("depth", "unique", "bili"),
                value.name = "bili1", 
                variable.name = "site")

p2 <- ggplot(data = bili.d2, 
             mapping = aes(depth, bili1),) +
  geom_point(
    color = "grey",
    size = 4,
    na.rm = T) + 
  geom_smooth(method = 'lm', 
              color = "#ee6d70",
              fill = "#358bac",
              na.rm = T) +
  stat_cor(method = "pearson",size = 5,
           na.rm = T) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 14)) +
  labs(x = "depth (m)", y = "percentage of shared OTUs (%)") +
  scale_y_continuous(labels = function(x) x * 100)


p2 


####6. PCA####
pcares <- get_pca(obj=ps, method="hellinger")

pcares@sampleda$depth <- factor(pcares@sampleda$depth, levels = depth_level)

pcaplot <- ggordpoint(obj=pcares, 
                      pc = c(1, 2), 
                      factorNames=c("depth"),
                      #showsample = TRUE, 
                      poinsize = 6, 
                      stroke = 0.1, 
                      biplot=F, 
                      speciesannot=TRUE,
                      # linesize = 0.5, 
                      # arrowsize = 1.5,
                      # arrowlinecolour = "grey",
                      # textlinesize = 0.02,
                      ellipse=TRUE,
                      ellipse_alpha = 1, 
                      ellipse_pro = 0.9,
                      #ellipse_linewd = 1.5, 
                      #ellipse_lty = 2, 
                      fontface = "bold.italic",
                      fontfamily = "sans"
)  + 
  guides(fill=guide_legend(override.aes = list(size=8)))+
  theme(plot.title = element_blank(),
        legend.position = "right",
        #strip.background = element_rect(fill = "#FFFFF0", colour = "black"),
        #strip.text = element_text(size = 14,colour = "black"),
        #axis.text.x = element_text(angle = -90, hjust = 1, size = 12,colour = "black"),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 14,colour = "black"),
        legend.text = element_text(size = 14,colour = "black"),
        legend.title = element_text(size = 14,colour = "black", face = "bold")
  ) +
  scale_fill_paletteer_d("ggsci::category10_d3") +
  scale_color_paletteer_d("ggsci::category10_d3")

pcaplot 




####7. PCOA####
pcoares <- get_pcoa(obj=ps, 
                    distmethod="bray", 
                    method="hellinger")

pcoares@sampleda$depth <- factor(pcoares@sampleda$depth, levels = depth_level)

pcoaplot <- ggordpoint(obj=pcoares, 
                       pc = c(1, 2), 
                       factorNames=c("depth"),
                       #showsample = TRUE,
                       poinsize = 5, 
                       stroke = 0.1, 
                       biplot=TRUE, 
                       speciesannot=TRUE,
                       # linesize = 0.5, 
                       # arrowsize = 1.5,
                       # arrowlinecolour = "grey",
                       # textlinesize = 0.02,
                       ellipse=TRUE,
                       ellipse_alpha = 1, 
                       ellipse_pro = 0.9, 
                       #ellipse_linewd = 1.5, 
                       #ellipse_lty = 2, 
                       fontface = "bold.italic",
                       fontfamily = "sans"
                       ) +
  scale_fill_paletteer_d("ggthemes::Tableau_10")
pcoaplot


####8. UPGMA ####   

hcsample <- get_clust(obj=ps, distmethod="bray",
                      method="hellinger", hclustmethod="average")
hcsample@sampleda$depth <- factor(hcsample@sampleda$depth, levels = depth_level)

cplot1 <- ggclust(obj=hcsample,
                  layout = "rectangular",
                  pointsize=4,
                  fontsize=3,
                  factorNames=c("depth")
                  ) +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  theme_tree2(legend.position="right",
              plot.title = element_text(face="bold", lineheight=25,hjust=0.5))

cplot1 

phybar1
phybars <- phybar1 + 
  #scale_y_continuous(expand=c(0,5)) +
  coord_flip() + 
  theme_void() 
phybars 

julei <- phybars %>% insert_left(cplot1+ xlim(0, 0.48), width = 0.3)
julei



####9. Community Assembly ####
dataset <- qiime2meco(ASV_data = "rare-feature-table97.qza",
                      taxonomy_data = "taxonomy.qza", 
                      phylo_tree = "rooted-tree.qza",  
                      #rep_fasta = "rare-rep-seqs-dn-97.qza",
                      sample_data = "sample-metadata2.tsv",
                      auto_tidy = TRUE)

dataset$tidy_dataset()
dataset

env_data_16S <- read.table("yingyangyan3.txt", header = T, sep = "\t", quote = "F")
env_data_16S <- column_to_rownames(env_data_16S, "X")


t1 <- trans_nullmodel$new(dataset, filter_thres = 0.0001, add_data = env_data_16S)
t1$cal_ses_betamntd(runs=500, abundance.weighted = TRUE)
dataset$beta_diversity[["betaNTI"]] <- t1$res_ses_betamntd

t2 <- trans_beta$new(dataset = dataset, group = "group", measure = "betaNTI")
t2$cal_group_distance()
g2 <- t2$plot_group_distance(distance_pair_stat = F,
                             #plot_group_order = depth_level,
                             ) + 
  geom_hline(yintercept = -2, linetype = 2) + 
  geom_hline(yintercept = 2, linetype = 2) +
  scale_color_paletteer_d("ggsci::default_igv")
g2



##RCbray （Bray-Curtis-based Raup-Crick）
t1$cal_rcbray(runs = 1000)
t1$cal_process(use_betamntd = TRUE)
t1$res_process

