######## Packages required + custom functions
Packages <- c("phyloseq", "data.table", "ggplot2", "plyr","dplyr","reshape2","grid",
              "gridExtra","scales","dplyr", "ggpubr","vegan","multcompView","rcompanion","betapart")

lapply(Packages, library, character.only = TRUE)

########################
# Dataset constitution #
########################
load("Fungi_Data.RData")
data_Fung$Month.Inoc <- factor(data_Fung$Month.Inoc, levels = c("Jul.DR", "Jul.DS","Oct.DR", "Oct.DS", "Dec.DR","Dec.DS", "Feb.DR","Feb.DS"))
data_Fung$Month.ent<-factor(data_Fung$Month.Ent, levels=c("Jul.","Oct.","Dec.","Feb."))

###
phyloRunFungi <- phyloseq(otu_table((haploITS), taxa_are_rows=TRUE),
                          tax_table(taxITS),
                          sample_data(data_Fung))
##

phyloMock<-subset_samples(phyloRunFungi, cond == "mock") #see Mock Analysis section
phyloRunFungi<-subset_samples(phyloRunFungi, cond != "mock")

condition <- function(x) { x > 0 } 
taxaToKeep <- genefilter_sample(phyloRunFungi, condition, 1)
phyloRunFungi<-prune_taxa(taxaToKeep, phyloRunFungi)

taxa_names(phyloRunFungi) <- paste0("F", seq(ntaxa(phyloRunFungi)))
tax_table(phyloRunFungi)[,1] <- gsub("k__", "", tax_table(phyloRunFungi)[,1]);tax_table(phyloRunFungi)[,2] <- gsub("p__", "", tax_table(phyloRunFungi)[,2]);tax_table(phyloRunFungi)[,3] <- gsub("c__", "", tax_table(phyloRunFungi)[,3]);tax_table(phyloRunFungi)[,4] <- gsub("o__", "", tax_table(phyloRunFungi)[,4]);tax_table(phyloRunFungi)[,5] <- gsub("f__", "", tax_table(phyloRunFungi)[,5]);tax_table(phyloRunFungi)[,6] <- gsub("g__", "", tax_table(phyloRunFungi)[,6])


##
phyloRunFungifilt <- transform_sample_counts(phyloRunFungi,function(x) ifelse(x>=0.003*sum(x),x,0)) #see Mock Analysis section

condition <- function(x) { x > 0 } 
taxaToKeep <- genefilter_sample(phyloRunFungifilt, condition, 1)
phyloRunFungifilt<-prune_taxa(taxaToKeep, phyloRunFungifilt)

table(tax_table(phyloRunFungifilt)[, "Phylum"])
phyloRunFungifilt <- subset_taxa(phyloRunFungifilt, !Phylum %in% c("", "Cyanobacteria/Chloroplast", "Unclassified"))
taxa_names(phyloRunFungifilt)
phyloFungnorm <- transform_sample_counts(phyloRunFungifilt, function(x) round(x/sum(x) *100000 ))


###########################
####--- Alpha - div ---####
###########################

measures = c("Shannon","Observed")
Rich<-estimate_richness(phyloFungnorm, measures = measures)
Rich2<-cbind(as(sample_data(phyloFungnorm),"matrix"),Rich)


my_comparisons <- list( c("Jul.","Oct."),c("Oct." , "Dec."),c("Dec." , "Feb."))

Rich2$Month.ent <- factor(Rich2$Month.ent, levels = c("Jul.", "Oct.", "Dec.", "Feb."))


RichPlotITS<-ggplot(Rich2, aes(x=Month.ent, y=Shannon,fill=cond)) +
  stat_compare_means(aes(group = cond), label = "p.signif",label.y=2.8)+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", label.y=3)+ # Add pairwise comparisons p-value
  geom_boxplot(lwd=0.15, width = 0.5)+scale_fill_manual(values=c("#99CCFF","#FFCC66"))+
  ylab(paste("Shannon index (Fungal diversity)")) + xlab("")+ylim(0,5)
RichPlotITS+theme_bw()



##########################
####--- Beta - div ---####
##########################

OTU.ordITS<-ordinate(phyloFungnorm,"MDS","bray")
#pdf("MDS_bray_Bact_all.pdf",width=7,height=3)
MdsITS=plot_ordination(phyloFungnorm,OTU.ordITS,type="sample",color="Month.Inoc",axes = 2:1)+
  scale_shape_manual(values = c(16,17)) + geom_point(size=1)+
  theme(panel.spacing = unit(1, "lines"))+ 
  scale_color_manual(values = c("#009900","#66ff66","#cc0000","#ff9999","#0033cc","#99b3ff","#0d0d0d","#b3b3b3"))
MdsITS+theme_bw()+facet_grid(~Month.ent)+
  stat_ellipse(type = "t")+
  theme_bw()+theme(text = element_text(size=16), axis.text=element_text(size=9), legend.position = "none")
#dev.off()

set.seed(1)
adonis2(distance(phyloFungnorm, "bray") ~ cond+Month.ent,by="margin", data = as(sample_data(phyloFungnorm), "data.frame"))


source(file="~/R/parwise.adonis.r")
pair<-cbind(t(otu_table(phyloFungnorm)),sample_data(phyloFungnorm))
pairwise.adonis(pair[,1:335],pair$Month.Inoc,sim.method = 'bray')


#########################
####---  Heatmap  ---####
#########################

Funglom<-tax_glom(phyloFungnorm,"Genus")
glom   <- subset_taxa(Funglom, Genus !="Unclassified")
glom   <- subset_taxa(glom, Genus !="unclassified_Phaeosphaeriaceae")
glom   <- subset_taxa(glom, Genus !="unclassified_Chaetomiaceae")
glom   <- subset_taxa(glom, Genus !="unclassified_Holtermanniales")
TopNOTUs <- names(sort(taxa_sums(glom), TRUE)[1:25])
Firsts25   <- prune_taxa(TopNOTUs, Funglom)


nameX<- rev(tax_table(ent10)[,"Genus"][order(tax_table(Firsts25)[,"Genus"]),])


plot_heatmap(Firsts25,  method=NULL, taxa.label = "Genus", taxa.order = taxa_names(nameX), low="#F9F8F8", high="#000000", na.value = "#FFFFFF")+
  facet_grid(~date+cond, scales="free")+
  theme(panel.spacing = unit(0, "lines"))+theme_bw()+
  theme(axis.text.y = element_text(face="italic",size=10))



################################
####---  Abundance plot  ---####
################################

GenusGlom<-tax_glom(phyloFungnorm,"Genus")


require(microbiome)

glomtax <- transform_sample_counts(GenusGlom, function(x) round(x/sum(x) *100 ))


Leptoplot<-boxplot_abundance(glomtax, x="Month.ent", y="F3", violin = FALSE, na.rm = FALSE,show.points = FALSE)+
  aes(color=cond)+ geom_boxplot(aes(fill=cond))+
  scale_fill_manual(values=c("#99CCFF","#FFCC66"))+
  scale_color_manual(values=c("black", "black"))+
  labs(x=NULL, y="Reads percentage")+
  theme(axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12),
        axis.title.y = element_text(size=12),strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))+stat_compare_means(aes(group = cond), label = "p.signif",label.y=90)+
  ylim(0,100)+theme_bw()+
  ggtitle("Leptosphaeria")

Leptoplot

Plenoplot<-boxplot_abundance(glomtax, x="Month.ent", y="F1", violin = FALSE, na.rm = FALSE,show.points = FALSE)+
  aes(color=cond)+ geom_boxplot(aes(fill=cond))+
  scale_fill_manual(values=c("#99CCFF","#FFCC66"))+
  scale_color_manual(values=c("black", "black"))+
  labs(x=NULL, y="Reads percentage")+
  theme(axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12),
        axis.title.y = element_text(size=12),strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))+stat_compare_means(aes(group = cond), label = "p.signif",label.y=90)+
  ylim(0,100)+theme_bw()+
  ggtitle("Plenodomus")

Plenoplot

Alterplot<-boxplot_abundance(glomtax, x="Month.ent", y="F8", violin = FALSE, na.rm = FALSE,show.points = FALSE)+
  aes(color=cond)+ geom_boxplot(aes(fill=cond))+
  scale_fill_manual(values=c("#99CCFF","#FFCC66"))+
  scale_color_manual(values=c("black", "black"))+
  labs(x=NULL, y="Reads percentage")+
  theme(axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12),
        axis.title.y = element_text(size=12),strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))+stat_compare_means(aes(group = cond), label = "p.signif",label.y=90)+
  ylim(0,100)+theme_bw()+
  ggtitle("Alternaria")

Alterplot

Hydroplot<-boxplot_abundance(glomtax, x="Month.ent", y="F6", violin = FALSE, na.rm = FALSE,show.points = FALSE)+
  aes(color=cond)+ geom_boxplot(aes(fill=cond))+
  scale_fill_manual(values=c("#99CCFF","#FFCC66"))+
  scale_color_manual(values=c("black", "black"))+
  labs(x=NULL, y="Reads percentage")+
  theme(axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12),
        axis.title.y = element_text(size=12),strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))+stat_compare_means(aes(group = cond), label = "p.signif",label.y=90)+
  ylim(0,100)+theme_bw()+
  ggtitle("Hydropisphaerella")
Hydroplot

Monoplot<-boxplot_abundance(glomtax, x="Month.ent", y="F2", violin = FALSE, na.rm = FALSE,show.points = FALSE)+
  aes(color=cond)+ geom_boxplot(aes(fill=cond))+
  scale_fill_manual(values=c("#99CCFF","#FFCC66"))+
  scale_color_manual(values=c("black", "black"))+
  labs(x=NULL, y="Reads percentage")+
  theme(axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12),
        axis.title.y = element_text(size=12),strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))+stat_compare_means(aes(group = cond), label = "p.signif",label.y=90)+
  ylim(0,100)+theme_bw()+
  ggtitle("Hydropisphaerella")
Monoplot

library(cowplot)
plot_grid(Leptoplot, Plenoplot, Alterplot, Hydroplot)
#save 5.3*3.5


####################
## Mock Analysis ##
####################

sample_data(phyloMock)
taxa_names(phyloMock) <- paste0("MF", seq(ntaxa(phyloMock)))


Mocknorm <- transform_sample_counts(phyloMock, function(x) round(x/sum(x) *100000 ))

condition <- function(x) { x > 0 } 
taxaToKeep <- genefilter_sample(Mocknorm, condition, 1)
Mocknorm<-prune_taxa(taxaToKeep, Mocknorm)
sample_names(Mocknorm)<-c("Mock Run1","Mock Run3")
sample_data(Mocknorm)


#Filtering Threshold
Mocknorm <- transform_sample_counts(phyloMock, function(x) round(x/sum(x) *100000 ))

tdt = data.table(tax_table(Mocknorm),
                 TotalCounts = taxa_sums(Mocknorm),
                 OTU = taxa_names(Mocknorm))
ggplot(tdt, aes(TotalCounts)) + 
  geom_histogram() + 
  ggtitle("Histogram of Total Counts")


taxcumsum = tdt[, .N, by = TotalCounts]
setkey(taxcumsum, TotalCounts)
taxcumsum[, CumSum := cumsum(N)]
# Define the plot
pCumSum = ggplot(taxcumsum, aes(TotalCounts, CumSum)) + 
  geom_point() +
  xlab("Filtering Threshold, Minimum Total Counts") +
  ylab("OTUs Filtered") +
  ggtitle("OTUs that would be filtered vs. the minimum count threshold")
pCumSum+ geom_vline(xintercept = 300, color="red",size=1.5)+ scale_x_continuous(trans='log10')



#Heatmap
TopNOTUs <- names(sort(taxa_sums(Mocknorm), TRUE)[1:40])
Top40   <- prune_taxa(TopNOTUs, Mocknorm)

plot_heatmap(Top40, method = NULL,taxa.label = "Genus" )+facet_grid(~Sample, scales="free")

