######## Packages required + custom functions
Packages <- c("phyloseq", "data.table", "ggplot2", "plyr","dplyr","reshape2","grid",
              "gridExtra","scales","dplyr", "ggpubr","vegan","multcompView","rcompanion","betapart")

lapply(Packages, library, character.only = TRUE)



########################
# Dataset constitution #
########################

tax16S <- replace(tax16S, is.na(tax16S), "Unclassified")
data_16S$Month.ent <- factor(data_16S$Month.ent, levels = c("Jul.", "Oct.", "Dec.", "Feb."))
data_16S$Month.Inoc <- factor(data_16S$Month.Inoc, levels = c("Jul.DR", "Jul.DS","Oct.DR", "Oct.DS", "Dec.DR","Dec.DS", "Feb.DR","Feb.DS"))

###
phyloRunBact <- phyloseq(otu_table((haplo16S), taxa_are_rows=TRUE),
                          tax_table((tax16S)),
                          sample_data(data_16S))

taxa_names(phyloRunBact) <- paste0("B", seq(ntaxa(phyloRunBact)))
tax_table(phyloRunBact)[,1] <- gsub("k__", "", tax_table(phyloRunBact)[,1]);tax_table(phyloRunBact)[,2] <- gsub("p__", "", tax_table(phyloRunBact)[,2]);tax_table(phyloRunBact)[,3] <- gsub("c__", "", tax_table(phyloRunBact)[,3]);tax_table(phyloRunBact)[,4] <- gsub("o__", "", tax_table(phyloRunBact)[,4]);tax_table(phyloRunBact)[,5] <- gsub("f__", "", tax_table(phyloRunBact)[,5]);tax_table(phyloRunBact)[,6] <- gsub("g__", "", tax_table(phyloRunBact)[,6])

##
phyloMock<-subset_samples(phyloRunBact, cond == "mock")  #see Mock Analysis section
phyloRunBact<-subset_samples(phyloRunBact, cond != "mock")

condition <- function(x) { x > 0 } 
taxaToKeep <- genefilter_sample(phyloRunBact, condition, 1)
phyloRunBact<-prune_taxa(taxaToKeep, phyloRunBact)


##
phyloRunBactfilt <- transform_sample_counts(phyloRunBact,function(x) ifelse(x>=0.003*sum(x),x,0))  #see Mock Analysis section

condition <- function(x) { x > 0 } 
taxaToKeep <- genefilter_sample(phyloRunBactfilt, condition, 1)
phyloRunBactfilt<-prune_taxa(taxaToKeep, phyloRunBactfilt)


table(tax_table(phyloRunBactfilt)[, "Phylum"])
phyloRunBactfilt <- subset_taxa(phyloRunBactfilt, !Phylum %in% c("", "Cyanobacteria/Chloroplast", "Unclassified"))

phyloBactnorm <- transform_sample_counts(phyloRunBactfilt, function(x) round(x/sum(x) *100000 ))
phyloBactnorm <-subset_samples(phyloBactnorm, sample_names(phyloBactnorm) !="DR.2017.07.5")


###########################
####--- Alpha - div ---####
###########################


measures=c("Observed", "Shannon")
#various possible indexes : c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
Rich<-estimate_richness(phyloBactnorm, measures = measures)
data2<-subset(data_16S, cond != "mock")
data2<-subset(data2, row.names(data2) !="DR.2017.07.5")
row.names(data2)
Rich<-cbind(Rich,data2)


my_comparisons <- list( c("Jul.","Oct."),c("Oct." , "Dec."),c("Dec." , "Feb."))


compare_means(Shannon ~ cond,  data = Rich, method = "kruskal.test")
compare_means(Observed ~ cond,  data = Rich, method = "kruskal.test")

RichPlotITS<-ggplot(Rich, aes(x=Month.ent, y=Shannon,fill=cond)) +
  stat_compare_means(aes(group = cond), label = "p.signif",label.y=4.5)+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", label.y=5)+ # Add pairwise comparisons p-value
  geom_boxplot(lwd=0.15, width = 0.5)+scale_fill_manual(values=c("#99CCFF","#FFCC66"))+ ylim(0,5)+
  ylab(paste("Shannon index (Bacterial diversity)")) + xlab("")


RichPlotITS+theme_bw()


##########################
####--- Beta - div ---####
##########################

OTU.ordITS<-ordinate(phyloBactnorm,"MDS","bray")

MdsITS=plot_ordination(phyloBactnorm,OTU.ordITS,type="sample",color="Month.Inoc",axes = 2:1)+
  scale_shape_manual(values = c(16,17)) + geom_point(size=1)+
  theme(panel.spacing = unit(1, "lines"))+ 
  scale_color_manual(values = c("#009900","#66ff66","#cc0000","#ff9999","#0033cc","#99b3ff","#0d0d0d","#b3b3b3"))

MdsITS+theme_bw()+facet_grid(~Month.ent)+
  stat_ellipse(type = "t")+
  theme_bw()+theme(text = element_text(size=16), axis.text=element_text(size=9), legend.position = "none")


##
set.seed(1)
adonis2(distance(phyloBactnorm, "bray") ~ cond+Month.ent,by="margin", data = as(sample_data(phyloBactnorm), "data.frame"))

source(file="~/R/parwise.adonis.r")
pair<-cbind(t(otu_table(phyloBactnorm)),sample_data(phyloBactnorm))
pairwise.adonis(pair[,1:610],pair$Month.Inoc,sim.method = 'bray')


#########################
####---  Heatmap  ---####
#########################

Bactglom<-tax_glom(phyloBactnorm,"Genus")
test   <- subset_taxa(Bactglom, Genus !="Unclassified")
test   <- subset_taxa(test, Genus !="unclassified_Phaeosphaeriaceae")
TopNOTUs <- names(sort(taxa_sums(test), TRUE)[1:35])
ent10   <- prune_taxa(TopNOTUs, Bactglom)


nameX<- rev(tax_table(ent10)[,"Genus"][order(tax_table(ent10)[,"Genus"]),])

plot_heatmap(ent10,  method=NULL, taxa.label = "Genus", taxa.order = taxa_names(nameX), low="#F9F8F8", high="#000000", na.value = "#FFFFFF")+
  facet_grid(~date+cond, scales="free")+  
  theme_bw()+
  theme(axis.text.y = element_text(face="italic",size=10))
