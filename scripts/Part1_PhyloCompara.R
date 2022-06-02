

##########################################################################################
############################## 1. Building Phylogenetic Tree ############################
##########################################################################################

TreeInputMatrix = PruneTree(unique(InputMatrix$Species), treeVascularPlants)

extra_spp = setdiff(unique(InputMatrix$Species), TreeInputMatrix$tip.label)
extra_spp

## After PruneTree, there are two species are not on the tree: extra_spp
## These two species are: Callirhoe_leiocarpa and Digitaria_californic
## According two their taxonomic relationship assume their sister species: sister_spp
## Sister species are: Callirhoe_involucrata and Digitaria_ciliaris

sister_spp = c("Callirhoe_involucrata", "Digitaria_ciliaris")

## Adding "Callirhoe_leiocarpa"
tip1 <- extra_spp[1]
sister1 <- sister_spp[1]
tree1 = bind_sister(tip1,sister1, TreeInputMatrix)

## Adding "Digitaria_californic"
tip2 <- extra_spp[2]
sister2 <- sister_spp[2]
TreeAllMatrix = bind_sister(tip2, sister2, tree1)


##########################################################################################
############################## 2. Calculate Phylogenetic Positiion #######################
##########################################################################################
## Phylogenetic Position
par(mar=c(5,5,5,1))
dist <- cophenetic.phylo(TreeAllMatrix)
phyloposi <- MASS::isoMDS(dist) %>% as.data.frame()

phyloposi_species = phyloposi %>% 
  mutate(Species = row.names(phyloposi)) %>% 
  mutate(phy1 = round(scale(points.1),digits = 1)) %>% 
  mutate(phy2 = round(scale(points.2), digits = 1))

AllMatrix = merge (dat, phyloposi_species)

phyloposi_family = AllMatrix %>% group_by(Family) %>% 
  summarise(Family, Classification, phy1, phy2) %>%
  distinct() 

p_family = ggplot(phyloposi_family, aes(x=phy1, y=phy2, color=Classification)) +
  geom_point() + labs(x="phy1", y="phy2") + 
  geom_text_repel(aes(label = Family), size =3.5) + 
  scale_color_manual(values=c("#FFC20A", "#64D294", "#C76BA2")) +
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "white", color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'))

## Figure 2: Phylogenetic Position
## save phylogenetic position figure into results

pdf(file = file.path("results", paste0("Figure2_Phylogenetic_Position.pdf")))
print(p_family)
graphics.off()


##########################################################################################
################### 3. Array Trait Value with Phylogenetic Tree ##########################
##########################################################################################
# Phylogenetic Tree: TreeAllMatrix



# Trait Data: dat
A = TreeAllMatrix$tip.label

dat_G = merge(dat, Germination)

B1 = dat_G %>% mutate (Species1 = factor(Species, level = rev(A)),
                       Classification1 = factor(Classification, levels = c("Asteracea", "Other_Dicot", "Monoct"), 
                                                ordered = TRUE)) 

G1 = B1[with(B1, order(Classification1, Species1)),] 
G2 = G1 %>% mutate(Species1 = NULL, Classification1 = NULL, SpeciesNum = 1:length(TreeAllMatrix$tip.label))

#3.1 Trait Value Figures (Seed Mass and Seed Height, Colored)
Mass_Figure <- ggplot(G2, aes(y = mean_Mass, x = reorder(Species, -SpeciesNum))) + 
  geom_point(aes(shape = Classification, color = Classification)) + 
  geom_errorbar(aes(ymin=mean_Mass - se_Mass, ymax=mean_Mass + se_Mass, color = Classification))+
  labs(y = "Mass(log)") +
  scale_color_manual(values=c("#FFC20A", "#64D294", "#C76BA2")) +
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "white", color = 'black'),
        panel.grid.major = element_line(color = 'grey93', linetype = 'solid'),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip()
#Mass_Figure

Height_Figure <- ggplot(G2, aes(y = mean_Height, x = reorder(Species, -SpeciesNum))) + 
  geom_point(aes(shape = Classification, color = Classification)) + 
  geom_errorbar(aes(ymin=mean_Height - se_Height, ymax=mean_Height + se_Height, color = Classification))+
  labs(y = "Height (mm)") +
  scale_color_manual(values=c("#FFC20A", "#64D294", "#C76BA2")) +
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "white", color = 'black'),
        panel.grid.major = element_line(color = 'grey93', linetype = 'solid'),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip()
#Height_Figure

Area_Figure <- ggplot(G2, aes(y = mean_Area, x = reorder(Species, -SpeciesNum))) + 
  geom_point(aes(shape = Classification, color = Classification)) + 
  geom_errorbar(aes(ymin=mean_Area - se_Area, ymax=mean_Area + se_Area, color = Classification))+
  labs(y = "Area (sq mm)") +
  scale_color_manual(values=c("#FFC20A", "#64D294", "#C76BA2")) +
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "white", color = 'black'),
        panel.grid.major = element_line(color = 'grey93', linetype = 'solid'),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip()
#Area_Figure

Germination_Figure <- ggplot(G2, aes(y = mean_Germination, x = reorder(Species, -SpeciesNum))) + 
  geom_point(aes(shape = Classification, color = Classification)) + 
  geom_errorbar(aes(ymin=mean_Germination - se_Germination, ymax=mean_Germination + se_Germination, color = Classification))+
  labs(y = "Germination (%)") +
  scale_color_manual(values=c("#FFC20A", "#64D294", "#C76BA2")) +
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "white", color = 'black'),
        panel.grid.major = element_line(color = 'grey93', linetype = 'solid'),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip()
#Germination_Figure


#3.2 Plotting/Coloring Phylogenetic Tree, Colored

classification = list(asteracea = G2$Species[1:8],
                      other_dicoct = G2$Species[9:30],
                      monocot = G2$Species[31:45])

t = ggtree(TreeAllMatrix, branch.length ='none')


t2 = groupOTU(t, classification, 'Classification') + aes(color = Classification) + geom_treescale(x=20, y =1, offset = 2, color = "white") + geom_tiplab(size = 3.8, fontface = 'italic') + theme(legend.position = "none") + scale_color_manual(values=c("#FFC20A", "#64D294", "#C76BA2"))


#3.3 Figure 1
## save patch figure into results

pdf(file = file.path("results", paste0("Figure1_Phylogenetic_compara.pdf")), height = 8, width = 12)

(t2 + guides(colour = "none")) + Mass_Figure + Height_Figure + Area_Figure + Germination_Figure + patchwork::plot_layout(widths = c(2.5,1,1,1,1), guides = "collect") & theme(legend.position = 'bottom')

graphics.off()


##########################################################################################
############################## 4. Calculate Phylogenetic Signal ##########################
##########################################################################################



phyloSig_all = PhyloSig(TreeAllMatrix, dat_G)

## save table into results 

write.csv(phyloSig_all, file = "results/table2_phylogenetic_signal.csv")


##########################################################################################
############################## 4.Model Selection  ########################################
##########################################################################################

Model_AICc_all = Model_AICc(TreeAllMatrix, dat_G)

## Appendix I: AICc values of model selection
## save table into results 

write.csv(Model_AICc_all, file = "results/appendix_model_selection.csv")



##########################################################################################
###################### 5. Correlation between traits #####################################
##########################################################################################


# Correlation between mass and height 
cor.test(dat_G_scaled$scaled_Height, dat_G_scaled$scaled_Mass, 
         method = "pearson")
cor.test(dat_G_scaled$scaled_Mass, dat_G_scaled$scaled_Area, 
         method = "pearson")
cor.test(dat_G_scaled$scaled_Height, dat_G_scaled$scaled_Area, 
         method = "pearson")



