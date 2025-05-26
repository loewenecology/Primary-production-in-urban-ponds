
##### Road salt pollution elevates risk of cyanobacteria #####

##################
# Clear any variables from the R environment
##################

rm(list=ls())


##################
# Load R packages
##################

library(tidyverse)
library(vegan)
library(ggeffects)

library(ggplot2)
library(ggrepel)

library(lubridate)
library(grid)

library(sf)
library(ggspatial)


##################
# Load data in R environment
##################

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) #set working directory

pigments = read.table("Pond_pigments_2022_ms.csv", header = T, sep = ",") #read in sample pigment data
sites = read.table("Pond_sites_2022_ms.csv", header = T, sep = ",") #read in pond site data


##################
# Create theme for plotting
##################

workingtheme_bars <- theme(strip.background = element_blank(),
                           
                           panel.border = element_rect(colour = "black", fill = NA),
                           panel.grid.major = element_line(colour = "grey"),
                           panel.grid.minor = element_blank(),
                           panel.background = element_rect(colour = "black", fill = "white", linewidth = 0.75),
                           panel.spacing = unit(0.12, "cm"),
                           
                           axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 12),
                           axis.text.y = element_text(colour = "black", size = 12),
                           axis.ticks = element_line(colour = "black", linewidth = 0.4),
                           axis.title.y = element_text(colour = "black", size = 16),
                           axis.title.x = element_text(colour = "black", size = 16))

workingtheme_pca <- theme(strip.background = element_blank(),
                          
                          panel.border = element_rect(colour = "black", fill = NA),
                          panel.grid.major = element_line(colour = "grey"),
                          panel.grid.minor = element_blank(),
                          panel.background = element_rect(colour = "black", fill = "white", linewidth = 0.75),
                          panel.spacing = unit(0.12, "cm"),
                          
                          axis.text.x = element_text(colour = "black", size = 12),
                          axis.text.y = element_text(colour = "black", size = 12),
                          axis.ticks = element_line(colour = "black", linewidth = 0.4),
                          axis.title.y = element_text(colour = "black", size = 16),
                          axis.title.x = element_text(colour = "black", size = 16))

workingtheme_regressions <- theme(strip.background = element_blank(),
                                  
                                  panel.border = element_rect(colour = "black", fill = NA),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  panel.background = element_rect(colour = "black", fill = "white", linewidth = 0.75),
                                  panel.spacing = unit(0.12, "cm"),
                                  
                                  axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 12),
                                  axis.text.y = element_text(colour = "black", size = 12),
                                  axis.ticks = element_line(colour = "black", linewidth = 0.4),
                                  axis.title.y = element_text(colour = "black", size = 16),
                                  axis.title.x = element_text(colour = "black", size = 16),
                                  legend.position.inside = c(0.3, 0.7))


##################
# Prepare pigment data
##################

# Run the following two lines to exclude C56 with elevated chloride concentrations for sensitivity analysis
#pigments <- subset(pigments, Site != "C56")
#sites <- subset(sites, Site != "C56")

# Create pigment data frame for filter data (phytoplankton) dropping Chlorophyll and Pheophytin
top_filter <- subset(pigments, Type == "Top") %>%
  select(!Type & !Chlorophyll.a & !Pheophytin.a) %>%
  remove_rownames() %>%
  column_to_rownames(var = "Site")

# Standardize data frame (divide by margin total to calculate relative concentrations) 
top_filter_total <- decostand(top_filter, "total")

# Create data frame for slide data (periphyton) dropping Chlorophyll and Pheophytin
top_slide <- subset(pigments, Type == "Slides") %>%
  select(!Type & !Chlorophyll.a & !Pheophytin.a) %>%
  remove_rownames() %>%
  column_to_rownames(var = "Site")

# Standardize data frame (divide by margin total to calculate relative concentrations) 
top_slide_total <- decostand(top_slide, "total")

# Create column for site name to join filter and slide data
top_filter_total <- rownames_to_column(top_filter_total, var = "Site")
top_slide_total <- rownames_to_column(top_slide_total, var = "Site")
top_filterslide_merge <- left_join(top_filter_total, top_slide_total, by = c("Site"))

# Reset row names
top_filter_total <- column_to_rownames(top_filter_total, "Site")
top_slide_total <- column_to_rownames(top_slide_total, "Site")

# Create data frame of taxonomic groups for filter data (phytoplankton)
top_filter_groups <- top_filter %>%
  mutate(Chromophytes = Fucoxanthin + Diadinoxanthin + Diatoxanthin) %>%
  mutate(Chlorophytes = Chlorophyll.b + Violaxanthin + Neoxanthin + Lutein) %>%
  mutate(Cryptophytes = Alloxanthin) %>%
  mutate(Cyanobacteria = Zeaxanthin + Canthaxanthin + Myxoxanthophyll + Lutein.Zeaxanthin) %>%
  select(c("Chromophytes", "Chlorophytes", "Cryptophytes", "Cyanobacteria"))

# Standardize data frame (divide by margin total to calculate relative concentrations) 
top_filter_groups_total <- decostand(top_filter_groups, "total")

# Create data frame of taxonomic groups for slide data (periphyton)
top_slide_groups <- top_slide %>%
  mutate(Chromophytes = Fucoxanthin + Diadinoxanthin + Diatoxanthin) %>%
  mutate(Chlorophytes = Chlorophyll.b + Violaxanthin + Neoxanthin + Lutein) %>%
  mutate(Cryptophytes = Alloxanthin) %>%
  mutate(Cyanobacteria = Zeaxanthin + Canthaxanthin + Myxoxanthophyll + Lutein.Zeaxanthin) %>%
  select(c("Chromophytes", "Chlorophytes", "Cryptophytes", "Cyanobacteria"))

# Standardize data frame (divide by margin total to calculate relative concentrations) 
top_slide_groups_total <- decostand(top_slide_groups, "total")

# Create column for site name to join filter and slide data
top_filter_groups_total <- rownames_to_column(top_filter_groups_total, var = "Site")
top_slide_groups_total <- rownames_to_column(top_slide_groups_total, var = "Site")
top_slide_group_join <- left_join(top_filter_groups_total, top_slide_groups_total, by = c("Site"))

# Reset row names
top_filter_groups_total <- column_to_rownames(top_filter_groups_total, var = "Site")
top_slide_groups_total <- column_to_rownames(top_slide_groups_total, var = "Site")

# Rename variables
names(sites)[names(sites) == 'Total.Phosphorus..as.P.'] <- 'Total Phosphorus (as P)'
names(sites)[names(sites) == 'NO3.NO2..as.N.'] <- 'NO3 + NO2 (as N)'
names(sites)[names(sites) == 'Total.Suspended.Solids'] <- 'Total Suspended Solids'
names(sites)[names(sites) == 'Field.pH'] <- 'Field pH'
names(sites)[names(sites) == 'Field.Temp'] <- 'Field Temp'
names(sites)[names(sites) == 'Total.Copper'] <- 'Total Copper'
names(sites)[names(sites) == 'Total.Sodium'] <- 'Total Sodium'

# Create new data frame of selected parameters
site_select <- select(sites, !any_of(c("Latitude",
                                       "Longitude",
                                       "Total Sodium",
                                       "Conductance")))


##################
# Create variable correlation heatmaps
##################

# Create data frame of log-transformed variables
site_select_log <- site_select
site_select_log$`Total Phosphorus (as P)` <- log(site_select$`Total Phosphorus (as P)`)
site_select_log$`NO3 + NO2 (as N)` <- log(site_select$`NO3 + NO2 (as N)`)
site_select_log$`Total Suspended Solids` <- log(site_select$`Total Suspended Solids`)
site_select_log$`Field Temp` <- log(site_select$`Field Temp`)
site_select_log$`Chloride` <- log(site_select$`Chloride`)
site_select_log$`Total Copper` <- log(site_select$`Total Copper`)

# Set row names
site_select_log <- remove_rownames(site_select_log) %>% column_to_rownames(var = "Site")

# Calculate Pearson Correlation Coefficients
env.cor <- as.data.frame(round(cor(site_select_log) , 2))
env.cor <- rownames_to_column(env.cor, var = "var1")

env.cor.long <- as.data.frame(env.cor) %>%
  pivot_longer(where(is.numeric), names_to = "var2", values_to = "cor")

# Plot Pearson Correlation Coefficient Heatmap (6 x 6)
ggplot(data = env.cor.long, aes(x = var1, y = var2, fill = cor)) +
  geom_tile() +
  scale_fill_gradient(name = "Pearson Correlation Coefficient", low = "lightblue", high = "tomato1") +
  geom_text(aes(label = cor),
    color = "black",
    size = 3.5) +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_blank(),
        legend.position = "bottom") +
  coord_fixed()


##################
# Conduct PCAs of pigment data
##################

## PCA of combined pigments
# Create combined data frame for filter (phytoplankton) and slide (periphyton) data dropping Chlorophyll and Pheophytin
combined_pigments <- subset(pigments, Type == "Top" | Type == "Slides") %>%
  select(!Chlorophyll.a & !Pheophytin.a) %>%
  remove_rownames()

combined_pigments$Type <- sub('Top', 'phyto', combined_pigments$Type)
combined_pigments$Type <- sub('Slides', 'peri', combined_pigments$Type)

row.names(combined_pigments) <- paste(combined_pigments$Site, combined_pigments$Type)
combined_pigments <- select(combined_pigments, -Site, -Type)

# Hellinger transform combined (phyto/peri) data frame
combined_pigments_hell <- decostand(combined_pigments, "hellinger")

# Conduct PCA of Hellinger transformed pigment data
combined_pca <- rda(combined_pigments_hell)
summary(combined_pca)

# Obtain PCA scores
pig_scores <- scores(combined_pca, display = "species", choices = c(1,2), scaling = 2)
site_scores <- scores(combined_pca, display = "sites", choices = c(1,2), scaling = 2)
site_scores <- as.data.frame(site_scores)

# Create column for site name
site_scores <- rownames_to_column(site_scores, var = "Data")

# Create Site and Type columns
site_scores_delim <-  separate_wider_delim(data = site_scores, cols = Data, delim = " ", names = c("site", "type"))

# Plot PCA of pigments (5 x 5)
(plot.PCA.pigs <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "lightgrey", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey", linewidth = 0.5) +
    geom_line(data = site_scores_delim, aes(x = PC1, y = PC2, group = site), color = "grey", linewidth = 1) +
    geom_point(data = site_scores_delim, aes(x = PC1, y = PC2, color = type), size = 2) +
    geom_point(data = pig_scores, aes(x = PC1, y = PC2), size = 3, color = "black") +
    geom_label_repel(data = pig_scores, aes(x = PC1, y = PC2, label = rownames(pig_scores)), size = 3) +
    scale_x_continuous("PC1 (29%)", limits = c(-0.8, 0.6), breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
    scale_y_continuous("PC2 (24%)", limits = c(-0.9, 0.75), breaks = c(-0.8, -0.4, 0, 0.4, 0.8)) +
    workingtheme_pca +
    theme(legend.position = "none",
          panel.grid.major = element_blank()))

# Create matching site data
site_select_phyto <- site_select
site_select_peri <- site_select

site_select_phyto$Type <- "phyto"
site_select_peri$Type <- "peri"

combined_site_select <- bind_rows(site_select_phyto, site_select_peri)

# Fit environmental vectors to PCA of pigments
combined_pca_env <- envfit(combined_pca ~ 
                             log(`Total Phosphorus (as P)`) + 
                             log(`NO3 + NO2 (as N)`) + 
                             log(`Total Suspended Solids`) + 
                             `Field pH` + 
                             log(`Field Temp`) + 
                             log(Chloride) + 
                             log(`Total Copper`),
                           data = combined_site_select)

rownames(combined_pca_env$vectors$arrows) <- c("Total Phosphorus (as P)",
                                               "NO3 + NO2 (as N)",
                                               "Total Suspended Solids",
                                               "Field pH",
                                               "Field Temp",
                                               "Chloride",
                                               "Total Copper")

# Obtain environmental vectors
env_scores <- as.data.frame(scores(combined_pca_env, display = "vectors"))

# Plot PCA of pigments with environmental vectors (5 x 5)
(plot.PCA.pigs.env <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "lightgrey", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey", linewidth = 0.5) +
    geom_line(data = site_scores_delim, aes(x = PC1, y = PC2, group = site), color = "grey", linewidth = 1) +
    geom_point(data = site_scores_delim, aes(x = PC1, y = PC2, color = type), size = 2) +
    geom_segment(data = env_scores, aes(xend = PC1 * ordiArrowMul(combined_pca_env), yend = PC2 * ordiArrowMul(combined_pca_env)), x = 0, y = 0,
                 color = "black", linewidth = 1,
                 arrow = arrow(type = "closed", length = unit(0.3,"cm"))) +
    geom_label_repel(data = env_scores, aes(x = PC1 * ordiArrowMul(combined_pca_env), y = PC2 * ordiArrowMul(combined_pca_env), label = rownames(env_scores)), size = 3) +
    scale_x_continuous("PC1 (29%)", limits = c(-0.8, 0.6), breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
    scale_y_continuous("PC2 (24%)", limits = c(-0.9, 0.75), breaks = c(-0.8, -0.4, 0, 0.4, 0.8)) +
    workingtheme_pca +
    theme(legend.position = "none",
          panel.grid.major = element_blank()))

## PCA of combined taxonomic groups
# Create combined data frame of taxonomic groups for filter (phytoplankton) and slide (periphyton) data
top_filter_groups$Type <- "phyto"
top_slide_groups$Type <- "peri"

top_filter_groups <- rownames_to_column(top_filter_groups, var = "Site")
top_slide_groups <- rownames_to_column(top_slide_groups, var = "Site")

combined_groups <- bind_rows(top_filter_groups, top_slide_groups)

row.names(combined_groups) <- paste(combined_groups$Site, combined_groups$Type)
combined_groups <- select(combined_groups, -Site, -Type)

# Hellinger transform combined (phyto/peri) data frame
combined_groups_trans <- decostand(combined_groups, "hellinger")

# Conduct PCA of Hellinger transformed taxonomic group data
combined_groups_pca <- rda(combined_groups_trans)
summary(combined_groups_pca)

# Obtain PCA scores
pig_scores <- scores(combined_groups_pca, display = "species", choices = c(1,2), scaling = 2)
site_scores <- scores(combined_groups_pca, display = "sites", choices = c(1,2), scaling = 2)
site_scores <- as.data.frame(site_scores)

# Create column for site name
site_scores <- rownames_to_column(site_scores, var = "Data")

# Create Site and Type columns
site_score_delim <-  separate_wider_delim(data = site_scores, cols = Data, delim = " ", names = c("site", "type"))

# Plot PCA of taxonomic groups (5 x 5)
(plot.PCA.group <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "lightgrey", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey", linewidth = 0.5) +
    geom_line(data = site_score_delim, aes(x = PC1, y = PC2, group = site), color = "grey", linewidth = 1) +
    geom_point(data = site_score_delim, aes(x = PC1, y = PC2, color = type), size = 2) +
    geom_point(data = pig_scores, aes(x = PC1, y = PC2), color = "black", size = 3) +
    geom_label_repel(data = pig_scores, aes(x = PC1, y = PC2, label = rownames(pig_scores)), size = 3) +
    scale_x_continuous("PC1 (43%)", limits = c(-0.7, 1.4), breaks = c(-0.8, -0.4, 0, 0.4, 0.8, 1.2)) +
    scale_y_continuous("PC2 (36%)", limits = c(-0.7, 0.9), breaks = c(-0.6, -0.3, 0, 0.3, 0.6, 0.6, 0.9)) +
    workingtheme_pca +
    theme(legend.position = "none",
          panel.grid.major = element_blank()))

# Fit environmental vectors to PCA of taxonomic groups
combined_groups_pca_env <- envfit(combined_groups_pca ~ 
                                    log(`Total Phosphorus (as P)`) + 
                                    log(`NO3 + NO2 (as N)`) + 
                                    log(`Total Suspended Solids`) + 
                                    `Field pH` + 
                                    log(`Field Temp`) + 
                                    log(Chloride) + 
                                    log(`Total Copper`), 
                                  data = combined_site_select)

rownames(combined_groups_pca_env$vectors$arrows) <- c("Total Phosphorus (as P)",
                                                      "NO3 + NO2 (as N)",
                                                      "Total Suspended Solids",
                                                      "Field pH",
                                                      "Field Temp",
                                                      "Chloride",
                                                      "Total Copper") 

# Obtain environmental vectors
env_scores <- as.data.frame(scores(combined_groups_pca_env, display = "vectors"))

# Plot PCA of taxonomic groups with environmental vectors (5 x 5)
(plot.PCA.group.env <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "lightgrey", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey", linewidth = 0.5) +
    geom_line(data = site_score_delim, aes(x = PC1, y = PC2, group = site), color = "grey", linewidth = 1) +
    geom_point(data = site_score_delim, aes(x = PC1, y = PC2, color = type), size = 2) +
    geom_segment(data = env_scores, aes(xend = PC1 * ordiArrowMul(combined_groups_pca_env), yend = PC2 * ordiArrowMul(combined_groups_pca_env)), x = 0, y = 0,
                 color = "black", linewidth = 1,
                 arrow = arrow(type = "closed", length = unit(0.3,"cm"))) +
    geom_label_repel(data = env_scores, aes(x = PC1 * ordiArrowMul(combined_groups_pca_env), y = PC2 * ordiArrowMul(combined_groups_pca_env), label = rownames(env_scores)), 
                     max.overlaps = 11, size = 3) +
    #scale_x_continuous("PC1 (43%)", limits = c(-0.7, 1.4), breaks = c(-0.8, -0.4, 0, 0.4, 0.8, 1.2)) +
    scale_x_continuous("PC1 (43%)", limits = c(-0.8, 1.4), breaks = c(-0.8, -0.4, 0, 0.4, 0.8, 1.2)) + #use for sensitivity analysis plot
    #scale_y_continuous("PC2 (36%)", limits = c(-0.7, 0.9), breaks = c(-0.6, -0.3, 0, 0.3, 0.6, 0.6, 0.9)) +
    scale_y_continuous("PC2 (36%)", limits = c(-0.7, 1.0), breaks = c(-0.6, -0.3, 0, 0.3, 0.6, 0.6, 0.9)) + #use for sensitivity analysis plot
    workingtheme_pca +
    theme(legend.position = "none",
          panel.grid.major = element_blank()))


##################
# Create a study area map
##################

# Set coordinate reference system
samplingLocs <- sites %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

# Plot local study area map (6 x 6)
(Sampling_Mapa <- ggplot() +
    annotation_map_tile(zoom = 12, "cartolight") +
    geom_sf(data = samplingLocs, size = 3,
            aes(color = Chloride)) +
    scale_colour_gradient(name = Chloride~(mg/L), low = "blue4", high = "red") +
    annotation_scale(location = "bl", width_hint = 0.5) +
    annotation_north_arrow(location = "bl", which_north = "true",
                           pad_x = unit(0.25, "in"), pad_y = unit(0.25, "in"),
                           style = north_arrow_fancy_orienteering) +
    xlab("Longitude") + ylab("Latitude") +
    coord_sf(xlim = c(-79.65205, -79.84537),
             ylim = c(43.615, 43.79379),
             crs = 4326) +
    theme(legend.position = "bottom"))

# Plot regional study area map (5 x 5)
(Sampling_Mapb <- ggplot() +
    annotation_map_tile(zoom = 5, "cartolight") +
    geom_sf(data = samplingLocs, size = 2,
            aes(color = Chloride)) +
    scale_colour_gradient(name = Chloride~(mg/L), low = "blue4", high = "red") +
    xlab("Longitude") + ylab("Latitude") +
    coord_sf(xlim = c(-68.65205, -90.84537),
             ylim = c(32.615, 54.79379),
             crs = 4326))


##################
# Contrast pigment composition across sites
##################

# Pivot pigment data longer
top_filterslide_merge_long <- top_filterslide_merge %>%
  pivot_longer(!Site, names_to = "pigment", values_to = "value")

top_filterslide_merge_long <- top_filterslide_merge_long %>%
  mutate(type = ifelse(grepl("x$", pigment), "Phyto", "Peri"))

top_filterslide_merge_long$pigment = gsub("[.]x", "", top_filterslide_merge_long$pigment)
top_filterslide_merge_long$pigment = gsub("[.]y", "", top_filterslide_merge_long$pigment)

top_filterslide_merge_long <- left_join(top_filterslide_merge_long, sites, by = c("Site"))

# Reorder sites based on chloride concentrations
top_filterslide_merge_long$Site <- reorder(top_filterslide_merge_long$Site, top_filterslide_merge_long$Chloride)

# Reorder types
top_filterslide_merge_long$type <- ordered(top_filterslide_merge_long$type, levels = c("Peri", "Phyto"))

# Plot pigment composition across sites (8 x 6 portrait)
(plot.bar.pig.facet <- ggplot() +
    geom_col(data = top_filterslide_merge_long,
             aes(x = type, y = value, fill = pigment)) +
    ylab("Proportion") +
    xlab("Community") +
    facet_wrap(~Site, nrow = 6) +
    workingtheme_bars +
    guides(fill = guide_legend(title = "Pigment")) +
    theme(panel.grid.major = element_blank(),
          legend.title=element_blank(),
          legend.position = "bottom"))


##################
# Contrast taxonomic group composition across sites
##################

# Pivot taxonomic group data longer
top_slide_group_join_long <- top_slide_group_join %>%
  pivot_longer(!Site, names_to = "group", values_to = "value")

top_slide_group_join_long <- top_slide_group_join_long %>%
  mutate(type = ifelse(grepl("x$", group), "Phyto", "Peri"))

top_slide_group_join_long$group = gsub("[.]x", "", top_slide_group_join_long$group)
top_slide_group_join_long$group = gsub("[.]y", "", top_slide_group_join_long$group)

top_slide_group_join_long <- left_join(top_slide_group_join_long, sites, by = c("Site"))
sites <- remove_rownames(sites) %>% column_to_rownames("Site")

# Reorder sites based on chloride concentrations
top_slide_group_join_long$Site <- reorder(top_slide_group_join_long$Site, top_slide_group_join_long$Chloride)

# Reorder types
top_slide_group_join_long$type <- ordered(top_slide_group_join_long$type, levels = c("Peri", "Phyto"))

# Plot taxonomic group composition across sites (8 x 6 portrait)
(plot.bar.group.facet <- ggplot() +
    geom_col(data = top_slide_group_join_long, aes(x = type, y = value, fill = group)) +
    ylab("Proportion") +
    xlab("Community") +
    facet_wrap(~Site, nrow = 6) +
    workingtheme_bars +
    guides(fill = guide_legend(title = "Group")) +
    theme(panel.grid.major = element_blank(),
          legend.title=element_blank(),
          legend.position = "bottom"))


##################
# Run regression models
##################

### Filter pigments (phytoplankton)
## Chlorophyll.a
# Create data frame with chlorophyll a data from filters (phytoplankton)
top_filter_chla <- subset(pigments, Type == "Top") %>%
  select(!Type) %>%
  remove_rownames() %>%
  column_to_rownames(var = "Site")

top_filter_chla <- rownames_to_column(top_filter_chla, var = "Site")

top_filter_chla <- left_join(top_filter_chla, site_select, by = c("Site"))

top_filter_chla <- column_to_rownames(top_filter_chla, "Site")
site_select <- remove_rownames(site_select) %>% column_to_rownames("Site")

# Run Chlorophyll.a GLM
fit.chla.filter <- glm(Chlorophyll.a ~
                         scale(log(Chloride)) + 
                         scale(log(`Total Phosphorus (as P)`)) +
                         scale(log(`NO3 + NO2 (as N)`)) +
                         scale(log(`Total Suspended Solids`)) +
                         scale(`Field pH`) +
                         scale(log(`Field Temp`)) +
                         scale(log(`Total Copper`)),
                       data = top_filter_chla, family = Gamma(link = "log"))

summary(fit.chla.filter)

# Run P-value adjustments
p.adjust(coef(summary(fit.chla.filter))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.chla.filter, terms = "Chloride [all]")

# Plot Chlorophyll.a against chloride (4 x 4)
(plot.filter.Chlorophyll.a.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_filter_chla, alpha = 0.8, size = 4, aes(x = Chloride, y = Chlorophyll.a)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Chlorophyll a (μg/L)") +
    coord_cartesian(ylim = c(0, 114.9816)) +
    workingtheme_regressions)

##Secondary pigments
# Create site column from row names
top_filter_total <- rownames_to_column(top_filter_total, var = "Site")
site_select <- rownames_to_column(site_select, var = "Site")

# Join phytoplankton pigment data with site data based on site
top_filter_select <- left_join(top_filter_total, site_select, by = c("Site"))

# Reset row names
top_filter_total <- column_to_rownames(top_filter_total, "Site")
site_select <- column_to_rownames(site_select, "Site")

## Fucoxanthin
# Run Fucoxanthin GLM
fit.fuco.filter <- glm(Fucoxanthin ~
                         scale(log(Chloride)) + 
                         scale(log(`Total Phosphorus (as P)`)) +
                         scale(log(`NO3 + NO2 (as N)`)) +
                         scale(log(`Total Suspended Solids`)) +
                         scale(`Field pH`) +
                         scale(log(`Field Temp`)) +
                         scale(log(`Total Copper`)),
                       data = top_filter_select, family = quasibinomial(link = "logit"))

summary(fit.fuco.filter)

# Run P-value adjustments
p.adjust(coef(summary(fit.fuco.filter))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.fuco.filter, terms = "Chloride [all]")

# Plot Fucoxanthin against chloride (4 x 4)
(plot.filter.Fucoxanthin.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_filter_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Fucoxanthin)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Fucoxanthin (proportion)") +
    coord_cartesian(ylim = c(0, 0.4068113)) +
    workingtheme_regressions)

## Neoxanthin
# Run Neoxanthin GLM
fit.neox.filter <- glm(Neoxanthin ~
                         scale(log(Chloride)) + 
                         scale(log(`Total Phosphorus (as P)`)) +
                         scale(log(`NO3 + NO2 (as N)`)) +
                         scale(log(`Total Suspended Solids`)) +
                         scale(`Field pH`) +
                         scale(log(`Field Temp`)) +
                         scale(log(`Total Copper`)),
                       data = top_filter_select, family = quasibinomial(link = "logit"))

summary(fit.neox.filter)

# Run P-value adjustments
p.adjust(coef(summary(fit.neox.filter))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.neox.filter, terms = "Chloride [all]")

# Plot Neoxanthin against chloride (4 x 4)
(plot.filter.Neoxanthin.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_filter_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Neoxanthin)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Neoxanthin (proportion)") +
    coord_cartesian(ylim = c(0, 0.02785201)) +
    workingtheme_regressions)

## Violaxanthin
# Run Violaxanthin GLM
fit.viol.filter <- glm(Violaxanthin ~
                         scale(log(Chloride)) + 
                         scale(log(`Total Phosphorus (as P)`)) +
                         scale(log(`NO3 + NO2 (as N)`)) +
                         scale(log(`Total Suspended Solids`)) +
                         scale(`Field pH`) +
                         scale(log(`Field Temp`)) +
                         scale(log(`Total Copper`)),
                       data = top_filter_select, family = quasibinomial(link = "logit"))

summary(fit.viol.filter)

# Run P-value adjustments
p.adjust(coef(summary(fit.viol.filter))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.viol.filter, terms = "Chloride [all]")

# Plot Violaxanthin against chloride (4 x 4)
(plot.filter.Violaxanthin.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_filter_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Violaxanthin)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Violaxanthin (proportion)") +
    coord_cartesian(ylim = c(0, 0.05061069)) +
    workingtheme_regressions)

### Diadinoxanthin
# Run Diadinoxanthin GLM
fit.diad.filter <- glm(Diadinoxanthin ~
                         scale(log(Chloride)) + 
                         scale(log(`Total Phosphorus (as P)`)) +
                         scale(log(`NO3 + NO2 (as N)`)) +
                         scale(log(`Total Suspended Solids`)) +
                         scale(`Field pH`) +
                         scale(log(`Field Temp`)) +
                         scale(log(`Total Copper`)),
                       data = top_filter_select, family = quasibinomial(link = "logit"))

summary(fit.diad.filter)

# Run P-value adjustments
p.adjust(coef(summary(fit.diad.filter))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.diad.filter, terms = "Chloride [all]")

# Plot Diadinoxanthin against chloride (4 x 4)
(plot.filter.Diadinoxanthin.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_filter_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Diadinoxanthin)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Diadinoxanthin (proportion)") +
    coord_cartesian(ylim = c(0, 0.9809331)) +
    workingtheme_regressions)

## Myxoxanthophyll
# Run Myxoxanthophyll GLM
fit.myxo.filter <- glm(Myxoxanthophyll ~
                         scale(log(Chloride)) + 
                         scale(log(`Total Phosphorus (as P)`)) +
                         scale(log(`NO3 + NO2 (as N)`)) +
                         scale(log(`Total Suspended Solids`)) +
                         scale(`Field pH`) +
                         scale(log(`Field Temp`)) +
                         scale(log(`Total Copper`)),
                       data = top_filter_select, family = quasibinomial(link = "logit"))

summary(fit.myxo.filter)

# Run P-value adjustments
p.adjust(coef(summary(fit.myxo.filter))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.myxo.filter, terms = "Chloride [all]")

# Plot Myxoxanthophyll against chloride (4 x 4)
(plot.filter.Myxoxanthophyll.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_filter_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Myxoxanthophyll)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Myxoxanthophyll (proportion)") +
    coord_cartesian(ylim = c(0, 0.8458461)) +
    workingtheme_regressions)

## Alloxanthin
# Run Alloxanthin GLM
fit.allo.filter <- glm(Alloxanthin ~
                         scale(log(Chloride)) + 
                         scale(log(`Total Phosphorus (as P)`)) +
                         scale(log(`NO3 + NO2 (as N)`)) +
                         scale(log(`Total Suspended Solids`)) +
                         scale(`Field pH`) +
                         scale(log(`Field Temp`)) +
                         scale(log(`Total Copper`)),
                       data = top_filter_select, family = quasibinomial(link = "logit"))

summary(fit.allo.filter)

# Run P-value adjustments
p.adjust(coef(summary(fit.allo.filter))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.allo.filter, terms = "Chloride [all]")

# Plot Alloxanthin against chloride (4 x 4)
(plot.filter.Alloxanthin.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_filter_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Alloxanthin)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Alloxanthin (proportion)") +
    coord_cartesian(ylim = c(0, 0.7309081)) +
    workingtheme_regressions)

## Diatoxanthin
# Run Diatoxanthin GLM
fit.diat.filter <- glm(Diatoxanthin ~
                         scale(log(Chloride)) + 
                         scale(log(`Total Phosphorus (as P)`)) +
                         scale(log(`NO3 + NO2 (as N)`)) +
                         scale(log(`Total Suspended Solids`)) +
                         scale(`Field pH`) +
                         scale(log(`Field Temp`)) +
                         scale(log(`Total Copper`)),
                       data = top_filter_select, family = quasibinomial(link = "logit"))

summary(fit.diat.filter)

# Run P-value adjustments
p.adjust(coef(summary(fit.diat.filter))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.diat.filter, terms = "Chloride [all]")

# Plot Diatoxanthin against chloride (4 x 4)
(plot.filter.Diatoxanthin.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_filter_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Diatoxanthin)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Diatoxanthin (proportion)") +
    coord_cartesian(ylim = c(0, 0.1060211)) +
    workingtheme_regressions)

## Lutein
# Run Lutein GLM
fit.lute.filter <- glm(Lutein ~
                         scale(log(Chloride)) + 
                         scale(log(`Total Phosphorus (as P)`)) +
                         scale(log(`NO3 + NO2 (as N)`)) +
                         scale(log(`Total Suspended Solids`)) +
                         scale(`Field pH`) +
                         scale(log(`Field Temp`)) +
                         scale(log(`Total Copper`)),
                       data = top_filter_select, family = quasibinomial(link = "logit"))

summary(fit.lute.filter)

# Run P-value adjustments
p.adjust(coef(summary(fit.lute.filter))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.lute.filter, terms = "Chloride [all]")

# Plot Lutein against chloride (4 x 4)
(plot.filter.Lutein.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_filter_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Lutein)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Lutein (proportion)") +
    coord_cartesian(ylim = c(0, 0.2460803)) +
    workingtheme_regressions)

## Zeaxanthin
# Run Zeaxanthin GLM
fit.zeax.filter <- glm(Zeaxanthin ~
                         scale(log(Chloride)) + 
                         scale(log(`Total Phosphorus (as P)`)) +
                         scale(log(`NO3 + NO2 (as N)`)) +
                         scale(log(`Total Suspended Solids`)) +
                         scale(`Field pH`) +
                         scale(log(`Field Temp`)) +
                         scale(log(`Total Copper`)),
                       data = top_filter_select, family = quasibinomial(link = "logit"))

summary(fit.zeax.filter)

# Run P-value adjustments
p.adjust(coef(summary(fit.zeax.filter))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.zeax.filter, terms = "Chloride [all]")

# Plot Zeaxanthin against chloride (4 x 4)
(plot.filter.Zeaxanthin.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_filter_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Zeaxanthin)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Zeaxanthin (proportion)") +
    coord_cartesian(ylim = c(0, 0.6560729)) +
    workingtheme_regressions)

## Lutein.Zeaxanthin
# Run Lutein.Zeaxanthin GLM
fit.luze.filter <- glm(Lutein.Zeaxanthin ~
                         scale(log(Chloride)) + 
                         scale(log(`Total Phosphorus (as P)`)) +
                         scale(log(`NO3 + NO2 (as N)`)) +
                         scale(log(`Total Suspended Solids`)) +
                         scale(`Field pH`) +
                         scale(log(`Field Temp`)) +
                         scale(log(`Total Copper`)),
                       data = top_filter_select, family = quasibinomial(link = "logit"))

summary(fit.luze.filter)

# Run P-value adjustments
p.adjust(coef(summary(fit.luze.filter))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.luze.filter, terms = "Chloride [all]")

# Plot Lutein.Zeaxanthin against chloride (4 x 4)
(plot.filter.Lutein.Zeaxanthin.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_filter_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Lutein.Zeaxanthin)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)", limits = c(0, 1600), breaks = c(0, 400, 800, 1200, 1600)) +
    scale_y_continuous("Lutein.Zeaxanthin (proportion)",limits = c(0, 0), breaks = c(0, 0.3, 0.6, 0.9)) +
    workingtheme_regressions)

## Canthaxanthin
# Run Canthaxanthin GLM
fit.cant.filter <- glm(Canthaxanthin ~
                         scale(log(Chloride)) + 
                         scale(log(`Total Phosphorus (as P)`)) +
                         scale(log(`NO3 + NO2 (as N)`)) +
                         scale(log(`Total Suspended Solids`)) +
                         scale(`Field pH`) +
                         scale(log(`Field Temp`)) +
                         scale(log(`Total Copper`)),
                       data = top_filter_select, family = quasibinomial(link = "logit"))

summary(fit.cant.filter)

# Run P-value adjustments
p.adjust(coef(summary(fit.cant.filter))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.cant.filter, terms = "Chloride [all]")

# Plot Canthaxanthin against chloride (4 x 4)
(plot.filter.Canthaxanthin.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_filter_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Canthaxanthin)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Canthaxanthin (proportion)") +
    coord_cartesian(ylim = c(0, 0.1735692)) +
    workingtheme_regressions)

## Chlorophyll.b
# Run Chlorophyll.b GLM
fit.chlb.filter <- glm(Chlorophyll.b ~
                         scale(log(Chloride)) + 
                         scale(log(`Total Phosphorus (as P)`)) +
                         scale(log(`NO3 + NO2 (as N)`)) +
                         scale(log(`Total Suspended Solids`)) +
                         scale(`Field pH`) +
                         scale(log(`Field Temp`)) +
                         scale(log(`Total Copper`)),
                       data = top_filter_select, family = quasibinomial(link = "logit"))

summary(fit.chlb.filter)

# Run P-value adjustments
p.adjust(coef(summary(fit.chlb.filter))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.chlb.filter, terms = "Chloride [all]")

# Plot Chlorophyll.b against chloride (4 x 4)
(plot.filter.Chlorophyll.b.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_filter_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Chlorophyll.b)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Chlorophyll b (proportion)") +
    coord_cartesian(ylim = c(0, 0.7939238)) +
    workingtheme_regressions)


### Filter taxonomic groups (phytoplankton)
# Create site column from row names
top_filter_groups_total <- rownames_to_column(top_filter_groups_total, var = "Site")
site_select <- rownames_to_column(site_select, var = "Site")

# Join phytoplankton pigment data with site data based on site
top_filter_groups_select <- left_join(top_filter_groups_total, site_select, by = c("Site"))

# Reset row names
top_filter_groups_select <- column_to_rownames(top_filter_groups_select, "Site")
site_select <- column_to_rownames(site_select, "Site")

## Chromophytes
# Run Chromophytes GLM
fit.chrom.filter <- glm(Chromophytes ~
                         scale(log(Chloride)) + 
                         scale(log(`Total Phosphorus (as P)`)) +
                         scale(log(`NO3 + NO2 (as N)`)) +
                         scale(log(`Total Suspended Solids`)) +
                         scale(`Field pH`) +
                         scale(log(`Field Temp`)) +
                          scale(log(`Total Copper`)),
                        data = top_filter_groups_select, family = quasibinomial(link = "logit"))

summary(fit.chrom.filter)

# Run P-value adjustments
p.adjust(coef(summary(fit.chrom.filter))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.chrom.filter, terms = "Chloride [all]")

# Plot Chromophytes against chloride (4 x 4)
(plot.filter.Chromophytes.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_filter_groups_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Chromophytes)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Chromophytes (proportion)") +
    coord_cartesian(ylim = c(0, 0.9809331)) +
    workingtheme_regressions)

## Chlorophytes 
# Run Chlorophytes GLM
fit.chlo.filter <- glm(Chlorophytes ~
                          scale(log(Chloride)) + 
                          scale(log(`Total Phosphorus (as P)`)) +
                          scale(log(`NO3 + NO2 (as N)`)) +
                          scale(log(`Total Suspended Solids`)) +
                          scale(`Field pH`) +
                          scale(log(`Field Temp`)) +
                         scale(log(`Total Copper`)),
                       data = top_filter_groups_select, family = quasibinomial(link = "logit"))

summary(fit.chlo.filter)

# Run P-value adjustments
p.adjust(coef(summary(fit.chlo.filter))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.chlo.filter, terms = "Chloride [all]")

# Plot Chlorophytes against chloride (4 x 4)
(plot.filter.Chlorophytes.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_filter_groups_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Chlorophytes)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Chlorophytes (proportion)") +
    coord_cartesian(ylim = c(0, 0.9609725)) +
    workingtheme_regressions)

## Cryptophytes 
# Run Cryptophytes GLM
fit.cryp.filter <- glm(Cryptophytes ~
                         scale(log(Chloride)) + 
                         scale(log(`Total Phosphorus (as P)`)) +
                         scale(log(`NO3 + NO2 (as N)`)) +
                         scale(log(`Total Suspended Solids`)) +
                         scale(`Field pH`) +
                         scale(log(`Field Temp`)) +
                         scale(log(`Total Copper`)),
                       data = top_filter_groups_select, family = quasibinomial(link = "logit"))

summary(fit.cryp.filter)

# Run P-value adjustments
p.adjust(coef(summary(fit.cryp.filter))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.cryp.filter, terms = "Chloride [all]")

# Plot Cryptophytes against chloride (4 x 4)
(plot.filter.Cryptophytes.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_filter_groups_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Cryptophytes)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Cryptophytes (proportion)") +
    coord_cartesian(ylim = c(0, 0.7309081)) +
    workingtheme_regressions)

## Cyanobacteria
# Run Cyanobacteria GLM
fit.cyan.filter <- glm(Cyanobacteria ~
                         scale(log(Chloride)) + 
                         scale(log(`Total Phosphorus (as P)`)) +
                         scale(log(`NO3 + NO2 (as N)`)) +
                         scale(log(`Total Suspended Solids`)) +
                         scale(`Field pH`) +
                         scale(log(`Field Temp`)) +
                         scale(log(`Total Copper`)),
                       data = top_filter_groups_select, family = quasibinomial(link = "logit"))

summary(fit.cyan.filter)

# Run P-value adjustments
p.adjust(coef(summary(fit.cyan.filter))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.cyan.filter, terms = "Chloride [all]")

# Plot Cyanobacteria against chloride (4 x 4)
(plot.filter.Cyanobacteria.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_filter_groups_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Cyanobacteria)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Cyanobacteria (proportion)") +
    coord_cartesian(ylim = c(0, 0.8949164)) +
    workingtheme_regressions)

### Slides (periphyton)
## Chlorophyll.a
# Create data frame with chlorophyll a data from slides (periphyton)
top_slide_chla <- subset(pigments, Type == "Slides") %>%
  select(!Type) %>%
  remove_rownames() %>%
  column_to_rownames(var = "Site")

top_slide_chla <- rownames_to_column(top_slide_chla, var = "Site")
site_select <- rownames_to_column(site_select, var = "Site")

top_slide_chla <- left_join(top_slide_chla, site_select, by = c("Site"))

top_slide_chla <- column_to_rownames(top_slide_chla, "Site")
site_select <- column_to_rownames(site_select, "Site")

# Run Chlorophyll.a GLM
fit.chla.slide <- glm(Chlorophyll.a ~
                         scale(log(Chloride)) + 
                         scale(log(`Total Phosphorus (as P)`)) +
                         scale(log(`NO3 + NO2 (as N)`)) +
                         scale(log(`Total Suspended Solids`)) +
                         scale(`Field pH`) +
                         scale(log(`Field Temp`)) +
                        scale(log(`Total Copper`)),
                      data = top_slide_chla, family = Gamma(link = "log"))

summary(fit.chla.slide)

# Run P-value adjustments
p.adjust(coef(summary(fit.chla.slide))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.chla.slide, terms = "Chloride [all]")

# Plot Chlorophyll.a against chloride (4 x 4)
(plot.slide.Chlorophyll.a.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_slide_chla, alpha = 0.8, size = 4, aes(x = Chloride, y = Chlorophyll.a)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous(bquote('Chlorophyll a ' (μg/cm^2/day))) +
    coord_cartesian(ylim = c(0, 50.32194)) +
    workingtheme_regressions)

##Secondary pigments
# Create site column from row names
top_slide_total <- rownames_to_column(top_slide_total, var = "Site")
site_select <- rownames_to_column(site_select, var = "Site")

# Join periphyton pigment data with site data based on site
top_slide_select <- left_join(top_slide_total, site_select, by = c("Site"))

# Reset row names
top_slide_total <- column_to_rownames(top_slide_total, "Site")
site_select <- column_to_rownames(site_select, "Site")

## Fucoxanthin
# Run Fucoxanthin GLM
fit.fuco.slide <- glm(Fucoxanthin ~
                        scale(log(Chloride)) + 
                        scale(log(`Total Phosphorus (as P)`)) +
                        scale(log(`NO3 + NO2 (as N)`)) +
                        scale(log(`Total Suspended Solids`)) +
                        scale(`Field pH`) +
                        scale(log(`Field Temp`)) +
                        scale(log(`Total Copper`)),
                      data = top_slide_select, family = quasibinomial(link = "logit"))

summary(fit.fuco.slide)

# Run P-value adjustments
p.adjust(coef(summary(fit.fuco.slide))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.fuco.slide, terms = "Chloride [all]")

# Plot Fucoxanthin against chloride (4 x 4)
(plot.slide.Fucoxanthin.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_slide_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Fucoxanthin)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Fucoxanthin (proportion)") +
    coord_cartesian(ylim = c(0, 0.9403374)) +
    workingtheme_regressions)

## Neoxanthin
# Run Neoxanthin GLM
fit.neox.slide <- glm(Neoxanthin ~
                        scale(log(Chloride)) + 
                        scale(log(`Total Phosphorus (as P)`)) +
                        scale(log(`NO3 + NO2 (as N)`)) +
                        scale(log(`Total Suspended Solids`)) +
                        scale(`Field pH`) +
                        scale(log(`Field Temp`)) +
                        scale(log(`Total Copper`)),
                      data = top_slide_select, family = quasibinomial(link = "logit"))

summary(fit.neox.slide)

# Run P-value adjustments
p.adjust(coef(summary(fit.neox.slide))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.neox.slide, terms = "Chloride [all]")

# Plot Neoxanthin against chloride (4 x 4)
(plot.slide.Neoxanthin.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_slide_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Neoxanthin)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Neoxanthin (proportion)") +
    coord_cartesian(ylim = c(0, 0.02157141)) +
    workingtheme_regressions)

## Violaxanthin
# Run Violaxanthin GLM
fit.viol.slide <- glm(Violaxanthin ~
                        scale(log(Chloride)) + 
                        scale(log(`Total Phosphorus (as P)`)) +
                        scale(log(`NO3 + NO2 (as N)`)) +
                        scale(log(`Total Suspended Solids`)) +
                        scale(`Field pH`) +
                        scale(log(`Field Temp`)) +
                        scale(log(`Total Copper`)),
                      data = top_slide_select, family = quasibinomial(link = "logit"))

summary(fit.viol.slide)

# Run P-value adjustments
p.adjust(coef(summary(fit.viol.slide))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.viol.slide, terms = "Chloride [all]")

# Plot Violaxanthin against chloride (4 x 4)
(plot.slide.Violaxanthin.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_slide_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Violaxanthin)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Violaxanthin (proportion)") +
    coord_cartesian(ylim = c(0, 0.04812117)) +
    workingtheme_regressions)

## Diadinoxanthin
# Run Diadinoxanthin GLM
fit.diad.slide <- glm(Diadinoxanthin ~
                        scale(log(Chloride)) + 
                        scale(log(`Total Phosphorus (as P)`)) +
                        scale(log(`NO3 + NO2 (as N)`)) +
                        scale(log(`Total Suspended Solids`)) +
                        scale(`Field pH`) +
                        scale(log(`Field Temp`)) +
                        scale(log(`Total Copper`)),
                      data = top_slide_select, family = quasibinomial(link = "logit"))

summary(fit.diad.slide)

# Run P-value adjustments
p.adjust(coef(summary(fit.diad.slide))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.diad.slide, terms = "Chloride [all]")

# Plot Diadinoxanthin against chloride (4 x 4)
(plot.slide.Diadinoxanthin.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_slide_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Diadinoxanthin)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Diadinoxanthin (proportion)") +
    coord_cartesian(ylim = c(0, 0.1400746)) +
    workingtheme_regressions)

## Myxoxanthophyll
# Run Myxoxanthophyll GLM
fit.myxo.slide <- glm(Myxoxanthophyll ~
                        scale(log(Chloride)) + 
                        scale(log(`Total Phosphorus (as P)`)) +
                        scale(log(`NO3 + NO2 (as N)`)) +
                        scale(log(`Total Suspended Solids`)) +
                        scale(`Field pH`) +
                        scale(log(`Field Temp`)) +
                        scale(log(`Total Copper`)),
                      data = top_slide_select, family = quasibinomial(link = "logit"))

summary(fit.myxo.slide)

# Run P-value adjustments
p.adjust(coef(summary(fit.myxo.slide))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.myxo.slide, terms = "Chloride [all]")

# Plot Myxoxanthophyll against chloride (4 x 4)
(plot.slide.Myxoxanthophyll.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_slide_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Myxoxanthophyll)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Myxoxanthophyll (proportion)") +
    coord_cartesian(ylim = c(0, 0.5381339)) +
    workingtheme_regressions)

## Alloxanthin
# Run Alloxanthin GLM
fit.allo.slide <- glm(Alloxanthin ~
                        scale(log(Chloride)) + 
                        scale(log(`Total Phosphorus (as P)`)) +
                        scale(log(`NO3 + NO2 (as N)`)) +
                        scale(log(`Total Suspended Solids`)) +
                        scale(`Field pH`) +
                        scale(log(`Field Temp`)) +
                        scale(log(`Total Copper`)),
                      data = top_slide_select, family = quasibinomial(link = "logit"))

summary(fit.allo.slide)

# Run P-value adjustments
p.adjust(coef(summary(fit.allo.slide))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.allo.slide, terms = "Chloride [all]")

# Plot Alloxanthin against chloride (4 x 4)
(plot.slide.Alloxanthin.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_slide_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Alloxanthin)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Alloxanthin (proportion)") +
    coord_cartesian(ylim = c(0, 0.03183886)) +
    workingtheme_regressions)

## Diatoxanthin
# Run Diatoxanthin GLM
fit.diat.slide <- glm(Diatoxanthin ~
                        scale(log(Chloride)) + 
                        scale(log(`Total Phosphorus (as P)`)) +
                        scale(log(`NO3 + NO2 (as N)`)) +
                        scale(log(`Total Suspended Solids`)) +
                        scale(`Field pH`) +
                        scale(log(`Field Temp`)) +
                        scale(log(`Total Copper`)),
                      data = top_slide_select, family = quasibinomial(link = "logit"))

summary(fit.diat.slide)

# Run P-value adjustments
p.adjust(coef(summary(fit.diat.slide))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.diat.slide, terms = "Chloride [all]")

# Plot Diatoxanthin against chloride (4 x 4)
(plot.slide.Diatoxanthin.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_slide_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Diatoxanthin)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Diatoxanthin (proportion)") +
    coord_cartesian(ylim = c(0, 0.06932443)) +
    workingtheme_regressions)

## Lutein
# Run Lutein GLM
fit.lute.slide <- glm(Lutein ~
                        scale(log(Chloride)) + 
                        scale(log(`Total Phosphorus (as P)`)) +
                        scale(log(`NO3 + NO2 (as N)`)) +
                        scale(log(`Total Suspended Solids`)) +
                        scale(`Field pH`) +
                        scale(log(`Field Temp`)) +
                        scale(log(`Total Copper`)),
                      data = top_slide_select, family = quasibinomial(link = "logit"))

summary(fit.lute.slide)

# Run P-value adjustments
p.adjust(coef(summary(fit.lute.slide))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.lute.slide, terms = "Chloride [all]")

# Plot Lutein against chloride (4 x 4)
(plot.slide.Lutein.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_slide_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Lutein)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Lutein (proportion)") +
    coord_cartesian(ylim = c(0, 0.2456123)) +
    workingtheme_regressions)

## Zeaxanthin
# Run Zeaxanthin GLM
fit.zeax.slide <- glm(Zeaxanthin ~
                        scale(log(Chloride)) + 
                        scale(log(`Total Phosphorus (as P)`)) +
                        scale(log(`NO3 + NO2 (as N)`)) +
                        scale(log(`Total Suspended Solids`)) +
                        scale(`Field pH`) +
                        scale(log(`Field Temp`)) +
                        scale(log(`Total Copper`)),
                      data = top_slide_select, family = quasibinomial(link = "logit"))

summary(fit.zeax.slide)

# Run P-value adjustments
p.adjust(coef(summary(fit.zeax.slide))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.zeax.slide, terms = "Chloride [all]")

# Plot Zeaxanthin against chloride (4 x 4)
(plot.slide.Zeaxanthin.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_slide_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Zeaxanthin)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Zeaxanthin (proportion)") +
    coord_cartesian(ylim = c(0, 0.2021819)) +
    workingtheme_regressions)

## Lutein.Zeaxanthin
# Run Lutein.Zeaxanthin GLM
fit.luze.slide <- glm(Lutein.Zeaxanthin ~
                        scale(log(Chloride)) + 
                        scale(log(`Total Phosphorus (as P)`)) +
                        scale(log(`NO3 + NO2 (as N)`)) +
                        scale(log(`Total Suspended Solids`)) +
                        scale(`Field pH`) +
                        scale(log(`Field Temp`)) +
                        scale(log(`Total Copper`)),
                      data = top_slide_select, family = quasibinomial(link = "logit"))

summary(fit.luze.slide)

# Run P-value adjustments
p.adjust(coef(summary(fit.luze.slide))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.luze.slide, terms = "Chloride [all]")

# Plot Lutein.Zeaxanthin against chloride (4 x 4)
(plot.slide.Lutein.Zeaxanthin.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_slide_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Lutein.Zeaxanthin)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Lutein.Zeaxanthin (proportion)") +
    coord_cartesian(ylim = c(0, 0.1578048)) +
    workingtheme_regressions)

## Canthaxanthin
# Run Canthaxanthin GLM
fit.cant.slide <- glm(Canthaxanthin ~
                        scale(log(Chloride)) + 
                        scale(log(`Total Phosphorus (as P)`)) +
                        scale(log(`NO3 + NO2 (as N)`)) +
                        scale(log(`Total Suspended Solids`)) +
                        scale(`Field pH`) +
                        scale(log(`Field Temp`)) +
                        scale(log(`Total Copper`)),
                      data = top_slide_select, family = quasibinomial(link = "logit"))

summary(fit.cant.slide)

# Run P-value adjustments
p.adjust(coef(summary(fit.cant.slide))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.cant.slide, terms = "Chloride [all]")

# Plot Canthaxanthin against chloride (4 x 4)
(plot.slide.Canthaxanthin.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_slide_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Canthaxanthin)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Canthaxanthin (proportion)") +
    coord_cartesian(ylim = c(0, 0.1291782)) +
    workingtheme_regressions)

## Chlorophyll.b
# Run Chlorophyll.b GLM
fit.chlb.slide <- glm(Chlorophyll.b ~
                        scale(log(Chloride)) + 
                        scale(log(`Total Phosphorus (as P)`)) +
                        scale(log(`NO3 + NO2 (as N)`)) +
                        scale(log(`Total Suspended Solids`)) +
                        scale(`Field pH`) +
                        scale(log(`Field Temp`)) +
                        scale(log(`Total Copper`)),
                      data = top_slide_select, family = quasibinomial(link = "logit"))

summary(fit.chlb.slide)

# Run P-value adjustments
p.adjust(coef(summary(fit.chlb.slide))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.chlb.slide, terms = "Chloride [all]")

# Plot Chlorophyll.b against chloride (4 x 4)
(plot.slide.Chlorophyllb.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_slide_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Chlorophyll.b)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Chlorophyll b (proportion)") +
    coord_cartesian(ylim = c(0, 0.6302423)) +
    workingtheme_regressions)


### Slide taxonomic groups (periphyton)
# Create site column from row names
top_slide_groups_total <- rownames_to_column(top_slide_groups_total, var = "Site")
site_select <- rownames_to_column(site_select, var = "Site")

# Join phytoplankton pigment data with site data based on site
top_slide_groups_select <- left_join(top_slide_groups_total, site_select, by = c("Site"))

# Reset row names
top_slide_groups_select <- column_to_rownames(top_slide_groups_select, "Site")
site_select <- column_to_rownames(site_select, "Site")

## Chromophytes
# Run Chromophytes GLM
fit.chrom.slide <- glm(Chromophytes ~
                         scale(log(Chloride)) + 
                         scale(log(`Total Phosphorus (as P)`)) +
                         scale(log(`NO3 + NO2 (as N)`)) +
                         scale(log(`Total Suspended Solids`)) +
                         scale(`Field pH`) +
                         scale(log(`Field Temp`)) +
                         scale(log(`Total Copper`)),
                       data = top_slide_groups_select, family = quasibinomial(link = "logit"))

summary(fit.chrom.slide)

# Run P-value adjustments
p.adjust(coef(summary(fit.chrom.slide))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.chrom.slide, terms = "Chloride [all]")

# Plot Chromophytes against chloride (4 x 4)
(plot.slide.Chromophytes.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_slide_groups_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Chromophytes)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Chromophytes (proportion)") +
    coord_cartesian(ylim = c(0, 1)) +
    workingtheme_regressions)

## Chlorophytes 
# Run Chlorophytes GLM
fit.chlo.slide <- glm(Chlorophytes ~
                        scale(log(Chloride)) + 
                        scale(log(`Total Phosphorus (as P)`)) +
                        scale(log(`NO3 + NO2 (as N)`)) +
                        scale(log(`Total Suspended Solids`)) +
                        scale(`Field pH`) +
                        scale(log(`Field Temp`)) +
                        scale(log(`Total Copper`)),
                      data = top_slide_groups_select, family = quasibinomial(link = "logit"))

summary(fit.chlo.slide)

# Run P-value adjustments
p.adjust(coef(summary(fit.chlo.slide))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.chlo.slide, terms = "Chloride [all]")

# Plot Chlorophytes against chloride (4 x 4)
(plot.slide.Chlorophytes.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_slide_groups_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Chlorophytes)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Chlorophytes (proportion)") +
    coord_cartesian(ylim = c(0, 0.8159064)) +
    workingtheme_regressions)

## Cryptophytes 
# Run Cryptophytes GLM
fit.cryp.slide <- glm(Cryptophytes ~
                        scale(log(Chloride)) + 
                        scale(log(`Total Phosphorus (as P)`)) +
                        scale(log(`NO3 + NO2 (as N)`)) +
                        scale(log(`Total Suspended Solids`)) +
                        scale(`Field pH`) +
                        scale(log(`Field Temp`)) +
                        scale(log(`Total Copper`)),
                      data = top_slide_groups_select, family = quasibinomial(link = "logit"))

summary(fit.cryp.slide)

# Run P-value adjustments
p.adjust(coef(summary(fit.cryp.slide))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.cryp.slide, terms = "Chloride [all]")

# Plot Cryptophytes against chloride (4 x 4)
(plot.slide.Cryptophytes.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_slide_groups_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Cryptophytes)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Cryptophytes (proportion)") +
    coord_cartesian(ylim = c(0, 0.03183886)) +
    workingtheme_regressions)

## Cyanobacteria 
# Run Cyanobacteria GLM
fit.cyan.slide <- glm(Cyanobacteria ~
                        scale(log(Chloride)) + 
                        scale(log(`Total Phosphorus (as P)`)) +
                        scale(log(`NO3 + NO2 (as N)`)) +
                        scale(log(`Total Suspended Solids`)) +
                        scale(`Field pH`) +
                        scale(log(`Field Temp`)) +
                        scale(log(`Total Copper`)),
                      data = top_slide_groups_select, family = quasibinomial(link = "logit"))

summary(fit.cyan.slide)

# Run P-value adjustments
p.adjust(coef(summary(fit.cyan.slide))[,4], "fdr")

# Predict model fit
fit.pred <- ggpredict(fit.cyan.slide, terms = "Chloride [all]")

# Plot Cyanobacteria against chloride (4 x 4)
(plot.slide.Cyanobacteria.Chloride <- ggplot() +
    geom_ribbon(data = fit.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = top_slide_groups_select, alpha = 0.8, size = 4, aes(x = Chloride, y = Cyanobacteria)) +
    geom_line(data = fit.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    scale_x_continuous("Chloride (mg/L)") +
    scale_y_continuous("Cyanobacteria (proportion)") +
    coord_cartesian(ylim = c(0, 0.6518175)) +
    workingtheme_regressions)


##################
# Plot effect sizes
##################

# Create vector of pigments
pigmentord <- as.factor(c("Diadinoxanthin",
                          "Diatoxanthin",
                          "Fucoxanthin",
                          "Alloxanthin",
                          "Chlorophyll.b",
                          "Lutein",
                          "Neoxanthin",
                          "Violaxanthin",
                          "Canthaxanthin",
                          "Lutein.Zeaxanthin",
                          "Myxoxanthophyll",
                          "Zeaxanthin"))

# Obtain filter pigment regression coefficients 
filter.pig.coef <- c(summary(fit.diad.filter)$coefficients[2, 1],
                     summary(fit.diat.filter)$coefficients[2, 1],
                     summary(fit.fuco.filter)$coefficients[2, 1],
                     summary(fit.allo.filter)$coefficients[2, 1],
                     summary(fit.chlb.filter)$coefficients[2, 1],
                     summary(fit.lute.filter)$coefficients[2, 1],
                     summary(fit.neox.filter)$coefficients[2, 1],
                     summary(fit.viol.filter)$coefficients[2, 1],
                     summary(fit.cant.filter)$coefficients[2, 1],
                     summary(fit.luze.filter)$coefficients[2, 1],
                     summary(fit.myxo.filter)$coefficients[2, 1],
                     summary(fit.zeax.filter)$coefficients[2, 1])

# Obtain filter pigment regression error bars 
filter.pig.se <- c(summary(fit.diad.filter)$coefficients[2, 2],
                   summary(fit.diat.filter)$coefficients[2, 2],
                   summary(fit.fuco.filter)$coefficients[2, 2],
                   summary(fit.allo.filter)$coefficients[2, 2],
                   summary(fit.chlb.filter)$coefficients[2, 2],
                   summary(fit.lute.filter)$coefficients[2, 2],
                   summary(fit.neox.filter)$coefficients[2, 2],
                   summary(fit.viol.filter)$coefficients[2, 2],
                   summary(fit.cant.filter)$coefficients[2, 2],
                   summary(fit.luze.filter)$coefficients[2, 2],
                   summary(fit.myxo.filter)$coefficients[2, 2],
                   summary(fit.zeax.filter)$coefficients[2, 2])

# Indicate significance of filter pigment regression coefficients
filter.cl.allpigs <-  as.data.frame(pigmentord)
filter.cl.allpigs$coef <-  filter.pig.coef
filter.cl.allpigs$se <-  filter.pig.se
filter.cl.allpigs$sig <- as.factor(c("P>0.10", #fit.diad.filter
                                     "P>0.10", #fit.diat.filter
                                     "P>0.10", #fit.fuco.filter
                                     "P>0.10", #fit.allo.filter
                                     "P>0.10", #fit.chlb.filter
                                     "P>0.10", #fit.lute.filter
                                     "P>0.10", #fit.neox.filter
                                     "P>0.10", #fit.viol.filter
                                     "P>0.10", #fit.cant.filter
                                     "P>0.10", #fit.luze.filter
                                     "P>0.10", #fit.myxo.filter
                                     "P<0.05")) #fit.zeax.filter

# Obtain slide pigment regression coefficients 
slide.pig.coef <- c(summary(fit.diad.slide)$coefficients[2, 1],
                    summary(fit.diat.slide)$coefficients[2, 1],
                    summary(fit.fuco.slide)$coefficients[2, 1],
                    summary(fit.allo.slide)$coefficients[2, 1],
                    summary(fit.chlb.slide)$coefficients[2, 1],
                    summary(fit.lute.slide)$coefficients[2, 1],
                    summary(fit.neox.slide)$coefficients[2, 1],
                    summary(fit.viol.slide)$coefficients[2, 1],
                    summary(fit.cant.slide)$coefficients[2, 1],
                    summary(fit.luze.slide)$coefficients[2, 1],
                    summary(fit.myxo.slide)$coefficients[2, 1],
                    summary(fit.zeax.slide)$coefficients[2, 1])

# Obtain slide pigment regression error bars 
slide.pig.se <- c(summary(fit.diad.slide)$coefficients[2, 2],
                  summary(fit.diat.slide)$coefficients[2, 2],
                  summary(fit.fuco.slide)$coefficients[2, 2],
                  summary(fit.allo.slide)$coefficients[2, 2],
                  summary(fit.chlb.slide)$coefficients[2, 2],
                  summary(fit.lute.slide)$coefficients[2, 2],
                  summary(fit.neox.slide)$coefficients[2, 2],
                  summary(fit.viol.slide)$coefficients[2, 2],
                  summary(fit.cant.slide)$coefficients[2, 2],
                  summary(fit.luze.slide)$coefficients[2, 2],
                  summary(fit.myxo.slide)$coefficients[2, 2],
                  summary(fit.zeax.slide)$coefficients[2, 2])

# Indicate significance of slide pigment regression coefficients
slide.cl.allpigs <-  as.data.frame(pigmentord)
slide.cl.allpigs$coef <- slide.pig.coef
slide.cl.allpigs$se <- slide.pig.se
slide.cl.allpigs$sig <- as.factor(c("P>0.10", #fit.diad.slide
                                    "P>0.10", #fit.diat.slide
                                    "P<0.05", #fit.fuco.slide
                                    "P>0.10", #fit.allo.slide
                                    "P>0.10", #fit.chlb.slide
                                    "P>0.10", #fit.lute.slide
                                    "P<0.10", #fit.neox.slide
                                    "P>0.10", #fit.viol.slide
                                    "P<0.05", #fit.cant.slide
                                    "P>0.10", #fit.luze.slide
                                    "P>0.10", #fit.myxo.slide
                                    "P>0.10")) #fit.zeax.slide

# Reorder pigment factor
slide.cl.allpigs$pigmentord <- ordered(slide.cl.allpigs$pigmentord, levels = c("Diadinoxanthin",
                                                                               "Diatoxanthin",
                                                                               "Fucoxanthin",
                                                                               "Alloxanthin",
                                                                               "Chlorophyll.b",
                                                                               "Lutein",
                                                                               "Neoxanthin",
                                                                               "Violaxanthin",
                                                                               "Canthaxanthin",
                                                                               "Lutein.Zeaxanthin",
                                                                               "Myxoxanthophyll",
                                                                               "Zeaxanthin"))

# Reorder significance factor
slide.cl.allpigs$sig <- ordered(slide.cl.allpigs$sig, levels = c("P>0.10", "P<0.10", "P<0.05"))

# Plot standardized chloride slope coefficients for pigments (8 x 6 landscape)
(plot.cl.allpigs <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
    geom_errorbar(data = slide.cl.allpigs, width = 0.2, color = "grey",
                  aes(x = pigmentord, ymin = coef - se, ymax = coef + se)) +
    geom_errorbar(data = filter.cl.allpigs, width = 0.2, color = "grey",
                  aes(x = pigmentord, ymin = coef - se, ymax = coef + se)) +
    geom_point(data = slide.cl.allpigs, size = 4, shape = 15, aes(x = pigmentord, y = coef, color = sig)) +
    geom_point(data = filter.cl.allpigs, size = 4, shape = 16, aes(x = pigmentord, y = coef, color = sig)) +
    scale_color_manual(values = c("black", "orange", "red")) +
    scale_x_discrete("Pigments") +
    scale_y_continuous("Standardized Partial Chloride Slope Coefficient") +
    coord_cartesian(ylim = c(-1, 3.1)) +
    guides(color = guide_legend(title = "Significance")) +
    workingtheme_regressions)

# Create vector of taxonomic groups
taxaord <- c("Chromophytes",
             "Cryptophytes",
             "Chlorophytes",
             "Cyanobacteria")

# Obtain filter taxonomic group regression coefficients 
filter.taxa.coef <- c(summary(fit.chrom.filter)$coefficients[2, 1],
                      summary(fit.cryp.filter)$coefficients[2, 1],
                      summary(fit.chlo.filter)$coefficients[2, 1],
                      summary(fit.cyan.filter)$coefficients[2, 1])

# Obtain filter taxonomic group regression error bars 
filter.taxa.se <- c(summary(fit.chrom.filter)$coefficients[2, 2],
                    summary(fit.cryp.filter)$coefficients[2, 2],
                    summary(fit.chlo.filter)$coefficients[2, 1],
                    summary(fit.cyan.filter)$coefficients[2, 2])

# Indicate significance of filter taxonomic group regression coefficients
filter.cl.alltaxa <-  as.data.frame(taxaord)
filter.cl.alltaxa$coef <- filter.taxa.coef
filter.cl.alltaxa$se <- filter.taxa.se
filter.cl.alltaxa$sig <- as.factor(c("P>0.10", #fit.chrom.filter
                                     "P>0.10", #fit.cryp.filter
                                     "P<0.10", #fit.chlo.filter
                                     "P<0.05")) #fit.cyan.filter

# Obtain slide taxonomic group regression coefficients 
slide.taxa.coef <- c(summary(fit.chrom.slide)$coefficients[2, 1],
                     summary(fit.cryp.slide)$coefficients[2, 1],
                     summary(fit.chlo.slide)$coefficients[2, 1],
                     summary(fit.cyan.slide)$coefficients[2, 1])

# Obtain slide taxonomic group regression error bars 
slide.taxa.se <- c(summary(fit.chrom.slide)$coefficients[2, 2],
                   summary(fit.cryp.slide)$coefficients[2, 2],
                   summary(fit.chlo.slide)$coefficients[2, 2],
                   summary(fit.cyan.slide)$coefficients[2, 2])

# Indicate significance of slide taxonomic group regression coefficients
slide.cl.alltaxa <-  as.data.frame(taxaord)
slide.cl.alltaxa$coef <-  slide.taxa.coef
slide.cl.alltaxa$se <-  slide.taxa.se
slide.cl.alltaxa$sig <- as.factor(c("P<0.05",
                                    "P>0.10",
                                    "P>0.10",
                                    "P<0.10"))

# Reorder taxonomic group factor
slide.cl.alltaxa$taxaord <- ordered(slide.cl.alltaxa$taxaord, levels = c("Chromophytes",
                                                                         "Cryptophytes",
                                                                         "Chlorophytes",
                                                                         "Cyanobacteria"))

# Reorder significance factor
slide.cl.alltaxa$sig <- ordered(slide.cl.alltaxa$sig, levels = c("P>0.10", "P<0.10", "P<0.05"))

# Plot standardized chloride slope coefficients for taxonomic groups (6 x 4 portrait)
(plot.cl.alltaxa <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
    geom_errorbar(data = slide.cl.alltaxa, width = 0.2, color = "grey",
                  aes(x = taxaord, ymin = coef - se, ymax = coef + se)) +
    geom_errorbar(data = filter.cl.alltaxa, width = 0.2, color = "grey",
                  aes(x = taxaord, ymin = coef - se, ymax = coef + se)) +
    geom_point(data = slide.cl.alltaxa, size = 4, shape = 15, aes(x = taxaord, y = coef, color = sig)) +
    geom_point(data = filter.cl.alltaxa, size = 4, shape = 16, aes(x = taxaord, y = coef, color = sig)) +
    scale_color_manual(values = c("black", "orange", "red")) +
    scale_x_discrete("Group") +
    scale_y_continuous("Standardized Partial Chloride Slope Coefficient") +
    coord_cartesian(ylim = c(-1, 3.1)) +
    guides(color = guide_legend(title = "Significance")) +
    workingtheme_regressions)

# Create chlorophyll a label
chlaord <- "Chlorophyll.a"

# Obtain filter chlorophyll a regression coefficient
filter.chla.coef <- summary(fit.chla.filter)$coefficients[2, 1]

# Obtain filter chlorophyll a regression error bars 
filter.chla.se <- summary(fit.chla.filter)$coefficients[2, 2]

# Indicate significance of filter chlorophyll a regression coefficients
filter.cl.chla <-  as.data.frame(chlaord)
filter.cl.chla$coef <- filter.chla.coef
filter.cl.chla$se <- filter.chla.se
filter.cl.chla$sig <- as.factor("P>0.10") #fit.chla.filter

# Obtain slide chlorophyll a regression coefficients 
slide.chla.coef <- summary(fit.chla.slide)$coefficients[2, 1]

# Obtain slide chlorophyll a regression error bars 
slide.chla.se <- summary(fit.chla.slide)$coefficients[2, 2]

# Indicate significance of slide chlorophyll a regression coefficients
slide.cl.chla <-  as.data.frame(chlaord)
slide.cl.chla$coef <-  slide.chla.coef
slide.cl.chla$se <- slide.chla.se
slide.cl.chla$sig <- as.factor("P<0.05") #fit.chla.slide

# Reorder significance factor
slide.cl.chla$sig <- ordered(slide.cl.chla$sig, levels = c("P>0.10", "P<0.10", "P<0.05"))

# Plot standardized chloride slope coefficients for chlorophyll a (6 x 3 portrait)
(plot.cl.chla <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
    geom_errorbar(data = slide.cl.chla, width = 0.2, color = "grey",
                  aes(x = chlaord, ymin = coef - se, ymax = coef + se)) +
    geom_errorbar(data = filter.cl.chla, width = 0.2, color = "grey",
                  aes(x = chlaord, ymin = coef - se, ymax = coef + se)) +
    geom_point(data = slide.cl.chla, size = 4, shape = 15, aes(x = chlaord, y = coef, color = sig)) +
    geom_point(data = filter.cl.chla, size = 4, shape = 16, aes(x = chlaord, y = coef, color = sig)) +
    scale_color_manual(values = c("red", "black")) +
    scale_x_discrete("Total") +
    scale_y_continuous("Standardized Partial Chloride Slope Coefficient") +
    coord_cartesian(ylim = c(-1, 3.1)) +
    guides(color = guide_legend(title = "Significance")) +
    workingtheme_regressions)

# Create vector of covariates
covariateord <- c("Chloride",
                "Total Phosphorus (as P)",
                "NO3 + NO2 (as N)",
                "Total Suspended Solids",
                "Field pH",
                "Field Temp",
                "Total Copper")

# Obtain filter cyanobacteria chloride regression coefficients 
filter.cyan.coef <- c(summary(fit.cyan.filter)$coefficients[2, 1],
                      summary(fit.cyan.filter)$coefficients[3, 1],
                      summary(fit.cyan.filter)$coefficients[4, 1],
                      summary(fit.cyan.filter)$coefficients[5, 1],
                      summary(fit.cyan.filter)$coefficients[6, 1],
                      summary(fit.cyan.filter)$coefficients[7, 1],
                      summary(fit.cyan.filter)$coefficients[8, 1])

# Obtain filter cyanobacteria chloride regression error bars 
filter.cyan.se <- c(summary(fit.cyan.filter)$coefficients[2, 2],
                    summary(fit.cyan.filter)$coefficients[3, 2],
                    summary(fit.cyan.filter)$coefficients[4, 2],
                    summary(fit.cyan.filter)$coefficients[5, 2],
                    summary(fit.cyan.filter)$coefficients[6, 2],
                    summary(fit.cyan.filter)$coefficients[7, 2],
                    summary(fit.cyan.filter)$coefficients[8, 2])

# Indicate significance of filter cyanobacteria chloride regression coefficients
filter.cyano.all <-  as.data.frame(covariateord)
filter.cyano.all$coef <- filter.cyan.coef
filter.cyano.all$se <- filter.cyan.se
filter.cyano.all$sig <- as.factor(c("P<0.05", 
                                    "P>0.10",
                                    "P>0.10",
                                    "P>0.10",
                                    "P>0.10",
                                    "P>0.10",
                                    "P>0.10"))

# Obtain slide cyanobacteria chloride regression coefficients 
slide.cyan.coef <- c(summary(fit.cyan.slide)$coefficients[2, 1],
                     summary(fit.cyan.slide)$coefficients[3, 1],
                     summary(fit.cyan.slide)$coefficients[4, 1],
                     summary(fit.cyan.slide)$coefficients[5, 1],
                     summary(fit.cyan.slide)$coefficients[6, 1],
                     summary(fit.cyan.slide)$coefficients[7, 1],
                     summary(fit.cyan.slide)$coefficients[8, 1])

# Obtain slide cyanobacteria chloride regression error bars 
slide.cyan.se <- c(summary(fit.cyan.slide)$coefficients[2, 2],
                   summary(fit.cyan.slide)$coefficients[3, 2],
                   summary(fit.cyan.slide)$coefficients[4, 2],
                   summary(fit.cyan.slide)$coefficients[5, 2],
                   summary(fit.cyan.slide)$coefficients[6, 2],
                   summary(fit.cyan.slide)$coefficients[7, 2],
                   summary(fit.cyan.slide)$coefficients[8, 2])

# Indicate significance of slide cyanobacteria chloride regression coefficients
slide.cyano.all <-  as.data.frame(covariateord)
slide.cyano.all$coef <- slide.cyan.coef
slide.cyano.all$se <- slide.cyan.se
slide.cyano.all$sig <- as.factor(c("P<0.10", 
                                   "P>0.10", 
                                   "P>0.10", 
                                   "P>0.10", 
                                   "P>0.10", 
                                   "P>0.10", 
                                   "P>0.10"))

# Reorder covariates
filter.cyano.all$covariateord = with(filter.cyano.all, reorder(covariateord, coef, FUN = max))

# Reorder significance factor
slide.cyano.all$sig <- ordered(slide.cyano.all$sig, levels = c("P>0.10", "P<0.10", "P<0.05"))

# Plot standardized chloride slope coefficients for taxonomic groups (9 x 6 landscape)
(plot.cyano.all <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
    geom_errorbar(data = slide.cyano.all, width = 0.2, color = "grey",
                  aes(x = covariateord, ymin = coef - se, ymax = coef + se)) +
    geom_errorbar(data = filter.cyano.all, width = 0.2, color = "grey",
                  aes(x = covariateord, ymin = coef - se, ymax = coef + se)) +
    geom_point(data = slide.cyano.all, size = 4, shape = 15, aes(x = covariateord, y = coef, color = sig)) +
    geom_point(data = filter.cyano.all, size = 4, shape = 16, aes(x = covariateord, y = coef, color = sig)) +
    scale_color_manual(values = c("black", "orange", "red")) +
    scale_x_discrete("Predictors of Cyanobacteria Relative Concentration") +
    scale_y_continuous("Standardized Parital Slope Coefficient", limits = c(-0.9, 1.8), breaks = c(-0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6)) +
    guides(color = guide_legend(title = "Significance")) +
    workingtheme_regressions)


##################
# Plot effect sizes
##################

# Model and plot specific conductance against dissolved chloride
spc.cl <- glm(Conductance ~ Chloride, data = sites, family = gaussian(link = "identity"))
summary(spc.cl)

spc.cl.pred <- ggpredict(spc.cl, terms = "Chloride [all]")

(plot.spc.cl.pred <- ggplot() +
    geom_ribbon(data = spc.cl.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = sites, alpha = 0.8, size = 4, aes(x = Chloride, y = Conductance)) +
    geom_line(data = spc.cl.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    labs(x = "Dissolved chloride (mg/L)",
         y = Specific~conductance~(mu*S/cm)) +
    workingtheme_regressions)

# Model and plot specific conductance against total sodium
spc.na <- glm(Conductance ~ `Total Sodium`, data = sites, family = gaussian(link = "identity"))
summary(spc.na)

spc.na.pred <- ggpredict(spc.na, terms = "Total Sodium [all]")

(plot.spc.na.pred <- ggplot() +
    geom_ribbon(data = spc.na.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = sites, alpha = 0.8, size = 4, aes(x = `Total Sodium`, y = Conductance)) +
    geom_line(data = spc.na.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    labs(x = "Total Sodium (mg/L)",
         y = Specific~conductance~(mu*S/cm)) +
    workingtheme_regressions)

# Model and plot dissolved chloride against total sodium
cl.na <- glm(Chloride ~ `Total Sodium`, data = sites, family = gaussian(link = "identity"))
summary(cl.na)

cl.na.pred <- ggpredict(cl.na, terms = "Total Sodium [all]")

(plot.cl.na.pred <- ggplot() +
    geom_ribbon(data = cl.na.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = sites, alpha = 0.8, size = 4, aes(x = `Total Sodium`, y = Chloride)) +
    geom_line(data = cl.na.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    labs(x = "Total Sodium (mg/L)",
         y = "Dissolved Chloride (mg/L)") +
    workingtheme_regressions)

################## END


