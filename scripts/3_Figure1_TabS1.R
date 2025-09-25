# -------------------------------------------------------------------------
#Created 1.13.25 by Emily Wilson
#last update 09.22.25 by Sawyer Balint
#Title: Plant Species Drive Global Coastal Wetland Methane Fluxes Meta-analysis
#Figure 1

rm(list=ls()) #clear the environment

# install packages --------------------------------------------------------
library(tidyverse) #for data manipulation
library(scales) #for y-axis
library(multcompView) #for multiple comparison tests and displaying letters
library(FSA) #for non-parametric pairwise comparisons
library(cowplot) #for combining subplots
library(patchwork) #for combining subplots
library(ggforce) #for sina
library(ggsci) #for color scales
library(magick) #for combining images

# configure graphing ------------------------------------------------------

#change default aesthetics for point - used for geom_sina
update_geom_defaults("point", list(shape = 21,  size=1, stroke=0.3))

#function for plotting
asinh_trans <- trans_new(
  name = "asinh",
  transform = asinh,
  inverse = sinh
)

#define the dot and bars
median_IQR <- function(x) {
  data.frame(y = median(x), # Median
             ymin = quantile(x)[2], # 1st quartile
             ymax = quantile(x)[4])  # 3rd quartile
}

#graphing theme to reduce repetition
sina_theme <- list(
  theme_classic(),
  geom_sina(
   shape = 21,
   size = 1.2,
   stroke = .3,
   alpha = 0.3),
  scale_y_continuous(
    trans = asinh_trans,
    limits=c(-2000,10000),
    breaks = c(-1000,-100, -10, 0, 10, 100, 1000,10000)),
  labs(y = expression(paste("CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")"))),
  stat_summary(geom = "linerange",
              fun.data = median_IQR, 
              position = position_dodge(width=0.2), 
              linewidth=0.5,
              col = "black"), 
  stat_summary(geom = "point", 
              fun = "median", 
              position = position_dodge(width=0.2),
              shape=19,
              size = 2, 
              col = "black"),
  geom_hline(yintercept = 0, linetype = "dashed"),
  theme(legend.position="none",
        axis.text.x = element_text(size = 7, colour = "black", angle = 30, hjust = .90, family = "Arial"), 
        axis.text.y = element_text(size = 7, colour = "black", family = "Arial"),
        axis.title.x = element_text(size = 7, color = "black", family = "Arial", face="bold"),
        axis.title.y = element_text(size = 7, color = "black", family = "Arial"))
  )

# import data -------------------------------------------------------------

df <- read.csv("raw/ch4_full_dataset.csv")

#figures and kruskall-wallis are on untransformed data
#sinh(yi) puts data back on the orignal scale not asinh scale

# Rename plants for graph -------------------------------------------------

df <- df %>%
  mutate(plant_species = recode(plant_species,
                                "Aeluropus littoralis" = "A. littoralis",
                                "Carex sp." = "Carex. sp.",
                                "Cladium jamaicense" = "C. jamaicense",
                                "Elytrigia repens" = "E. repens",
                                "Juncus sp." = "Juncus sp.",
                                "Plantago maritima" = "Pl. maritima",
                                "Puccinellia maritima" = "P. maritima",
                                "Schoenoplectus americanus" = "S. americanus",
                                "Schoenoplectus triqueter" = "S. triqueter",
                                "Spartina angelica" = "S. angelica",
                                "Spartina foliosa" = "S. foliosa",
                                "Sporobolus virginicus" = "S. virginicus", 
                                "Spartina patens" = "S. patens",
                                "Spartina alterniflora" = "S. alterniflora",
                                "Distichlis Spicata" = "D. spicata",
                                "Phragmites australis" = "P. australis",
                                "Cyperus malaccensis" = "C. malaccensis",
                                "Scirpus mariqueter" = "S. mariqueter",
                                "Suaeda salsa" = "S. salsa",
                                "Tamarix chinensis" = "T. chinensis"
  ))

#only include plants that are from more than one paper
plot1.df <- df %>%
  group_by(plant_species) %>%
  filter(n_distinct(p_num) > 1) %>%
  ungroup()


#Run Kruskal-Wallis and Dunn --------------------------------------------

kruskal <- kruskal.test(sinh(yi) ~ plant_species, data = plot1.df)
print(kruskal)
dunn <- dunnTest(sinh(yi) ~ plant_species, data = plot1.df)
print(dunn)

#extract the p-values
dunn_results <- dunn$res

#create comparison letters
pvals <- dunn_results %>%
  mutate(comparison = gsub(" - ", "-", Comparison)) %>%
  select(comparison, P.adj)

names(pvals) <- c("Comparison", "p.value")

#create a named vector of p-values with group comparisons as names
pval_vector <- setNames(pvals$p.value, pvals$Comparison)

#get compact letter display
letters <- multcompLetters(pval_vector, threshold = 0.01)

cld <- as.data.frame(letters$Letters)
cld <- rownames_to_column(cld)
colnames(cld) <- c("plant_species", "letter")


# Plants Figure ----------------------------------------------------------------

#organize by median value
medians <- plot1.df %>%
  group_by(plant_species) %>%
  summarize(median_value = median(yi, na.rm = TRUE)) %>%
  arrange(desc(median_value)) #order plants by median

plot1.df <- plot1.df %>%
  mutate(plant_species = factor(plant_species, levels = medians$plant_species))

p5 <- ggplot(plot1.df, aes(x=plant_species, y=sinh(yi), color=plant_species, fill=plant_species, group=plant_species))+
  sina_theme+
  scale_fill_manual(values = c("C. malaccensis"="#025766", "S. patens" = "#5c3972", 
                               "S. alterniflora" = "#2FC45D", "P. australis" ="#20bcff",
                               "D. spicata" = "#47724c", "S. salsa" = "#c6b822",
                               "C. jamaicense" =  "#564e95", "Juncus sp." = "#86a291", "Juncae sp" = "lightblue", 
                               "S. portulacastrum" = "#90E669", "P. maritima" = "#30BCAD",
                               "S. mariqueter" = "#60563f", "T. chinensis"="#FF46A2",
                               "Salicornia sp."="black", "Suaeda sp." = "black",
                               "Carex. sp."="#EEAAEE", "S. americanus"="#FF9E4A" ,
                               "S. angelica"="#D62728"))+
  scale_color_manual(values = c("C. malaccensis"="#025766", "S. patens" = "#5c3972", 
                                "S. alterniflora" = "#2FC45D", "P. australis" ="#20bcff",
                                "D. spicata" = "#47724c", "S. salsa" = "#c6b822",
                                "C. jamaicense" =  "#564e95", "Juncus sp." = "#86a291", "Juncae sp" = "lightblue", 
                                "S. portulacastrum" = "#90E669", "P. maritima" = "#30BCAD",
                                "S. mariqueter" = "#60563f", "T. chinensis"="#FF46A2",
                                "Salicornia sp."="black", "Suaeda sp." = "black",
                                "Carex. sp."="#EEAAEE", "S. americanus"="#FF9E4A" ,
                                "S. angelica"="#D62728"))+
  theme(
  axis.text.x = element_text(size = 7, colour = "black", face="italic", angle = 30, hjust = .90, family = "Arial"))+
  labs(x = "Plant Species")
  
p5

# Salinity ----------------------------------------------------------------

plot2.df <- df

#Run Kruskal-Wallis and Dunn --------------------------------------------

kruskal <- kruskal.test(sinh(yi) ~ salinity, data = plot2.df)
print(kruskal)
dunn <- dunnTest(sinh(yi) ~ salinity, data = plot2.df)
print(dunn)
dunn_results <- dunn$res

# Create comparison letters
# dunnTest uses comparisons like "A - B", so we reformat for multcompLetters
pvals <- dunn_results %>%
  mutate(comparison = gsub(" - ", "-", Comparison)) %>%
  select(comparison, P.adj)

names(pvals) <- c("Comparison", "p.value")

# Create a named vector of p-values with group comparisons as names
pval_vector <- setNames(pvals$p.value, pvals$Comparison)

# Get compact letter display
letters <- multcompLetters(pval_vector, threshold = 0.01)

cld <- as.data.frame(letters$Letters)
cld <- rownames_to_column(cld)
colnames(cld) <- c("salinity", "letter")

medians <- plot2.df %>%
  group_by(salinity) %>%
  summarize(median_value = median(yi, na.rm = TRUE)) %>%
  arrange(desc(median_value))

plot2.df <- plot2.df %>%
  mutate(salinity = factor(salinity, levels = medians$salinity))

p4 <- ggplot(plot2.df, aes(x=salinity, y=sinh(yi), color=plant_species, fill=plant_species, group=salinity))+
  sina_theme+
  scale_fill_manual(values = c("C. malaccensis"="#025766", "S. patens" = "#5c3972", 
                               "S. alterniflora" = "#2FC45D", "P. australis" ="#20bcff",
                               "D. spicata" = "#47724c", "S. salsa" = "#c6b822",
                               "C. jamaicense" =  "#564e95", "Juncus sp." = "#86a291", "Juncae sp" = "lightblue", 
                               "S. portulacastrum" = "#90E669", "P. maritima" = "#30BCAD",
                               "S. mariqueter" = "#60563f", "T. chinensis"="#FF46A2",
                               "Salicornia sp."="black", "Suaeda sp." = "black",
                               "Carex. sp."="#EEAAEE", "S. americanus"="#FF9E4A" ,
                               "S. angelica"="#D62728"))+
  scale_color_manual(values = c("C. malaccensis"="#025766", "S. patens" = "#5c3972", 
                                "S. alterniflora" = "#2FC45D", "P. australis" ="#20bcff",
                                "D. spicata" = "#47724c", "S. salsa" = "#c6b822",
                                "C. jamaicense" =  "#564e95", "Juncus sp." = "#86a291", "Juncae sp" = "lightblue", 
                                "S. portulacastrum" = "#90E669", "P. maritima" = "#30BCAD",
                                "S. mariqueter" = "#60563f", "T. chinensis"="#FF46A2",
                                "Salicornia sp."="black", "Suaeda sp." = "black",
                                "Carex. sp."="#EEAAEE", "S. americanus"="#FF9E4A" ,
                                "S. angelica"="#D62728"))+
  labs(x = "Salinity")

p4


# Season ------------------------------------------------------------------

# Clean filtering and overwriting
plot3.df<-plot2.df%>%
  mutate(season=as.factor(season)) %>%
  filter(!is.na(season), season != "annual") %>%
  mutate(season = droplevels(season))


# Kruskal Wallis and Dunn ----------------------------------------------------------------
kruskal <- kruskal.test(sinh(yi) ~ season, data = plot3.df)
print(kruskal)
dunn <- dunnTest(sinh(yi) ~ season, data = plot3.df)
print(dunn)

# Extract the p-values for compact letter display
dunn_results <- dunn$res

# Create comparison letters
pvals <- dunn_results %>%
  mutate(comparison = gsub(" - ", "-", Comparison)) %>%
  select(comparison, P.adj)

names(pvals) <- c("Comparison", "p.value")

# Create a named vector of p-values with group comparisons as names
pval_vector <- setNames(pvals$p.value, pvals$Comparison)

# Get compact letter display
letters <- multcompLetters(pval_vector, threshold = 0.01)

cld <- as.data.frame(letters$Letters)
cld <- rownames_to_column(cld)
colnames(cld) <- c("season", "letter")

medians <- plot3.df %>%
  group_by(season) %>%
  summarize(median_value = median(yi, na.rm = TRUE)) %>%
  arrange(desc(median_value))

plot3.df <- plot3.df %>%
  mutate(season = factor(season, levels = medians$season))

p2 <- ggplot(plot3.df, aes(x=season, y=sinh(yi), color=plant_species, fill=plant_species, group=season))+
  sina_theme+
  scale_fill_manual(values = c("C. malaccensis"="#025766", "S. patens" = "#5c3972", 
                               "S. alterniflora" = "#2FC45D", "P. australis" ="#20bcff",
                               "D. spicata" = "#47724c", "S. salsa" = "#c6b822",
                               "C. jamaicense" =  "#564e95", "Juncus sp." = "#86a291", "Juncae sp" = "lightblue", 
                               "S. portulacastrum" = "#90E669", "P. maritima" = "#30BCAD",
                               "S. mariqueter" = "#60563f", "T. chinensis"="#FF46A2",
                               "Salicornia sp."="black", "Suaeda sp." = "black",
                               "Carex. sp."="#EEAAEE", "S. americanus"="#FF9E4A" ,
                               "S. angelica"="#D62728"))+
  scale_color_manual(values = c("C. malaccensis"="#025766", "S. patens" = "#5c3972", 
                                "S. alterniflora" = "#2FC45D", "P. australis" ="#20bcff",
                                "D. spicata" = "#47724c", "S. salsa" = "#c6b822",
                                "C. jamaicense" =  "#564e95", "Juncus sp." = "#86a291", "Juncae sp" = "lightblue", 
                                "S. portulacastrum" = "#90E669", "P. maritima" = "#30BCAD",
                                "S. mariqueter" = "#60563f", "T. chinensis"="#FF46A2",
                                "Salicornia sp."="black", "Suaeda sp." = "black",
                                "Carex. sp."="#EEAAEE", "S. americanus"="#FF9E4A" ,
                                "S. angelica"="#D62728"))+
  labs(x = "Season")

p2


# Tide ---------------------------------------------------

plot4.df <- df %>%
  mutate(tide = factor(tide, 
                       levels = c("macrotidal", "mesotidal", "microtidal")))

# Kruskal Wallis ----------------------------------------------------------------

kruskal <- kruskal.test(sinh(yi) ~ tide, data = plot4.df)
print(kruskal)
dunn <- dunnTest(sinh(yi) ~ tide, data = plot4.df)
print(dunn)

# Extract the p-values for compact letter display
dunn_results <- dunn$res

# Create comparison letters
# dunnTest uses comparisons like "A - B", so we reformat for multcompLetters
pvals <- dunn_results %>%
  mutate(comparison = gsub(" - ", "-", Comparison)) %>%
  select(comparison, P.adj)

names(pvals) <- c("Comparison", "p.value")

# Create a named vector of p-values with group comparisons as names
pval_vector <- setNames(pvals$p.value, pvals$Comparison)

# Get compact letter display
letters <- multcompLetters(pval_vector, threshold = 0.01)

cld <- as.data.frame(letters$Letters)
cld <- rownames_to_column(cld)
colnames(cld) <- c("tide", "letter")

medians <- plot4.df %>%
  group_by(tide) %>%
  summarize(median_value = median(yi, na.rm = TRUE)) %>%
  arrange(desc(median_value))

plot4.df <- plot4.df %>%
  mutate(tide = factor(tide, levels = medians$tide))

p3 <- ggplot(plot4.df, aes(x=tide, y=sinh(yi), color=plant_species, fill=plant_species, group=tide))+
  sina_theme+
  scale_fill_manual(values = c("C. malaccensis"="#025766", "S. patens" = "#5c3972", 
                               "S. alterniflora" = "#2FC45D", "P. australis" ="#20bcff",
                               "D. spicata" = "#47724c", "S. salsa" = "#c6b822",
                               "C. jamaicense" =  "#564e95", "Juncus sp." = "#86a291", "Juncae sp" = "lightblue", 
                               "S. portulacastrum" = "#90E669", "P. maritima" = "#30BCAD",
                               "S. mariqueter" = "#60563f", "T. chinensis"="#FF46A2",
                               "Salicornia sp."="black", "Suaeda sp." = "black",
                               "Carex. sp."="#EEAAEE", "S. americanus"="#FF9E4A" ,
                               "S. angelica"="#D62728"))+
  scale_color_manual(values = c("C. malaccensis"="#025766", "S. patens" = "#5c3972", 
                                "S. alterniflora" = "#2FC45D", "P. australis" ="#20bcff",
                                "D. spicata" = "#47724c", "S. salsa" = "#c6b822",
                                "C. jamaicense" =  "#564e95", "Juncus sp." = "#86a291", "Juncae sp" = "lightblue", 
                                "S. portulacastrum" = "#90E669", "P. maritima" = "#30BCAD",
                                "S. mariqueter" = "#60563f", "T. chinensis"="#FF46A2",
                                "Salicornia sp."="black", "Suaeda sp." = "black",
                                "Carex. sp."="#EEAAEE", "S. americanus"="#FF9E4A" ,
                                "S. angelica"="#D62728"))+
  labs(x = "Tidal Range")

p3

# Climate Region ----------------------------------------------------------

plot5.df <- df %>%
  mutate(climate_region = recode(climate_region,
                                "subtropical/warm temperate"="warm temperate"))
# Kruskal Wallis ----------------------------------------------------------------
kruskal <- kruskal.test(sinh(yi) ~ climate_region, data = plot5.df)
print(kruskal)

# Run Dunn's test
dunn <- dunnTest(sinh(yi) ~ climate_region, data = plot5.df)
print(dunn)

# Extract the p-values for compact letter display
dunn_results <- dunn$res

# Create comparison letters
# dunnTest uses comparisons like "A - B", so we reformat for multcompLetters
pvals <- dunn_results %>%
  mutate(comparison = gsub(" - ", "-", Comparison)) %>%
  select(comparison, P.adj)

names(pvals) <- c("Comparison", "p.value")

# Create a named vector of p-values with group comparisons as names
pval_vector <- setNames(pvals$p.value, pvals$Comparison)

# Get compact letter display
letters <- multcompLetters(pval_vector, threshold = 0.01)

cld <- as.data.frame(letters$Letters)
cld <- rownames_to_column(cld)
colnames(cld) <- c("Climate Region", "letter")

medians <- plot5.df %>%
  group_by(climate_region) %>%
  summarize(median_value = median(yi, na.rm = TRUE)) %>%
  arrange(desc(median_value))

plot5.df <- plot5.df %>%
  mutate(climate_region = factor(climate_region, levels = medians$climate_region))

p1 <- ggplot(plot5.df, aes(x=climate_region, y=sinh(yi), color=plant_species, fill=plant_species, group=climate_region))+
  sina_theme+
  scale_fill_manual(values = c("C. malaccensis"="#025766", "S. patens" = "#5c3972", 
                               "S. alterniflora" = "#2FC45D", "P. australis" ="#20bcff",
                               "D. spicata" = "#47724c", "S. salsa" = "#c6b822",
                               "C. jamaicense" =  "#564e95", "Juncus sp." = "#86a291", "Juncae sp" = "lightblue", 
                               "S. portulacastrum" = "#90E669", "P. maritima" = "#30BCAD",
                               "S. mariqueter" = "#60563f", "T. chinensis"="#FF46A2",
                               "Salicornia sp."="black", "Suaeda sp." = "black",
                               "Carex. sp."="#EEAAEE", "S. americanus"="#FF9E4A" ,
                               "S. angelica"="#D62728"))+
  scale_color_manual(values = c("C. malaccensis"="#025766", "S. patens" = "#5c3972", 
                                "S. alterniflora" = "#2FC45D", "P. australis" ="#20bcff",
                                "D. spicata" = "#47724c", "S. salsa" = "#c6b822",
                                "C. jamaicense" =  "#564e95", "Juncus sp." = "#86a291", "Juncae sp" = "lightblue", 
                                "S. portulacastrum" = "#90E669", "P. maritima" = "#30BCAD",
                                "S. mariqueter" = "#60563f", "T. chinensis"="#FF46A2",
                                "Salicornia sp."="black", "Suaeda sp." = "black",
                                "Carex. sp."="#EEAAEE", "S. americanus"="#FF9E4A" ,
                                "S. angelica"="#D62728"))+
  labs(x = "Climate Region")

p1

# combine plots -----------------------------------------------------------

p1 <- p1 + labs(tag = "b.")
p2 <- p2 + labs(tag = "c.")
p3 <- p3 + labs(tag = "d.")
p4 <- p4 + labs(tag = "e.")
p5 <- p5 + labs(tag = "f.")

combined_plot <- 
  (p1 + p2 + p3 + plot_layout(widths = c(2,2,2))) /
  (p4 + p5 + plot_layout(widths = c(1, 2))) /
  plot_layout(guides = "collect") +
  plot_annotation() & 
  theme(plot.tag = element_text( size = 12, family="Arial Black"))

combined_plot
ggsave("resources/enviro_variables.png", combined_plot, width = 160, height = 120, dpi=300, units = "mm")

# Load images
img1 <- image_read("resources/ghg_meta_coords_fullmap.png")
img2 <- image_read("resources/enviro_variables.png")

# Get width of img2
img2_info <- image_info(img2)
target_width <- img2_info$width*.95

# Resize img1 to match width, keeping aspect ratio
img1_resized <- image_scale(img1, paste0(target_width))

img1_annotated <- image_annotate(
  img1_resized,
  text = "a.",
  font = "Arial Black",  # Use the bold variant directly
  size = 50,
  gravity = "northwest",
  location = "+10+10",
  color = "black",
  boxcolor = "none"
)
print(img1_annotated)

# Combine top-to-bottom
combined_v <- image_append(c(img1_annotated, img2), stack = TRUE)

image_write(combined_v, "figures/figure1.png")

