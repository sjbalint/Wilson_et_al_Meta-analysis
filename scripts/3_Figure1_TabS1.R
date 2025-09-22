# -------------------------------------------------------------------------
#Created 1.13.25 by Emily Wilson, last update 09.22.25 by Emily Wilson
#Title: Plant Species Drive Global Coastal Wetland Methane Fluxes Meta-analysis
#Figure 1
rm(list=ls()) #clear the environment
#dev.off()



# install packages --------------------------------------------------------
library(tidyverse) #for data manipulation
library(scales) #for y-axis
library(multcompView) #for multiple comparison tests and displaying letters
library(FSA) #for non-parametric pairwise comparisons
library(cowplot) #for combining subplots
library(patchwork) #for combining subplots
library(ggforce) #for sina
library(ggsci) #for color scales



# Get the GHG flux data ---------------------------------------------------
#meta_data<-readRDS("full_data.RDS")
meta_data<-read.csv("raw/ch4_full_dataset.csv") #if using csv file

#figures and kruskall-wallis are on untransformed data
#sinh(yi) puts data back on the orignal scale not asinh scale


# Rename plants for graph -------------------------------------------------
meta_data <- meta_data %>%
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
meta_data <- meta_data %>%
  group_by(plant_species) %>%
  filter(n_distinct(p_num) > 1) %>%
  ungroup()



#Run Kruskal-Wallis and Dunn --------------------------------------------
kruskal <- kruskal.test(sinh(yi) ~ plant_species, data = meta_data)
print(kruskal)
dunn <- dunnTest(sinh(yi) ~ plant_species, data = meta_data)
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
#change default aesthetics for point - used for geom_sina
update_geom_defaults("point", list(shape = 21,  size=1, stroke=0.3))

#define the dot and bars
median_IQR <- function(x) {
  data.frame(y = median(x), # Median
             ymin = quantile(x)[2], # 1st quartile
             ymax = quantile(x)[4])  # 3rd quartile
}

medians <- meta_data %>%
  group_by(plant_species) %>%
  summarize(median_value = median(yi, na.rm = TRUE)) %>%
  arrange(desc(median_value)) #order plants by median

meta_data <- meta_data %>%
  mutate(plant_species = factor(plant_species, levels = medians$plant_species))

asinh_trans <- trans_new( #function for plotting, made by Sawyer Balint
  name = "asinh",
  transform = asinh,
  inverse = sinh
)

p5 <- ggplot(meta_data, aes(x=plant_species, y=sinh(yi), color=plant_species, fill=plant_species, group=plant_species))+
  theme_classic()+
  geom_sina(
    aes(fill = plant_species, color = plant_species),  # group aesthetics
    shape = 21,                        # filled circle
    size = 1.2,
    stroke = .3,
    alpha = 0.3)+
  scale_y_continuous(
    trans = asinh_trans,
    limits=c(-2000,10000),
    breaks = c(-1000,-100, -10, 0, 10, 100, 1000,10000))+
  scale_fill_manual(values = c("C. malaccensis"="#025766", "S. patens" = "#5c3972", 
                               "S. alterniflora" = "#2FC45D", "P. australis" ="#20bcff",
                               "D. spicata" = "#47724c", "S. salsa" = "#c6b822",
                               "C. jamaicense" =  "#564e95", "Juncus sp." = "#86a291", "Juncae sp" = "lightblue", 
                               "S. portulacastrum" = "#90E669", "P. maritima" = "#30BCAD",
                               "S. mariqueter" = "#60563f", "T. chinensis"="#FF46A2",
                               "Salicornia sp."="black", "Suaeda sp." = "black",
                               "Carex. sp."="#EEAAEE", "S. americanus"="#FF9E4A" ,
                               "S. angelica"="#D62728")) +
  scale_color_manual(values = c("C. malaccensis"="#025766", "S. patens" = "#5c3972", 
                               "S. alterniflora" = "#2FC45D", "P. australis" ="#20bcff",
                               "D. spicata" = "#47724c", "S. salsa" = "#c6b822",
                               "C. jamaicense" =  "#564e95", "Juncus sp." = "#86a291", "Juncae sp" = "lightblue", 
                               "S. portulacastrum" = "#90E669", "P. maritima" = "#30BCAD",
                               "S. mariqueter" = "#60563f", "T. chinensis"="#FF46A2",
                               "Salicornia sp."="black", "Suaeda sp." = "black",
                               "Carex. sp."="#EEAAEE", "S. americanus"="#FF9E4A" ,
                               "S. angelica"="#D62728")) +
  theme(legend.position="none")+
  
  labs(y = expression(paste("CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), 
       x = "Plant Species") +
  stat_summary(geom = "linerange",
               fun.data = median_IQR, 
               position = position_dodge(width=0.2), 
               linewidth=0.5,
               col = "black") + 
  stat_summary(geom = "point", 
               fun = "median", 
               position = position_dodge(width=0.2),
               shape=19,
               size = 2, 
               col = "black")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7, colour = "black", face="italic", angle = 30, hjust = .90, family = "Arial"), 
        axis.text.y = element_text(size = 7, colour = "black", family = "Arial"),
        axis.title.x = element_text(size = 7, color = "black", family = "Arial", face="bold"),
        axis.title.y = element_text(size = 7, color = "black", family = "Arial")) +
  geom_hline(yintercept = 0, linetype = "dashed")
  
p5



# Salinity ----------------------------------------------------------------
#get the GHG flux data 
meta_data<-read.csv("raw/ch4_full_dataset.csv") #if using csv file, reload bc of changes
#meta_data<-readRDS("full_data.RDS")

#change default aesthetics for point - used for geom_sina
update_geom_defaults("point", list(shape = 21,  size=1, stroke=0.3))

#define the dot and bars
median_IQR <- function(x) {
  data.frame(y = median(x), # Median
             ymin = quantile(x)[2], # 1st quartile
             ymax = quantile(x)[4])  # 3rd quartile
}

asinh_trans <- trans_new( #function for plotting
  name = "asinh",
  transform = asinh,
  inverse = sinh
)

medians <- meta_data %>%
  group_by(salinity) %>%
  summarize(median_value = median(yi, na.rm = TRUE)) %>%
  arrange(desc(median_value))

meta_data <- meta_data %>%
  mutate(salinity = factor(salinity, levels = medians$salinity))

#Run Kruskal-Wallis and Dunn --------------------------------------------
kruskal <- kruskal.test(sinh(yi) ~ salinity, data = meta_data)
print(kruskal)
dunn <- dunnTest(sinh(yi) ~ salinity, data = meta_data)
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

medians <- meta_data %>%
  group_by(salinity) %>%
  summarize(median_value = median(yi, na.rm = TRUE)) %>%
  arrange(desc(median_value))

meta_data <- meta_data %>%
  mutate(salinity = factor(salinity, levels = medians$salinity))

p4 <- ggplot(meta_data, aes(x=salinity, y=sinh(yi), color=salinity, fill=salinity, group=salinity))+
  #basetheme+
  theme_classic()+
  geom_sina(
    aes(fill = salinity, color = salinity),  # group aesthetics
    shape = 21,                        # filled circle
    size = 1.2,
    stroke = .3,
    alpha = 0.3)+
  scale_y_continuous(
    trans = asinh_trans,
    limits=c(-2000,13000),
    breaks = c(-1000,-100, -10, 0, 10, 100, 1000,10000))+
  scale_fill_manual(values = c("oligohaline" = "lightblue", 
                                 "mesohaline" = "#1F78B4",
                                "polyhaline" = "darkblue",
                               "euhaline" = "darkgreen")) +
  scale_color_manual(values = c("oligohaline" = "lightblue", 
                               "mesohaline" = "#1F78B4",
                               "polyhaline" = "darkblue",
                               "euhaline" = "darkgreen")) +
  labs(y = expression(paste("CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), 
       x = "Salinity") +
  stat_summary(geom = "linerange",
               fun.data = median_IQR, 
               position = position_dodge(width=0.2), 
               linewidth=0.5,
               col = "black") + 
  stat_summary(geom = "point", 
               fun = "median", 
               position = position_dodge(width=0.2),
               shape=19,
               size = 2, 
               col = "black")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7, colour = "black", angle = 30, hjust = .90, family = "Arial"), 
        axis.text.y = element_text(size = 7, colour = "black", family = "Arial"),
        axis.title.x = element_text(size = 7, color = "black", family = "Arial", face="bold"),
        axis.title.y = element_text(size = 7, color = "black", family = "Arial")) +
  geom_hline(yintercept = 0, linetype = "dashed")

p4



# Season ------------------------------------------------------------------
# Clean filtering and overwriting
meta_data<-meta_data%>%
  mutate(season=as.factor(season))
meta_data_season <- meta_data %>%
  filter(!is.na(season), season != "annual") %>%
  mutate(season = droplevels(season))

medians <- meta_data_season %>%
  group_by(season) %>%
  summarize(median_value = median(yi, na.rm = TRUE)) %>%
  arrange(desc(median_value))

meta_data_season <- meta_data_season %>%
  mutate(season = factor(season, levels = medians$season))



# Kruskal Wallis and Dunn ----------------------------------------------------------------
kruskal <- kruskal.test(sinh(yi) ~ season, data = meta_data_season)
print(kruskal)
dunn <- dunnTest(sinh(yi) ~ season, data = meta_data_season)
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

medians <- meta_data_season %>%
  group_by(season) %>%
  summarize(median_value = median(yi, na.rm = TRUE)) %>%
  arrange(desc(median_value))


meta_data_season <- meta_data_season %>%
  mutate(season = factor(season, levels = medians$season))


p2 <- ggplot(meta_data_season, aes(x=season, y=sinh(yi), color=season, fill=season, group=season))+
  #basetheme+
  theme_classic()+
  geom_sina(
    aes(fill = season, color = season),  # group aesthetics
    shape = 21,                        # filled circle
    size = 1.2,
    stroke = .3,
    alpha = 0.3)+
  scale_y_continuous(
    trans = asinh_trans,
    limits=c(-2000,14000),
    breaks = c(-1000,-100, -10, 0, 10, 100, 1000,10000))+
  scale_fill_manual(values = c("fall" = "#D95F02", 
                               "winter" = "#1F78B4", 
                               "spring" = "#33A02C", 
                               "summer" = "darkgreen"))+
  scale_color_manual(values = c("fall" = "#D95F02", 
                              "winter" = "#1F78B4", 
                              "spring" = "#33A02C", 
                              "summer" = "darkgreen"))+
  labs(y = expression(paste("CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), 
       x = "Season") +
  stat_summary(geom = "linerange",
               fun.data = median_IQR, 
               position = position_dodge(width=0.2), 
               linewidth=0.5,
               col = "black") + 
  stat_summary(geom = "point", 
               fun = "median", 
               position = position_dodge(width=0.2),
               shape=19,
               size = 2, 
               col = "black")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7, colour = "black", angle = 30, hjust = .90, family = "Arial"), 
        axis.text.y = element_text(size = 7, colour = "black", family = "Arial"),
        axis.title.x = element_text(size = 7, color = "black", family = "Arial", face="bold"),
        axis.title.y = element_text(size = 7, color = "black", family = "Arial")) +
   geom_hline(yintercept = 0, linetype = "dashed")

p2



# Tide ---------------------------------------------------
#meta_data<-readRDS("full_data.RDS")
meta_data<-read.csv("raw/ch4_full_dataset.csv") #if using csv file

meta_data <- meta_data %>%
  mutate(tide = factor(tide, 
                           levels = c("macrotidal", "mesotidal", "microtidal")))
meta_data$tide <- factor(meta_data$tide)

# Kruskal Wallis ----------------------------------------------------------------
kruskal <- kruskal.test(sinh(yi) ~ tide, data = meta_data)
print(kruskal)
dunn <- dunnTest(sinh(yi) ~ tide, data = meta_data)
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

medians <- meta_data %>%
  group_by(tide) %>%
  summarize(median_value = median(yi, na.rm = TRUE)) %>%
  arrange(desc(median_value))

meta_data <- meta_data %>%
  mutate(tide = factor(tide, levels = medians$tide))

p3 <- ggplot(meta_data, aes(x=tide, y=sinh(yi), color=tide, fill=tide, group=tide))+
  #basetheme+
  theme_classic()+
  geom_sina(
    aes(fill = tide, color = tide),  # group aesthetics
    shape = 21,                        # filled circle
    size = 1.2,
    stroke = .3,
    alpha = 0.3)+
  scale_y_continuous(
    trans = asinh_trans,
    limits=c(-2000,10000),
    breaks = c(-1000,-100, -10, 0, 10, 100, 1000,10000))+
  scale_fill_manual(values = c(
    "microtidal" = "#B57FA4", 
   "mesotidal" = "#998eff",
  "macrotidal"="#2C5985"))+
  scale_color_manual(values = c(
    "microtidal" = "#B57FA4", 
    "mesotidal" = "#998eff",
    "macrotidal"="#2C5985"))+
  labs(y = expression(paste("CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), 
       x = "Tidal Range") +
  stat_summary(geom = "linerange",
               fun.data = median_IQR, 
               position = position_dodge(width=0.2), 
               linewidth=0.5,
               col = "black") + 
  stat_summary(geom = "point", 
               fun = "median", 
               position = position_dodge(width=0.2),
               shape=19,
               size = 2, 
               col = "black")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7, colour = "black", angle = 30, hjust = .90, family = "Arial"), 
        axis.text.y = element_text(size = 7, colour = "black", family = "Arial"),
        axis.title.x = element_text(size = 7, color = "black", family = "Arial", face= "bold"),
        axis.title.y = element_text(size = 7, color = "black", family = "Arial")) +
  geom_hline(yintercept = 0, linetype = "dashed")

p3



# Climate Region ----------------------------------------------------------
#meta_data<-readRDS("full_data.RDS")
meta_data<-read.csv("raw/ch4_full_dataset.csv") #if using csv file
meta_data <- meta_data %>%
  mutate(climate_region = recode(climate_region,
                                "subtropical/warm temperate"="warm temperate"))
# Kruskal Wallis ----------------------------------------------------------------
kruskal <- kruskal.test(sinh(yi) ~ climate_region, data = meta_data)
print(kruskal)

# Run Dunn's test
dunn <- dunnTest(sinh(yi) ~ climate_region, data = meta_data_season)
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

medians <- meta_data %>%
  group_by(climate_region) %>%
  summarize(median_value = median(yi, na.rm = TRUE)) %>%
  arrange(desc(median_value))

meta_data <- meta_data %>%
  mutate(climate_region = factor(climate_region, levels = medians$climate_region))

p1 <- ggplot(meta_data, aes(x=climate_region, y=sinh(yi), color=climate_region, fill=climate_region, group=climate_region))+
  #basetheme+
  theme_classic()+
  geom_sina(
    aes(fill = climate_region, color = climate_region),  # group aesthetics
    shape = 21,                        # filled circle
    size = 1.2,
    stroke = .3,
    alpha = 0.3)+
  scale_y_continuous(
    trans = asinh_trans,
    limits=c(-2000,10000),
    breaks = c(-1000,-100, -10, 0, 10, 100, 1000,10000))+
  scale_fill_manual(values = c(
    "continental" = "#4EA6DC", 
    "warm temperate"="#E32D91",
    "mediterranean"="#F14124",
    "oceanic"="#6997AF",
    "polar"="#014D64"
   ))+
  scale_color_manual(values = c(
    "continental" = "#4EA6DC", 
    "warm temperate"="#E32D91",
    "mediterranean"="#F14124",
    "oceanic"="#6997AF",
    "polar"="#014D64"
  ))+
 
  labs(y = expression(paste("CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), 
       x = "Climate Region") +
  stat_summary(geom = "linerange",
               fun.data = median_IQR, 
               position = position_dodge(width=0.2), 
               linewidth=0.5,
               col = "black") + 
  stat_summary(geom = "point", 
               fun = "median", 
               position = position_dodge(width=0.2),
               shape=19,
               size = 2, 
               col = "black")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7, colour = "black", angle = 30, hjust = .90, family = "Arial"), 
        axis.text.y = element_text(size = 7, colour = "black", family = "Arial"),
        axis.title.x = element_text(size = 7, color = "black", family = "Arial", face="bold"),
        axis.title.y = element_text(size = 7, color = "black", family = "Arial")) +
  geom_hline(yintercept = 0, linetype = "dashed")

p1

library(patchwork)
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
ggsave("figures/enviro_variables.png", combined_plot, width = 160, height = 120, dpi=300, units = "mm")


