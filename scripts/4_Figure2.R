# -------------------------------------------------------------------------
#Created 1.13.25 by Emily Wilson
#last update 09.25.25 by Sawyer Balint
#Title: Plant Species Drive Global Coastal Wetland Methane Fluxes Meta-analysis
#Figure 2

rm(list=ls()) #clear the environment

# install packages --------------------------------------------------------
library(scales) #for stats
library(tidyverse) #for data manipulation
library(patchwork) #to put plots together
library(grid)
library(ggtext)
library(pheatmap)
library(viridis)

source("functions/utilities.r")

# import data ------------------------------------------------------------

df <- read.csv("raw/ch4_soilsalinity_dataset.csv") %>%
  mutate(
    plant_species=as.factor(plant_species),
    salinity=as.factor(salinity),
    tide=as.factor(tide),
    ch4_method_simple=as.factor(ch4_method_simple),
    season=as.factor(season),
    lat=abs(as.numeric(lat)))

full.df <- read.csv("raw/ch4_full_dataset.csv")

rf.model <- readRDS("Rdata/RandomForest.rds")
gam1.model <- readRDS("Rdata/GAM_SoilSalinity.rds")
gam2.model <- readRDS("Rdata/GAM_SalinityCategory.rds")

# configure graphing ------------------------------------------------------

theme <- list(
  theme_classic(),
  theme(
    axis.text.x = element_text(size = 5, colour = "black", family = "Arial"), 
    axis.text.y = element_text(size = 5, colour = "black", family = "Arial"),
    axis.title.x = element_text(size = 7, color = "black", family = "Arial"),
    axis.title.y = element_text(size = 7, color = "black", family = "Arial"))
)

# plot salinity and yi ----------------------------------------------------

meta_data <- df%>% #clean up the data
  filter(!is.na(salinity_conductivity)) %>%  #only include fluxes with soil salinity
  filter(plant_species!="Juncae sp")

# Plot salinity and yi -------------------------------------------------------

#Model Comparison -------------------------------------------------------

summary(meta_data$season)

#subset data to unique combinations of salinity and plant_species
#the model uses latitude and season, so assign values for each
new_data <- meta_data %>%
  group_by(salinity, plant_species) %>%
  summarise(
    yi = median(yi, na.rm = TRUE),
    n = n(),
    lat=median(lat)
  ) %>%
  ungroup() %>%
  mutate(season="fall")

new_data$predicted_plant_sal <- predict(gam2.model, newdata = new_data, type = "response", se.fit = FALSE)

plot_data <- data.frame(
  salinity = new_data$salinity,
  plant_species = new_data$plant_species,
  actual = new_data$yi,
  predicted_plant_sal = new_data$predicted_plant_sal) %>%
  mutate(plant_species = recode_factor(plant_species,
                                       "Cladium jamaicense" = "C. jamaicense",
                                       "Juncus sp." = "Juncus sp.",
                                       "Plantago maritima" = "Pl. maritima",
                                       "Spartina patens" = "S. patens",
                                       "Spartina alterniflora" = "S. alterniflora",
                                       "Distichlis Spicata" = "D. spicata",
                                       "Phragmites australis" = "P. australis",
                                       "Cyperus malaccensis" = "C. malaccensis",
                                       "Scirpus mariqueter" = "S. mariqueter",
                                       "Suaeda salsa" = "S. salsa",
                                       "Tamarix chinensis" = "T. chinensis"))

plant_sal_r2 <- cor(plot_data$actual, plot_data$predicted_plant_sal)^2

p6 <- ggplot(plot_data, aes(x = sinh(actual), y = sinh(predicted_plant_sal),
                          fill = plant_species, shape=salinity)) +
  theme+
  geom_point(size = 2, alpha = 0.8, stroke = 0.3) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.4, linetype = "dashed", color = "black")+
  annotate("text", x = -10, y = -15, label = "1:1 Line", hjust = 0, size = 2, family = "Arial") +
  scale_x_continuous(transform = asinh_trans,
                     limits=c(-100,8000),
                     breaks=c(-100,-10, 0, 10, 100, 1000))+
  scale_y_continuous(
    trans = asinh_trans,
    limits=c(-100,8000),
    breaks = c(-100, -10, 0, 10, 100, 1000))+
  annotate("text", x = 1000, y = 5000,  # Adjust x and y based on your plot range
           label =paste0("Plant Species by Salinity Category Model \n Predicted Fluxes ~ Actual Median Fluxes, R² = ", round(plant_sal_r2, 2)),
           size = 2, hjust = 1, color = "black", family = "Arial")+
  scale_shape_manual(values = c("oligohaline"=21, "mesohaline" = 24, "polyhaline"=22,
                                "euhaline" = 23))+
  scale_fill_manual(
    values = c(
      "C. malaccensis"  = "#025766",
      "S. patens"       = "#5c3972",
      "S. alterniflora" = "#2FC45D",
      "P. australis"    = "#20bcff",
      "D. spicata"      = "#47724c",
      "S. salsa"        = "#c6b822",
      "C. jamaicense"   = "#564e95",
      "Juncus sp."      = "#86a291",
      "Pl. maritima"    = "#30BCAD",
      "S. mariqueter"   = "#60563f",
      "T. chinensis"    = "#8B5F65",
      "Salicornia sp."  = "#CD8162",
      "Suaeda sp."      = "black"),
    labels = c(
      "<i>C. malaccensis</i>", "<i>S. patens</i>", "<i>S. alterniflora</i>", "<i>P. australis</i>",
      "<i>D. spicata</i>", "<i>S. salsa</i>", "<i>C. jamaicense</i>", "<i>Juncus sp.</i>",
      "<i>Pl. maritima</i>", "<i>S. mariqueter</i>", "<i>T. chinensis</i>",
      "<i>Salicornia sp.</i>", "<i>Suaeda sp.</i>"
    ))+
  guides(
    fill = guide_legend(ncol = 5, override.aes = list(shape = 21, face = "italics")),
    shape = guide_legend(ncol = 2, override.aes = list(fill = "grey80")) # keep shape legend clear
  )+
  labs(fill="Plant Species", shape = "Salinity Category",
       x = expression(paste("Actual Median CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), 
       y = expression(paste("Predicted CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), 
       family="Arial")+
  theme(legend.title = element_blank(),
          legend.key.height = unit(0.1, "lines"),  # height of each key; smaller = tighter rows
          legend.key.width = unit(0.1, "lines"),   # width of each key; smaller = more items per row
          legend.spacing.y = unit(0.1, "lines"),   # vertical spacing between rows
          legend.spacing.x = unit(0.1, "lines"),   # horizontal spacing between items
        legend.position="bottom",
        legend.text = element_markdown(size = 5, family = "Arial"))

p6


# model comparison --------------------------------------------------------

meta_data <- df

meta_data<-meta_data %>%
  filter(!is.na(salinity_conductivity),
         season!="winter")

#subset your data to unique combinations of salinity and plant_species
new_data <- meta_data %>%
  mutate(salinity_conductivity = round(salinity_conductivity, 0)) %>%
  mutate(lat = abs(as.numeric(round(lat, 1)))) %>%
  group_by(plant_species, salinity_conductivity, lat, season, salinity) %>%
  summarise(
    yi = mean(yi, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

new_data$predicted_plant_lat <- predict(gam1.model, newdata = new_data, type = "response", se.fit = FALSE)

sal_pw<-new_data$salinity_conductivity

# Define the conversion formula as a function
predict_poff <- function(sal) {
  (
    (10^((-0.056 * sal) + 1.38)) #untransform
    /0.34) * # divide by 0.34 puts it in mg/m2/d
    (1 / 16.04) * (1000 / 1) * (1/24) #convert to mmols, then umol, then h
}

# Apply function to each salinity value
predicted_sal_pw <- predict_poff(sal_pw)

plot_data <- data.frame(
  salinity = new_data$salinity,
  salinity_conductivity = new_data$salinity_conductivity,
  plant_species = new_data$plant_species,
  actual = new_data$yi,
  poffenbarger = asinh(predicted_sal_pw),
  predicted_plant_lat = new_data$predicted_plant_lat) %>%
  mutate(plant_species = recode_factor(plant_species,
                                       "Cladium jamaicense" = "C. jamaicense",
                                       "Juncus sp." = "Juncus sp.",
                                       "Plantago maritima" = "Pl. maritima",
                                       "Spartina patens" = "S. patens",
                                       "Spartina alterniflora" = "S. alterniflora",
                                       "Distichlis Spicata" = "D. spicata",
                                       "Phragmites australis" = "P. australis",
                                       "Cyperus malaccensis" = "C. malaccensis",
                                       "Scirpus mariqueter" = "S. mariqueter",
                                       "Suaeda salsa" = "S. salsa",
                                       "Tamarix chinensis" = "T. chinensis"))

#get R2 values
poffenbarger_r2 <- cor(plot_data$actual, plot_data$poffenbarger)^2
plant_lat_r2 <- cor(plot_data$actual, plot_data$predicted_plant_lat)^2


# Plot
p2<-ggplot(plot_data, aes(x = sinh(actual), y = sinh(poffenbarger),
                          fill = plant_species, shape=salinity)) +
  theme+
  geom_point(size = 2, alpha = 0.8, stroke = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")+
  annotate("text", x = -10, y = -15, label = "1:1 Line", hjust = 0, size = 2, family = "Arial") +
  scale_x_continuous(transform = asinh_trans,
                     limits=c(-100,10000),
                     breaks=c(-100,-100,-10, 0, 10,100, 1000))+
  scale_y_continuous(
    trans = asinh_trans,
    limits=c(-100,10000),
    breaks = c(-100, -100, -10, 0, 10, 100, 1000)
  )+
  annotate("text", x = 1300, y = 5000,
           label =paste0("Salinity-Based Proxy Predicted \n Fluxes ~ Actual Fluxes, R² = ", round(poffenbarger_r2, 2)),
           size = 2, hjust = 1, color = "darkred", family = "Arial")+
  scale_shape_manual(values = c("oligohaline"=21, "mesohaline" = 24, "polyhaline"=22,
                                "euhaline" = 23))+
  scale_fill_manual(
    values = c(
      "C. malaccensis"  = "#025766",
      "S. patens"       = "#5c3972",
      "S. alterniflora" = "#2FC45D",
      "P. australis"    = "#20bcff",
      "D. spicata"      = "#47724c",
      "S. salsa"        = "#c6b822",
      "C. jamaicense"   = "#564e95",
      "Juncus sp."      = "#86a291",
      "Pl. maritima"    = "#30BCAD",
      "S. mariqueter"   = "#60563f",
      "T. chinensis"    = "#8B5F65",
      "Salicornia sp."  = "#CD8162",
      "Suaeda sp."      = "black"))+
  guides(fill = "none", shape = "none", color="none")+
  labs(
    x = expression(paste("Actual CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), 
    y = expression(paste("Predicted CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), 
    family="Arial")

p2

p8 <-ggplot(plot_data, aes(x = sinh(actual), y = sinh(predicted_plant_lat),
                           fill = plant_species, shape=salinity)) +
  theme+
  geom_point(size = 2, alpha = 0.8, stroke = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")+
  annotate("text", x = -10, y = -15, label = "1:1 Line", hjust = 0, size = 2, family = "Arial") +
  scale_x_continuous(transform = asinh_trans,
                     limits=c(-100,10000),
                     breaks=c(-100,-10, 0, 10, 100, 1000))+
  scale_y_continuous(
    trans = asinh_trans,
    limits=c(-100,10000),
    breaks = c(-100, -10, 0, 10, 100, 1000))+
  annotate("text", x = 1000, y = 5000,  # Adjust x and y based on your plot range
           label =paste0("This Study's Full Model Predicted \n Fluxes ~ Actual Fluxes, R² = ", round(plant_lat_r2, 2)),
           size = 2, hjust = 1, color = "black", family = "Arial")+
  scale_shape_manual(values = c("oligohaline"=21, "mesohaline" = 24, "polyhaline"=22,
                                "euhaline" = 23))+
  scale_fill_manual(
    values = c(
      "C. malaccensis"  = "#025766",
      "S. patens"       = "#5c3972",
      "S. alterniflora" = "#2FC45D",
      "P. australis"    = "#20bcff",
      "D. spicata"      = "#47724c",
      "S. salsa"        = "#c6b822",
      "C. jamaicense"   = "#564e95",
      "Juncus sp."      = "#86a291",
      "Pl. maritima"    = "#30BCAD",
      "S. mariqueter"   = "#60563f",
      "T. chinensis"    = "#8B5F65",
      "Salicornia sp."  = "#CD8162",
      "Suaeda sp."      = "black")) +
  guides(fill = "none", shape = "none", color="none")+
  labs(x = expression(paste("Actual CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), 
       y = expression(paste("Predicted CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), 
       family="Arial")

# random forest -----------------------------------------------------------

var_imp <- data.frame(
  Variable = names(rf.model$importance),
  Importance = rf.model$importance)

var_imp <- var_imp %>%
  mutate(Variable = case_when(
    Variable == "tide" ~ "Tidal Range",
    Variable ==  "climate_region" ~ "Climate Region",
    Variable ==  "season" ~ "Season",
    Variable ==  "salinity" ~ "Salinity Category",
    Variable == "ch4_method_simple" ~ "Sampling Method",
    Variable == "plant_species" ~ "Plant Species",
    Variable == "lat" ~  "Absolute Latitude",
    Variable == "salinity_conductivity"~"Soil salinity",
    TRUE ~   Variable
  ))
# Arrange from highest to lowest importance
var_imp <- var_imp %>%
  arrange(desc(Importance)) %>%
  mutate(Variable=str_replace(Variable, " ","\n"))

# Plot using ggplot2
RF <- ggplot(var_imp, aes(x = reorder(Variable, Importance), y = Importance)) +
  theme+
  geom_col(fill = "darkgreen", color="black", linewidth=0.3) +
  coord_flip() +
  labs(
    x = "Variable",
    y = "Importance")+
  theme(axis.text.y = element_text(size = 6, colour = "black", family = "Arial"))

# heat map ----------------------------------------------------------------

plot.df <- full.df %>%
  group_by(plant_species) %>%
  filter(n_distinct(p_num) > 1) %>%
  ungroup() %>%
  mutate(plant_species = recode_factor(plant_species,
                                "Aeluropus littoralis" = "A. littoralis",
                                "Carex sp." = "Carex sp.",
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
  )) %>%
  mutate(salinity_category = factor(salinity,
                                    levels = rev(c("oligohaline", "mesohaline", "polyhaline", "euhaline"))))

#average yi for each species and salinity category
plot.df <- plot.df %>%
  group_by(plant_species, salinity_category) %>%
  summarise(median_yi = median(yi, na.rm = TRUE), .groups = "drop")

species_order <- plot.df %>%
  group_by(plant_species) %>%
  summarize(yi = mean(median_yi, na.rm=TRUE)) %>%
  arrange(yi) %>%
  pull(plant_species)

plot.df <- plot.df %>%
  mutate(plant_species=factor(plant_species, levels=rev(species_order)))%>%
  mutate(ch4_flux = sinh(median_yi)) %>%
  complete(plant_species, salinity_category)

map <- ggplot(plot.df, aes(plant_species, salinity_category, fill=ch4_flux))+
  theme_bw()+
  geom_tile(color="black")+
  scale_fill_viridis(na.value="white", 
                     trans = asinh_trans,
                     limits=c(-10,NA),
                     breaks = c(-10, 0, 10, 100, 1000, 10000),
                     name=expression(paste("CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")))+
  guides(fill = guide_colourbar(
    barwidth=0.5,
    harheight=10,
    title.position="right",
    theme = theme(
      legend.ticks = element_blank(),
      text = element_text(size=7),
      title = element_text(angle = 90, vjust=0.5, size=7, margin = margin(r = 0, l = 0, t = 0, b = 0)),
    )
  ))+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  labs(x=NULL, y=NULL)+
  theme(axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.spacing = unit(0, "pt"),
        legend.box.spacing = unit(5, "pt"),
        axis.text.x = element_text(size = 7, colour = "black", face="italic", angle = 45, hjust = 1, family = "Arial"),
        axis.text.y = element_text(size = 7, colour = "black", family = "Arial"),
        axis.title.x = element_text(size = 7, color = "black", family = "Arial"),
        axis.title.y = element_text(size = 7, color = "black", family = "Arial")
  )

#wrap the gtable directly for patchwork
map_wrapped <- wrap_elements(map)


# combine plots -----------------------------------------------------------

combined_plot <- (RF + p8 + p2 + plot_layout(widths = c(1.0,1.5,1.5))) / (p6 + map_wrapped + plot_layout(widths = c(1.6,2.4))) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = 'a', 
    tag_suffix = "."  # use lower-case letters
  ) &
  theme(
    legend.text = element_markdown(size = 7), 
    plot.tag = element_text(family = "Arial Black", size = 12),
    legend.position = "bottom"           # place at bottom
    )

combined_plot

ggsave("figures/figure2.png", combined_plot, width = 180, height = 175, dpi=300, units ="mm")

