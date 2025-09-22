# -------------------------------------------------------------------------
#Created 1.13.25 by Emily Wilson, last update 08.20.25 by Emily Wilson
#Title: Plant Species Drive Global Coastal Wetland Methane Fluxes Meta-analysis

#############################################################################
#############################################################################
#############################################################################
#Figure 2
rm(list=ls()) #clear the environment
#dev.off()



# install packages --------------------------------------------------------
library(scales) #for stats
library(tidyverse) #for data manipulation
library(patchwork) #to put plots together
library(mgcv) #for gam
library(grid)

#meta_data<-readRDS("pw_data.RDS")
meta_data<-read.csv("raw/ch4_soilsalinity_dataset.csv")

meta_data <- meta_data %>%
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

str(meta_data)
meta_data <- meta_data %>% #clean up the data
  filter(!is.na(salinity_conductivity)) %>%  #only include fluxes with soil salinity
  mutate(lat=abs(as.numeric(lat)), #make lat absolute for analysis
         plant_species=as.factor(plant_species),
         season=as.factor(season),
         salinity_conductivity=as.numeric(salinity_conductivity),
         ch4_method_simple=as.factor(ch4_method_simple),
         tide=as.factor(tide))%>%
  filter(plant_species!="Juncae sp")



# Plot salinity and yi -------------------------------------------------------
#custom asinh transform for ggplot
asinh_trans <- trans_new(
  name = "asinh",
  transform = asinh,
  inverse = sinh
)

#Model Comparison -------------------------------------------------------

meta_data$plant_species<-as.factor(meta_data$plant_species)

gam_model6 <- gam(yi ~ plant_species+salinity, data = meta_data, weights = vi)
summary(gam_model6) #R2=0.62

#subset data to unique combinations of salinity and plant_species
new_data <- meta_data %>%
  group_by(salinity, plant_species) %>%
  summarise(
    yi = median(yi, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
new_data$predicted_plant_sal <- predict(gam_model6, newdata = new_data, type = "response", se.fit = FALSE)

plot_data <- data.frame(
  salinity = new_data$salinity,
  plant_species = new_data$plant_species,
  actual = new_data$yi,
  predicted_plant_sal = new_data$predicted_plant_sal)

plant_sal_r2 <- cor(plot_data$actual, plot_data$predicted_plant_sal)^2

library(ggtext)

p6<-ggplot(plot_data, aes(x = sinh(actual), y = sinh(predicted_plant_sal),
                          fill = plant_species, shape=salinity)) +
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
  annotate("text", x = 900, y = 5000,  # Adjust x and y based on your plot range
           label =paste0("Plant Species by Salinity Category Model \n Predicted Fluxes ~ Actual Median Fluxes, R² = ", round(plant_sal_r2, 2)),
           size = 2, hjust = 1, color = "black", family = "Arial")+
  scale_shape_manual(values = c("oligohaline"=21, "mesohaline" = 24, "polyhaline"=22,
                                "euhaline" = 23))+
  scale_fill_manual(values = c("C. malaccensis"="#025766", "S. patens" = "#5c3972", 
                               "S. alterniflora" = "#2FC45D", "P. australis" ="#20bcff",
                               "D. spicata" = "#47724c", "S. salsa" = "#c6b822",
                               "C. jamaicense" =  "#564e95", "Juncus sp." = "#86a291", "Juncae sp" = "lightblue", 
                               "S. portulacastrum" = "#90E669", "Pl. maritima" = "#30BCAD",
                               "S. mariqueter" = "#60563f", "T. chinensis"="#FF46A2",
                               "Salicornia sp."="black", "Suaeda sp." = "black",
                               "Carex. sp."="#EEAAEE", "S. americanus"="#FF9E4A" ,
                               "S. angelica"="#D62728"),
                    labels = c("<i>C. malaccensis</i>", "<i>S. patens</i>", "<i>S. alterniflora</i>", "<i>P. australis</i>",
                               "<i>D. spicata</i>", "<i>S. salsa</i>", "<i>C. jamaicense</i>", "<i>Juncus sp.</i>", "<i>Juncae sp</i>",
                               "<i>S. portulacastrum</i>", "<i>Pl. maritima</i>", "<i>S. mariqueter</i>", "<i>T. chinensis</i>",
                               "<i>Salicornia sp.</i>", "<i>Suaeda sp.</i>", "<i>Carex. sp.</i>", "<i>S. americanus</i>", "<i>S. angelica</i>")) +
  guides(
    fill = guide_legend(ncol = 5, override.aes = list(shape = 21, face = "italics")),   # force filled legend
    shape = guide_legend(ncol = 2, override.aes = list(fill = "grey80")) # keep shape legend clear
  )+
  labs(fill="Plant Species", shape = "Salinity Category",
       x = expression(paste("Actual Median CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), 
       y = expression(paste("Predicted CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), 
       family="Arial")+
  theme_classic()+
  theme(legend.title = element_blank(),
          legend.key.height = unit(0.1, "lines"),  # height of each key; smaller = tighter rows
          legend.key.width = unit(0.1, "lines"),   # width of each key; smaller = more items per row
          legend.spacing.y = unit(0.1, "lines"),   # vertical spacing between rows
          legend.spacing.x = unit(0.1, "lines"),   # horizontal spacing between items
        legend.position="bottom",
        legend.text = element_markdown(size = 5, family = "Arial"),
        axis.text.x = element_text(size = 5, colour = "black", family = "Arial"), 
        axis.text.y = element_text(size = 5, colour = "black", family = "Arial"),
        axis.title.x = element_text(size = 7, color = "black", family = "Arial"),
        axis.title.y = element_text(size = 7, color = "black", family = "Arial"))
p6
ggsave("figures/model_plants_lat.png",width=8,height=6)



#Model Comparison
#meta_data<-readRDS("pw_data.RDS")
meta_data<-read.csv("raw/ch4_soilsalinity_dataset.csv")

meta_data <- meta_data %>%
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
meta_data$plant_species<-as.factor(meta_data$plant_species)
meta_data<-meta_data %>%
  filter(!is.na(salinity_conductivity),
         season!="winter")
gam_model4 <- gam(yi ~ plant_species+s(lat)+season+s(salinity_conductivity), data = meta_data, weights = vi)
summary(gam_model4) #R2=0.69


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
new_data$predicted_plant_lat <- predict(gam_model4, newdata = new_data, type = "response", se.fit = FALSE)

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
  predicted_plant_lat = new_data$predicted_plant_lat)

poffenbarger_r2 <- cor(plot_data$actual, plot_data$poffenbarger)^2
plant_lat_r2 <- cor(plot_data$actual, plot_data$predicted_plant_lat)^2


# Plot
p2<-ggplot(plot_data, aes(x = sinh(actual), y = sinh(poffenbarger),
                          fill = plant_species, shape=salinity)) +
  geom_point(size = 2, alpha = 0.8, stroke = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")+
  annotate("text", x = -10, y = -15, label = "1:1 Line", hjust = 0, size = 2, family = "Arial") +
  scale_x_continuous(transform = asinh_trans,
                     limits=c(-1000,8000),
                     breaks=c(-100,-100,-10, 0, 10,100, 1000))+
  scale_y_continuous(
    trans = asinh_trans,
    limits=c(-1000,8000),
    breaks = c(-100, -100, -10, 0, 10, 100, 1000)
  )+
  annotate("text", x = 900, y = 5000,  # Adjust x and y based on your plot range
           label =paste0("Salinity-Based Proxy Equation \n Predicted Fluxes ~ Actual Fluxes, R² = ", round(poffenbarger_r2, 2)),
           size = 2, hjust = 1, color = "darkred", family = "Arial")+
  scale_shape_manual(values = c("oligohaline"=21, "mesohaline" = 24, "polyhaline"=22,
                                "euhaline" = 23))+
  scale_fill_manual(values = c("C. malaccensis"="#025766", "S. patens" = "#5c3972", 
                               "S. alterniflora" = "#2FC45D", "P. australis" ="#20bcff",
                               "D. spicata" = "#47724c", "S. salsa" = "#c6b822",
                               "C. jamaicense" =  "#564e95", "Juncus sp." = "#86a291", "Juncae sp" = "lightblue", 
                               "S. portulacastrum" = "#90E669", "Pl. maritima" = "#30BCAD",
                               "S. mariqueter" = "#60563f", "T. chinensis"="#FF46A2",
                               "Salicornia sp."="black", "Suaeda sp." = "black",
                               "Carex. sp."="#EEAAEE", "S. americanus"="#FF9E4A" ,
                               "S. angelica"="#D62728")) +
  theme(legend.position="none")+
  guides(fill = "none", shape = "none", color="none")+
  labs(
    x = expression(paste("Actual CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), 
    y = expression(paste("Predicted CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), 
    family="Arial")+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 5, colour = "black", family = "Arial"), 
    axis.text.y = element_text(size = 5, colour = "black", family = "Arial"),
    axis.title.x = element_text(size = 7, color = "black", family = "Arial"),
    axis.title.y = element_text(size = 7, color = "black", family = "Arial"))

p2
ggsave("figures/poff_actual.png",width=8,height=6)



p8 <-ggplot(plot_data, aes(x = sinh(actual), y = sinh(predicted_plant_lat),
                           fill = plant_species, shape=salinity)) +
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
  annotate("text", x = 900, y = 5000,  # Adjust x and y based on your plot range
           label =paste0("This Study's Full Model \n Predicted Fluxes ~ Actual Fluxes, R² = ", round(plant_lat_r2, 2)),
           size = 2, hjust = 1, color = "black", family = "Arial")+
  scale_shape_manual(values = c("oligohaline"=21, "mesohaline" = 24, "polyhaline"=22,
                                "euhaline" = 23))+
  scale_fill_manual(values = c("C. malaccensis"="#025766", "S. patens" = "#5c3972", 
                               "S. alterniflora" = "#2FC45D", "P. australis" ="#20bcff",
                               "D. spicata" = "#47724c", "S. salsa" = "#c6b822",
                               "C. jamaicense" =  "#564e95", "Juncus sp." = "#86a291", "Juncae sp" = "lightblue", 
                               "S. portulacastrum" = "#90E669", "Pl. maritima" = "#30BCAD",
                               "S. mariqueter" = "#60563f", "T. chinensis"="#FF46A2",
                               "Salicornia sp."="black", "Suaeda sp." = "black",
                               "Carex. sp."="#EEAAEE", "S. americanus"="#FF9E4A" ,
                               "S. angelica"="#D62728")) +
  theme(legend.position="none")+
  guides(fill = "none", shape = "none", color="none")+
  labs(x = expression(paste("Actual CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), 
       y = expression(paste("Predicted CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), 
       family="Arial")+
  theme_classic()+
  theme(
    axis.text.x.top = element_text(size = 7, family = "Arial"),
    axis.title.x.top = element_text(size = 7, family = "Arial"),
    axis.text.x = element_text(size = 5, colour = "black", family = "Arial"), 
    axis.text.y = element_text(size = 5, colour = "black", family = "Arial"),
    axis.title.x = element_text(size = 7, color = "black", family = "Arial"),
    axis.title.y = element_text(size = 7, color = "black", family = "Arial"))

p8
ggsave("figures/model_plants_lat_season_sal.png",width=8,height=6)


library(patchwork)
library(tidyverse) #for data manipulation
library(randomForestSRC) #parallel random forest
library(scales) #for transformation function
# RANDOM FOREST
# With PW conductivity ----------------------------------------------------
#meta_data <- readRDS("pw_data.RDS") #get data
meta_data<-read.csv("raw/ch4_soilsalinity_dataset.csv")

# Prepare final meta_data -------------------------------------------------------
meta_data <- meta_data %>%
  select(yi, lat, plant_species, season, salinity, salinity_conductivity, ch4_method_simple, tide)

meta_data<-data.frame(meta_data)
meta_data<-meta_data%>%
  mutate(
    plant_species=as.factor(plant_species),
    salinity=as.factor(salinity),
    tide=as.factor(tide),
    ch4_method_simple=as.factor(ch4_method_simple),
    season=as.factor(season),
    lat=abs(as.numeric(lat)))
str(meta_data)
#tune(yi~., data=meta_data, trace=TRUE)
set.seed(39) #42
model <- rfsrc(yi~.,
               ntree=5000, 
               data=meta_data, 
               nodesize = 12, #lowest error, highest r2
               mtry=5, #suggest by tune
               case.wt = meta_data$vi,
               forest=TRUE,
               importance = TRUE)
#looks good
plot(model)
print(model)

var_imp <- data.frame(
  Variable = names(model$importance),
  Importance = model$importance)

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
  arrange(desc(Importance))

# Plot using ggplot2
RF<-ggplot(var_imp, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_col(fill = "darkgreen", color="black", linewidth=0.3) +
  coord_flip() +
  theme_classic()+
  labs(
    x = "Variable",
    y = "Importance")+ 
  theme(
    axis.text.x = element_text(size = 5, colour = "black", angle = 20, hjust = 0.9, family = "Arial"), 
    axis.text.y = element_text(size = 5, colour = "black", family = "Arial"),
    axis.title.x = element_text(size = 7, color = "black", family = "Arial"),
    axis.title.y = element_text(size = 7, color = "black", family = "Arial"))
RF
ggsave("figures/CH4_forest_soil.png",width=6,height=5)




#HEAT MAP
# Get the GHG flux data ---------------------------------------------------
#call up the dataset
# Load libraries
library(tidyverse)
library(pheatmap)
library(viridis)

#meta_data<-readRDS("full_data.RDS")
meta_data<-read.csv("raw/ch4_full_dataset.csv")

meta_data <- meta_data %>%
  group_by(plant_species) %>%
  filter(n_distinct(p_num) > 1) %>%
  ungroup()

meta_data <- meta_data %>%
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
  ))

meta_data <- meta_data %>%
  mutate(salinity_category = factor(salinity,
                                    levels = c("oligohaline", "mesohaline", "polyhaline", "euhaline")))

#average yi for each species and salinity category
species_salinity_avg <- meta_data %>%
  group_by(plant_species, salinity_category) %>%
  summarise(median_yi = median(yi, na.rm = TRUE), .groups = "drop")

species_matrix <- species_salinity_avg %>%
  pivot_wider(
    names_from = plant_species,
    values_from = median_yi  # NAs remain
  ) %>%
  column_to_rownames("salinity_category") %>%
  as.matrix()

#reorder rows based on salinity
species_matrix <- species_matrix[levels(meta_data$salinity_category), ]

#impute row means for clustering, bc they cannot be blank
species_matrix_imputed <- species_matrix
for (i in 1:nrow(species_matrix_imputed)) {
  row_mean <- mean(species_matrix_imputed[i, ], na.rm = TRUE)
  species_matrix_imputed[i, is.na(species_matrix_imputed[i, ])] <- row_mean
}

#euclidean distance and clustering on species
dist_mat <- dist(t(species_matrix_imputed), method = "euclidean")  # transpose to cluster columns
clust <- hclust(dist_mat, method = "complete")

#sort columns by overall flux
species_order <- order(colMeans(species_matrix, na.rm = TRUE), decreasing = TRUE)
species_matrix <- species_matrix[, species_order]

rownames(species_matrix) <- as.character(rownames(species_matrix))

asinh_trans <- trans_new( #function for plotting, made by Sawyer Balint
  name = "asinh",
  transform = asinh,
  inverse = sinh
)

color_breaks <- seq(min(species_matrix, na.rm=TRUE),
                    max(species_matrix, na.rm=TRUE),
                    length.out = 100)

legend_labels <- round(sinh(color_breaks), 0)
pretty_vals <- c(-10, 0, 10, 100, 1000)
legend_breaks <- asinh(pretty_vals)
library(pheatmap)
library(viridis)

map<-pheatmap(
  species_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = viridis(99),
  breaks = seq(min(species_matrix, na.rm=TRUE),
               max(species_matrix, na.rm=TRUE),
               length.out = 100),
  legend_breaks = legend_breaks,   # transformed positions
  legend_labels = pretty_vals,     # raw values to display
  fontsize_row = 6,
  fontsize_col = 6,
  angle_col = 45,
  na_col = "white",
  border_color = "black",
  fontsize = 6,
  silent = TRUE
)

map$gtable$grobs[[which(map$gtable$layout$name == "col_names")]]$gp <- gpar(fontface = "italic")
#wrap the gtable directly for patchwork
map_wrapped <- wrap_elements(map$gtable)
map_wrapped

combined_plot <- (RF + p8 + p2 + plot_layout(widths = c(1,1.5,1.5))) / (p6 + map_wrapped + plot_layout(widths = c(1.6,2.4))) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = 'a', 
    tag_suffix = "."  # use lower-case letters
  ) &
  theme(
    legend.text = element_markdown(size = 7), 
    plot.tag = element_text(family = "Arial Black", size = 12.5),
    legend.position = "bottom"           # place at bottom
    )

combined_plot

ggsave("figures/plot_test.png", combined_plot, width = 200, height = 175, dpi=300, units ="mm")

