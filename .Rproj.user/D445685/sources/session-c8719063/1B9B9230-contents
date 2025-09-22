# -------------------------------------------------------------------------
#Created 1.13.25 by Emily Wilson, last update 08.20.25 by Emily Wilson
#Title: Plant Species Drive Global Coastal Wetland Methane Fluxes Meta-analysis

#############################################################################
#############################################################################
#############################################################################
#Created 9.08.25 by Emily Wilson
#Supplemental Figure
rm(list=ls()) #clear the environment
#dev.off()



# install packages --------------------------------------------------------
library(scales) #for stats
library(tidyverse) #for data manipulation
library(patchwork) #to put plots together
library(mgcv) #for gam



#meta_data<-readRDS("pw_data.RDS")
meta_data<-read.csv("raw/ch4_soilsalinity_dataset.csv")

meta_data<-meta_data %>%
  filter(!is.na(salinity_conductivity), #no winter bc Poffenbarger et al. is growing season
         season!="winter")

str(meta_data)
meta_data <- meta_data %>% #clean up the data
  filter(!is.na(salinity_conductivity),
         salinity_conductivity < 60) %>%  #only include fluxes with soil salinity
  mutate(lat=abs(as.numeric(lat)), #make lat absolute for analysis
         plant_species=as.factor(plant_species),
         season=as.factor(season),
         salinity_conductivity=round(as.numeric(salinity_conductivity), 2),
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

# Step 1: Create a sequence of x values
line_data <- data.frame(x = seq(0,
                                max(meta_data$salinity_conductivity, na.rm = TRUE),
                                length.out = 100))

# Step 2: Add y values based on your equation
line_data <- line_data %>%
  mutate(y = (((-0.056 * x) + 1.38))/0.34)

line_data <- line_data %>%
  mutate(y = 10^y)

line_data <- line_data %>%
  mutate(y=y*(1/24)*(1/16.04)*(1000/1))


# Predict values for porewater -------------------------
# No winter -------------------------------------------------------------

# Get R²
model1<-lm(yi~salinity_conductivity, data=meta_data)
summary(model1)
r2 <- summary(model1)$r.squared
r2_label <- expression("_ This Study's Relationship between Salinity and CH"[4] * " flux, R"^2 * " = 0.08")

# Get predictions with prediction interval
preds <- predict(model1, newdata = meta_data, interval = "prediction", se.fit = TRUE)
meta_data <- meta_data %>%
  mutate(
    fit = exp(preds$fit[, "fit"]),
    fit.low = exp(preds$fit[, "lwr"]),
    fit.high = exp(preds$fit[, "upr"]))

p1<-ggplot(meta_data, aes(x = salinity_conductivity)) + 
  geom_point(aes(y = sinh(yi), shape=plant_species, fill=plant_species), size=2, alpha = 0.6) +
  geom_ribbon(aes(ymin = fit.low, ymax = fit.high), alpha = 0.1, fill = "black") +
  geom_line(aes(y = fit), color = "black", linewidth = 1) +
  geom_vline(xintercept=18, linetype="dashed", color="darkblue", linewidth=1)+
  geom_hline(yintercept=0, linetype="solid", color="black")+
  geom_line(data = line_data, aes(x = x, y =y), color = "darkred", linetype="dotdash", linewidth=1) +
  scale_y_continuous(transform = asinh_trans,
                     breaks=c(-1000,-100,-10,0,10,100,1000,10000))+
  theme_classic()+
  scale_shape_manual(values = c("Cyperus malaccensis"=21, "Spartina alterniflora"=22,
                                "Distichlis Spicata" = 8, "Spartina patens" = 24, 
                                "Suaeda salsa" = 25, "Phragmites australis" = 23,
                                "Cladium jamaicense" = 13, "Juncus sp." = 21,
                                "Juncae sp" = 21,
                                "Sesuvium portulacastrum" = 18, "Plantago maritima" = 10,
                                "Scirpus mariqueter" = 11, "Tamarix chinensis"=14,
                                "Salicornia sp."=12, "Suaeda sp." = 3))+
  scale_fill_manual(values = c("Cyperus malaccensis"="#025766", "Spartina patens" = "#5c3972", 
                               "Spartina alterniflora" = "#2FC45D", "Phragmites australis" ="#20bcff",
                               "Distichlis Spicata" = "#47724c", "Suaeda salsa" = "#c6b822",
                               "Cladium jamaicense" =  "#564e95", "Juncus sp." = "#86a291", "Juncae sp" = "lightblue", 
                               "Sesuvium portulacastrum" = "#90E669", "Plantago maritima" = "#30BCAD",
                               "Scirpus mariqueter" = "#60563f", "Tamarix chinensis"="#FF46A2",
                               "Salicornia sp."="black", "Suaeda sp." = "black")) +
  theme(legend.position="right")+
  theme(legend.title=element_blank())+
  annotate("text", x = 31, y = 2200,  # Adjust x and y based on your plot range
           label = expression(paste(". _ . _ Salinity-Based Proxy Equation Line")),
           size = 2.3, color = "darkred", family = "Arial")+
  annotate("text", x = 31, y = 3500,  # Adjust x and y based on your plot range
           label = expression(paste("---- 18 ppt Salinity Threshold")),
           size = 2.3, color = "darkblue", family = "Arial")+
  annotate("text", x = 30, y = 6000, label = r2_label, 
           , size = 2.3, family = "Arial" ) +
  theme(
    legend.text =element_text(size = 7, colour = "black", family = "Arial"),
    axis.text.x = element_text(size = 7, colour = "black", family = "Arial"), 
    axis.text.y = element_text(size = 7, colour = "black", family = "Arial"),
    axis.title.x = element_text(size = 7, color = "black", family = "Arial"),
    axis.title.y = element_text(size = 7, color = "black", family = "Arial")) +
  labs(x = "Salinity (ppt)", fill = "Plant Species", shape="Plant Species", y = expression(paste("CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), family="Arial")
p1

#meta_data<-readRDS("pw_data.RDS")
meta_data<-read.csv("raw/ch4_soilsalinity_dataset.csv")

meta_data<-meta_data %>%
  filter(!is.na(salinity_conductivity), #no winter bc Poffenbarger et al. is growing season
         season!="winter")

str(meta_data)
meta_data <- meta_data %>% #clean up the data
  filter(!is.na(salinity_conductivity),
         salinity_conductivity < 60) %>%  #only include fluxes with soil salinity
  mutate(lat=abs(as.numeric(lat)), #make lat absolute for analysis
         plant_species=as.factor(plant_species),
         season=as.factor(season),
         salinity_conductivity=round(as.numeric(salinity_conductivity), 2),
         ch4_method_simple=as.factor(ch4_method_simple),
         tide=as.factor(tide))%>%
  filter(plant_species!="Juncae sp")

summary(lm(yi~salinity_conductivity, data=meta_data))

meta_data <- meta_data %>%
  mutate(repeat_flag = ifelse("repeat" == "Y", TRUE, FALSE)) %>%
  mutate(group_key = ifelse(!repeat_flag,
                            paste(p_num, lat, ch4_chamber, estuary, marsh, plant_species, salinity, specific_site, section, sep = "_"),
                            NA)) %>%
  # Create numeric unique IDs for each group_key
  mutate(unique_id = ifelse(!repeat_flag,
                            as.integer(factor(group_key)),
                            NA)) %>%
  relocate(unique_id, .after = p_num)

meta_data <- meta_data %>%
  group_by(unique_id) %>%
  summarise(
    # Calculate mean of yi
    yi = mean(yi, na.rm = TRUE),
    
    # Calculate sum of vi (n = number of rows per group)
    vi = as.numeric(sum(vi)), #used to weight gams
    
    ch4_n=as.integer(sum(ch4_n)),
    
    salinity_conductivity = round(as.numeric(mean(salinity_conductivity), 2)),
    
    #Aggregate the columns for data analysis
    p_num = as.factor(paste(unique(p_num), collapse = ", ")),
    plant_species = as.factor(paste(unique(plant_species), collapse = ", ")),
    salinity = as.factor(paste(unique(salinity), collapse = ", ")),
    marsh = as.factor(paste(unique(marsh), collapse = ", ")),
    season = as.factor(paste(unique(season), collapse = ", ")),
    tide = as.factor(paste(unique(tide), collapse = ", ")),
    lat = paste(unique(lat), collapse = ", "),  
    tide = as.factor(paste(unique(tide), collapse = ", ")),
    ch4_method_simple = as.factor(paste(unique(ch4_method_simple), collapse = ", ")),
    climate_region = as.factor(paste(unique(climate_region), collapse = ", ")),
    month = as.factor(paste(unique(month), collapse = ", ")),
    #sampling_year = as.numeric(paste(unique(sampling_year), collapse = ", "))
    
    # Add other columns as needed here...
  )

# Fit model
model1 <- lm(yi ~ salinity_conductivity, data = meta_data)
summary(model1)
plot(model1)

# Get R²
r2 <- summary(model1)$r.squared
r2_label <- expression("_ This Study's Relationship between Salinity and CH"[4] * " flux, R"^2 * " = 0.05")


# Get predictions with prediction interval
preds <- predict(model1, newdata = meta_data, interval = "prediction", se.fit = TRUE)
meta_data <- meta_data %>%
  mutate(
    fit = exp(preds$fit[, "fit"]),
    fit.low = exp(preds$fit[, "lwr"]),
    fit.high = exp(preds$fit[, "upr"]))

p2<-ggplot(meta_data, aes(x = salinity_conductivity)) + 
  geom_point(aes(y = sinh(yi), shape=plant_species, fill=plant_species), size=2, alpha = 0.6) +
  geom_ribbon(aes(ymin = fit.low, ymax = fit.high), alpha = 0.1, fill = "black") +
  geom_line(aes(y = fit), color = "black", linewidth = 1) +
  geom_vline(xintercept=18, linetype="dashed", color="darkblue", linewidth=1)+
  geom_hline(yintercept=0, linetype="solid", color="black")+
  geom_line(data = line_data, aes(x = x, y =y), color = "darkred", linetype="dotdash", linewidth=1) +
  scale_y_continuous(transform = asinh_trans,
                     breaks=c(-1000,-100,-10,0,10,100,1000,10000))+
  theme_classic()+
  scale_shape_manual(values = c("Cyperus malaccensis"=21, "Spartina alterniflora"=22,
                                "Distichlis Spicata" = 8, "Spartina patens" = 24, 
                                "Suaeda salsa" = 25, "Phragmites australis" = 23,
                                "Cladium jamaicense" = 13, "Juncus sp." = 21,
                                "Juncae sp" = 21,
                                "Sesuvium portulacastrum" = 18, "Plantago maritima" = 10,
                                "Scirpus mariqueter" = 11, "Tamarix chinensis"=14,
                                "Salicornia sp."=12, "Suaeda sp." = 3))+
  scale_fill_manual(values = c("Cyperus malaccensis"="#025766", "Spartina patens" = "#5c3972", 
                               "Spartina alterniflora" = "#2FC45D", "Phragmites australis" ="#20bcff",
                               "Distichlis Spicata" = "#47724c", "Suaeda salsa" = "#c6b822",
                               "Cladium jamaicense" =  "#564e95", "Juncus sp." = "#86a291", "Juncae sp" = "lightblue", 
                               "Sesuvium portulacastrum" = "#90E669", "Plantago maritima" = "#30BCAD",
                               "Scirpus mariqueter" = "#60563f", "Tamarix chinensis"="#FF46A2",
                               "Salicornia sp."="black", "Suaeda sp." = "black")) +
  theme(legend.position="none")+
  annotate("text", x = 31, y = 2200,  # Adjust x and y based on your plot range
           label = expression(paste(". _ . _ Salinity-Based Proxy Equation Line")),
           size = 2.3, color = "darkred", family = "Arial")+
  annotate("text", x = 31, y = 3500,  # Adjust x and y based on your plot range
           label = expression(paste("---- 18 ppt Salinity Threshold")),
           size = 2.3, color = "darkblue", family = "Arial")+
  annotate("text", x = 30, y = 6000, label = r2_label, 
           , size = 2.3, family = "Arial" ) +
  theme(
    axis.text.x = element_text(size = 7, colour = "black", family = "Arial"), 
    axis.text.y = element_text(size = 7, colour = "black", family = "Arial"),
    axis.title.x = element_text(size = 7, color = "black", family = "Arial"),
    axis.title.y = element_text(size = 7, color = "black", family = "Arial")) +
  labs(x = "Salinity (ppt)", fill = "Plant Species", shape="Plant Species", y = expression(paste("CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), family="Arial")
p2
p2 <- p2 + theme(legend.position = "none")
p2 <- p2 + guides(fill = "none", color = "none", shape = "none")
library(patchwork)

combined_plot <- (p1/p2) +
  plot_layout(nrow = 2,guides = "collect") +
  plot_annotation(
    tag_levels = 'a', 
    tag_suffix = "."  # use lower-case letters
  ) &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6),
    plot.tag = element_text(family = "Arial Black", size = 20)
  )
combined_plot
ggsave("figures/supp_fig.png", combined_plot, width = 170, height = 260, dpi=300, units = "mm")



#meta_data<-readRDS("pw_data.RDS")
meta_data<-read.csv("raw/ch4_soilsalinity_dataset.csv")
library(purrr)
library(broom)

# 1. Filter species with >10 complete observations
species_to_keep <- meta_data %>%
  filter(!is.na(yi) & !is.na(salinity_conductivity)) %>%
  group_by(plant_species) %>%
  filter(n() > 15) %>%
  ungroup()

# 2. Safe GAM function in case a species still causes errors
safe_gam <- possibly(
  ~ gam(yi ~ s(salinity_conductivity, k=3), data = ., weights = .$vi),
  otherwise = NULL
)

#Fit models and generate predictions
models <- species_to_keep %>%
  group_by(plant_species) %>%
  nest() %>%
  mutate(
    model = map(data, safe_gam),
    r2 = map_dbl(model, ~ if (!is.null(.x)) summary(.x)$r.sq else NA),
    pred_grid = map2(data, model, ~ {
      if (is.null(.y)) return(NULL)
      newdata <- data.frame(salinity_conductivity = seq(min(.x$salinity_conductivity, na.rm = TRUE),
                                                        max(.x$salinity_conductivity, na.rm = TRUE),
                                                        length.out = 100))
      newdata$yi_pred <- predict(.y, newdata = newdata)
      newdata
    })) %>%
  filter(!map_lgl(model, is.null))  # Drop failed fits


#Unnest prediction and observation data separately
pred_data <- models %>%
  select(plant_species, pred_grid, r2) %>%
  unnest(pred_grid) %>%
  mutate(label = paste0("R² = ", round(r2, 2)))

models <- species_to_keep %>%
  group_by(plant_species) %>%
  nest() %>%
  mutate(
    model = map(data, safe_gam),
    r2 = map_dbl(model, ~ if (!is.null(.x)) summary(.x)$r.sq else NA),
    pval = map_dbl(model, ~ {
      if (is.null(.x)) return(NA_real_)
      s <- summary(.x)
      # Example: get the p-value for the first smooth term
      if (!is.null(s$s.table)) {
        s$s.table[1, "p-value"]
      } else if (!is.null(s$p.table)) {
        s$p.table[1, "Pr(>|t|)"]
      } else {
        NA_real_
      }
    }),
    pred_grid = map2(data, model, ~ {
      if (is.null(.y)) return(NULL)
      newdata <- data.frame(
        salinity_conductivity = seq(min(.x$salinity_conductivity, na.rm = TRUE),
                                    max(.x$salinity_conductivity, na.rm = TRUE),
                                    length.out = 100)
      )
      newdata$yi_pred <- predict(.y, newdata = newdata)
      newdata
    })
  ) %>%
  filter(!map_lgl(model, is.null))

# Example: unnesting predictions with R² and p-values
pred_data <- models %>%
  select(plant_species, pred_grid, r2, pval) %>%
  unnest(pred_grid) %>%
  mutate(label = paste0("R² = ", round(r2, 2), 
                        ", p = ", signif(pval, 2)))

obs_data <- models %>%
  select(plant_species, data) %>%
  unnest(data)

# 5. Plot: raw points, GAM smooths, R² labels
ggplot() +
  geom_point(data = obs_data, aes(x = salinity_conductivity, y = sinh(yi), fill=plant_species, shape=plant_species), alpha = 0.5, size=2.5) +
  geom_line(data = pred_data, aes(x = salinity_conductivity, y = sinh(yi_pred)), color = "#025766", size = 1) +
  scale_shape_manual(values = c("Cyperus malaccensis"=21, "Spartina alterniflora"=22,
                                "Distichlis Spicata" = 8, "Spartina patens" = 24, 
                                "Suaeda salsa" = 25, "Phragmites australis" = 23,
                                "Cladium jamaicense" = 13, "Juncae sp" = 20,
                                "Sesuvium portulacastrum" = 18, "Plantago maritima" = 10,
                                "Scirpus mariqueter" = 11, "Tamarix chinensis"=14,
                                "Salicornia sp."=12, "Suaeda sp." = 3))+
  scale_fill_manual(values = c("Cyperus malaccensis"="#025766", "Spartina patens" = "#5c3972", 
                               "Spartina alterniflora" = "#2FC45D", "Phragmites australis" ="#20bcff",
                               "Distichlis Spicata" = "#47724c", "Suaeda salsa" = "#c6b822",
                               "Cladium jamaicense" =  "#564e95", "Juncae sp" = "#86a291",
                               "Sesuvium portulacastrum" = "#90E669", "Plantago maritima" = "#47724c",
                               "Scirpus mariqueter" = "#60563f", "Tamarix chinensis"="#FF46A2",
                               "Salicornia sp."="black", "Suaeda sp." = "black")) +
  scale_y_continuous(transform = asinh_trans,
                     breaks=c(-1000,-100,-10,0,10,100,1000,10000))+
  geom_hline(yintercept=0, linetype="dashed", color="black")+
  geom_text(
    data = pred_data %>% group_by(plant_species) %>% slice(1),
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.1, vjust = 1.2, inherit.aes = FALSE, size = 3,
    family="Helvetica") +
  facet_wrap(~ plant_species) +
  theme_classic() +
  theme(
    legend.position="none",
    axis.text.x = element_text(size = 11, colour = "black", family = "Arial"), 
    axis.text.y = element_text(size = 11, colour = "black", family = "Arial"),
    axis.title.x = element_text(size = 11, color = "black", family = "Arial"),
    axis.title.y = element_text(size = 1, color = "black", family = "Arial")) +
  labs(x = "Salinity", y =  expression(paste("CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")))
ggsave("figures/ch4_salinity_gam_facet.png",width=9,height=6)


