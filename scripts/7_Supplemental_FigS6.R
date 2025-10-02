#last update 09.29.25 by Sawyer Balint
#Title: Plant Species Drive Global Coastal Wetland Methane Fluxes Meta-analysis
#SMD analysis

# load packages -----------------------------------------------------------

rm(list = ls()) #clear environment

library(tidyverse)
library(metafor)
library(patchwork)

# import data -------------------------------------------------------------

raw.df <- read.csv("raw/SMD_data.csv")

# graphing theme ----------------------------------------------------------

#custom graphing theme to reduce repetition
theme <- list(
  theme_linedraw(),
  geom_point(shape=21, size = 3),
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40"),
  scale_color_manual(values = c("Cyperus malaccensis"="#025766", "Spartina patens" = "#5c3972", 
                                "Spartina alterniflora" = "#2FC45D", "Phragmites australis" ="#20bcff",
                                "Distichlis Spicata" = "#47724c", "Suaeda salsa" = "#c6b822",
                                "Cladium jamaicense" =  "#564e95", "Juncae sp" = "lightblue",
                                "Juncus sp." = "#86a291",
                                "Sesuvium portulacastrum" = "#90E669", "Plantago maritima" = "#47724c",
                                "Scirpus mariqueter" = "#60563f", "Tamarix chinensis"="#FF46A2",
                                "Salicornia sp."="black", "Suaeda sp." = "black")),
    labs(
      x = "Standardized Mean Difference (SMD)",
      y = NULL),
    theme(
      legend.text =element_text(size = 12, colour = "black", family = "Helvetica"),
      axis.text.x = element_text(size = 10, colour = "black", family = "Helvetica"), 
      axis.text.y = element_text(size = 11, colour = "black", family = "Helvetica"),
      axis.title.x = element_text(size = 11, color = "black", family = "Helvetica"),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.position="none",
      panel.grid.major.y = element_blank()
    )
  
)

#function to calculate effect size
calc_es <- function(data, species){
  
  #separate out Spartina alterniflora rows
  df1 <- data %>%
    filter(plant_species == species) %>%
    select(p_num, month, sampling_year,
           ch4_var_t = ch4_var, 
           ch4_n_t = ch4_n, 
           ch4_flux_t = ch4_flux)
  
  #filter non-Spartina rows
  df2 <- data %>%
    filter(plant_species != species)
  
  #join non-Spartina with matching Spartina rows by p_num and season
  paired.df <- df2 %>%
    left_join(df1, by = c("p_num", "sampling_year", "month"), relationship = "many-to-many") %>%
    filter(!is.na(ch4_flux_t))  # keep only pairs where Spartina data exists
  
  meta.df <- escalc(measure = "SMD", 
                    m1i = ch4_flux_t, 
                    sd1i = ch4_var_t, 
                    n1i = ch4_n_t,
                    m2i = ch4_flux, 
                    sd2i = ch4_var, 
                    n2i = ch4_n, 
                    data = paired.df)
  
  return(meta.df)
  
}

#function to calculate average yi and vi
calc_yi <- function(data){
  
  data %>%
    summarise(
      avg_yi = mean(yi, na.rm = TRUE),
      avg_vi = sum(vi, na.rm = TRUE) / (n()^2),
      ch4_n = sum(ch4_n),
      p_num = paste(unique(p_num), collapse = ", "),
      salinity = paste(unique(salinity), collapse = ", "),  
    )
}

#function to calculate average CH4 flux
calc_mean <- function(data){
  
  data %>%
    summarise(
      ch4_flux = mean(ch4_flux, na.rm = TRUE),
      ch4_n = sum(ch4_n),
      ch4_var=mean(ch4_var),
      
      plant_species = as.factor(paste(unique(plant_species), collapse = ", ")),
      salinity = as.factor(paste(unique(salinity), collapse = ", ")),
      month = as.factor(paste(unique(month), collapse = ", ")),
      sampling_year = as.numeric(paste(unique(sampling_year), collapse = ", ")),
      p_num = as.numeric(paste(unique(p_num), collapse = ", "))
    )
}

#function to calculate confidence intrval
calc_ci <- function(data, model){
  
    data %>%
      mutate(
        yi = model$yi,
        vi = model$vi,
        se = sqrt(vi),
        ci_lower = yi - 1.96 * se,
        ci_upper = yi + 1.96 * se
      )
  
}
  
# perform SMD analysis ----------------------------------------------------

df <- raw.df %>% #clean up the data
  drop_na(ch4_flux, ch4_n, ch4_var, plant_species) %>%
  filter(plant_species != "") %>%
  group_by(plant_species) %>%
  filter(n_distinct(p_num) > 1) %>% #make sure there is more than one study that measures the plant species
  ungroup()
  
#Data from the same study, latitude, estuary, marsh, year, month, plant species, and salinity are assign same ID.
df <- df %>%
  mutate(repeat_flag = ifelse("repeat" == "Y", TRUE, FALSE)) %>%
  mutate(group_key = ifelse(!repeat_flag,
                            paste(p_num, lat, estuary, marsh, sampling_year, month, plant_species, salinity, specific_site, section, sep = "_"),
                            NA)) %>%
  # Create numeric unique IDs for each group_key
  mutate(unique_id = ifelse(!repeat_flag,
                            as.integer(factor(group_key)),
                            NA)) %>%
  relocate(unique_id, .after = p_num)

meta.df <- calc_es(df, "Spartina alterniflora")

# graph by plant species --------------------------------------------------

#Count number of unique studies per plant_species
count.df <- meta.df %>%
  group_by(plant_species) %>%
  summarise(n_studies = n_distinct(p_num), .groups = "drop")

#Aggregate the data
aggregated.df <- meta.df %>%
  group_by(plant_species) %>%
  calc_yi() %>%
  left_join(count.df, by = "plant_species") %>%
  mutate(slab_label = paste0(plant_species, " | S=", n_studies, " | n=", ch4_n))

model <- rma(avg_yi, avg_vi, data = aggregated.df, method = "REML", knha=TRUE)

# Extract model results
plot.df <- aggregated.df %>%
  calc_ci(model)

p1<- ggplot(plot.df, aes(x = yi, y = reorder(slab_label, yi), color=plant_species)) +
  theme+
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color=plant_species), height = 0.2) +
  scale_x_continuous(limits = c(-max(plot.df$ci_upper), max(plot.df$ci_upper) + 0.2)) +
  labs(title = "Spartina alterniflora")

# by salinity -------------------------------------------------------------

salinity.df <- df %>%
  group_by(unique_id) %>%
  calc_mean()

meta.df <- calc_es(salinity.df, "Spartina alterniflora")

# Step 1: Count number of unique studies per plant_species
count.df <- meta.df %>%
  group_by(salinity) %>%
  summarise(n_studies = n_distinct(p_num), .groups = "drop")

aggregated.df <- meta.df %>%
  group_by(salinity) %>%
  calc_yi()%>%
  left_join(count.df, by = "salinity") %>%
  mutate(slab_label = paste0(salinity, " | S=", n_studies, " | n=", ch4_n)) %>%
  mutate(salinity = factor(salinity, 
                           levels = c("oligohaline", "mesohaline", "polyhaline", "euhaline"), 
                           ordered = TRUE))

model <- rma(avg_yi, avg_vi, data = aggregated.df, method = "REML", knha=TRUE)
summary(model)

# Extract model results
plot.df <- aggregated.df %>%
  calc_ci(model)%>%
  mutate(
    salinity = factor(salinity, 
                      levels = c("oligohaline", "mesohaline", "polyhaline", "euhaline")),
    slab_label = fct_reorder(slab_label, as.numeric(salinity))  # order by salinity
  )


p2 <- ggplot(plot.df, aes(x = yi, y = slab_label, color=salinity)) +
  theme+
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color=salinity), height = 0.2) +
  scale_color_manual(values = c("oligohaline" = "lightblue", 
                                "mesohaline" = "#1F78B4",
                                "polyhaline" = "darkblue",
                                "euhaline" = "darkgreen")) +
  scale_x_continuous(limits = c(-max(plot.df$ci_upper), max(plot.df$ci_upper) + 0.2)) +
  labs(title = "Spartina alterniflora"
  )

# Phragmites --------------------------------------------------------------

meta.df <- calc_es(df, "Phragmites australis")

#count number of unique studies per plant_species
count.df <- meta.df %>%
  group_by(plant_species) %>%
  summarise(n_studies = n_distinct(p_num), .groups = "drop")

#aggregate the data
aggregated.df <- meta.df %>%
  group_by(plant_species) %>%
  calc_yi() %>%
  left_join(count.df, by = "plant_species") %>%  
  mutate(slab_label = paste0(plant_species, " | S=", n_studies, " | n=", ch4_n))

model <- rma(avg_yi, avg_vi, data = aggregated.df, method = "REML", knha=TRUE)
summary(model)

# Extract model results
plot.df <- aggregated.df %>%
  calc_ci(model)

p3<- ggplot(plot.df, aes(x = yi, y = reorder(slab_label, yi), color=plant_species)) +
  theme+
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color=plant_species), height = 0.2) +
  scale_x_continuous(limits = c(-max(plot.df$ci_upper), max(plot.df$ci_upper) + 0.2)) +
  labs(title = "Phragmites australis"
  )

# by salinity -------------------------------------------------------------

meta.df <- df %>%
  group_by(unique_id) %>%
  calc_mean() %>%
  calc_es("Phragmites australis")

count.df <- meta.df %>%
  group_by(salinity) %>%
  summarise(n_studies = n_distinct(p_num), .groups = "drop")

aggregated.df <- meta.df %>%
  group_by(salinity) %>%
  calc_yi() %>%
  left_join(count.df, by = "salinity") %>%
  mutate(slab_label = paste0(salinity, " | S=", n_studies, " | n=", ch4_n)) %>%
  mutate(salinity = factor(salinity, 
                           levels = c("oligohaline", "mesohaline", "polyhaline", "euhaline"), 
                           ordered = TRUE))

model <- rma(avg_yi, avg_vi, data = aggregated.df, method = "REML", knha=TRUE)
summary(model)

# Extract model results
plot.df <- aggregated.df %>%
  calc_ci(model) %>%
  mutate(salinity = factor(salinity, 
                           levels = c("oligohaline", "mesohaline", "polyhaline", "euhaline"), 
                           ordered = TRUE),
         slab_label = fct_reorder(slab_label, as.numeric(salinity))
  )

p4 <- ggplot(plot.df, aes(x = yi, y = slab_label, color=salinity)) +
  theme+
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color=salinity), height = 0.2) +
  scale_color_manual(values = c("oligohaline" = "lightblue", 
                                "mesohaline" = "#1F78B4",
                                "polyhaline" = "darkblue",
                                "euhaline" = "darkgreen")) +
  scale_x_continuous(limits = c(-max(plot.df$ci_upper), max(plot.df$ci_upper) + 0.2)) +
  labs(title = "Phragmites australis"
  )

# combine plots -----------------------------------------------------------

combined_plot <- wrap_plots(p1 + p2+ p3+p4) + 
  plot_annotation(tag_levels = 'a') &
  theme(
    plot.tag = element_text(face = "bold"))

combined_plot

ggsave("figures/Figure_S6.png",width=12,height=8)

