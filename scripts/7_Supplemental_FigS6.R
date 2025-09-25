rm(list = ls()) #clear environment
library(tidyverse)
library(metafor)
library(dplyr)

#ghgflux <- readRDS("ghgflux.rds") #get data
ghgflux<-read.csv("raw/SMD_data.csv")
data_ch4 <- ghgflux %>% #clean up the data
  filter(
    !is.na(ch4_flux),
    !is.na(ch4_n),
    !is.na(ch4_var),
    !is.na(plant_species),
    plant_species != "",
  )

data_ch4 <- data_ch4 %>%
  group_by(plant_species) %>%
  filter(n_distinct(p_num) > 1) %>% #make sure there is more than one study that measures the plant species
  ungroup()

# Assign unique ID to aggregate data --------------------------------------
#Data from the same study, latitude, estuary, marsh, year, month, plant species, and salinity are assign same ID.
data_ch4 <- data_ch4 %>%
  mutate(repeat_flag = ifelse("repeat" == "Y", TRUE, FALSE)) %>%
  mutate(group_key = ifelse(!repeat_flag,
                            paste(p_num, lat, estuary, marsh, sampling_year, month, plant_species, salinity, specific_site, section, sep = "_"),
                            NA)) %>%
  # Create numeric unique IDs for each group_key
  mutate(unique_id = ifelse(!repeat_flag,
                            as.integer(factor(group_key)),
                            NA)) %>%
  relocate(unique_id, .after = p_num)

meta_data <- data_ch4 %>%
  dplyr::select(ch4_flux, ch4_var, ch4_n, sampling_year, unique_id, p_num, lat, month, ch4_method_simple, season, plant_species, salinity, tide)

#separate out Spartina alterniflora rows
spartina <- meta_data %>%
  filter(plant_species == "Spartina alterniflora") %>%
  dplyr::select(p_num, month, sampling_year,
         ch4_var_t = ch4_var, 
         ch4_n_t = ch4_n, 
         ch4_flux_t = ch4_flux, 
        )

#filter non-Spartina rows
non_spartina <- meta_data %>%
  filter(plant_species != "Spartina alterniflora")

#join non-Spartina with matching Spartina rows by p_num and season
paired_data <- non_spartina %>%
  left_join(spartina, by = c("p_num", "sampling_year", "month"), relationship = "many-to-many") %>%
  filter(!is.na(ch4_flux_t))  # keep only pairs where Spartina data exists

ch4_meta_data <- escalc(measure = "SMD", 
                        m1i = ch4_flux_t, 
                        sd1i = ch4_var_t, 
                        n1i = ch4_n_t,
                        m2i = ch4_flux, 
                        sd2i = ch4_var, 
                        n2i = ch4_n, 
                        data = paired_data)

aggregated_data <- ch4_meta_data %>%
  group_by(unique_id) %>%
  summarise(
    # Calculate mean of yi
    avg_yi = mean(yi, na.rm = TRUE),
    
    # Calculate sum of vi / n^2 (n = number of rows per group)
    avg_vi = sum(vi, na.rm = TRUE) / (n()^2),
    
    ch4_n=sum(ch4_n),
    
    # Example aggregation for other columns
    sampling_year = paste(unique(sampling_year), collapse = ", "),
    p_num = paste(unique(p_num), collapse = ", "),
    season = paste(unique(season), collapse = ", "),
    salinity = paste(unique(salinity), collapse = ", "),  
    plant_species = paste(unique(plant_species), collapse = ", "),
    tide = paste(unique(tide), collapse = ", ")

    # Add other columns as needed here...
  )

#Count number of unique studies per plant_species
study_counts <- ch4_meta_data %>%
  group_by(plant_species) %>%
  summarise(n_studies = n_distinct(p_num), .groups = "drop")

#Aggregate the data
plants_aggregated <- ch4_meta_data %>%
  group_by(plant_species) %>%
  summarise(
    avg_yi = mean(yi, na.rm = TRUE),
    avg_vi = sum(vi, na.rm = TRUE) / (n()^2),
    ch4_n = sum(ch4_n),
    p_num = paste(unique(p_num), collapse = ", "),
    salinity = paste(unique(salinity), collapse = ", "),  
    tide = paste(unique(tide), collapse = ", "),
    ch4_method_simple = paste(unique(ch4_method_simple), collapse = ", ")
  ) %>%
  left_join(study_counts, by = "plant_species") %>%  # ← This is key!
  mutate(slab_label = paste0(plant_species, " | S=", n_studies, " | n=", ch4_n))

print(plants_aggregated)

meta_model <- rma(avg_yi, avg_vi, data = plants_aggregated, method = "REML", knha=TRUE)
summary(meta_model)

# Extract model results
plot_data <- plants_aggregated %>%
  mutate(
    yi = meta_model$yi,
    vi = meta_model$vi,
    se = sqrt(vi),
    ci_lower = yi - 1.96 * se,
    ci_upper = yi + 1.96 * se
  )


p1<- ggplot(plot_data, aes(x = yi, y = reorder(slab_label, yi), color=plant_species)) +
  geom_point(aes(color = plant_species), size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color=plant_species), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("Cyperus malaccensis"="#025766", "Spartina patens" = "#5c3972", 
                               "Spartina alterniflora" = "#2FC45D", "Phragmites australis" ="#20bcff",
                               "Distichlis Spicata" = "#47724c", "Suaeda salsa" = "#c6b822",
                               "Cladium jamaicense" =  "#564e95", "Juncae sp" = "lightblue",
                               "Juncus sp." = "#86a291",
                               "Sesuvium portulacastrum" = "#90E669", "Plantago maritima" = "#47724c",
                               "Scirpus mariqueter" = "#60563f", "Tamarix chinensis"="#FF46A2",
                               "Salicornia sp."="black", "Suaeda sp." = "black")) +
  scale_x_continuous(limits = c(-max(plot_data$ci_upper), max(plot_data$ci_upper) + 0.2)) +
  labs(
    x = "Standardized Mean Difference (SMD)",
    y = NULL,
    title = "Spartina alterniflora"
  ) +
  theme(
  legend.text =element_text(size = 12, colour = "black", family = "Helvetica"),
    axis.text.x = element_text(size = 12, colour = "black", family = "Helvetica"), 
    axis.text.y = element_text(size = 12, colour = "black", family = "Helvetica"),
    axis.title.x = element_text(size = 12, color = "black", family = "Helvetica"),
    axis.title.y = element_text(size = 12, color = "black", family = "Helvetica")) +
  theme(plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA))+
  theme_linedraw()+
  theme(legend.position="none",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 11)
  )
p1

# by salinity -------------------------------------------------------------
data_ch4 <- data_ch4 %>%
  mutate(repeat_flag = ifelse("repeat" == "Y", TRUE, FALSE)) %>%
  mutate(group_key = ifelse(!repeat_flag,
                            paste(p_num, lat, estuary, marsh, sampling_year, month, plant_species, salinity, specific_site, section, sep = "_"),
                            NA)) %>%
  # Create numeric unique IDs for each group_key
  mutate(unique_id = ifelse(!repeat_flag,
                            as.integer(factor(group_key)),
                            NA)) %>%
  relocate(unique_id, .after = p_num)


ch4_meta_data <- data_ch4 %>%
  group_by(unique_id) %>%
  summarise(
    # Calculate mean of yi
    ch4_flux = mean(ch4_flux, na.rm = TRUE),
    
    # Calculate sum of vi / n^2 (n = number of rows per group)
    ch4_n = sum(ch4_n),
    
    ch4_var=mean(ch4_var),
    
    # Example aggregation for other columns
    month = as.factor(paste(unique(month), collapse = ", ")),
    p_num = as.numeric(paste(unique(p_num), collapse = ", ")),
    plant_species = as.factor(paste(unique(plant_species), collapse = ", ")),
    salinity = as.factor(paste(unique(salinity), collapse = ", ")),
    season = as.factor(paste(unique(season), collapse = ", ")),
    tide = as.factor(paste(unique(tide), collapse = ", ")),
    lat = paste(unique(lat), collapse = ", "),  
    tide = as.factor(paste(unique(tide), collapse = ", ")),
    ch4_method_simple = as.factor(paste(unique(ch4_method_simple), collapse = ", ")),
    climate_region = as.factor(paste(unique(climate_region), collapse = ", ")),
    month = as.factor(paste(unique(month), collapse = ", ")),
    sampling_year = as.numeric(paste(unique(sampling_year), collapse = ", "))
    
    # Add other columns as needed here...
  )

meta_data <- ch4_meta_data %>%
  dplyr::select(ch4_flux, ch4_var, ch4_n, p_num, sampling_year, lat, month, ch4_method_simple, season, plant_species, salinity, tide)

# Step 1: Separate out Spartina alterniflora rows
spartina <- meta_data %>%
  filter(plant_species == "Spartina alterniflora") %>%
  dplyr::select(p_num, month, sampling_year,
                ch4_var_t = ch4_var, 
                ch4_n_t = ch4_n, 
                ch4_flux_t = ch4_flux, 
  )

# Step 2: Filter non-Spartina rows
non_spartina <- meta_data %>%
  filter(plant_species != "Spartina alterniflora")

# Step 3: Join non-Spartina with matching Spartina rows by p_num and season
paired_data <- non_spartina %>%
  left_join(spartina, by = c("p_num", "month", "sampling_year"), relationship = "many-to-many") %>%
  filter(!is.na(ch4_flux_t))  # keep only pairs where Spartina data exists

ch4_meta_data <- escalc(measure = "SMD", 
                        m1i = ch4_flux_t, 
                        sd1i = ch4_var_t, 
                        n1i = ch4_n_t,
                        m2i = ch4_flux, 
                        sd2i = ch4_var, 
                        n2i = ch4_n, 
                        data = paired_data)

# Step 1: Count number of unique studies per plant_species
study_counts <- ch4_meta_data %>%
  group_by(salinity) %>%
  summarise(n_studies = n_distinct(p_num), .groups = "drop")

sal_aggregated <- ch4_meta_data %>%
  group_by(salinity) %>%
  summarise(
    avg_yi = mean(yi, na.rm = TRUE),
    avg_vi = sum(vi, na.rm = TRUE) / (n()^2),
    ch4_n = sum(ch4_n),
    p_num = paste(unique(p_num), collapse = ", "),
    tide = paste(unique(tide), collapse = ", "),
    ch4_method_simple = paste(unique(ch4_method_simple), collapse = ", ")
  ) %>%
  left_join(study_counts, by = "salinity") %>%  # ← This is key!
  mutate(slab_label = paste0(salinity, " | S=", n_studies, " | n=", ch4_n))

sal_aggregated <- sal_aggregated %>%
  mutate(salinity = factor(salinity, 
                           levels = c("oligohaline", "mesohaline", "polyhaline", "euhaline"), 
                           ordered = TRUE))

print(sal_aggregated)

meta_model <- rma(avg_yi, avg_vi, data = sal_aggregated, method = "REML", knha=TRUE)
summary(meta_model)


levels(sal_aggregated$salinity)

# Extract model results
plot_data <- sal_aggregated %>%
  mutate(
    yi = meta_model$yi,
    vi = meta_model$vi,
    se = sqrt(vi),
    ci_lower = yi - 1.96 * se,
    ci_upper = yi + 1.96 * se
  )

plot_data <- plot_data %>%
  mutate(
    salinity = factor(salinity, 
                      levels = c("oligohaline", "mesohaline", "polyhaline", "euhaline")),
    slab_label = fct_reorder(slab_label, as.numeric(salinity))  # order by salinity
  )


p2 <- ggplot(plot_data, aes(x = yi, y = slab_label, color=salinity)) +
  geom_point(aes(color = salinity), size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color=salinity), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("oligohaline" = "lightblue", 
                                "mesohaline" = "#1F78B4",
                                "polyhaline" = "darkblue",
                                "euhaline" = "darkgreen")) +
  scale_x_continuous(limits = c(-max(plot_data$ci_upper), max(plot_data$ci_upper) + 0.2)) +
  labs(
    x = "Standardized Mean Difference (SMD)",
    y = NULL,
    title = "Spartina alterniflora"
  ) +
  theme(
    legend.text =element_text(size = 12, colour = "black", family = "Helvetica"),
    axis.text.x = element_text(size = 12, colour = "black", family = "Helvetica"), 
    axis.text.y = element_text(size = 12, colour = "black", family = "Helvetica"),
    axis.title.x = element_text(size = 12, color = "black", family = "Helvetica"),
    axis.title.y = element_text(size = 12, color = "black", family = "Helvetica")) +
  theme(plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA))+
  theme_linedraw()+
  theme(legend.position="none",
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 11)
  )
p2

# Phragmites --------------------------------------------------------------
data_ch4 <- ghgflux %>% #clean up the data
  filter(
    !is.na(ch4_flux),
    !is.na(ch4_n),
    !is.na(plant_species),
    plant_species != "",
  )

data_ch4 <- data_ch4 %>%
  group_by(plant_species) %>%
  filter(n_distinct(p_num) > 1) %>%
  ungroup()

# Assign unique ID to aggregate data --------------------------------------
#Data from the same study, latitude, estuary, marsh, year, month, plant species, and salinity are assign same ID.
data_ch4 <- data_ch4 %>%
  mutate(repeat_flag = ifelse("repeat" == "Y", TRUE, FALSE)) %>%
  mutate(group_key = ifelse(!repeat_flag,
                            paste(p_num, lat, estuary, marsh, sampling_year, month, plant_species, salinity, specific_site, section, sep = "_"),
                            NA)) %>%
  # Create numeric unique IDs for each group_key
  mutate(unique_id = ifelse(!repeat_flag,
                            as.integer(factor(group_key)),
                            NA)) %>%
  relocate(unique_id, .after = p_num)

meta_data <- data_ch4 %>%
  dplyr::select(ch4_flux, ch4_var, ch4_n, sampling_year, unique_id, p_num, lat, month, ch4_method_simple, season, plant_species, salinity, tide)

#separate out phrag rows
phrag <- meta_data %>%
  filter(plant_species == "Phragmites australis") %>%
  dplyr::select(p_num, month, sampling_year,
                ch4_var_t = ch4_var, 
                ch4_n_t = ch4_n, 
                ch4_flux_t = ch4_flux, 
  )

#filter non-phrag rows
non_phrag <- meta_data %>%
  filter(plant_species != "Phragmites australis")

#join non-phrag with matching phrag rows by p_num and season
paired_data <- non_phrag %>%
  left_join(phrag, by = c("p_num", "month", "sampling_year"), relationship = "many-to-many") %>%
  filter(!is.na(ch4_flux_t))  # keep only pairs where Spartina data exists

ch4_meta_data <- escalc(measure = "SMD", 
                        m1i = ch4_flux_t, 
                        sd1i = ch4_var_t, 
                        n1i = ch4_n_t,
                        m2i = ch4_flux, 
                        sd2i = ch4_var, 
                        n2i = ch4_n, 
                        data = paired_data)

aggregated_data <- ch4_meta_data %>%
  group_by(unique_id) %>%
  summarise(
    # Calculate mean of yi
    avg_yi = mean(yi, na.rm = TRUE),
    
    # Calculate sum of vi / n^2 (n = number of rows per group)
    avg_vi = sum(vi, na.rm = TRUE) / (n()^2),
    
    ch4_n=sum(ch4_n),
    
    # Example aggregation for other columns
    p_num = paste(unique(p_num), collapse = ", "),
    season = paste(unique(season), collapse = ", "),
    salinity = paste(unique(salinity), collapse = ", "),  
    plant_species = paste(unique(plant_species), collapse = ", "),
    tide = paste(unique(tide), collapse = ", ")
    
    # Add other columns as needed here...
  )

#count number of unique studies per plant_species
study_counts <- ch4_meta_data %>%
  group_by(plant_species) %>%
  summarise(n_studies = n_distinct(p_num), .groups = "drop")

#aggregate the data
plants_aggregated <- ch4_meta_data %>%
  group_by(plant_species) %>%
  summarise(
    avg_yi = mean(yi, na.rm = TRUE),
    avg_vi = sum(vi, na.rm = TRUE) / (n()^2),
    ch4_n = sum(ch4_n),
    p_num = paste(unique(p_num), collapse = ", "),
    salinity = paste(unique(salinity), collapse = ", "),  
    tide = paste(unique(tide), collapse = ", "),
    ch4_method_simple = paste(unique(ch4_method_simple), collapse = ", ")
  ) %>%
  left_join(study_counts, by = "plant_species") %>%  # ← This is key!
  mutate(slab_label = paste0(plant_species, " | S=", n_studies, " | n=", ch4_n))

print(plants_aggregated)

meta_model <- rma(avg_yi, avg_vi, data = plants_aggregated, method = "REML", knha=TRUE)
summary(meta_model)

# Extract model results
plot_data <- plants_aggregated %>%
  mutate(
    yi = meta_model$yi,
    vi = meta_model$vi,
    se = sqrt(vi),
    ci_lower = yi - 1.96 * se,
    ci_upper = yi + 1.96 * se
  )

p3<- ggplot(plot_data, aes(x = yi, y = reorder(slab_label, yi), color=plant_species)) +
  geom_point(aes(color = plant_species), size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color=plant_species), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("Cyperus malaccensis"="#025766", "Spartina patens" = "#5c3972", 
                                "Spartina alterniflora" = "#2FC45D", "Phragmites australis" ="#20bcff",
                                "Distichlis Spicata" = "#47724c", "Suaeda salsa" = "#c6b822",
                                "Cladium jamaicense" =  "#564e95", "Juncus sp." = "#86a291",
                                "Sesuvium portulacastrum" = "#90E669", "Plantago maritima" = "#47724c",
                                "Scirpus mariqueter" = "#60563f", "Tamarix chinensis"="#FF46A2",
                                "Salicornia sp."="black", "Suaeda sp." = "black")) +
  scale_x_continuous(limits = c(-max(plot_data$ci_upper), max(plot_data$ci_upper) + 0.2)) +
  labs(
    x = "Standardized Mean Difference (SMD)",
    y = NULL,
    title = "Phragmites australis"
  ) +
  theme(
    legend.text =element_text(size = 12, colour = "black", family = "Helvetica"),
    axis.text.x = element_text(size = 12, colour = "black", family = "Helvetica"), 
    axis.text.y = element_text(size = 12, colour = "black", family = "Helvetica"),
    axis.title.x = element_text(size = 12, color = "black", family = "Helvetica"),
    axis.title.y = element_text(size = 12, color = "black", family = "Helvetica")) +
  theme_linedraw()+
  theme(legend.position="none",
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 11))+
          theme(plot.background = element_rect(fill = "transparent", colour = NA),
                panel.background = element_rect(fill = "transparent", colour = NA),
                legend.background = element_rect(fill = "transparent", colour = NA))
p3


# by salinity -------------------------------------------------------------
data_ch4 <- data_ch4 %>%
  mutate(repeat_flag = ifelse("repeat" == "Y", TRUE, FALSE)) %>%
  mutate(group_key = ifelse(!repeat_flag,
                            paste(p_num, lat, estuary, marsh, sampling_year, month, plant_species, salinity, specific_site, section, sep = "_"),
                            NA)) %>%
  # Create numeric unique IDs for each group_key
  mutate(unique_id = ifelse(!repeat_flag,
                            as.integer(factor(group_key)),
                            NA)) %>%
  relocate(unique_id, .after = p_num)


ch4_meta_data <- data_ch4 %>%
  group_by(unique_id) %>%
  summarise(
    # Calculate mean of yi
    ch4_flux = mean(ch4_flux, na.rm = TRUE),
    
    # Calculate sum of vi / n^2 (n = number of rows per group)
    ch4_n = sum(ch4_n),
    
    ch4_var=mean(ch4_var),
    
    # Example aggregation for other columns
    month = as.factor(paste(unique(month), collapse = ", ")),
    p_num = as.numeric(paste(unique(p_num), collapse = ", ")),
    plant_species = as.factor(paste(unique(plant_species), collapse = ", ")),
    salinity = as.factor(paste(unique(salinity), collapse = ", ")),
    season = as.factor(paste(unique(season), collapse = ", ")),
    tide = as.factor(paste(unique(tide), collapse = ", ")),
    lat = paste(unique(lat), collapse = ", "),  
    tide = as.factor(paste(unique(tide), collapse = ", ")),
    ch4_method_simple = as.factor(paste(unique(ch4_method_simple), collapse = ", ")),
    climate_region = as.factor(paste(unique(climate_region), collapse = ", ")),
    month = as.factor(paste(unique(month), collapse = ", ")),
    sampling_year = as.numeric(paste(unique(sampling_year), collapse = ", "))
    
    # Add other columns as needed here...
  )

meta_data <- ch4_meta_data %>%
  dplyr::select(ch4_flux, ch4_var, ch4_n, sampling_year, p_num, lat, month, ch4_method_simple, season, plant_species, salinity, tide)

#separate out Phragmites rows
phrag <- meta_data %>%
  filter(plant_species == "Phragmites australis") %>%
  dplyr::select(p_num, month, sampling_year,
                ch4_var_t = ch4_var, 
                ch4_n_t = ch4_n, 
                ch4_flux_t = ch4_flux, 
  )

#filter non-Phragmites rows
non_phrag <- meta_data %>%
  filter(plant_species != "Phragmites australis")

#join non-Phragmites with matching Phragmites rows by p_num and season
paired_data <- non_phrag %>%
  left_join(phrag, by = c("p_num", "month", "sampling_year"), relationship = "many-to-many") %>%
  filter(!is.na(ch4_flux_t))  # keep only pairs where Spartina data exists

ch4_meta_data <- escalc(measure = "SMD", 
                        m1i = ch4_flux_t, 
                        sd1i = ch4_var_t, 
                        n1i = ch4_n_t,
                        m2i = ch4_flux, 
                        sd2i = ch4_var, 
                        n2i = ch4_n, 
                        data = paired_data)

study_counts <- ch4_meta_data %>%
  group_by(salinity) %>%
  summarise(n_studies = n_distinct(p_num), .groups = "drop")


sal_aggregated <- ch4_meta_data %>%
  group_by(salinity) %>%
  summarise(
    avg_yi = mean(yi, na.rm = TRUE),
    avg_vi = sum(vi, na.rm = TRUE) / (n()^2),
    ch4_n = sum(ch4_n),
    p_num = paste(unique(p_num), collapse = ", "),
    tide = paste(unique(tide), collapse = ", "),
    ch4_method_simple = paste(unique(ch4_method_simple), collapse = ", ")
  ) %>%
  left_join(study_counts, by = "salinity") %>%  # ← This is key!
  mutate(slab_label = paste0(salinity, " | S=", n_studies, " | n=", ch4_n))

sal_aggregated <- sal_aggregated %>%
  mutate(salinity = factor(salinity, 
                           levels = c("oligohaline", "mesohaline", "polyhaline", "euhaline"), 
                           ordered = TRUE))

print(sal_aggregated, n=1000)

meta_model <- rma(avg_yi, avg_vi, data = sal_aggregated, method = "REML", knha=TRUE)
summary(meta_model)

levels(sal_aggregated$salinity)

sal_aggregated <- sal_aggregated %>%
  mutate(salinity = factor(salinity, 
                           levels = c("oligohaline", "mesohaline", "polyhaline", "euhaline"), 
                           ordered = TRUE))

# Extract model results
plot_data <- sal_aggregated %>%
  mutate(
    yi = meta_model$yi,
    vi = meta_model$vi,
    se = sqrt(vi),
    ci_lower = yi - 1.96 * se,
    ci_upper = yi + 1.96 * se
  )

plot_data <- plot_data %>%
  mutate(salinity = factor(salinity, 
                           levels = c("oligohaline", "mesohaline", "polyhaline", "euhaline"), 
                           ordered = TRUE))

plot_data <- plot_data %>%
  mutate(
    salinity = factor(salinity, 
                      levels = c("oligohaline", "mesohaline", "polyhaline", "euhaline")),
    slab_label = fct_reorder(slab_label, as.numeric(salinity))  # order by salinity
  )

p4 <- ggplot(plot_data, aes(x = yi, y = slab_label, color=salinity)) +
  geom_point(aes(color = salinity), size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color=salinity), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("oligohaline" = "lightblue", 
                                "mesohaline" = "#1F78B4",
                                "polyhaline" = "darkblue",
                                "euhaline" = "darkgreen")) +
  scale_x_continuous(limits = c(-max(plot_data$ci_upper), max(plot_data$ci_upper) + 0.2)) +
  labs(
    x = "Standardized Mean Difference (SMD)",
    y = NULL,
    title = "Phragmites australis"
  ) +
  theme(
    legend.text =element_text(size = 12, colour = "black", family = "Helvetica"),
    axis.text.x = element_text(size = 12, colour = "black", family = "Helvetica"), 
    axis.text.y = element_text(size = 12, colour = "black", family = "Helvetica"),
    axis.title.x = element_text(size = 12, color = "black", family = "Helvetica"),
    axis.title.y = element_text(size = 12, color = "black", family = "Helvetica")) +
  theme(plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA))+
  theme_linedraw()+
  theme(legend.position="none",
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 11)
  )

p4


library(patchwork)
combined_plot <- wrap_plots(p1 + p2+ p3+p4) + 
  plot_annotation(tag_levels = 'a') &
  theme(
    plot.tag = element_text(face = "bold"),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA)
  )
combined_plot

ggsave("figures/Figure_S6.png",width=12,height=8)

