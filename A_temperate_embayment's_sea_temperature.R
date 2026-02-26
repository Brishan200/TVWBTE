# Code for A temperate embayment's sea temperature
# Libraries & constants ----------------------------------------------------------
month_names <- c("January", "February", "March","April","May","June",
                 "July", "August","September","October","November","December")

month_names2 <- c("Jan", "Feb", "Mar","Apr","May","Jun",
                 "Jul", "Aug","Sep","Oct","Nov","Dec")

coloursforgraphs <- c("#8C613C","#0072B2", "#009E73", "#D55E00", "#9440DD", 
                      "#C51B8A","#56B4E9", "#FF988B", "#55183F", "darkgoldenrod",
                      "olivedrab","seagreen2","orchid2","violetred2","deepskyblue2" )



sitenames_abb <- c("Alexandria Dune Fields Inshore 30m" = "ADF30m",
                   "Algoa Bay Central 60m" = "ABC60m",
                   "Algoa Bay Mouth 80m" = "ABM80m",
                   "Bird Island Inshore 30m" = "BII30m",
                   "Bird Island Offshore 30m" = "BIO30m",
                   "Bird Island Offshore 80m" = "BIO80m",
                   "Cape Recife Inshore 30m" = "CRI30m",
                   "St Croix Island Inshore 30m" = "SCI30m",
                   "Sundays River Inshore 30m" = "SRI30m",
                   "Woody Cape Inshore 30m" = "WCI30m"
)


acronym_map <- c(
  "Alexandria Dune Fields Inshore 30m" = "ADF",
  "St Croix Island Inshore 30m" = "SCI",
  "Cape Recife Inshore 30m" = "CRI",
  "Sundays River Inshore 30m" = "SRI",
  "Woody Cape Inshore 30m" = "WCI",
  "Algoa Bay Central 60m" = "ABC",
  "Algoa Bay Mouth 80m" = "ABM",
  "Bird Island Inshore 30m" = "BII",
  "Bird Island Offshore 30m" = "BIO30",
  "Bird Island Offshore 80m" = "BIO80"
  
)

sites <- c("St Croix Island Inshore 30m",
           "Sundays River Inshore 30m",
           "Alexandria Dune Fields Inshore 30m",
           "Bird Island Inshore 30m",
           
           "Algoa Bay Central 60m",
           "Bird Island Offshore 30m",
           
           "Cape Recife Inshore 30m",
           "Algoa Bay Mouth 80m",
           "Bird Island Offshore 80m",
           "Woody Cape Inshore 30m")



site_coordinates <- data.frame(
  siteName = c("Alexandria Dune Fields Inshore 30m","Algoa Bay Central 60m",
               "Algoa Bay Mouth 80m","Bird Island Inshore 30m",
               "Bird Island Offshore 30m","Bird Island Offshore 80m",
               "Cape Recife Inshore 30m", "St Croix Island Inshore 30m", 
               "Sundays River Inshore 30m","Woody Cape Inshore 30m"),
  lat = c(-33.73579, -33.88195, 
          -33.94816, -33.81273, 
          -33.86849, -33.90385,
          -34.03411, -33.82751, 
          -33.76546,-33.75607),
  lon = c(26.06014, 25.98585, 
          26.11775, 26.31346, 
          26.29145, 26.30076,
          25.72771, 25.75143, 
          25.89722, 26.22909)
)

library(metR)
library(cowplot)
library(emmeans)
library(tidyverse)
library(scales)
library(reshape2)
library(gridExtra)
library(heatwaveR)
library(lubridate)
library(MBA)
library(mgcv)
library(lme4)
library(broom)
library(gratia)
library(MuMIn)
library(ggpubr)
library(sjPlot)
library(purrr)
library(htmlTable)
library(forecast)
library(grateful)
library(packrat)
library(funchir)
library(kableExtra)
library(cluster)
library(pheatmap)
library(factoextra)
library(TSclust)
library(ggdendro)
library(patchwork)
library(trend)


###


# Functions ---------------------------------------------------------------

total_days_per_month <- function(start_year, end_year) {
  # Generate sequence of dates
  dates <- seq(as.Date(paste0(start_year, "-01-01")), 
               as.Date(paste0(end_year, "-12-31")), 
               by = "day")
  
  # Extract month and year
  data <- data.frame(date = dates,
                     caldate_month = format(dates, "%m"),
                     year = format(dates, "%Y"))
  
  
  # Calculate cumulative sum of days
  data <- data %>%
    group_by(caldate_month) %>%
    summarise(cumulative_days = n())
  
  return(data)
}


total_days_per_year <- function(start_year, end_year) {
  # Generate sequence of dates
  dates <- seq(as.Date(paste0(start_year, "-01-01")), 
               as.Date(paste0(end_year, "-12-31")), 
               by = "day")
  
  # Extract month and year
  data <- data.frame(date = dates,
                     caldate_year = format(dates, "%Y"))
  
  
  # Calculate cumulative sum of days
  data <- data %>%
    group_by(caldate_year) %>%
    summarise(cumulative_days = n())
  
  return(data)
}

##

# data loading and cleaning -----------------------------------------------
##Download data from the http://googlie.shinyapps.io/SAEON_data_download
##and used to analysis the temperature of Algoa Bay from thermistor strings.

##
UTR_data <- read_csv("C:/Resources/Campus Work/PhD/Writing/data doi/ABSS PELTER 2008-2024 - ref st phyto/abss_utr_allstations.csv")
##

str(UTR_data)

UTR_data %>%
  mutate(caldate = as.Date(caldate)) %>%
  group_by(siteName) %>%
  summarise(
    start_date = min(caldate, na.rm = TRUE),
    end_date = max(caldate, na.rm = TRUE),
    .groups = "drop")

ggplot(data = UTR_data %>% mutate(caldate = as_date(caldate)) %>% 
         filter(siteName == "St Croix Island Inshore 30m")) +
  geom_point(aes(
    x = caldate,
    y = temp, 
    color = depth)) +
  labs(x = "",
       y = "") +
  scale_x_date(date_labels = "%Y", breaks = "year") +
  facet_wrap(~siteName,
             labeller = labeller(siteName = function(value) {
               wrapped_name <- str_wrap(value, width = 20)  # Adjust width as needed
               return(wrapped_name)
             })) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey70"),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "transparent")) 

unique(UTR_data$caldate_day)
unique(UTR_data$timeutc)
unique(UTR_data$hour)
unique(UTR_data$caldate_year)
###
# Satellite data  ---------------------------------------------------------
sst_data <- read.table("data/nasa_sst_point_series.txt",sep = ",", header = TRUE)

sst_data <- sst_data %>% left_join(site_coordinates, by = join_by(lon,lat))

sst_data <- sst_data %>% 
  mutate(caldate = as_date(time_start)) %>%
  filter(pixel_value > 5) %>% 
  select(!c("file","lat","lon","time_start","time_end",
            "pixel_count","valid","invalid","min","max","mean","median","std")) %>% 
  rename("d_avg_temp" ="pixel_value",
         "depth" = "variable") %>% 
  mutate(caldate_year = str_split_i(caldate,"-",1),
         caldate_month = str_split_i(caldate,"-",2),
         caldate_day = str_split_i(caldate,"-",3)) 

ggplot(data = sst_data) +
  geom_point(aes(
    x = caldate,
    y = d_avg_temp)) +
  labs(x = "",
       y = "") +
  scale_x_date(date_labels = "%Y", breaks = "year") +
  facet_wrap(~siteName,
             labeller = labeller(siteName = function(value) {
               wrapped_name <- str_wrap(value, width = 20)  # Adjust width as needed
               return(wrapped_name)
             })) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey70"),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "transparent")) 

##
# Diurnal variation  ------------------------------------------------------
diurnaldata_with_zero_range <- UTR_data %>% 
  group_by(siteName,depth,caldate) %>% 
  summarise(count = n(),
            diurnal_range = max(temp)-min(temp)) %>% 
  ungroup() %>% 
  filter(diurnal_range == 0, 
         count != 1) #Exclude records with only 1 temperature reading in a 24 hour cycle

# Filter UTR_data to include only dates, sites, and depths that have diurnal range == 0
matching_data <- UTR_data %>%
  inner_join(diurnaldata_with_zero_range, by = c("siteName", "depth", "caldate"))

# showing the records that have a 0 diurnal temeprature range
# Plot temperature over time (24-hours) for each site and depth
ggplot(matching_data %>% 
         mutate(date = as_date(date))) +
  geom_point(aes(x = date, y = temp, color = siteName,shape = depth),
             position = position_jitter()) +
  labs(title = "24-hour Temperature Profile for Dates with Zero Diurnal Range",
       x = "Time (UTC)",
       y = "Temperature (°C)") +
  scale_x_date(date_breaks = "year") +
theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###Remove records from UTR data set that have 1 record for a day in a 24h cycle
UTR_data %>% 
  group_by(siteName,depth,caldate) %>% 
  summarise(count = n()) %>% 
  filter(count == 1)
#Total of 20 records

UTR_data <- anti_join(UTR_data, UTR_data %>% 
                        group_by(siteName,depth,caldate) %>% 
                        summarise(count = n()) %>% 
                        filter(count == 1),
                      by = c("siteName", "depth", "caldate"))
str(UTR_data)
#Calculate the diurnal range for the dataset
diurnaldata <- UTR_data %>% 
  group_by(siteName,depth,caldate) %>% 
  summarise(numrecords = n(),
            diurnal_range = max(temp)-min(temp)) %>% 
  ungroup() 

#Exploring diurnal data
min(diurnaldata$diurnal_range)
max(diurnaldata$diurnal_range)
median(diurnaldata$diurnal_range)
mean(diurnaldata$diurnal_range)
sd(diurnaldata$diurnal_range)
range(diurnaldata$caldate)

hist(diurnaldata$diurnal_range)
quantile(diurnaldata$diurnal_range, probs = 0.95)
quantile(diurnaldata$diurnal_range, probs = 0.95)+sd(diurnaldata$diurnal_range) 

# save(diurnaldata, file = "C:/Code/AB-UTR-Temp/SAEON-UTR-static-db/diurnaldata.RData")


diurnaldata %>% 
  filter(caldate == "2014-01-31",
         siteName == "Cape Recife Inshore 30m",
         depth == "20m")

UTR_data %>% 
  filter(caldate == "2014-01-31" | caldate == "2014-02-01" | caldate == "2014-01-30",
         siteName == "Cape Recife Inshore 30m",
         depth == "20m")


ggplot(data = diurnaldata) +
  geom_violin(aes(x = depth, 
                  y = diurnal_range,
                  fill = depth),
              scale = "count",
              # trim = FALSE,
              # draw_quantiles = c(0.10, 0.90)
  ) +
  scale_fill_brewer(palette = "Spectral") +
  labs(title = "Daily diurnal temperature range",
       x = "",
       y = "") +
  facet_wrap(~siteName,
             labeller = labeller(siteName = function(value) {
               # Remove numeric values and units from siteName
               #formatted_name <- str_remove_all(value, "\\s*\\d+m$")
               # Wrap the facet label to ensure it displays fully
               wrapped_name <- str_wrap(value, width = 16)  # Adjust width as needed
               return(wrapped_name)
             }), scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey70"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "transparent")) 


##
#  daily diurnal range temp Plot
plots <- list()

# Loop through each site and create individual plots
for (site in sites) {
  site_data <- diurnaldata %>% filter(siteName == site)
  
  p <- ggplot(site_data) + #change between hourly data and the diurnal range data
    geom_tile(aes(x = as.Date(caldate),
                  y = fct_relevel(depth,rev),
                  fill = diurnal_range,
                  group = depth)) +
    scale_fill_viridis_c(name = "Temp.(°C)",
                         option = "plasma",
                         direction = -1,
                         aesthetics = "fill",
                         breaks =  seq(0,14, by = 2),
                         trans = "reverse") +
    scale_x_date(date_labels = "%Y",
                 date_breaks = "year") +
    labs(title = "",
         x = "",
         y = "") +
    facet_wrap(~siteName, scale = "free_y",
               labeller = labeller(siteName = function(value) {
                 # Remove numeric values and units from siteName
                 #formatted_name <- str_remove_all(value, "\\s*\\d+m$")
                 # Wrap the facet label to ensure it displays fully
                 wrapped_name <- str_wrap(value, width = 15)  # Adjust width as needed
                 return(wrapped_name)
               })) +
    theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.2),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "grey70"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(fill = "transparent"),
          legend.position = "none",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          strip.text = element_text(size = 12),
          text = element_text(family = "Arial", colour = "black"),
          strip.background = element_rect(fill = "transparent"))
  
  plots[[site]] <- p  # Store the plot in the list
}

# Arrange all plots in a grid
grid_plot <- grid.arrange(grobs = plots, ncol = 4, 
                          layout_matrix = rbind(
                            c(1, 2, 3, 4),  
                            c(NA, 5, 6, NA),
                            c(8, 9, 10, 11) 
                          ))



legend <- cowplot::get_legend(plots[[1]] + theme(legend.key = element_rect(fill = "transparent"),
                                                 legend.position="right"))

# Add legend to the combined plot
cplot <- cowplot::plot_grid(grid_plot, legend, ncol = 2, rel_widths = c(4, 0.5))

cplot


ggsave(filename ="plotsPaper/diurnalrange_temp.png", width = 18, height = 14, dpi = 600)


##
##

# --- tweak these ---
gap_days <- 2L        # allowed gap between days to still be same event
# -------------------

# 1) Tag each record with q95 + sd (per siteName x depth) and assign range_category
diurnal_tagged <- diurnaldata %>%
  mutate(caldate = as_date(caldate)) %>%
  group_by(siteName, depth) %>%
  mutate(
    q95    = quantile(diurnal_range, probs = 0.95, na.rm = TRUE),
    sd_all = sd(diurnal_range, na.rm = TRUE),
    range_category = case_when(
      is.na(diurnal_range) | is.na(q95) ~ NA_character_,
      diurnal_range < q95 ~ "Low",
      diurnal_range >= q95 & diurnal_range < (q95 + sd_all) ~ "Medium",
      diurnal_range >= (q95 + sd_all) ~ "High",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(range_category))

# 2) Build site-day-category event ids using ONLY the time-gap rule (no tol_abs/tol_rel)
site_day_events <- diurnal_tagged %>%
  group_by(siteName, caldate, range_category) %>%
  summarise(
    dr_med   = median(diurnal_range, na.rm = TRUE),
    n_depths = n_distinct(depth),
    .groups = "drop"
  ) %>%
  arrange(siteName, range_category, caldate) %>%
  group_by(siteName, range_category) %>%
  mutate(
    date_gap  = as.integer(caldate - lag(caldate)),
    new_group = if_else(is.na(date_gap) | date_gap > gap_days, 1L, 0L),
    event_id  = cumsum(new_group)
  ) %>%
  ungroup()

# 3) Attach the site-level event_id back to each depth record
diurnal_events <- diurnal_tagged %>%
  left_join(
    site_day_events %>% select(siteName, caldate, range_category, event_id),
    by = c("siteName", "caldate", "range_category")
  ) %>%
  arrange(siteName, depth, range_category, caldate)

# 4) Event counts per depth (unique events) + wide by category
event_counts_per_depth <- diurnal_events %>%
  distinct(siteName, depth, range_category, event_id) %>%
  count(siteName, depth, range_category, name = "n_events") %>%
  pivot_wider(
    names_from = range_category,
    values_from = n_events,
    values_fill = 0
  ) %>%
  arrange(siteName, depth)

print(event_counts_per_depth, n = Inf)

# 5) Summarise each event run per category per depth
diurnal_runs_summary <- diurnal_events %>%
  group_by(siteName, depth, range_category, event_id) %>%
  summarise(
    start_date        = min(caldate, na.rm = TRUE),
    end_date          = max(caldate, na.rm = TRUE),
    ave_diurnal_range = mean(diurnal_range, na.rm = TRUE),
    n_days            = n_distinct(caldate),
    total_records     = sum(numrecords, na.rm = TRUE),
    q95               = first(q95),
    sd_all            = first(sd_all),
    .groups = "drop"
  ) %>%
  arrange(siteName, depth, range_category, start_date)

diurnal_runs_summary
##
# Create the summary table of diurnal range 
# indicating low med high events using a 95th %tile and the sd

quantile(diurnaldata$diurnal_range, probs = 0.95) 
quantile(diurnaldata$diurnal_range, probs = 0.95) + sd(diurnaldata$diurnal_range)

diurnal_summary <- 
  # diurnaldata %>%
  # mutate(
  #   range_category = case_when(
  #     diurnal_range < quantile(diurnal_range, probs = 0.95) ~ "Low",
  #     diurnal_range >= quantile(diurnal_range, probs = 0.95) & 
  #       diurnal_range < quantile(diurnal_range, probs = 0.95) + sd(diurnal_range) ~ "Medium",
  #     diurnal_range >= quantile(diurnal_range, probs = 0.95) + sd(diurnal_range) ~ "High"),
  #   caldate = as_date(caldate)) %>% 
  diurnal_runs_summary %>% 
  group_by(siteName, start_date, depth, range_category) %>%
  summarise(num_records = n(), .groups = "drop") %>%
  pivot_wider(names_from = range_category, values_from = num_records, values_fill = 0) 

diurnal_summary %>% 
  mutate(year = format(start_date, "%Y")) %>%
  group_by(siteName, year, depth) %>%
  summarise(
    Low = sum(Low),
    Medium = sum(Medium),
    High = sum(High),
    .groups = "drop"
  ) %>%
  arrange(siteName, year, depth) %>%
  htmlTable(
    header = c("Site Name", "Year", "Depth", "Low", "Medium", "High"),
    align = "lccc",
    ctable = TRUE,
    rowlabel = "Depth"
  )

diurnal_summary %>% 
  select(siteName,start_date,depth,Medium,High) %>% 
  mutate(depth = as.character(depth),
         year = year(start_date)) %>%
  group_by(siteName,year,depth) %>%
  reframe(Medium = sum(Medium),
          High = sum(High)) %>% 
  mutate(summ = paste0(Medium,"; ", High)) %>% 
  pivot_wider(
    names_from = depth, 
    values_from = summ, 
    values_fill = "NA"
  ) %>% 
  select(siteName,year,"10m", "15m", "20m", "30m", "40m", "50m", "60m", "70m") %>% 
  drop_na() %>% 
  htmlTable()

diurnaldata %>% filter(diurnal_range == 0)

diurnal_summary %>% 
  pivot_longer(cols = c("Low", "Medium", "High"), 
               names_to = "range", 
               values_to = "Count") %>%
  mutate(year = as.character(year(start_date)),
         range = factor(range, levels = c("Low","Medium","High"))) %>% 
  filter(range != "Low") %>% 
  group_by(range) %>%
  # group_by(siteName) %>%
  summarize(total_count = sum(Count, na.rm = TRUE), .groups = 'drop')

ggplot(diurnaldata, aes(x = diurnal_range)) +
  geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black") +
  labs(title = "", x = "Diurnal Range", y = "Num. occurrences") +
  geom_vline(xintercept = quantile(diurnaldata$diurnal_range, 0.95), 
             color = "brown", linetype = "dashed", size = 1.2) +
  geom_vline(xintercept = quantile(diurnaldata$diurnal_range, 0.95) + sd(diurnaldata$diurnal_range, na.rm = TRUE), 
             color = "black", linetype = "dashed", size = 1.2) +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    # strip.text.x = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.position = "right",
    axis.title = element_text(size = 12),
    legend.key = element_rect(fill = "transparent"),
    panel.grid.major.x = element_line(color = "grey75"),
    panel.grid.major.y = element_line(color = "grey75"), 
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(size = 12,colour = "black"),
    strip.placement = "outside",
    text = element_text(family = "Arial", colour = "black")) 

ggsave(filename = "plotsPaper/hist_diurnalrange.jpg", width = 16, height = 9, dpi = 600)

range_labels <- c("Medium" = "a", "High" = "b")

ggplot(diurnal_summary %>% 
         rename("caldate" = "start_date") %>% 
         pivot_longer(cols = c("Low", "Medium", "High"), 
                      names_to = "range", 
                      values_to = "Count") %>%
         mutate(year = as.character(year(caldate)),
                range = factor(range, levels = c("Low","Medium","High"))) %>% 
         filter(range == "Medium") %>% 
         group_by(siteName, year, depth, range) %>%
         summarize(total_count = sum(Count, na.rm = TRUE), .groups = 'drop')) +
  geom_bar(stat = "identity", position = "stack",
           aes(x = year ,
               y = total_count,
               fill = depth,
               group = range)) +
  facet_grid(range ~ siteName,
             labeller = labeller(siteName = function(value) {
               wrapped_name <- str_wrap(value, width = 16)  # Adjust width as needed
               return(wrapped_name)
             })) +
  scale_fill_manual(values = coloursforgraphs) +
  labs(
    title = "",
    x = "",
    y = "Number of Records",
    fill = "")  +
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    # strip.text.x = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.position = "right",
    axis.title = element_text(size = 12),
    legend.key = element_rect(fill = "transparent"),
    panel.grid.major.x = element_line(color = "grey75"),
    panel.grid.major.y = element_line(color = "grey75"), 
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(size = 12,colour = "black"),
    strip.placement = "outside")+
  ggtitle("a") -> plt_diurnal_medium

ggplot(diurnal_summary %>% 
         rename("caldate" = "start_date") %>% 
         pivot_longer(cols = c("Low", "Medium", "High"), 
                      names_to = "range", 
                      values_to = "Count") %>%
         mutate(year = as.character(year(caldate)),
                range = factor(range, levels = c("Low","Medium","High"))) %>% 
         filter(range == "High") %>% 
         group_by(siteName, year, depth, range) %>%
         summarize(total_count = sum(Count, na.rm = TRUE), .groups = 'drop')) +
  geom_bar(stat = "identity", position = "stack",
           aes(x = year ,
               y = total_count,
               fill = depth,
               group = range)) +
  facet_grid(range ~ siteName,
             labeller = labeller(siteName = function(value) {
               wrapped_name <- str_wrap(value, width = 16)  # Adjust width as needed
               return(wrapped_name)
             })) +
  scale_fill_manual(values = coloursforgraphs) +
  labs(
    title = "",
    x = "",
    y = "Number of Records",
    fill = "")  +
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    # strip.text.x = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.position = "right",
    axis.title = element_text(size = 12),
    legend.key = element_rect(fill = "transparent"),
    panel.grid.major.x = element_line(color = "grey75"),
    panel.grid.major.y = element_line(color = "grey75"), 
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(size = 12,colour = "black"),
    strip.placement = "outside") +
  ggtitle("b") -> plt_diurnal_high

plt_diurnal_medium / plt_diurnal_high +
  plot_layout(guides = "collect", axis_titles = "collect",axes = "collect") &
  theme(legend.position = "bottom") -> diurnal_event_plt

ggsave(diurnal_event_plt,
       filename = "plotsPaper/DiurnalHighMedChanges_year.jpg", width = 16, height = 9, dpi = 600)


# Daily averaged temperature ----------------------------------------------

#Averaging the raw UTR hourly data in daily average temperatures

utr_daily_avg_temp <- UTR_data %>% 
  select(!c("station","sensor","date",
            "latitude","longitude","elevation")) %>% 
  mutate(caldate = as_date(caldate)) %>% 
  group_by(siteName,depth,caldate) %>% 
  summarise(d_avg_temp = mean(temp),
            d_sd_temp = sd(temp),
            d_cv = (d_sd_temp/d_avg_temp)*100) %>% 
  ungroup() %>% 
  mutate(caldate_year = str_split_i(caldate,"-",1),
         caldate_month = str_split_i(caldate,"-",2),
         caldate_day = str_split_i(caldate,"-",3))


utr_daily_avg_temp %>% 
  group_by(siteName,depth) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  summarise(ave = mean(count))

##Combining utr and sst data

daily_avg_temp <-  full_join(utr_daily_avg_temp,sst_data)

daily_avg_temp <- daily_avg_temp %>% 
  mutate(season = as.factor(case_when(
    month(caldate) %in% c(9, 10, 11) ~ "spring",
    month(caldate) %in% c(12, 1, 2) ~ "summer",
    month(caldate) %in% c(3, 4, 5) ~ "autumn",
    month(caldate) %in% c(6, 7, 8) ~ "winter"
  )),
  sector = as.factor(case_when(
    siteName %in% c("Alexandria Dune Fields Inshore 30m", "Bird Island Inshore 30m", "Woody Cape Inshore 30m") ~ "eastern",
    siteName %in% c("Algoa Bay Central 60m", "Algoa Bay Mouth 80m", "Bird Island Offshore 30m","Bird Island Offshore 80m") ~ "offshore",
    siteName %in% c("Cape Recife Inshore 30m", "St Croix Island Inshore 30m", "Sundays River Inshore 30m") ~ "western",
  )),
  depth = factor(depth, labels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m"),
                 levels = c("sst","10m","15m","20m","30m","40m","50m","60m","70m")))


# save(daily_avg_temp, file = "C:/Code/AB-UTR-Temp/SAEON-UTR-static-db/daily_avg_temp.RData")


temp_summary <- daily_avg_temp %>%
  group_by(depth, season) %>%
  summarise(
    mean_temp = round(mean(d_avg_temp, na.rm = TRUE), 2),
    min_temp  = round(min(d_avg_temp, na.rm = TRUE), 2),
    max_temp  = round(max(d_avg_temp, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(
    value = paste0(mean_temp, " (", min_temp, "–", max_temp, ")")
  ) %>%
  select(depth, season, value) %>% 
  mutate(season = factor(season, levels = c("spring", "summer", "autumn", "winter"))) %>% 
  arrange(depth, season) %>%         
  pivot_wider(
    names_from = season,
    values_from = value)

temp_summary %>%
  kable("html", escape = FALSE, align = "c") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))

ggplot(data = daily_avg_temp) +
  geom_line(aes(x = caldate, 
                y = d_avg_temp , 
                color = d_avg_temp )) +
  scale_color_viridis_c(option = "B") +
  facet_grid(depth~siteName, scales = "free")



utr_daily_avg_temp %>% 
  # filter(depth != "SST") %>% 
  group_by(siteName) %>% 
  reframe(mind = min(caldate), 
          maxd = max(caldate)) %>% 
  kable("html", escape = FALSE, align = "c") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))


##

# monthly climatology (average temperate) ----------------------------------------------------
monthly_clim <- daily_avg_temp %>% 
  group_by(siteName,depth, caldate_month) %>% 
  summarise(m_avgclim_temp = mean(d_avg_temp),
            m_sdclim_temp = sd(d_avg_temp),
            m_cv_clim_temp = ((m_sdclim_temp/m_avgclim_temp)*100),
            m_10clim_perc = quantile(d_avg_temp, probs = 0.10),
            m_90clim_perc = quantile(d_avg_temp, probs = 0.90)) %>% 
  ungroup() %>% 
  mutate(caldate = paste("2000-",caldate_month,"-15", sep = ""),
         depth = factor(depth, labels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m"),
                        levels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m")),
         sector = as.factor(case_when(
           siteName %in% c("Alexandria Dune Fields Inshore 30m", "Bird Island Inshore 30m", "Woody Cape Inshore 30m") ~ "eastern",
           siteName %in% c("Algoa Bay Central 60m", "Algoa Bay Mouth 80m", "Bird Island Offshore 30m","Bird Island Offshore 80m") ~ "offshore",
           siteName %in% c("Cape Recife Inshore 30m", "St Croix Island Inshore 30m", "Sundays River Inshore 30m") ~ "western",
         )),
         season = as.factor(case_when(
           month(caldate) %in% c(9, 10, 11) ~ "spring",
           month(caldate) %in% c(12, 1, 2) ~ "summer",
           month(caldate) %in% c(3, 4, 5) ~ "autumn",
           month(caldate) %in% c(6, 7, 8) ~ "winter"
         ))) 

# save(monthly_clim, file = "C:/Code/AB-UTR-Temp/SAEON-UTR-static-db/monthly_clim.RData")

monthly_clim %>% 
  group_by(siteName, depth) %>%
  summarize(
    month_max = month.name[as.integer(caldate_month[which.max(m_avgclim_temp)])],
    month_min = month.name[as.integer(caldate_month[which.min(m_avgclim_temp)])],
    "max_temp(°C)" = round(max(m_avgclim_temp, na.rm = TRUE),2),
    "min_temp(°C)" = round(min(m_avgclim_temp, na.rm = TRUE),2)
  ) %>%
  ungroup() %>% 
  htmlTable(rownames = FALSE)

##Get missing data from specified years

# Step 1: Calculate the total number of expected records for each siteName, depth, and caldate_month
expected_records <- total_days_per_month(2008,2024)

# Step 2: Calculate the total number of actual records for each siteName, depth, and caldate_month
actual_records <- daily_avg_temp %>%
  filter(caldate_year %in% 2008:2024) %>%
  group_by(siteName, depth, caldate_month) %>%
  summarise(actual_rec = n()) %>%
  ungroup()

# Step 3: Calculate the percentage of missing data
missing_percentage <- expected_records %>%
  left_join(actual_records, by = c("caldate_month")) %>%
  mutate(percentage_missing = (100-(actual_rec / cumulative_days) * 100))

# summary(missing_percentage)

monthly_clim <- monthly_clim %>%
  left_join(missing_percentage, by = c("siteName", "depth", "caldate_month"))

monthly_clim <- monthly_clim %>% 
  mutate(depth = factor(depth, labels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m"),
                        levels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m")))
##
# Plots
#Missing data plot 

ggplot(data = monthly_clim %>% filter(depth != "SST") %>% 
         mutate(siteName = factor(siteName, levels = sites))) +
  geom_tile(aes(x = caldate_month,
                y = fct_relevel(depth, rev),
                fill = percentage_missing ,
                group = depth)) +
  scale_fill_viridis_c(name = "% missing",
                       option = "viridis",
                       direction = 1,
                       aesthetics = "fill",
                       breaks =  seq(0,100, by = 10),
                       trans = "reverse") +
  scale_x_discrete( labels= month_names2) +
  # scale_x_date(date_labels = "%B",
  #              limits = c(as.Date("2000-01-01"), as.Date("2000-12-31")),
  #              breaks = seq(as.Date("2000-01-01"), as.Date("2000-12-31"), by = "1 month")) +
  labs(title = "",
       x = "",
       y = "") +
  facet_wrap(~siteName,
             labeller = labeller(siteName = function(value) {
               #if (grepl("emp", value)) return("")  # Remove label for empty spaces
               wrapped_name <- str_wrap(value, width = 20)  
               return(wrapped_name)
             }),
             scales = "free_y") +
  theme(axis.text.x = element_text(size = 12, hjust = 1, angle = 90, vjust = 0.2),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey70"),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        text = element_text(family = "Arial", colour = "black"),
        strip.background = element_rect(fill = "transparent"))


ggsave(filename ="plotsPaper/percmissing-monthly.png",
       width = 16, height = 9, dpi = 600)


#Temperature

##Interpolation 
plots <- list()

daily_avg_plot <- monthly_clim %>%  
  mutate(depth = factor(depth, labels = c("sst","10","15","20","30","40","50","60","70"),
                        levels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m"))) %>% 
  rename(date = caldate,
         temp = m_avgclim_temp) %>% 
  na.omit() %>%
  select(siteName, date, depth, temp)

daily_avg_plot$depth <- as.numeric(gsub("sst", "0", as.character(daily_avg_plot$depth)))
daily_avg_plot$date <- as.Date(daily_avg_plot$date)
daily_avg_plot$date <- decimal_date(daily_avg_plot$date)

interpolate_and_plot <- function(site_name) {
  # Filter data for the current site
  daily_avg_plot_site <- daily_avg_plot %>% 
    filter(siteName == site_name) %>% 
    select(date, depth, temp)
  
  # Interpolate the data
  ctd_mba <- mba.surf(daily_avg_plot_site, no.X = 300, no.Y = 300, extend = TRUE)
  dimnames(ctd_mba$xyz.est$z) <- list(ctd_mba$xyz.est$x, ctd_mba$xyz.est$y)
  ctd_mba <- melt(ctd_mba$xyz.est$z, varnames = c('date', 'depth'), value.name = 'temp') %>% 
    mutate(temp = round(temp, 1))
  
  # Plot the interpolated data
  ggplot(data = ctd_mba,
         aes(x = date_decimal(date), y = depth)) +
    geom_raster(aes(fill = temp)) +
    scale_fill_viridis_c(option = "A", limits = c(9,25))  +
    scale_y_reverse() +
    scale_x_datetime(breaks = "1 month", date_labels = "%b") +
    geom_contour2(aes(z = round(temp,0)), binwidth = 1,
                  colour = "black", alpha = 1) +
    metR::geom_text_contour(aes(z = round(temp,0)),stroke = 0.2, size = 4,
                            check_overlap = TRUE,
                            label.placer = label_placer_n(n = 1,
                                                          rot_adjuster = isoband::angle_halfcircle_bottom())) +
    coord_cartesian(expand = 0) +
    labs(x = "",
         y = "",
         fill = "Temp. (°C)",
         title = site_name) +
    theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.2),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "grey70"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(fill = "transparent"),
          legend.position = "none",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          strip.text = element_text(size = 12),
          text = element_text(family = "Arial", colour = "black"),
          strip.background = element_rect(fill = "transparent"))
}


# Create and save plots for each site
plots <- lapply(sites, interpolate_and_plot)

# plots[[1]]
grid_plot <- grid.arrange(grobs = plots, ncol = 4, 
                          layout_matrix = rbind(
                            c(1, 2, 3, 4),  
                            c(NA, 5, 6, NA),
                            c(8, 9, 10, 11) 
                          ))

legend <- cowplot::get_legend(plots[[1]] + theme(legend.key = element_rect(fill = "transparent"),
                                                 legend.position="right"))

# Add legend to the combined plot
cplot <- cowplot::plot_grid(grid_plot, legend, ncol = 2, rel_widths = c(4, 0.5))

cplot

ggsave(cplot, filename ="plotsPaper/Interpolated-Monthly-climatology.png",
       width = 16, height = 9, dpi = 600)


#Coefficient of variation 

##Set plot layout

plots <- list()

# Loop through each site and create individual plots
for (site in sites) {
  site_data <- monthly_clim %>% filter(siteName == site) %>%
    mutate(depth = factor(depth, levels = rev(unique(depth))))  # Adjust depth order if needed
  
  p <- ggplot(data = site_data) +
    geom_tile(aes(x = caldate_month,
                  y = depth,
                  fill = m_cv_clim_temp)) +
    scale_fill_viridis_c(name = "CV",
                         option = "plasma",
                         direction = -1,
                         breaks = seq(0, 21, by = 3),
                         trans = "reverse") +
    scale_x_discrete( labels= month_names2) +
    labs(title = site,  # Set the title as the site name
         x = "",
         y = "") +
    theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
          axis.text.y = element_text(size = 12),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "grey70"),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          legend.key = element_rect(fill = "transparent"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          text = element_text(family = "Arial", colour = "black"))
  
  plots[[site]] <- p  # Store the plot in the list
}

# Arrange all plots in a grid
grid_plot <- grid.arrange(grobs = plots, ncol = 4, 
                          layout_matrix = rbind(
                            c(1, 2, 3, 4),  
                            c(NA, 5, 6, NA),
                            c(8, 9, 10, 11) 
                          ))


legend <- cowplot::get_legend(plots[[1]] + theme(legend.key = element_rect(fill = "transparent"),
                                                 legend.position="right"))

# Add legend to the combined plot
cplot <- cowplot::plot_grid(grid_plot, legend, ncol = 2, rel_widths = c(4, 0.5))

cplot


ggsave(filename ="plotsPaper/Monthly-climatology-cv.png", width = 16, height = 9, dpi = 600)

##
# yearly climatology ------------------------------------------------------

yearly_clim <-  daily_avg_temp %>% 
  #filter(caldate_year %in% (2010:2023)) %>%  # for convience sake 
  #can include season to get yearly seasonal climatology 
  group_by(siteName,depth, caldate_year) %>% 
  summarise(y_avgclim_temp = median(d_avg_temp),
            y_sdclim_temp = sd(d_avg_temp),
            y_10clim_perc = quantile(d_avg_temp, probs = 0.10),
            y_90clim_perc = quantile(d_avg_temp, probs = 0.90),
            y_cv = (y_sdclim_temp / y_avgclim_temp) * 100) %>% 
  ungroup() %>% 
  mutate(depth = factor(depth, labels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m"),
                        levels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m")),
         sector = as.factor(case_when(
           siteName %in% c("Alexandria Dune Fields Inshore 30m", "Bird Island Inshore 30m", "Woody Cape Inshore 30m") ~ "eastern",
           siteName %in% c("Algoa Bay Central 60m", "Algoa Bay Mouth 80m", "Bird Island Offshore 30m","Bird Island Offshore 80m") ~ "offshore",
           siteName %in% c("Cape Recife Inshore 30m", "St Croix Island Inshore 30m", "Sundays River Inshore 30m") ~ "western",
         )))

 # save(yearly_clim, file = "C:/Code/AB-UTR-Temp/SAEON-UTR-static-db/yearly_clim.RData")

yearly_clim %>% 
  filter(caldate_year %in% (2008:2024),
         depth != "SST") %>%
  group_by(siteName) %>% 
  summarise(t2022_temp = round(mean(y_avgclim_temp[caldate_year %in% c(2008:2022)]),2),
            t2023_temp = round(mean(y_avgclim_temp[caldate_year == "2023"]),2)) %>%
  mutate(diff = round(t2023_temp- t2022_temp,2)) %>%
  ungroup() %>% 
  select(siteName, diff) %>% 
  htmlTable(rownames = FALSE)
# filter(percentage_missing < 20) %>% 
# filter(y_avgclim_temp %in% c(min(y_avgclim_temp), max(y_avgclim_temp)))

yearly_clim %>% 
  group_by(siteName) %>% 
  reframe(round(min(y_avgclim_temp),2),
          round(max(y_avgclim_temp),2)) %>% 
  htmlTable(rownames = FALSE)


##Get missing data 
# Step 1: Calculate the total number of expected records for each siteName, depth, and caldate_month
expected_records <- total_days_per_year(2000,2024)

# Step 2: Calculate the total number of actual records for each siteName, depth, and caldate_month
actual_records <- daily_avg_temp %>%
  # filter(caldate_year %in% 2010:2023) %>%
  group_by(siteName, depth, caldate_year) %>%
  summarise(actual_rec = n()) %>%
  ungroup()

# Step 3: Calculate the percentage of missing data
missing_percentage <- expected_records %>%
  left_join(actual_records, by = c("caldate_year")) %>%
  mutate(percentage_missing = (100-(actual_rec / cumulative_days) * 100))


yearly_clim <- yearly_clim %>%
  left_join(missing_percentage, by = c("siteName", "depth", "caldate_year"))

yearly_clim <- yearly_clim %>% 
  mutate(depth = factor(depth, labels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m"),
                        levels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m")))
#
# Plots
#Missing data plot 

ggplot(data = yearly_clim %>% filter(caldate_year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", 
                                                         "2015", "2016", "2017", "2018", "2019", "2020", "2021", 
                                                         "2022", "2023", "2024")) %>% 
         mutate(siteName = factor(siteName, levels = sites))) +
  geom_tile(aes(x = caldate_year,
                y = fct_relevel(depth, rev),
                fill = percentage_missing ,
                group = depth)) +
  scale_fill_viridis_c(name = "% missing",
                       option = "viridis",
                       direction = 1,
                       aesthetics = "fill",
                       trans = "reverse") +
  labs(title = "",
       x = "",
       y = "") +
  facet_wrap(~siteName,
             labeller = labeller(siteName = function(value) {
               #if (grepl("emp", value)) return("")  # Remove label for empty spaces
               wrapped_name <- str_wrap(value, width = 20)  
               return(wrapped_name)
             }),
             scales = "free_y") +
  theme(axis.text.x = element_text(size = 12, hjust = 1, angle = 90, vjust = 0.2),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey70"),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        text = element_text(family = "Arial", colour = "black"),
        strip.background = element_rect(fill = "transparent"))


ggsave(filename ="plotsPaper/percmissing-yearly.png",
       width = 16, height = 9, dpi = 600)



#Temperature

#Interpolation
plots <- list()

daily_avg_plot <- yearly_clim %>%  
  mutate(depth = factor(depth, labels = c("sst","10","15","20","30","40","50","60","70"),
                        levels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m")),
         caldate = paste(caldate_year,"-01","-01", sep = "")) %>% 
  rename(date = caldate,
         temp = y_avgclim_temp ) %>% 
  filter(caldate_year>2007) %>% 
  # na.omit() %>%
  select(siteName, date, depth, temp)

daily_avg_plot$depth <- as.numeric(gsub("sst", "0", as.character(daily_avg_plot$depth)))
daily_avg_plot$date <- as.Date(daily_avg_plot$date)
daily_avg_plot$date <- decimal_date(daily_avg_plot$date)

interpolate_and_plot <- function(site_name) {
  # Filter data for the current site
  daily_avg_plot_site <- daily_avg_plot %>% 
    filter(siteName == site_name) %>% 
    select(date, depth, temp)
  
  # Interpolate the data
  ctd_mba <- mba.surf(daily_avg_plot_site, no.X = 300, no.Y = 300, extend = TRUE)
  dimnames(ctd_mba$xyz.est$z) <- list(ctd_mba$xyz.est$x, ctd_mba$xyz.est$y)
  ctd_mba <- melt(ctd_mba$xyz.est$z, varnames = c('date', 'depth'), value.name = 'temp') %>% 
    mutate(temp = round(temp, 1))
  
  # Plot the interpolated data
  ggplot(data = ctd_mba,
         aes(x = date_decimal(date), y = depth)) +
    geom_raster(aes(fill = temp)) +
    scale_fill_viridis_c(option = "A", limits = c(7,25)) +
    scale_y_reverse() +
    scale_x_datetime(breaks = "year", date_labels = "%Y") +
    geom_contour2(aes(z = round(temp,0)), binwidth = 1,
                  colour = "black", alpha = 1) +
    metR::geom_text_contour(aes(z = round(temp,0)),stroke = 0.2, size = 4,
                            check_overlap = TRUE,
                            label.placer = label_placer_n(n = 1,
                                                          rot_adjuster = isoband::angle_halfcircle_bottom())) +
    coord_cartesian(expand = 0) +
    labs(x = "",
         y = "",
         fill = "Temp. (°C)",
         title = site_name) +
    theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.2),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "grey70"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(fill = "transparent"),
          legend.position = "none",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          strip.text = element_text(size = 12),
          text = element_text(family = "Arial", colour = "black"),
          strip.background = element_rect(fill = "transparent"))
}


# Create and save plots for each site
plots <- lapply(sites, interpolate_and_plot)

grid_plot <- grid.arrange(grobs = plots, ncol = 4, 
                          layout_matrix = rbind(
                            c(1, 2, 3, 4),  
                            c(NA, 5, 6, NA),
                            c(8, 9, 10, 11) 
                          ))

legend <- cowplot::get_legend(plots[[1]] + theme(legend.key = element_rect(fill = "transparent"),
                                                 legend.position="right"))

# Add legend to the combined plot
cplot <- cowplot::plot_grid(grid_plot, legend, ncol = 2, rel_widths = c(4, 0.5))

cplot

ggsave(cplot, filename ="plotsPaper/Interpolated-Yearly-climatology.png", 
       width = 16, height = 9, dpi = 600)

#
#Coefficient of variation 

#Set plot layout
plots <- list()

# Loop through each site and create individual plots
for (site in sites) {
  site_data <- yearly_clim %>% filter(siteName == site) %>%
    mutate(depth = factor(depth, levels = rev(unique(depth))))  # Adjust depth order if needed
  
  p <- ggplot(data = site_data %>%  
                filter(caldate_year > 2007,
                       y_cv < 25)) +
    geom_tile(aes(x = caldate_year ,
                  y = depth,
                  fill = y_cv ,
                  group = depth)) +
    scale_fill_viridis_c(name = "CV", #"Std.Dev.(σ)",
                         option = "plasma",
                         direction = -1,
                         aesthetics = "fill",
                         breaks =  seq(0,21, by = 3),
                         trans = "reverse") +
    labs(title = site,
         x = "",
         y = "") +
    theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.2),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "grey70"),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          legend.key = element_rect(fill = "transparent"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          text = element_text(family = "Arial", colour = "black"))
  
  plots[[site]] <- p  # Store the plot in the list
}

# Arrange all plots in a grid
grid_plot <- grid.arrange(grobs = plots, ncol = 4, 
                          layout_matrix = rbind(
                            c(1, 2, 3, 4),  
                            c(NA, 5, 6, NA),
                            c(8, 9, 10, 11) 
                          ))


legend <- cowplot::get_legend(plots[[1]] + theme(legend.key = element_rect(fill = "transparent"),
                                                 legend.position="right"))

# Add legend to the combined plot
cplot <- cowplot::plot_grid(grid_plot, legend, ncol = 2, rel_widths = c(4, 0.5))

cplot

ggsave(filename ="plotsPaper/Yearly-climatology-CV.png", width = 16, height = 9, dpi = 600)

##

# LMM on temperature ossiclations over time  ------------------------------
#### GLMM Model
##Create a daily climatology from the dataset to see changes of each year relative this dataset
daily_clim <-  daily_avg_temp %>% 
  group_by(siteName,depth, caldate_month, caldate_day) %>% 
  summarise(d_avgclim_temp = mean(d_avg_temp),
            d_medclim_temp = median(d_avg_temp),
            d_sdclim_temp = sd(d_avg_temp),
            d_cv_clim = (d_sdclim_temp/d_avgclim_temp)*100, 
            d_10clim_perc = quantile(d_avg_temp, probs = 0.10),
            d_90clim_perc = quantile(d_avg_temp, probs = 0.90),
            d_maxclim = max(d_avg_temp),
            d_minclim = min (d_avg_temp)) %>% 
  ungroup() %>% 
  mutate(caldate = as_date(paste("2000-",caldate_month,"-",caldate_day, sep="")),
         depth = factor(depth, labels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m"),
                        levels = c("sst","10m","15m","20m","30m","40m","50m","60m","70m")),
         season = case_when(
           caldate_month %in% c("09","10","11") ~ "Spring",
           caldate_month %in% c("12","01","02") ~ "Summer",
           caldate_month %in% c("03","04","05") ~ "Autumn",
           caldate_month %in% c("06","07","08") ~ "Winter",
           TRUE ~ "Unknown"),
         sector = as.factor(case_when(
           siteName %in% c("Alexandria Dune Fields Inshore 30m", "Bird Island Inshore 30m", "Woody Cape Inshore 30m") ~ "eastern",
           siteName %in% c("Algoa Bay Central 60m", "Algoa Bay Mouth 80m", "Bird Island Offshore 30m","Bird Island Offshore 80m") ~ "offshore",
           siteName %in% c("Cape Recife Inshore 30m", "St Croix Island Inshore 30m", "Sundays River Inshore 30m") ~ "western",
         ))) 


glmdata <-  daily_avg_temp %>% 
  mutate(d_cv = (d_sd_temp/d_avg_temp)*100,
         depth = factor(depth, labels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m"),
                        levels = c("sst","10m","15m","20m","30m","40m","50m","60m","70m"))) %>% 
  select(siteName,depth,caldate,d_avg_temp,d_sd_temp,d_cv) %>% 
  # filter(year(caldate) %in% c("2009":"2023")) %>%
  distinct() %>% 
  full_join(daily_clim %>% 
              select(siteName,depth,caldate,d_avgclim_temp ,d_sdclim_temp, d_cv_clim ) %>% 
              distinct(), 
            by = join_by(siteName, depth, caldate,
                         d_avg_temp == d_avgclim_temp,
                         d_sd_temp == d_sdclim_temp,
                         d_cv == d_cv_clim)) %>% 
  mutate(caldate_year = as.factor(str_split_i(caldate,"-",1)),
         season = as.factor(case_when(
           month(caldate) %in% c(9, 10, 11) ~ "spring",
           month(caldate) %in% c(12, 1, 2) ~ "summer",
           month(caldate) %in% c(3, 4, 5) ~ "autumn",
           month(caldate) %in% c(6, 7, 8) ~ "winter"
         )),
         siteName = as.factor(siteName))

model <- lmer(d_avg_temp  ~ caldate_year + (1|siteName) +(1|depth) + (1|season) + (1|caldate),
              data = glmdata)

r.squaredGLMM(model)

summary(model)

anova(model)
qqnorm(resid(model))
qqline(resid(model))
plot(predict(model), resid(model))
coef(model)
extractAIC(model)

hist(resid(model))
acf(resid(model))

tab_model(model, dv.labels = "")

plot_model(model, type = "est", show.values = TRUE, title = "",value.offset = 0.40,
           value.size = 6,line.size = 1.5,
           colors = c("steelblue", "firebrick"),
           axis.labels = c(2024,2023,2022,2021,2020,2019,2018,2017,2016,2015,2014,2013,2012,2011,2010,2009,2008)) +
  scale_y_continuous(limits = c(-1.1, 1)) + 
  theme(axis.text.x = element_text(size = 12, vjust = 0.1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90",size = 1),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        strip.text = element_text(size = 12),
        text = element_text(colour = "black", size = 12),
        strip.background = element_rect(fill = "transparent"))


ggsave(filename = "plotsPaper/Estimate-lme-plot-temperature-vs-time.png", width = 16, height = 9, dpi = 600)

###
# Cluster analysis --------------------------------------------------------
seasonal_data <- UTR_data %>%
  mutate(
    caldate = as.Date(caldate),
    month = month(caldate),
    season = case_when(
      month %in% c(12, 1, 2) ~ "Summer",
      month %in% c(3, 4, 5) ~ "Autumn",
      month %in% c(6, 7, 8) ~ "Winter",
      month %in% c(9, 10, 11) ~ "Spring"
    )
  ) %>% 
  group_by(siteName, depth, season) %>%
  summarize(
    mean_temp = mean(temp, na.rm = TRUE),
    sd_temp = sd(temp, na.rm = TRUE)
  ) %>% 
  ungroup()

seasonal_wide <- seasonal_data %>%
  group_by(siteName,depth) %>%
  pivot_wider(names_from = season, values_from = mean_temp) %>% 
  summarise(
    autumn = mean(Autumn, na.rm = TRUE),
    spring = mean(Spring, na.rm = TRUE),
    summer = mean(Summer, na.rm = TRUE),
    winter = mean(Winter, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  mutate(
    site_depth = str_replace_all(siteName, acronym_map),
    depth = paste(depth, "depth", sep = " ")
  )


temp_matrix <- as.matrix(seasonal_wide %>% select(-siteName, -depth, -site_depth))
rownames(temp_matrix) <- paste(seasonal_wide$site_depth, seasonal_wide$depth, sep = " - ")

row_clusters <- hclust(diss(temp_matrix, METHOD = "DTWARP"))
col_clusters <- hclust(diss(t(temp_matrix), METHOD = "DTWARP"))

dendrogram_plot <- function(temp_matrix, row_clusters, col_clusters) {
  par(mar = c(5, 5, 4, 8)) # Adjust margins for labels
  plot(row_clusters, horiz = TRUE, xlab = "Distance", main = "Dendrogram for Sites", cex.main=1.5)
  
  par(xpd = TRUE)
  plot(col_clusters, main = "Dendrogram for Seasons", xlab = "", ylim = c(0, max(col_clusters$height)), cex.main=1.5)
  
  # Adding heatmap
  heatmap(temp_matrix, Rowv = as.dendrogram(row_clusters), Colv = as.dendrogram(col_clusters),
          scale = "none", col = heat.colors(256), margins = c(5, 10), 
          main = "Seasonal Temperature Patterns by Site and Depth", xlab = "Seasons", ylab = "Sites and Depths")
}

dendrogram_plot(temp_matrix, row_clusters, col_clusters)


## Silouette plot
# Define a maximum number of clusters to check
max_k <- 10
sil_widths <- numeric(max_k)

# Distance matrix using DTW
distance_matrix <- diss(temp_matrix, METHOD = "DTWARP")

# Loop to calculate average silhouette width for each k
for (k in 2:max_k) {
  # Cut tree at k clusters
  clusters <- cutree(row_clusters, k = k)
  
  # Calculate silhouette width
  sil <- silhouette(clusters, distance_matrix)
  sil_widths[k] <- mean(sil[, 3])  # Average silhouette width for each k
}

sil_data <- data.frame(
  k = 1:max_k, 
  silhouette_width = sil_widths[1:max_k]
)

ggplot(sil_data, aes(x = k, y = silhouette_width)) +
  geom_line() +
  geom_point(shape = 19) +
  labs(
    x = "Number of Clusters (k)",
    y = "Average Silhouette Width",
    title = "Optimal k based on Silhouette Width"
  ) +
  theme(axis.text.x = element_text(size = 12, vjust = 0.2),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey70"),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        text = element_text(family = "Arial", colour = "black"),
        strip.background = element_rect(fill = "transparent"))

ggsave(filename ="plotsPaper/silhouette_plot.png",
       width = 16, height = 9, dpi = 600)

#Heatmap & dendogram
pheatmap(temp_matrix, 
         cluster_rows = row_clusters,
         cluster_cols = col_clusters, 
         angle_col = 0,
         main = "",
         fontsize_number = 12,                     # Adjust font size for all text elements
         show_rownames = TRUE,              # Show row names
         show_colnames = TRUE,              # Show column names
         legend = TRUE, 
         fontsize = 12, 
         display_numbers = TRUE,
         number_color = "black", 
         border_color = "black",
         legend_breaks = c(12,14,16,18,20),
         cutree_rows = 7,
         cutree_cols = 2
) -> heatmap

ggsave(heatmap,filename ="plotsPaper/dissimilarityheatmap.png", width = 16, height = 9, dpi = 600)

# some validation
cophenetic_dist <- cophenetic(row_clusters)
original_dist <- dist(temp_matrix)
cophenetic_correlation <- cor(original_dist, cophenetic_dist)
print(cophenetic_correlation)

##Creating cluster list to do statistics on for H2 
seasonal_wide$cluster <-  cutree(row_clusters,k = 7)
clustered_list <- split(seasonal_wide %>% select(siteName, depth), seasonal_wide$cluster)



###
# statistical analysis on clusters ----------------------------------------------------
#Liner regression of temp over time 
str(daily_avg_temp)
str(clustered_list)
##using combinations from cluster analysis 

clustered_reodlis <- c(clustered_list[7], clustered_list[-7])
names(clustered_reodlis) <- as.character(seq_along(clustered_reodlis))

cluster_df <- map_df(
  seq_along(clustered_reodlis), 
  ~ clustered_reodlis[[.x]] %>%
    mutate(cluster = .x) %>%  # Assign cluster number
    mutate(depth = gsub(" depth", "", depth),
           cluster = as.factor(cluster))  # Remove " depth" for matching
)

get_season_start_date <- function(year, season) {
  switch(as.character(season),
         spring = as.Date(paste(year, "03-01", sep = "-")),   # March 1
         summer = as.Date(paste(year, "06-01", sep = "-")),   # June 1
         autumn = as.Date(paste(year, "09-01", sep = "-")),   # September 1
         winter = as.Date(paste(year, "12-01", sep = "-")))   # December 1
}

ggplot(daily_avg_temp %>% 
         filter(caldate_year > 2008, depth != "SST") %>%
         left_join(cluster_df, by = c("siteName", "depth")) %>% 
         group_by(cluster, caldate_year, season) %>%
         summarize( d_avg = mean(d_avg_temp, na.rm =TRUE),
                    d_avg_cv = mean(d_cv, na.rm = TRUE), .groups = "drop") %>%
         rowwise() %>% 
         mutate(year_season = get_season_start_date(caldate_year, season),  # Convert to date
                year_season = as.Date(year_season)) ,
       aes(x = year_season, y = d_avg)) +
  geom_line(aes(group = cluster, color = cluster)) +  # Line plot for temperature over time by depth
  geom_smooth(aes(group = cluster, color = cluster),
              method = "lm", color = "black") + 
  stat_regline_equation(aes(label = after_stat(eq.label)), 
                        label.x.npc = "left", label.y.npc = "top", size = 3) +  # Equation and p-value
  facet_wrap(cluster~.) +  
  labs(title = "",
       x = "",
       y = "CV") +
  scale_x_date(date_labels = "%Y", date_breaks = "2 year") +
  scale_color_manual(values = coloursforgraphs) +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.2),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey70"),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        strip.text = element_text(size = 12),
        text = element_text(family = "Arial", colour = "black"),
        strip.background = element_rect(fill = "transparent"))

summary_data <- daily_avg_temp %>% 
  filter(caldate_year > 2008, depth != "SST") %>%
  left_join(cluster_df, by = c("siteName", "depth")) %>% 
  group_by(cluster, caldate_year, season) %>%
  summarize( d_avg = mean(d_avg_temp, na.rm =TRUE),
             d_avg_cv = mean(d_cv, na.rm = TRUE), .groups = "drop") %>%
  rowwise() %>% 
  mutate(year_season = get_season_start_date(caldate_year, season),  # Convert to date
         year_season = as.Date(year_season))


# Fit linear models and extract coefficients and R-squared
model_summary <- summary_data %>%
  group_by( cluster ) %>%
  do({
    model <- lm(d_avg ~ year_season, data = .) #Change between d_avg and d_avg_cv
    tidied <- tidy(model)
    r_squared <- summary(model)$r.squared
    
    data.frame(
      cluster  = unique(.$cluster ),
      intercept = tidied$estimate[1],
      slope = tidied$estimate[2],
      r_squared = r_squared,
      equation = paste("y =",tidied$estimate[1], "+",tidied$estimate[2], "x")
    )
  })

print(model_summary)

model_summary <- summary_data %>%
  group_by( cluster ) %>%
  do({
    model <- lm(d_avg_cv  ~ year_season, data = .) #Change between d_avg and d_avg_cv
    tidied <- tidy(model)
    r_squared <- summary(model)$r.squared
    
    data.frame(
      cluster  = unique(.$cluster ),
      intercept = tidied$estimate[1],
      slope = tidied$estimate[2],
      r_squared = r_squared,
      equation = paste("y =",tidied$estimate[1], "+",tidied$estimate[2], "x")
    )
  })

print(model_summary)

# write.csv(model_summary, "plotsPaper/temp_lm_model_summary.csv", row.names = FALSE)
##
# heatwaveR ---------------------------------------------------------------
##Only using UTR data not satellite data
grouped_list <- daily_avg_temp %>%
  # filter(#caldate_year %in% c(2010:2023),
  #        depth != "SST") %>%
  group_by(siteName, depth) %>%
  group_split()

# Apply your operation to each group
ts_list <- lapply(grouped_list, function(group) {
  # Extract siteName, depth, start and end date for the group
  site_name <- unique(group$siteName)
  depth_val <- unique(group$depth)
  startdate <- first(group$caldate)
  enddate <- last(group$caldate)
  
  if (length(site_name) == 1 && length(depth_val) == 1) {
    # Convert the group data to clm format
    ts_clm <- ts2clm(group %>%
                       rename(t = caldate, temp = d_avg_temp),
                     climatologyPeriod = c(startdate, enddate)#, pctile = 10,
    )
    
    # Attach siteName and depth to the resulting clm object
    attr(ts_clm, "siteName") <- site_name
    attr(ts_clm, "depth") <- depth_val
    
    return(ts_clm)
  } else {
    # Handle the case where there are multiple siteName or depth values in a group
    warning("Group contains multiple siteName or depth values. Skipping.")
    return(NULL)
  }
})

# Remove NULL elements
ts_list <- ts_list[!sapply(ts_list, is.null)]

# Set names for the list based on siteName and depth
ts_list <- setNames(ts_list, sapply(ts_list, function(ts_clm) {
  paste(attr(ts_clm, "siteName"), attr(ts_clm, "depth"), sep = "_")
}))

mhw_list <- map(ts_list, detect_event)

### IMPORTANT !!!!!!!!
##Before making mcs_list ensure that you select pctile = 10 for 
##ts2clm function above

mcs_list <- map(ts_list, ~detect_event(.x, coldSpells = TRUE))


#Heatwaves
hwR_c <- tibble()
hwR_e <- tibble()

for (i in seq_along(mhw_list)) {
  sNd <- names(mhw_list[i])
  
  hwR_clim <- map_dfr(mhw_list[i], ~tibble(
    clim  = .x[[1]]
  ))
  hwR_event <- map_dfr(mhw_list[i], ~tibble(
    event  = .x[[2]]
  ))
  
  hwR_clim <- hwR_clim %>% mutate(siteName = sNd)
  hwR_event <- hwR_event %>% mutate(siteName = sNd)
  
  hwR_c <- bind_rows(hwR_c, hwR_clim)  
  hwR_e <- bind_rows(hwR_e, hwR_event)
}

hwR_c <- hwR_c %>% unnest(cols = clim) 
hwR_e <- hwR_e %>% unnest(cols = event) 

hwR_c <- hwR_c %>% mutate(depth = sub('.*_', '\\1', siteName),
                          siteName = sub('_.*', '',siteName),
                          sector = as.factor(case_when(
                            siteName %in% c("Bird Island Inshore 30m","Woody Cape Inshore 30m", "Alexandria Dune Fields Inshore 30m" ) ~ "eastern",
                            siteName %in% c("Bird Island Offshore 30m","Bird Island Offshore 80m" ,"Algoa Bay Mouth 80m","Algoa Bay Central 60m") ~ "offshore",
                            siteName %in% c("Sundays River Inshore 30m", "St Croix Island Inshore 30m","Cape Recife Inshore 30m") ~ "western",
                          )),
                          c_year = year(t),
                          c_month = month(t))

hwR_e <- hwR_e %>% mutate(depth = sub('.*_', '\\1', siteName),
                          siteName = sub('_.*', '',siteName),
                          sector = as.factor(case_when(
                            siteName %in% c("Bird Island Inshore 30m","Woody Cape Inshore 30m", "Alexandria Dune Fields Inshore 30m" ) ~ "eastern",
                            siteName %in% c("Bird Island Offshore 30m","Bird Island Offshore 80m" ,"Algoa Bay Mouth 80m","Algoa Bay Central 60m") ~ "offshore",
                            siteName %in% c("Sundays River Inshore 30m", "St Croix Island Inshore 30m","Cape Recife Inshore 30m") ~ "western",
                          )),
                          e_year = year(date_start),
                          e_month = month(date_start))
hwR_c <- hwR_c %>% 
  mutate(depth = factor(depth, labels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m"),
                        levels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m")))

hwR_e <- hwR_e %>% 
  mutate(depth = factor(depth, labels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m"),
                        levels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m")))

# write_csv(hwR_c, "hwR_c.csv")
# write_csv(hwR_e, "hwR_e.csv")

###

### cold spells   ###

###
mcs_c <- tibble()
mcs_e <- tibble()

#ensure that the list was created with 10th percentile data
for (i in seq_along(mcs_list)) {
  sNd <- names(mcs_list[i])
  
  mcs_clim <- map_dfr(mcs_list[i], ~tibble(
    clim  = .x[[1]]
  ))
  mcs_event <- map_dfr(mcs_list[i], ~tibble(
    event  = .x[[2]]
  ))
  
  mcs_clim <- mcs_clim %>% mutate(siteName = sNd)
  mcs_event <- mcs_event %>% mutate(siteName = sNd)
  
  mcs_c <- bind_rows(mcs_c, mcs_clim)  
  mcs_e <- bind_rows(mcs_e, mcs_event)
}

mcs_c <- mcs_c %>% unnest(cols = clim) 
mcs_e <- mcs_e %>% unnest(cols = event) 

mcs_c <- mcs_c %>% mutate(depth = sub('.*_', '\\1', siteName),
                          siteName = sub('_.*', '',siteName),
                          sector = as.factor(case_when(
                            siteName %in% c("Alexandria Dune Fields Inshore 30m", "Bird Island Inshore 30m", "Woody Cape Inshore 30m") ~ "eastern",
                            siteName %in% c("Algoa Bay Central 60m", "Algoa Bay Mouth 80m", "Bird Island Offshore 30m","Bird Island Offshore 80m") ~ "offshore",
                            siteName %in% c("Cape Recife Inshore 30m", "St Croix Island Inshore 30m", "Sundays River Inshore 30m") ~ "western",
                          )),
                          c_year = year(t),
                          c_month = month(t))

mcs_e <- mcs_e %>% mutate(depth = sub('.*_', '\\1', siteName),
                          siteName = sub('_.*', '',siteName),
                          sector = as.factor(case_when(
                            siteName %in% c("Alexandria Dune Fields Inshore 30m", "Bird Island Inshore 30m", "Woody Cape Inshore 30m") ~ "eastern",
                            siteName %in% c("Algoa Bay Central 60m", "Algoa Bay Mouth 80m", "Bird Island Offshore 30m","Bird Island Offshore 80m") ~ "offshore",
                            siteName %in% c("Cape Recife Inshore 30m", "St Croix Island Inshore 30m", "Sundays River Inshore 30m") ~ "western",
                          )),
                          e_year = year(date_start),
                          e_month = month(date_start))

mcs_c <- mcs_c %>% 
  mutate(depth = factor(depth, labels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m"),
                        levels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m")))

mcs_e <- mcs_e %>% 
  mutate(depth = factor(depth, labels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m"),
                        levels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m")))


# write_csv(mcs_c, "mcs_c.csv")
# write_csv(mcs_e, "mcs_e.csv")

##
### 
####   HEatwaves and coldspells of algoa bay ###
###
##
hwR_e %>% 
  filter(year(date_start) > "2008") %>% 
  left_join(cluster_df) %>% 
  group_by(date_start, cluster) %>% 
  reframe(count = n()) %>% 
  filter(count > 0) -> testheatwavesab

hwe_ab <- hwR_e %>% 
  filter(date_start %in% testheatwavesab$date_start) %>% 
  filter(!is.na(date_start), !is.na(date_end)) %>%
  arrange(date_start) %>%
  mutate(
    # cumulative max of date_end (handle Date via numeric)
    prev_max_end = lag(
      as.Date(cummax(as.numeric(date_end)), origin = "1970-01-01")
    ),
    
    # start a new event only if the gap is > 5 days
    new_event = if_else(
      is.na(prev_max_end) | date_start > (prev_max_end + 5),
      1L, 0L
    ),
    
    event_group = cumsum(new_event),
    event_name  = paste0("e", event_group)
  ) %>%
  ungroup() %>%
  select(-prev_max_end, -new_event, -event_group) %>%
  mutate(event_name = factor(event_name, levels = unique(event_name)))



hwe_ab %>%
  group_by(event_name) %>%
  summarise(
    start = min(date_start),
    end   = max(date_end),
    n_events = n()
  )

hwe_ab %>% 
  select(event_name,date_start,date_end,duration,intensity_mean) %>% 
  group_by(event_name) %>% 
  reframe(event_name = event_name, 
          startdate = min(date_start),
          enddate = max(date_end),
          dur = enddate-startdate,
          insmean = round(mean(intensity_mean),2)) %>% 
  distinct() -> hw_ab_sumary

hw_ab_sumary %>% 
  reframe(durmean = mean(dur))

hw_ab_sumary %>% 
  reframe(durmax = max(dur))

hw_ab_sumary %>% 
  reframe(inmean = mean(insmean))

hw_ab_sumary %>%
  group_by(year(startdate)) %>% 
  reframe(count = n_distinct(event_name))


htmlTable(
  hw_ab_sumary,
  header = c(
    "Event",
    "Start date",
    "End date",
    "Duration (days)",
    "Ave. Intensity"
  ),
  rnames = FALSE,
  caption = "Summary of MHW events (merged events)",
  css.cell = "padding: 6px;",
  css.table = "width: 80%;"
)


#gets data and plots a profile of the heatwaves with 10 days beofre and after each MHW
hw_profile_events <- map_df(
  seq_len(nrow(hw_ab_sumary)),
  ~ {
    ev <- hw_ab_sumary[.x, ]
    
    hwR_c %>%
      filter(t >= ev$startdate - days(10), t <= ev$enddate + days(10)) %>%
      mutate(event_name  = ev$event_name,
             date_start = ev$startdate,
             date_end   = ev$enddate,
             duration   = ev$dur)
  }) %>% 
  mutate(depth_num = as.numeric(fct_rev(depth)))

temp_summary <- hw_profile_events %>%
  summarise(temp = mean(temp, na.rm = TRUE),
            .by = c(event_name, t, depth, depth_num))

event_levels <- levels(temp_summary$event_name)

for (ev in event_levels) {
  
  ev_dates <- hw_profile_events %>%
    filter(event_name == ev) %>%
    summarise(
      start_date = min(date_start),
      end_date   = max(date_end),
      .groups = "drop"
    )
  
  
  p <- temp_summary %>%
    filter(event_name == ev) %>%
    ggplot(aes(x = t, y = depth_num)) +
    geom_tile(aes(fill = temp)) +
    geom_contour(aes(z = temp), colour = "white") +
    geom_vline(
      xintercept = as.numeric(ev_dates$start_date),
      colour = "white",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    geom_vline(
      xintercept = as.numeric(ev_dates$end_date),
      colour = "white",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    scale_y_continuous(
      breaks = sort(unique(hw_profile_events$depth_num), decreasing = TRUE),
      labels = levels(hw_profile_events$depth)
    ) +
    scale_x_date(
      date_labels = "%d %b %n %Y",
      date_breaks = "7 days"
    ) +
    scale_fill_viridis_c(option = "H") +
    labs(x = "", y = "Depth (m)") +
    theme(
      axis.text.x = element_text(size = 12, vjust = 0.2),
      axis.text.y = element_text(size = 12),
      axis.title  = element_text(size = 12),
      panel.background = element_rect(fill = "gray90"),
      panel.grid.major = element_line(color = "grey70"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(fill = "transparent"),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.title = element_blank(),
      strip.text = element_text(size = 12),
      text = element_text(family = "Arial", colour = "black")
    )
  
  ggsave(
    filename = paste0("plotsPaper/MHWplots/MHW_", ev, ".png"),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
}

# save(hwe_ab, file = "C:/Code/AB-UTR-Temp/SAEON-UTR-static-db/hwe_ab.RData")
# save(hw_profile_events, file = "C:/Code/AB-UTR-Temp/SAEON-UTR-static-db/hwc_ab.RData")
# save(hwR_c, file = "C:/Code/AB-UTR-Temp/SAEON-UTR-static-db/hwc_full_ab.RData")

     ##             ###
     ## MHW classes ###
     ##             ###

#MHW MALAN Classes

depth_lookup <- tibble(
  depth = factor(c("SST","10m","15m","20m","30m","40m","50m","60m","70m"),
                 levels = c("SST","10m","15m","20m","30m","40m","50m","60m","70m")),
  depth_m = c(0,10,15,20,30,40,50,60,70)
)

hw_profiles <- hw_profile_events %>% 
  left_join(depth_lookup, by = "depth") %>%
  mutate(temp_anom = temp - thresh)

MLD_daily <- hw_profiles %>%
  group_by(siteName, t) %>%
  arrange(depth_m) %>%
  mutate(SST = dplyr::first(temp[depth_m == 0 & !is.na(temp)])) %>%
  filter(!is.na(SST), depth_m > 0) %>%
  summarise(
    MLD_m = {
      valid <- !is.na(temp)
      if (any(abs(temp[valid] - SST) >= 0.2, na.rm = TRUE)) {
        min(depth_m[abs(temp - SST) >= 0.2 & valid], na.rm = TRUE)
      } else if (any(valid)) {
        max(depth_m[valid], na.rm = TRUE)
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  )

thermo_daily <- hw_profiles %>%
  # filter(depth_m > 0) %>%
  arrange(siteName, t, depth_m) %>%
  group_by(siteName, t) %>%
  mutate(
    dTdz = (temp - dplyr::lag(temp)) /
      (depth_m - dplyr::lag(depth_m))
  ) %>%
  filter(!is.na(dTdz)) %>%
  summarise(
    thermo_m = depth_m[which.max(abs(dTdz))],
    thermo_strength = max(abs(dTdz), na.rm = TRUE),
    .groups = "drop"
  )

event_profiles <- hw_profiles %>% 
  filter(t >= date_start, t <= date_end) %>%
  group_by(event_name, depth_m) %>%
  summarise(
    temp_mean = mean(temp, na.rm = TRUE),
    anom_mean = mean(temp_anom, na.rm = TRUE),
    .groups = "drop"
  )

MLD_event <- MLD_daily %>%
  inner_join(hwe_ab, by = c("siteName")) %>%
  filter(t >= date_start, t <= date_end) %>%
  group_by(event_name) %>%
  summarise(
    MLD_m = median(MLD_m, na.rm = TRUE),
    .groups = "drop"
  )

thermo_event <- thermo_daily %>%
  inner_join(hwe_ab, by = c("siteName")) %>%
  filter(t >= date_start, t <= date_end) %>%
  group_by(event_name) %>%
  summarise(
    thermo_m = median(thermo_m, na.rm = TRUE),
    thermo_strength = median(thermo_strength, na.rm = TRUE),
    .groups = "drop"
  )

event_metrics <- event_profiles %>%
  # filter(depth_m > 0) %>%
  group_by(event_name) %>%
  summarise(
    max_anom = max(anom_mean, na.rm = TRUE),
    max_anom_depth = depth_m[which.max(anom_mean)],
    mean_anom = mean(anom_mean, na.rm = TRUE), # average anomaly value per event 
    surface_anom = anom_mean[depth_m == min(depth_m)],
    .groups = "drop"
  ) %>%
  left_join(MLD_event, by = "event_name") %>%
  left_join(thermo_event, by = "event_name")


event_classification <- event_metrics %>%
  mutate(
    MHW_class = case_when(
      
      max_anom_depth <= MLD_m ~
        "Mixed-layer MHW",
      
      max_anom_depth == thermo_m &
        max_anom > surface_anom ~
        "Thermocline MHW",
      
      surface_anom >= max_anom &
        max_anom_depth >= MLD_m &
        max_anom_depth >= thermo_m ~ 
        "Deep MHW",
      
      max_anom_depth == max(depth_lookup$depth_m) &
        surface_anom >= max_anom~
        "Full-depth MHW",
      
      max_anom_depth >= max(depth_lookup$depth_m) &
        max_anom > surface_anom ~
        "Benthic MHW",
      
      max_anom_depth >= thermo_m &
        surface_anom < max_anom ~
        "Submerged MHW",
      
      TRUE ~ "Unclassified"
    )
  )

event_classification %>%
  select(event_name, MHW_class)

unique(event_classification$MHW_class)

event_classification %>% 
  group_by(MHW_class) %>% 
  reframe(counts= n())


htmlTable(
  hw_ab_sumary %>% 
    left_join(event_classification %>% select(event_name, MHW_class)),
  header = c(
    "Event",
    "Start date",
    "End date",
    "Duration (days)",
    "Ave. Intensity (°C)",
    "Category"
  ),
  rnames = FALSE,
  caption = "Summary of MHW events (merged events)",
  css.cell = "padding: 6px;",
  css.table = "width: 80%;"
)
###
###
#                 ###
#    coldspells   ###
##                ###
mcs_e %>% 
  filter(year(date_start) > "2008") %>% 
  left_join(cluster_df) %>% 
  group_by(date_start, cluster) %>% 
  reframe(count = n()) %>% 
  filter(count > 0) -> testcoldspellsab

mcs_ab <- mcs_e %>% 
  filter(date_start %in% testcoldspellsab$date_start) %>% 
  filter(!is.na(date_start), !is.na(date_end)) %>%
  arrange(date_start) %>%
  mutate(
    # cumulative max of date_end (handle Date via numeric)
    prev_max_end = lag(
      as.Date(cummax(as.numeric(date_end)), origin = "1970-01-01")
    ),
    
    # start a new event only if the gap is > 5 days
    new_event = if_else(
      is.na(prev_max_end) | date_start > (prev_max_end + 5),
      1L, 0L
    ),
    
    event_group = cumsum(new_event),
    event_name  = paste0("e", event_group)
  ) %>%
  ungroup() %>%
  select(-prev_max_end, -new_event, -event_group) %>%
  mutate(event_name = factor(event_name, levels = unique(event_name)))


mcs_ab %>% 
  select(event_name,date_start,date_end,duration) %>% 
  group_by(event_name) %>% 
  reframe(event_name = event_name, 
          startdate = min(date_start),
          enddate = max(date_end),
          dur = enddate-startdate) %>% 
  distinct() -> mc_ab_summary


mc_ab_summary %>% 
  reframe(durmean = mean(dur))

mc_ab_summary %>% 
  reframe(durmax = max(dur))


mc_ab_summary %>%
  group_by(year(startdate)) %>% 
  reframe(count = n_distinct(event_name))


htmlTable(
  mc_ab_summary,
  header = c(
    "Event",
    "Start date",
    "End date",
    "Duration (days)"
  ),
  rnames = FALSE,
  caption = "Summary of cold-spell events (merged events)",
  css.cell = "padding: 6px;",
  css.table = "width: 80%;"
)



mc_profile_events <- map_df(
  seq_len(nrow(mcs_ab)),
  ~ {
    ev <- mcs_ab[.x, ]
    
    mcs_c %>%
      filter(t >= ev$date_start - days(2), t <= ev$date_end + days(2)) %>%
      mutate(event_no   = ev$event_no,
             event_name  = ev$event_name,
             date_start = ev$date_start,
             date_end   = ev$date_end,
             duration   = ev$duration)
  }
) %>% 
  mutate(depth_num = as.numeric(fct_rev(depth)))

mctemp_summary <- mc_profile_events %>%
  summarise(temp = mean(temp, na.rm = TRUE),
            .by = c(event_name, t, depth, depth_num))

event_levels <- levels(mctemp_summary$event_name)

for (ev in event_levels) {
  
  ev_dates <- mc_profile_events %>%
    filter(event_name == ev) %>%
    summarise(
      start_date = min(date_start),
      end_date   = max(date_end),
      .groups = "drop"
    )
  
  
  p <- mctemp_summary %>%
    filter(event_name == ev) %>%
    ggplot(aes(x = t, y = depth_num)) +
    geom_tile(aes(fill = temp)) +
    geom_contour(aes(z = temp), colour = "white") +
    geom_vline(
      xintercept = as.numeric(ev_dates$start_date),
      colour = "white",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    geom_vline(
      xintercept = as.numeric(ev_dates$end_date),
      colour = "white",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    scale_y_continuous(
      breaks = sort(unique(hw_profile_events$depth_num), decreasing = TRUE),
      labels = levels(hw_profile_events$depth)
    ) +
    scale_x_date(
      date_labels = "%d %b %n %Y",
      date_breaks = "7 days"
    ) +
    scale_fill_viridis_c(option = "H") +
    labs(x = "", y = "Depth (m)") +
    theme(
      axis.text.x = element_text(size = 12, vjust = 0.2),
      axis.text.y = element_text(size = 12),
      axis.title  = element_text(size = 12),
      panel.background = element_rect(fill = "gray90"),
      panel.grid.major = element_line(color = "grey70"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(fill = "transparent"),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.title = element_blank(),
      strip.text = element_text(size = 12),
      text = element_text(family = "Arial", colour = "black")
    )
  
  ggsave(
    filename = paste0("plotsPaper/MCSplots/MCS_", ev, ".png"),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
}


# save(mcs_ab, file = "C:/Code/AB-UTR-Temp/SAEON-UTR-static-db/mcse_ab.RData")
# save(mc_profile_events, file = "C:/Code/AB-UTR-Temp/SAEON-UTR-static-db/mcsc_ab.RData")
# save(mcs_c, file = "C:/Code/AB-UTR-Temp/SAEON-UTR-static-db/mcsc_full_ab.RData")

##                                                        ##
## Freq of MHW and MCS identified in Algoa Bay over time  ##
##                                                        ##
# --- Cold Spells ---
cold_counts <- mc_ab_summary %>%
  mutate(type = "mcs",
         year  = year(startdate),
         month = month(startdate),
         month_name = month(startdate, label = TRUE, abbr = TRUE)) %>%
  count(year, month, month_name, type, name = "n_events") %>%
  arrange(year, month)

# --- Heatwaves ---
heat_counts <- hw_ab_sumary %>%
  mutate(type = "mhw",
         year  = year(startdate),
         month = month(startdate),
         month_name = month(startdate, label = TRUE, abbr = TRUE)) %>%
  count(year, month, month_name, type,  name= "n_events") %>%
  arrange(year, month)

# --- Combine both and ensure order by year, month ---
event_counts <- bind_rows(cold_counts, heat_counts) %>%
  arrange(type, year, month)

##
#yearly not month and year 

event_counts_year <- event_counts %>%
  group_by(year, type) %>%
  summarise(
    n_events = sum(n_events, na.rm = TRUE),
    .groups = "drop")


ggplot(event_counts_year,
       aes(x = year,
           y = n_events,
           color = type,
           group = type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal(base_size = 14) +
  labs(title = "Yearly Frequency of Marine Heatwave & Cold Spell Events",
       x = "Year",
       y = "Number of Events",
       color = "Event Type") +
  scale_color_manual(values = c("mcs" = "steelblue",
                                "mhw" = "firebrick")) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank())



event_full <- event_counts_year %>%
  group_by(type) %>%
  complete(year = seq(min(year, na.rm = TRUE),
                      max(year, na.rm = TRUE),
                      by = 1),
           fill = list(n_events = 0)) %>%
  ungroup() %>%
  arrange(type, year) %>%
  group_by(type) %>%
  mutate(time_index = row_number()) %>%
  ungroup()

trend_lines <- event_full %>%
  group_by(type) %>%
  summarise(
    sen    = coef(lm(n_events ~ time_index))[2],
    sen_p  = summary(lm(n_events ~ time_index))$coefficients[2, 4],
    tau    = mk.test(n_events)$estimates["tau"],
    mk_p   = mk.test(n_events)$p.value,
    sen5 = sen*5,
    .groups = "drop"
  )

trend_lines

event_counts_plot <- event_full %>%
  left_join(trend_lines, by = "type") %>%
  group_by(type) %>%
  mutate(
    fitted = (sen * time_index) +
      (mean(n_events, na.rm = TRUE) -
         sen * mean(time_index, na.rm = TRUE))
  ) %>%
  ungroup()


ggplot(event_counts_plot,
       aes(x = year,
           y = n_events,
           color = type,
           group = type)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_line(aes(y = fitted),
            linewidth = 1.2) +
  scale_color_manual(values = c("mcs" = "steelblue",
                                "mhw" = "firebrick")) +
  theme_minimal(base_size = 14) +
  labs(title = "Yearly Frequency of Marine Heatwave & Cold Spell Events\nwith Linear Trend Lines",
       x = "Year",
       y = "Number of Events",
       color = "Event Type") +
  theme(legend.position = "top",
        panel.grid.minor = element_blank())




##
###Citation-----

library(NCmisc)

list.functions.in.file("A_temperate_embayment's_sea_temperature.R")

renv::dependencies("A_temperate_embayment's_sea_temperature.R")

get_citations("lubridate")
annotater()
cite_packages(out.dir = ".",
              pkgs = "All" ,
              omit= NULL,  include.RStudio = TRUE)

#




###individual MHW and MCS events per site / depth -----
## Plots of individual events for each site/ depth

hwR_e %>% 
  # filter(depth != "SST") %>% 
  group_by(depth) %>% 
  reframe(count = n(),
          range = range(intensity_mean_abs)) %>% 
  group_by(depth) %>% 
  reframe(min = min(range),
          max = max (range))

hwR_e %>% 
  filter(depth != "SST") %>%
  summarise(min = min(duration),
            mean = mean(duration),
            max = max(duration))

hwR_e %>% 
  filter(depth != "SST") %>% 
  group_by(e_year) %>% 
  reframe(count = n()) %>% 
  print(n = 112)


hwR_e %>% 
  filter(depth != "SST") %>% 
  group_by(e_year, siteName, depth) %>% 
  reframe(count = n()) %>% 
  group_by(depth) %>% 
  reframe(highest =max(count))




mcs_e %>% 
  #filter(depth != "SST") %>% 
  group_by(depth) %>% 
  reframe(count = n(),
          range = range(intensity_mean_abs)) %>% 
  group_by(depth) %>% 
  reframe(min = min(range),
          max = max (range))

mcs_e %>% 
  filter(depth != "SST") %>%
  summarise(min = min(duration),
            mean = mean(duration),
            max = max(duration))

mcs_e %>% 
  filter(depth != "SST") %>% 
  group_by(e_year) %>% 
  reframe(count = n()) %>% 
  print(n = 112)


mcs_e %>% 
  filter(depth != "SST") %>% 
  group_by(e_year, siteName, depth) %>% 
  reframe(count = n()) %>% 
  group_by(depth) %>% 
  reframe(highest =max(count))

## Number of events per year, site  and depth 
# Heatwaves

ggplot(hwR_e %>%
         mutate(n_dep =  as.numeric(depth),
                e_year = as.character(e_year),
                e_month = as.character(e_month)) %>%
         group_by(e_year, siteName, depth) %>%
         summarise(total_events = n()) %>%
         filter(e_year > 2008),
       aes(x = e_year, y = total_events, color = siteName, group = siteName )) +
  # geom_bar(stat = "identity", position = "dodge2") +
  geom_jitter() +
  labs(x = "", y = "Frequency of warming events", title = "") +
  facet_wrap(~depth) +
  # stat_cor() +
  # scale_x_date(date_labels = "%b %Y", date_breaks = "6 months") +
  scale_color_manual(values = coloursforgraphs) +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.2),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey70"),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        strip.text = element_text(size = 12),
        text = element_text(family = "Arial", colour = "black"),
        strip.background = element_rect(fill = "transparent"))

ggsave(filename ="plotsPaper/FreqWE-per-year-site-depth.png", width = 16, height = 9, dpi = 600)

##
#Coldspells
##

ggplot(mcs_e %>%
         mutate(n_dep =  as.numeric(depth),
                e_year = as.character(e_year),
                e_month = as.character(e_month)) %>%
         group_by(e_year, siteName, depth) %>%
         summarise(total_events = n()) %>%
         filter(e_year > 2008),
       aes(x = e_year, y = total_events, color = siteName, group = siteName )) +
  # geom_bar(stat = "identity", position = "dodge2") +
  geom_jitter() +
  labs(x = "", y = "Frequency of cooling events", title = "") +
  facet_wrap(~depth) +
  # stat_cor() +
  # scale_x_date(date_labels = "%b %Y", date_breaks = "6 months") +
  scale_color_manual(values = coloursforgraphs) +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.2),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.background = element_rect(fill = "gray90"),
        panel.grid.major = element_line(color = "grey70"),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        strip.text = element_text(size = 12),
        text = element_text(family = "Arial", colour = "black"),
        strip.background = element_rect(fill = "transparent"))

ggsave(filename ="plotsPaper/FreqCE-per-year-site-depth.png", width = 16, height = 9, dpi = 600)
##