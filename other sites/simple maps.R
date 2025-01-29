install.packages("leaflet")
install.packages("readr")

# Load required packages
library(tidyverse)
library(sf)
library(ggplot2)
library(ggthemes)
library(leaflet)
library(readr)

# Read in the CSV file
data <- read.csv("metadata_v2.csv")

# Create a base map centered on an average location
# Create a base map centered on an average location
center_lat <- mean(data$latitude)
center_lng <- mean(data$longitude)

map <- leaflet() %>%
  addTiles() %>%
  setView(lng = center_lng, lat = center_lat, zoom = 5)

# Add points to the map
map <- map %>% addCircleMarkers(
  lng = ~longitude, 
  lat = ~latitude, 
  radius = 5, 
  color = "blue", 
  fill = TRUE, 
  fillOpacity = 0.6
)

# Display the map
map

# Optionally, save the map as an HTML file
library(htmlwidgets)
saveWidget(map, "sites_map.html", selfcontained = TRUE)



