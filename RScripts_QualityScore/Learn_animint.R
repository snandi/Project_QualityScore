library("animint")
library("plyr")
library("maps")
library("ggplot2")
library('servr')
data(UStornadoes, package = "animint")

USpolygons <- ggplot2::map_data("state")
USpolygons$state <- state.abb[match(USpolygons$region, tolower(state.name))]

map <- ggplot() + 
  geom_polygon(aes(x = long, y = lat, group = group), 
               data = USpolygons, fill = "black", colour = "grey") +
  geom_segment(aes(x = startLong, y = startLat, xend = endLong, yend = endLat, showSelected = year), 
               colour = "#55B1F7", data = UStornadoes) +
  ggtitle("Tornadoes in the US")

ts <- ggplot() + 
  stat_summary(aes(x=year, y=year, clickSelects = year), 
               data = UStornadoes, fun.y = length, geom = "bar") + 
  ggtitle("Number of Recorded Tornadoes, 1950-2006") + 
  ylab("Number of Tornadoes") + 
  xlab("Year")

# specify map width to be 970px
# theme_animint() requires Toby's fork of ggplot2
# devtools::install_github("tdhock/ggplot2")
# (we're hopeful this fork will be merged with ggplot2 master)
map <- map + theme_animint(width = 970)

tornado.bar <- list(map = map, ts = ts) 

animint2dir(plot.list=tornado.bar, out.dir="tornado-bar")
## opening a web browser with a file:// URL; if the web page is blank, try running
servr::httd("/ua/snandi/Project_QualityScore/RScripts_QualityScore/tornado-bar")
