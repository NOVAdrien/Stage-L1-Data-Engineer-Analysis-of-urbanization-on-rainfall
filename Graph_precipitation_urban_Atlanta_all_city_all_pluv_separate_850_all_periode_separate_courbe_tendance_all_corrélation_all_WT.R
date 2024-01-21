# Rain data analysis

# charge libraries
library(ggplot2)
library(ggstatsplot)
library(ggside)

# SET WORKING DIRECTORY
DIR = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(DIR)

# all cities to run : "Atlanta", "Chicago", "Dallas", "Minneapolis", "New York", "Saint-Louis"

# set options to run
all_city <- c("Atlanta", "Chicago", "Dallas", "Minneapolis", "New York", "Saint-Louis")
all_cor_type <- c("parametric") # option to run : "nonparametric", "robust", "bayes"
Rmin <- 0
Rmax <- 20
Angle <- 45
Plevel <- 850
all_WT <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")

# import wind data
for (city in all_city) {
  all_pluv_nb <- substr(list.files(paste0("F:Nouvelles donnees/analyses/", city)), 1,6)
  for (wind_type in all_WT) {
    Table <- data.frame()
    Table2 <- data.frame()
    for (pluv_nb in all_pluv_nb) {
      table <- read.table(paste0(city, '/', pluv_nb, '_Rmin_', Rmin, '_Rmax_', Rmax, '_Angle_', Angle, '_Plevel_', Plevel, ".txt"),
                          sep = ";",
                          row.names = 1,
                          header = TRUE)
      
      
      # test data here
      if (length(table$ID) == 0) {
        print("no data during this period for that station")
      }
      else {
        Table <- rbind(Table, table)
        }
    }
    # create and save graph
    figure <- plot(Table$urban[Table$WT == wind_type], Table$precip[Table$WT == wind_type],
                   xlab = "urban percentage (%)",
                   ylab = "precipitaion (mm/year)",
                   xlim = NULL,
                   ylim = NULL)
    
    Figure <- ggplot(figure, aes(Table$urban[Table$WT == wind_type], Table$precip[Table$WT == wind_type])) +
      labs(title = paste0("Precipitation amount in ", city, " depending on urban surface percentage with wind from", wind_type)) +
      geom_point(colour = "blue") +
      geom_smooth(colour = "red")
    
    png(file=paste0(city, '/', 'Graph_2002_2017_', city, '_', wind_type, '.png'), width = 20, height = 20, units = "cm", res = 600)
    
    print(Figure)
    
    dev.off()
    
    # only keep lines of Table with correct WT
    Table <- Table[Table$WT == wind_type,]
    
    # make a graph for each correlation type
    for (cor_type in all_cor_type) {
      a <- ggscatterstats(data = Table,
                          x = urban,
                          y = precip,
                          type = cor_type)
      
      # save graph
      png(file=paste0(city, '/', 'Graph_2002_2017_', city, '_', wind_type, '_cbe_tdce_cor_', cor_type, '.png'), width = 20, height = 20, units = "cm", res = 600)
      
      print(a)
      
      dev.off()
    }
  }
}
