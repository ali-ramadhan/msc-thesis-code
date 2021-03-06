library(hrbrthemes)
library(extrafont)
import_roboto_condensed()
windowsFonts()
font_import()

# Momentum data
library(readr)
OCS7fs <-
  read_delim("D:/Google Drive/GHome/Lab/Geometry reconstruction/Momentum data/OCS 222 momenta/OCS_222_7fs_correct_order.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
OCS7fs <-
  read_delim("D:/Google Drive/GHome/Lab/Geometry reconstruction/Momentum data/OCS 222 momenta/OCS_222_30fs_correct_order.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
OCS7fs <-
  read_delim("D:/Google Drive/GHome/Lab/Geometry reconstruction/Momentum data/OCS 222 momenta/OCS_222_60fs_correct_order.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
OCS7fs <-
  read_delim("D:/Google Drive/GHome/Lab/Geometry reconstruction/Momentum data/OCS 222 momenta/OCS_222_100fs_correct_order.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
OCS7fs <-
  read_delim("D:/Google Drive/GHome/Lab/Geometry reconstruction/Momentum data/OCS 222 momenta/OCS_222_200fs_correct_order.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

drops <- c("X10")
OCS7fs <- OCS7fs[ , !(names(OCS7fs) %in% drops)]

proper_cols <- c(4, 5, 6, 1, 2, 3, 7, 8, 9)
OCS7fs <- OCS7fs[,proper_cols]

colnames(OCS7fs) <- c("O p_x", "O p_y", "O p_z", "C p_x", "C p_y", "C p_z",
                      "S p_x", "S p_y", "S p_z")

melted = melt(OCS7fs)

# Lookup table geometries
library(readr)
OCS7fsg <- read_csv("D:/Google Drive/GHome/Lab/MSc thesis figures/OCS_222_7fs_LT_geometries.csv", col_names = FALSE)
OCS7fsg <- read_csv("D:/Google Drive/GHome/Lab/MSc thesis figures/OCS_222_30fs_LT_geometries.csv", col_names = FALSE)
OCS7fsg <- read_csv("D:/Google Drive/GHome/Lab/MSc thesis figures/OCS_222_60fs_LT_geometries.csv", col_names = FALSE)
OCS7fsg <- read_csv("D:/Google Drive/GHome/Lab/MSc thesis figures/OCS_222_100fs_LT_geometries.csv", col_names = FALSE)
OCS7fsg <- read_csv("D:/Google Drive/GHome/Lab/MSc thesis figures/OCS_222_200fs_LT_geometries.csv", col_names = FALSE)

OCS7fsg <- OCS7fsg[ , 1:3]
colnames(OCS7fsg) <- c("rCO", "rCS", "theta")
OCS7fsg <- OCS7fsg[rowSums(OCS7fsg[, -1])>0, ]

# fmincon/MultiStart geometries
library(readr)
OCS7fsg <- read_csv("D:/Google Drive/GHome/Lab/MSc thesis figures/filtered_OCS_222_7fs_mostLikelyGeometries.csv", col_names = FALSE)
OCS7fsg <- read_csv("D:/Google Drive/GHome/Lab/MSc thesis figures/filtered_OCS_222_30fs_mostLikelyGeometries.csv", col_names = FALSE)
OCS7fsg <- read_csv("D:/Google Drive/GHome/Lab/MSc thesis figures/filtered_OCS_222_60fs_mostLikelyGeometries.csv", col_names = FALSE)
OCS7fsg <- read_csv("D:/Google Drive/GHome/Lab/MSc thesis figures/filtered_OCS_222_100fs_mostLikelyGeometries.csv", col_names = FALSE)
colnames(OCS7fsg) <- c("rCO", "rCS", "theta")

# OCS 222 7fs momentum plots (3x3 grid)
library(reshape2)
library(ggplot2)
library(ggthemes)

# pdf('rplot.pdf')

ggplot(melted, aes(x=value)) +
  facet_wrap(~variable, nrow=3) +
  geom_histogram(aes(y=..density..),
                 bins = 30,
                 color = "#32AAB5",
                 fill = "#32AAB5",
                 alpha = 0.8) +
  geom_line(stat="density", size = 1, color = "#0076A1") +
  # geom_density(size = 1, color = "#0076A1") +
  ylab("Counts (arbitrary units)") +
  xlab(bquote('Momentum ('*10^-22~'kg m/s)')) +
  scale_x_continuous(
    breaks = c(-5e-22, 0, 5e-22),
    label = c("-5", "0", "5")
  ) +
  theme_grey() +
  theme(text=element_text(family="Roboto Condensed"),
        strip.text = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=11),
        axis.title.y = element_text(size=14),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# dev.off()

# OCS 222 7fs momentum pair plots (9x9 pair plot)
library(ggplot2)
library(ggthemes)
library(GGally)

colnames(OCS7fs) <- c("Op_x", "Op_y", "Op_z", "Cp_x", "Cp_y", "Cp_z",
                      "Sp_x", "Sp_y", "Sp_z")

pdf('rplot.pdf')

ggpairs(OCS7fs,
        diag = list(continuous = wrap('densityDiag', color = "#0076A1", fill = "#32AAB5")),
        upper = list(continuous = wrap('density')),
        lower = list(continuous = wrap('points', size = 0.1))) +
  theme_few() +
  theme(text=element_text(family="Roboto Condensed"),
        axis.ticks=element_blank(),
        axis.line=element_blank(), 
        axis.text=element_blank(), 
        panel.grid.major= element_blank(),
        panel.spacing=grid::unit(0,"lines"))

dev.off()

# Molecular geometry parameter distributions (3x3 pair plot)
OCS7fsg[,1] <- OCS7fsg[,1]/1e-10
OCS7fsg[,2] <- OCS7fsg[,2]/1e-10

ggpairs(OCS7fsg, columns = c("rCO", "rCS", "theta"),
        diag = list(continuous = wrap('barDiag', color = "#0076A1", fill = "#32AAB5", bins = 20)),
        upper = list(continuous = wrap('density')),
        lower = list(continuous = wrap('points', size = 0.5))) +
  theme_grey() +
  theme(text=element_text(family="Roboto Condensed"))

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# CO2 calibration plots (from simplex simulations)
library(readr)
CO2_r12_simplex_simulations <- read_csv("D:/Google Drive/GHome/Lab/Geometry reconstruction/Lookup table/Calibration plots/CO2_r12_simplex_simulations.csv", col_names = FALSE)
CO2_r23_simplex_simulations <- read_csv("D:/Google Drive/GHome/Lab/Geometry reconstruction/Lookup table/Calibration plots/CO2_r23_simplex_simulations.csv", col_names = FALSE)
CO2_theta_simplex_simulations <- read_csv("D:/Google Drive/GHome/Lab/Geometry reconstruction/Lookup table/Calibration plots/CO2_theta_simplex_simulations.csv", col_names = FALSE)

# Must be the OCS calibration plots (from simplex simulations)
CO2_r12_simplex_simulations <- read_csv("D:/Google Drive/GHome/Lab/Geometry reconstruction/Lookup table/Calibration plots/varying_r12_full.csv", col_names = FALSE)
CO2_r23_simplex_simulations <- read_csv("D:/Google Drive/GHome/Lab/Geometry reconstruction/Lookup table/Calibration plots/varying_r23_full.csv", col_names = FALSE)
CO2_theta_simplex_simulations <- read_csv("D:/Google Drive/GHome/Lab/Geometry reconstruction/Lookup table/Calibration plots/varying_theta_full.csv", col_names = FALSE)

require(ggplot2)
require(ggthemes)
library(hrbrthemes)

g11 <- ggplot(CO2_r12_simplex_simulations, aes(x=X1, y=X4)) +
  geom_abline(intercept = 0, slope = 1, color="black") +
  geom_point(size=1.5, color="#0076A1") +
  xlab(expression(r[12]*" input ("*ring(A)*")")) +
  ylab(expression(r[12]*" output ("*ring(A)*")")) +
  theme_grey() +
  theme(text=element_text(family="Roboto Condensed"))

g12 <- ggplot(CO2_r12_simplex_simulations, aes(x=X1, y=X5)) +
  geom_hline(yintercept = 1.16, color="black") +
  # geom_hline(yintercept = 1.56, color="black") +
  geom_point(size=1.5, color="#0076A1") +
  xlab(expression(r[12]*" input ("*ring(A)*")")) +
  ylab(expression(r[23]*" output ("*ring(A)*")")) +
  theme_grey() +
  theme(text=element_text(family="Roboto Condensed"))

g13 <- ggplot(CO2_r12_simplex_simulations, aes(x=X1, y=X6)) +
  geom_hline(yintercept = 180, color="black") +
  geom_point(size=1.5, color="#0076A1") +
  xlab(expression(r[12]*" input ("*ring(A)*")")) +
  ylab(expression(theta*" output (degrees)")) +
  theme_grey() +
  theme(text=element_text(family="Roboto Condensed"))

g21 <- ggplot(CO2_r23_simplex_simulations, aes(x=X2, y=X4)) +
  geom_hline(yintercept = 1.16, color="black") +
  geom_point(size=1.5, color="#0076A1") +
  xlab(expression(r[23]*" input ("*ring(A)*")")) +
  ylab(expression(r[12]*" output ("*ring(A)*")")) +
  theme_grey() +
  theme(text=element_text(family="Roboto Condensed"))

g22 <- ggplot(CO2_r23_simplex_simulations, aes(x=X2, y=X5)) +
  geom_abline(intercept = 0, slope = 1, color="black") +
  geom_point(size=1.5, color="#0076A1") +
  xlab(expression(r[23]*" input ("*ring(A)*")")) +
  ylab(expression(r[23]*" output ("*ring(A)*")")) +
  theme_grey() +
  theme(text=element_text(family="Roboto Condensed"))

g23 <- ggplot(CO2_r23_simplex_simulations, aes(x=X2, y=X6)) +
  geom_hline(yintercept = 180, color="black") +
  geom_point(size=1.5, color="#0076A1") +
  xlab(expression(r[23]*" input ("*ring(A)*")")) +
  ylab(expression(theta*" output (degrees)")) +
  theme_grey() +
  theme(text=element_text(family="Roboto Condensed"))

g31 <- ggplot(CO2_theta_simplex_simulations, aes(x=X3, y=X4)) +
  geom_hline(yintercept = 1.16, color="black") +
  geom_point(size=1.5, color="#0076A1") +
  xlab(expression(theta*" intput (degrees)")) +
  ylab(expression(r[12]*" output ("*ring(A)*")")) +
  theme_grey() +
  theme(text=element_text(family="Roboto Condensed"))

g32 <- ggplot(CO2_theta_simplex_simulations, aes(x=X3, y=X5)) +
  geom_hline(yintercept = 1.16, color="black") +
  # geom_hline(yintercept = 1.56, color="black") +
  geom_point(size=1.5, color="#0076A1") +
  xlab(expression(theta*" output (degrees)")) +
  ylab(expression(r[23]*" output ("*ring(A)*")")) +
  theme_grey() +
  theme(text=element_text(family="Roboto Condensed"))

g33 <- ggplot(CO2_theta_simplex_simulations, aes(x=X3, y=X6)) +
  geom_abline(intercept = 0, slope = 1, color="black") +
  # geom_abline(intercept = 360, slope = -1, color="black") +
  geom_point(size=1.5, color="#0076A1") +
  xlab(expression(r[23]*" input ("*ring(A)*")")) +
  ylab(expression(theta*" output (degrees)")) +
  theme_grey() +
  theme(text=element_text(family="Roboto Condensed"))

# require(cowplot)
# plot_grid(g11, g12, g13, g21, g22, g23, g31, g32, g33)

multiplot(g11, g12, g13, g21, g22, g23, g31, g32, g33,
          layout=matrix(c(1,2,3,4,5,6,7,8,9), nrow=3, byrow=TRUE))

# Lookup table reconstruction pair plots
