#! /usr/bin/env Rscript

# generate box plots for paired data


#################################################################
library(ggplot2)
library(reshape2)
library(gridExtra)


#################################################################
args <- commandArgs(trailingOnly = TRUE)
infile = args[1]
outfile = args[2]


#################################################################
# Read the CSV file
data <- read.csv(infile, sep="\t", header=F)

V1 <- data[, 1]
V2 <- data[, 2]


#################################################################
# 2 graphs
#par(mfrow=c(2,2))
pdf(outfile)

# Create a data frame for plotting
plot_data <- data.frame(
  Values = c(V1, V2),
  Group = factor(rep(c("V1", "V2"), each = length(V1)))
)

# Create the paired box plot
ggplot(plot_data, aes(x = Group, y = Values)) +
  geom_boxplot(aes(fill = Group), width = 0.4, position = position_dodge(0.6)) +
  geom_line(aes(group = rep(1:length(V1), 2)), color = "gray", position = position_dodge(0.6)) +
  theme_minimal() +
  labs(x = "Group", y = "Value", title = "Paired Box Plots for V1 and V2") +
  theme(legend.position = "none")


#################################################################
#dev.new()
# Set layout to plot 2 graphs (one below the other)
#layout(matrix(c(1, 2), 2, 1, byrow = TRUE))
# Plot the first graph using ggplot2 (need to use grid.draw to work with layout)
#grid.draw(ggplotGrob(p1))


# density compare
plot(density(V1))
lines(density(V2), col="red")
#curve(density(V2), from=min(V2)*0.8, to=max(V2)*1.2, add=T, col="red")
dev.off()

