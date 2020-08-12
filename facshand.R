# This script is designed for FACS results handling derived from Flowjo
# Input files: data file exported from Flowjo, cell type file, cell number file, group file
# Output files: data file containing proportion and cell number columns, studnet t test results file.

# Path file
args = commandArgs(T)
path_args = args[1]
if (is.na(path_args)) {
  setwd(path)
  print("Please input path.")
} else {
  setwd(path_args)
}

library(ggplot2)
library(Hmisc)

# sample name should be formated as group_whatever.fcs
print("Read the input files...")
rawdata = read.table("rawdata.csv", header = T, sep = ",") 
colnames(rawdata) = c("Depth", "Name", "Statistic", "Cells")
# celltype should not be identical
celltype = read.table("celltype.csv", header = T, sep = ",")
# cell number should be handled as number/10^4 level
cellnumber = read.table("cellnumber.csv", header = T, sep = ",")
group = read.table("group.csv", header = T, sep = ",")

# for flowjo version <= 10.1, transform subgating name to full path name

print("Detect the flowjo version...")

if (length(grep("fcs", rawdata[rawdata$Depth == "> ", ]$Name, value = T)) == 0) { # full path detection
  print("Flowjo version <= 10.1, rawdata handling...")
  rawdata$Rank = 0
  for (i in 1:nrow(rawdata)) {
    n = nchar(as.character(rawdata[i, ]$Depth)) / 2
    rawdata[i, ]$Rank = n
  }
  rawdata$Name = as.character(rawdata$Name)
  for (i in 1:nrow(rawdata)){
    rank = rawdata[i, ]$Rank
    if (rank == 0) {
      name = rawdata[i, ]$Name
    }
    if (rank > 0) {
      name = as.character(paste(name, rawdata[i, ]$Name, sep = "/"))
      rawdata[i, ]$Name = name
    }
  }
  write.table(rawdata, "rawdata_full_path.csv", quote = F, sep = ",", col.names = T, row.names = F)
} else {
  print("Flowjo version > 10.1, continues...")
}

# fetch total cell number in fcs files
fcs = subset(rawdata, rawdata$Depth == "", select = c(Name, Cells))
fcsname = as.data.frame(fcs$Name)
fcstotal = as.data.frame(fcs$Cells)

# sort treated groups
rawdata1 = data.frame()
for (i in 1:nrow(group)){
  group_sort = paste("group", i, sep = "")
  group_sort = rawdata[grep(group[i, ], rawdata$Name, ignore.case = T), ]
  group_sort$group = rep(paste(group[i, ]), nrow(group_sort))
  rawdata1 = rbind(rawdata1, group_sort)
}

# sort celltype
data = data.frame()
for (i in 1:nrow(celltype)){
  cell = paste(celltype[i, ], "$", sep = "")
  data1 = rawdata1[grep(cell, rawdata1$Name, ignore.case = T),]
  data1$celltype = rep(celltype[i, ], nrow(data1))
  data = rbind(data, data1)
}

# calculate percentage of total cell number
data3 = data.frame()
for (j in 1:nrow(fcsname)){
  data2 = data[grep(fcsname[j, ], data$Name, ignore.case = T),]
  data2$total = rep(fcstotal[j, ], nrow(data2))
  data2$percentage = data2$Cells/data2$total
  data3 = rbind(data3, data2)
}

datacal = data.frame(data3$Name, data3$group, data3$celltype, data3$Statistic, data3$percentage)

# attach cell number data
tissue = as.data.frame(cellnumber$tissue)
number = as.data.frame(cellnumber$number)
datafinal = data.frame()

for (i in 1:nrow(tissue)){
  addition = datacal[grep(tissue[i, ], datacal$data3.Name, ignore.case = T), ]
  addition$number = rep(number[i, ], nrow(addition))
  datafinal = rbind(datafinal, addition)
}

colnames(datafinal) = c("name", "group", "celltype", "population", "percentage", "total")
datafinal$cellnumber = datafinal$percentage*datafinal$total


print("Exprot data files...")
write.table(datafinal, "data.csv", quote = F, row.names = F, col.names = T, sep = ",")

# Plot as dot-errorbar fig
print("Plotting...")
ylab = expression(paste("Cell number", " (×", 10^"4", ")", sep = ""))

# Output statistics file
sink("statistic.txt")

for (i in 1:nrow(celltype)){
  name = paste(celltype[i, ], "Population")
  name_num = paste(celltype[i, ], "Cell number")
  name_ttest = paste(celltype[i, ], "ttest.txt")
  type = celltype[i, ]
  plotdata = subset(datafinal, datafinal$celltype == type)
  
  if (nrow(group) == 2) {
    print(name)
    ttest = t.test(population ~ group, data = plotdata, paired = F)
    p_pop = signif(ttest$p.value, 2)
    print(ttest)
    
    print(name_num)
    ttest1 = t.test(cellnumber ~ group, data = plotdata, paired = F)
    p_num = signif(ttest1$p.value, 2)
    print(ttest1)
  } else {
    print(name)
    anova = aov(population ~ group, plotdata)
    p_pop = signif(summary(anova)[[1]][, "Pr(>F)"][1], 2)
    print(summary(anova))
    
    print(name_num)
    anova1 = aov(population ~ group, plotdata)
    p_num = signif(summary(anova1)[[1]][, "Pr(>F)"][1], 2)
    print(summary(anova1))
  }

  plot = ggplot(plotdata, aes(group, population)) + 
    geom_point(position= position_jitter(height = 0, width = 0.25), size = 6, alpha = 0.9) +
    stat_summary(fun = "mean", fun.max = "mean", fun.min = "mean", geom = "errorbar", color = "red", size = 2) +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", color = "red", width = 0.4, size = 1) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 1.5*max(plotdata$population))) +
    xlab(NULL) + ylab("Population %") + labs(title = paste("p = ", p_pop)) +
    theme(axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 30), 
          axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(size = 1, lineend = "butt"),
          axis.ticks = element_line(size = 1),
          plot.title = element_text(size = 20, hjust = 0.5))
  ggsave(paste(name, ".tiff", sep = ""),  plot, device = "tiff", height = 4.86, width = 1.5 * nrow(group))
  
  plot_num = ggplot(plotdata, aes(group, cellnumber)) + 
    geom_point(position= position_jitter(height = 0, width = 0.25, seed = NULL), size = 6, alpha = 0.9) +
    stat_summary(fun = "mean", fun.max = "mean", fun.min = "mean", geom = "errorbar", color = "blue", size = 2) +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", color = "blue", width = 0.4, size = 1) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 1.8*max(plotdata$cellnumber))) + 
    xlab(NULL) + ylab(ylab) + labs(title = paste("p = ", p_num)) +
    theme(axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 30), 
          axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30),
          axis.line = element_line(size = 1, lineend = "butt"),
          axis.ticks = element_line(size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size = 20, hjust = 0.5))
  ggsave(paste(name_num, ".tiff", sep = ""),  plot_num, device = "tiff", height = 4.86, width = 1.5 * nrow(group))
}
sink()
print("Completed :p")
