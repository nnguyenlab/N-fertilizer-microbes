library(ggplot2)
theme_set(theme_bw())
#Color Brewer Pallette
#RColorBrewer::display.brewer.all()

#Manual color scales
Brewer9 <- c("Acidobacteria" = "#a6cee3", "Actinobacteria" = "#1f78b4", "Bacteroidetes" = "#b2df8a", "Chloroflexi" = "#33a02c", "Crenarchaeota" = "#fb9a99", "Nitrospirae" = "#e31a1c", "Planctomycetes" = "#fdbf6f", "Proteobacteria" = "#ff7f00", "Verrucomicrobia" = "#cab2d6")
Rick <- c("Acidobacteria" = "red1", "Actinobacteria" = "tomato", "Armatimonadetes" = "darkorange", "Bacteroidetes" = "gold1", "Bdellovibrionota" = "bisque3", "Chloroflexi" = "green", "Deinococcota" = "green", "Firmicutes" = "seagreen", "Myxococcota" = "darkturquoise", "Nitrospirae" = "blue", "OP3" = "darkblue", "Planctomycetes" = "deepskyblue","Proteobacteria" = "orchid2", "Thermi" = "hotpink", "Verrucomicrobia" = "purple1", "WS3" = "magenta1")
Dark <- c("Acidobacteria" = "#1c9e77", "Actinobacteria" = "#d95f02", "Armatimonadetes" = "darkorange", "Bacteroidetes" = "#756fb4", "Bdellovibrionota" = "bisque3", "Crenarchaeota (Archaea)" = "#e6ab01", "Chlamydiae" = "#ffd92f" , "Chlorobi" = "#00665d", "Chloroflexi" = "#66a621", "Cyanobacteria"="#02818a", "Deinococcota" = "green", "Elusimicrobia" = "#6a51a3", "Firmicutes" = "#4536e4", "Gemmatimonadetes" = "#adc9c5", "Myxococcota" = "darkturquoise", "Nitrospirae" = "#f45c85", "OP3" = "black", "Planctomycetes" = "#a6761d","Proteobacteria" = "#0c7be3", "Thermi" = "hotpink", "Verrucomicrobia" = "#646782", "WS3" = "magenta1", "Other Bacteria" = "black")


data <- read.csv("taxa-bar-data-level-2-tidy.csv", header=TRUE)
data <- data.frame(data)

ggplot(data, aes(x=Treatment_Time, y=Abundance, fill=Phylum)) +
  geom_col(position="fill") +
  scale_fill_manual(values = Dark) +
  labs(x = "Treatment", y = "Abundance (%)") +
  facet_grid(cols = vars(Time), scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(labels=c("Compost_0"="Compost","Compost_1"="Compost","Compost_2"="Compost","Compost_3"="Compost","Compost_4"="Compost","Compost_5"="Compost","Compost_6"="Compost","Control_0"="Control","Control_1"="Control","Control_2"="Control","Control_3"="Control","Control_4"="Control","Control_5"="Control","Control_6"="Control","Urea_0"="Urea","Urea_1"="Urea","Urea_2"="Urea","Urea_3"="Urea","Urea_4"="Urea","Urea_5"="Urea","Urea_6"="Urea"))   
