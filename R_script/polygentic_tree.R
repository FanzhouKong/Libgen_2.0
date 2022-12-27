# 
# list.of.packages <- c("ape","ggtree","phangorn","ggplot2")
library(ggtree)
library(ggplot2)
library(phangorn)
library(ape)
data_temp<- read.csv('/Users/fanzhoukong/Documents/GitHub/Libgen_data/EAD/proxi_matrix_UVPD.csv',header=TRUE,row.names = 1,sep=",")
# data_temp<- read.csv('/Users/fanzhoukong/Documents/GitHub/Libgen_data/EAD/proxi_matrix_HCD.csv',header=TRUE,row.names = 1,sep=",")
# data_temp<- read.csv('/Users/fanzhoukong/Documents/GitHub/Libgen_data/EAD/proxi_matrix_ead.csv',header=TRUE,row.names = 1,sep=",")
# data_temp
colnames(data_temp)
proxi_matrix <-subset(data_temp, select = -c(Library_id,class, superclass,kingdom) )
rownames(proxi_matrix)<- data_temp$Library_id
colnames(proxi_matrix)<- data_temp$Library_id
table(data_temp$superclass)
# proxi_matrix
proxi_matrix_dist <- as.dist(proxi_matrix)
tupgma <- upgma(proxi_matrix_dist, method="average")
class <- data_temp$superclass
ids <- data_temp$Library_id
class_data <- cbind(as.data.frame(data_temp$Library_id),as.data.frame(class))
cols <- c( '#FF0001', '#00FF00',"#1400FF", '#FEFF01','#FF00FF','#05FFFF','#810000',
           '#008000','#050080','#818100','#800080','#028180','#E6E6E6','#9A99FF')
t4 <- ggtree(tupgma, layout="circular", size = 0.2) 
t4 <- t4 %<+% class_data + 
  # geom_tiplab(aes(color = class),hjust = - 0.2, geom = 'text',show.legend = FALSE)+
  # geom_text(aes(color = class, label = label))+
  geom_tippoint(aes(color = class), alpha = 1, size = 0.8) +
  scale_color_manual(values = cols) +
  theme(legend.position="right") +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  theme(plot.background = element_rect(fill = 'white', colour = 'white')) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(legend.background = element_rect(fill = 'white'),
        legend.title = element_text(size = 8,face = "bold")) + # LEGEND TITLE
  theme(legend.text = element_text(size = 8, face = "bold")) + # LEGEND TEXT
  theme(text = element_text(color = "black", size = 8)) +
  geom_treescale(x = 0.1,y = 0.1, width = 0.1, offset = NULL,
                 color = "white", linesize = 1E-100, fontsize = 1E-100)

ggsave(file="/Users/fanzhoukong/Documents/GitHub/Libgen_data/EAD/uvpd.png", plot=t4, width=10, height=8, units = 'in', dpi = 300)


