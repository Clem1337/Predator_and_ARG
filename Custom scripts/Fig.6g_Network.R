library(dplyr)
library(igraph)
library(ggraph)
library(ggplot2)

file_path <- "/path/to/your/file.txt"
data <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

edge_data <- data %>%
  filter(Sample_id %in% c("Non_p1", "Non_p2", "Non_p3")) %>%
  group_by(Resistance, Genus) %>%
  summarise(weight = n(), .groups = "drop") 

nodes_resistance <- data.frame(name = unique(edge_data$Resistance), type = "Resistance")
nodes_genus <- data.frame(name = unique(edge_data$Genus), type = "Genus")
node_data <- bind_rows(nodes_resistance, nodes_genus)
node_data <- node_data[!duplicated(node_data$name), ]

net_graph <- graph_from_data_frame(d = edge_data, vertices = node_data, directed = FALSE)

p <- ggraph(net_graph, layout = 'circle') +
  
  geom_edge_arc(aes(edge_width = weight), 
                color = "#CCCCCC", 
                strength = 0.2, 
                alpha = 0.8) +
  
  geom_node_point(aes(shape = type, fill = type), color = "black", size = 6, stroke = 1) +
  
  geom_node_text(aes(label = name), repel = TRUE, size = 3.5, fontface = "bold") +
  
  scale_edge_width_continuous(
    range = c(0.5, 4),  
    breaks = function(x) seq(ceiling(min(x)), floor(max(x))), 
    name = "Connection\nTimes"
  ) +
  
  scale_shape_manual(values = c("Genus" = 22, "Resistance" = 21)) +
  
  scale_fill_manual(values = c("Genus" = "#72B063", "Resistance" = "#E38965")) +
  
  theme_void() +
  theme(legend.position = "right")

print(p)