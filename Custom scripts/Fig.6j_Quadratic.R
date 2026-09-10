library(ggplot2)
library(dplyr)

file_path <- "/path/to/your/file.txt"
df <- read.table(file_path, header = TRUE, sep = "\t")

df <- na.omit(df)

p <- ggplot(df, aes(x = Pressure, y = ARG)) +
  
  geom_point(aes(color = Pressure), shape = 16, alpha = 0.5, size = 4) +
  scale_color_gradient(low = "#E5DCC5", high = "#C46673") +
  
  geom_smooth(method = "lm", 
              formula = y ~ poly(x, 2), 
              color = "grey30", 
              linetype = "dashed",
              se = TRUE,        
              fill = "grey80",   
              alpha = 0.5) +
  
  theme_classic() +
  labs(
    x = expression(bold("Predation Pressure" ~ Log[10] * "(Predator / Bacteria)")),
    y = expression(bold("ARG Abundance")),
    color = "Pressure"
  ) +
  
  scale_x_continuous(limits = c(-9, 0), breaks = seq(-9, 0, by = 3), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-4, 3), breaks = seq(-4, 3, by = 1), expand = c(0, 0)) +
  
  theme(
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(size = 13),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 10)
  )

poly_mod <- lm(ARG ~ poly(Pressure, 2, raw = TRUE), data = df)
coefs <- coef(poly_mod)

a <- coefs[3]
b <- coefs[2]

vertex_x <- -b / (2 * a)
cat("Calculated threshold:", round(vertex_x, 3), "\n")

p_final <- p + 
  geom_vline(xintercept = vertex_x, color = "red", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = vertex_x + 0.2, y = max(df$ARG) * 0.9, 
           label = paste("Threshold:", round(vertex_x, 2)), 
           color = "red", hjust = 0)

print(p_final)