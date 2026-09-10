if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(VennDiagram)
library(dplyr)
library(grid)

file_path <- "/path/to/your/file.txt"

df <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-16LE")

colnames(df)[1] <- "Mode"

# head(df)

sample_cols <- 7:ncol(df)

get_detected_samples <- function(mode_name) {
  sub_df <- df %>% filter(tolower(Mode) == tolower(mode_name))
  
  if (nrow(sub_df) == 0) return(character(0))
  
  sample_sums <- colSums(sub_df[, sample_cols, drop = FALSE], na.rm = TRUE)
  
  detected_samples <- names(sample_sums[sample_sums > 0])
  return(detected_samples)
}

set_phototroph  <- get_detected_samples("phototroph")
set_unknown     <- get_detected_samples("unknown")
set_symbiotroph <- get_detected_samples("symbiotroph")
set_predator    <- get_detected_samples("predator")

venn_list <- list(
  Phototroph  = set_phototroph,
  Unknown     = set_unknown,
  Symbiotroph = set_symbiotroph,
  Predator    = set_predator
)

grid.newpage()

fill_colors <- c("#364B9A", "#C2E4EF", "#FEDA8B", "#A50026")

venn_plot <- venn.diagram(
  x = venn_list,
  category.names = c("Phototroph", "Unknown", "Symbiotroph", "Predator"),
  filename = NULL,
  fill = fill_colors,
  alpha = 0.7,
  col = "black",
  lwd = 2,
  lty = "dashed",
  cex = 1.2,
  cat.cex = 1.3,
  cat.col = "black",
  cat.dist = c(0.2, 0.2, 0.1, 0.1),
  margin = 0.1
)

grid.draw(venn_plot)