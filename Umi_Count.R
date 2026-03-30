library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)

group_col <- "celltype"

celltype_order <- c(
  "Artery",
  "Apq7+ CapEC",
  "CapEC2",
  "CapEC3",
  "PCV",
  "HEV"
)

df <- seu@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    celltype = .data[[group_col]],
    umi = .data[["nCount_originalexp"]]
  ) %>%
  filter(!is.na(celltype), !is.na(umi))

df$celltype <- factor(df$celltype, levels = rev(celltype_order))

cols <- c(
  "Artery" = "#00BF7D",
  "Apq7+ CapEC" = "#F8766D",
  "CapEC2" = "#E6C400",
  "CapEC3" = "#A3A500",
  "PCV" = "#1F9CF0",
  "HEV" = "#D65FDE"
)


p <- ggplot(df, aes(x = umi, y = celltype, fill = celltype)) +
  geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.9) +
  geom_jitter(
    width = 0, height = 0.10,
    size = 0.35, alpha = 0.25, color = "black"
  ) +
  stat_summary(
    fun.data = function(x) {
      data.frame(x = median(x), xmin = quantile(x, 0.25), xmax = quantile(x, 0.75))
    },
    geom = "errorbarh",
    height = 0.18, color = "black", linewidth = 0.5
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 18, size = 2.2, color = "black"
  ) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(position = "top", expand = expansion(mult = c(0.01, 0.05))) +
  labs(x = "Number of UMI", y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p

