library(plotly)
library(Seurat)

emb3d <- as.data.frame(Embeddings(seu_ec_tspace, "tPC")[, 1:3])
colnames(emb3d) <- c("tPC_1", "tPC_2", "tPC_3")
emb3d$cell <- rownames(emb3d)
emb3d$celltype <- as.character(Idents(seu_ec_tspace))

plot_ly(
  data = emb3d,
  x = ~tPC_1,
  y = ~tPC_2,
  z = ~tPC_3,
  color = ~celltype,
  colors = "Set1",
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, opacity = 0.85),
  text = ~paste(
    "Cell:", cell,
    "<br>Cell type:", celltype,
    "<br>tPC_1:", round(tPC_1, 3),
    "<br>tPC_2:", round(tPC_2, 3),
    "<br>tPC_3:", round(tPC_3, 3)
  ),
  hoverinfo = "text"
) %>%
  layout(
    title = "tSpace embedding (tPC_1, tPC_2, tPC_3)",
    scene = list(
      xaxis = list(title = "tPC_1"),
      yaxis = list(title = "tPC_2"),
      zaxis = list(title = "tPC_3")
    )
  )
library(plotly)
library(Seurat)

emb3d <- as.data.frame(Embeddings(seu_ec_tspace, "tPC")[, 1:3])
colnames(emb3d) <- c("tPC_1", "tPC_2", "tPC_3")
emb3d$cell <- rownames(emb3d)
emb3d$pseudotime <- seu_ec_tspace$MLN_tSpace_pseudotime

plot_ly(
  data = emb3d,
  x = ~tPC_1,
  y = ~tPC_2,
  z = ~tPC_3,
  color = ~pseudotime,
  colors = "Viridis",
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, opacity = 0.85),
  text = ~paste(
    "Cell:", cell,
    "<br>Pseudotime:", round(pseudotime, 3),
    "<br>tPC_1:", round(tPC_1, 3),
    "<br>tPC_2:", round(tPC_2, 3),
    "<br>tPC_3:", round(tPC_3, 3)
  ),
  hoverinfo = "text"
) %>%
  layout(
    title = "tSpace pseudotime (tPC_1, tPC_2, tPC_3)",
    scene = list(
      xaxis = list(title = "tPC_1"),
      yaxis = list(title = "tPC_2"),
      zaxis = list(title = "tPC_3")
    )
  )

p3d <- plot_ly(
  data = emb3d,
  x = ~tPC_1, y = ~tPC_2, z = ~tPC_3,
  color = ~pseudotime,
  colors = "Viridis",
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, opacity = 0.85)
)

htmlwidgets::saveWidget(p3d, "tspace_3d_pseudotime.html", selfcontained = TRUE)

dim(Embeddings(seu_ec_tspace, "tPC"))