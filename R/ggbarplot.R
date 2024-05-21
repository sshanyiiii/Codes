ggbarplot(celltype_num, x = "Celltype", y = "Num",
          fill = "Celltype", 
          color = NULL,
          palette = cell_colors,
          label = T) +
  theme(axis.text.x= element_text(angle = 45, hjust=1, vjust=1))
  theme(aspect.ratio = 0.8)
