ggplot(stat_cell_num_long, aes(x = sample, y = cell_num, fill = state)) +
  geom_bar(stat="identity", width=0.65, position=position_dodge()) +
  scale_fill_manual(values = c("#999999", "#ef8a62")) +
  theme_classic() +
  labs(x = "",
       y = "Number of cell",
       title = "QC") +
  theme(axis.text.x= element_text(angle = 45, hjust=1, vjust=1)) +
  theme(aspect.ratio = 1)
