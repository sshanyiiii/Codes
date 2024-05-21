#                         prediction.score.max annotation_correct
# WT125_AAACAGCCAATAGCAA-1            1.0000000               TRUE
# WT125_AAACAGCCACCTAAGC-1            1.0000000               TRUE

ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) + 
  geom_density(alpha = 0.5) + 
  theme_cowplot() + 
  scale_fill_discrete(name = "Annotation Correct",   # 设置fill的图例
                     labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) +
  scale_color_discrete(name = "Annotation Correct",   # 设置color的图例
                       labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + 
  xlab("Prediction Score") +
  theme(aspect.ratio = 1)
