library(fgsea)
library(ggplot2)
library(gggsea)

# Step 1: Perform GSEA
gsea_results <- fgsea::fgsea(pathways = your_pathways, stats = your_stats)

# Step 2: Create the base enrichment plot
base_plot <- plotEnrichment(gsea_results[gsea_results$pathway == "your_pathway", ],
                            gseaRes = gsea_results)

# Step 3: Convert to ggplot object
base_plot <- base_plot + 
    ggplot2::theme_minimal() + 
    labs(title = "Enhanced GSEA Enrichment Plot", 
         x = "Ranked Genes", 
         y = "Enrichment Score")

# Step 4: Add the GSEA layer using geom_gsea()
gsea_plot <- base_plot + 
    gggsea::geom_gsea(data = gsea_results,
                      aes(x = NES, y = pathway, color = pathway),
                      size = 1) + 
    scale_color_manual(values = c("blue", "red"))  # Customize colors

# Print the final plot
print(gsea_plot)

# Optionally save the plot
ggsave("gsea_plot_with_geom_gsea.png", plot = gsea_plot, width = 8, height = 6, dpi = 300)
