library(igraph)
pairwise_data <- read.table('repeat.aba.blast4plot')
View(pairwise_data)
pairwise_data=pairwise_data[,1:2]
graph <- graph_from_data_frame(pairwise_data, directed = FALSE)
communities <- cluster_louvain(graph)
print(membership(communities))
community_colors <- rainbow(max(membership(communities)) + 1)
pdf('repeat_class.pdf',width = 8,height = 8)
plot(
    graph,
    layout = layout_with_fr(graph),  # You can use different layout algorithms
    vertex.color = community_colors[membership(communities)],
    vertex.size = 10,  # Adjust node size
    vertex.label.cex = 0.3,  # Adjust label size
    main = "Repeat class"
)

# Add a legend
legend(
    "bottomright",
    legend = unique(membership(communities)),
    col = community_colors,
    pch = 20,
)
dev.off()