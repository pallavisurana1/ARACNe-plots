
## Reference - based on https://christophergandrud.github.io/networkD3/

# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(BiocManager, reshape2, EnsDb.Hsapiens.v86, dplyr, tidyverse, igraph, networkD3, magrittr, linkcomm)

# Helper Functions
process_data_level <- function(data, new_data){
  new_data$Target = data$Regulator
  new_data$Regulator = new_data$MI = NULL
  new_data = merge(new_data, data, by.x="Target", by.y="Target")
  new_data = new_data[!duplicated(new_data[1:2]), ]
  return(new_data[, c(2, 1, 3)])
}

# Import Data
df <- read.csv("/path/to/network/file/network.txt", sep="")
gene <- read.csv("/path/to/network/file/gene_of_interest.txt", sep="")

# Annotation
Regulator_geneID <- ensembldb::select(EnsDb.Hsapiens.v86, keys= c(df$Regulator, df$Target), 
                                      keytype = "GENEID", columns = c("SYMBOL","GENEID"))

# Mapping IDs
map <- setNames(Regulator_geneID$SYMBOL, Regulator_geneID$GENEID)
df1 <- cbind(apply(df, 2, function(x) map[as.character(x)]), df)
df1 <- df1[, -c(3:6, 8)]

# Subset Genes of Interest
df_s_l1 <- subset(df1, df1$Target %in% gene$Symbol)
df_s_l1 <- df_s_l1[order(df_s_l1$Target), ]

# Data Manipulation
df_s_l2 <- process_data_level(df_s_l1, df_s_l1)
df_s_l3 <- process_data_level(df_s_l2, df_s_l2)
df_s_l4 <- process_data_level(df_s_l3, df_s_l3)

# Join all levels
df2 <- bind_rows(df_s_l1, df_s_l2, df_s_l3, df_s_l4)
df2 <- df2[!duplicated(df2[1:2]), ]

# Plot Levels 1 and 2
df.level2 <- bind_rows(df_s_l1, df_s_l2)
df.level2 <- df.level2[!duplicated(df.level2[1:2]), ]

# Plot Networks
lm <- getLinkCommunities(df.level2[1:2])
g <- plot(lm, type = "graph", layout = "spencer.circle")

pdf("my_plot.pdf")
print(g)
dev.off()

# Interactive Plot
df.level2.1 <- subset(df.level2, df.level2$MI > 0.5)
simpleNetwork(df.level2.1) %>% saveNetwork(file = '/path/to/network/file/network.html')
