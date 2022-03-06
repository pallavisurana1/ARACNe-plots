
#---https://christophergandrud.github.io/networkD3/

#---get packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(BiocManager, reshape2, EnsDb.Hsapiens.v86, dplyr, tidyverse, igraph, networkD3, magrittr, linkcomm)

#---get ARACNe data to make network
df = read.csv("/path/to/network/file/network.txt", sep="")

#---read genes of interest
gene = read.csv("/path/to/network/file/gene_of_interest.txt", sep="")

#---annotate the ARACNe output
Regulator_geneID <- ensembldb::select(EnsDb.Hsapiens.v86, keys= c(df$Regulator, df$Target), 
                                      keytype = "GENEID", columns = c("SYMBOL","GENEID"))

#---map ids
map <- setNames(Regulator_geneID$SYMBOL, Regulator_geneID$GENEID)
x = apply(df, 2, function(x) map[as.character(x)])
df1 = cbind(x, df)

#---subset data
df1 = df1[ -c(3:6, 8)]

#---subset genes of interest for level 1
df_s_l1 = subset(df1, df1$Target %in% spingo$Symbol)
df_s_l1 = df_s_l1[order(df_s_l1$Target), ]

#----- Level 2 ------# Regulator becomes Target
df_s_l2 = df_s_l1
df_s_l2$Target = df_s_l1$Regulator
df_s_l2$Regulator = df_s_l2$MI = NULL
df_s_l2$Target <- df_s_l2[order(df_s_l2$Target),]

#---find targets for the regulators
df_s_l2 = merge(df_s_l2, df1, by.x="Target", by.y="Target")
# t = subset(t, t$Target %in% spingo$Symbol)
df_s_l2 = df_s_l2[!duplicated(df_s_l2[1:2]), ]
df_s_l2 = df_s_l2[, c(2,1,3)]


#----- Level 3 ------# Regulator becomes Target
df_s_l3 = df_s_l2
df_s_l3$Target = df_s_l2$Regulator
df_s_l3$Regulator = df_s_l3$MI = NULL
df_s_l3$Target <- df_s_l3[order(df_s_l3$Target),]

#---find targets for the regulators
df_s_l3 = merge(df_s_l3, df1, by.x="Target", by.y="Target")
df_s_l3 = df_s_l3[!duplicated(df_s_l3[1:2]), ]
df_s_l3 = df_s_l3[, c(2,1,3)]


#----- Level 4 ------# Regulator becomes Target
df_s_l4 = df_s_l3
df_s_l4$Target = df_s_l3$Regulator
df_s_l4$Regulator = df_s_l4$MI = NULL
df_s_l4$Target <- df_s_l4[order(df_s_l4$Target),]

# find targets for the regulators
df_s_l4 = merge(df_s_l4, df1, by.x="Target", by.y="Target")
df_s_l4 = df_s_l4[!duplicated(df_s_l4[1:2]), ]
df_s_l4 = df_s_l4[, c(2,1,3)]


#---Join
#---join all 4 levels
df2 = bind_rows(df_s_l1, df_s_l2, df_s_l3, df_s_l4)
df2 = df2[!duplicated(df2[1:2]), ]

# write.csv2(df2, "path/to/processed/data/prcoessed_network.txt", row.names = FALSE)


#---plot only levels 1 and 2
#---join ------ level 1 and 2 
df.level2 = bind_rows(df_s_l1, df_s_l2)
df.level2 = df.level2[!duplicated(df.level2[1:2]), ]


#---plot - extracts link communities from networks
lm <- getLinkCommunities(df.level2[1:2])

g = plot(lm, type = "graph", layout = "spencer.circle")

pdf("my_plot.pdf")
g
dev.off() 

htmlwidgets::saveWidget(ggplotly(g), "/path/to/network/file/network.html")


#----plot interactive simple network
df.level2.1 = subset(df.level2, df.level2$MI>0.5)
simpleNetwork(df.level2.1) %>% saveNetwork(file = '/path/to/network/file/network.html')


#------------ END ---------------#



