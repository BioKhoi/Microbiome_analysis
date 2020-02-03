
library(ggplot2)
library (vegan)
library (gtools)


data <- as.data.frame(read.table("/Users/Nguyenphuong/Documents/GitHub/Microbiome_analysis/input_data.txt", header = TRUE, sep="\t"))


metadata <- as.data.frame(read.table("/Users/Nguyenphuong/Documents/GitHub/Microbiome_analysis/mapping_data.txt", header=TRUE, sep="\t"))


samples_in_both <- intersect(data$samples, metadata$samples)
data <- data[data$samples %in% samples_in_both, ]
metadata <- metadata[metadata$samples %in% samples_in_both, ]

metadata <- metadata[match(data$samples, metadata$samples),]
data$samples <- factor(data$samples)
metadata$samples <- factor(metadata$samples)


mds <- metaMDS(data[,-1], distance="euclidean")

mds$points


# plot NMDS using basic plot function and color points by "Country" from metadata
plot(mds$points, col=metadata$Treatment)

# draw dispersion ellipses around data points
ordiellipse(mds, metadata$Treatment, display = "sites", kind = "sd", label = T)


# Now a fancier plot
# Data frame df_ell contains values to show ellipses. 
# It is calculated with function veganCovEllipse which is hidden in vegan package. This function is applied to each level of NMDS (group) and it uses also function cov.wt to calculate covariance matrix.
# Ideas taken from http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo and from Meren's script.
NMDS = data.frame(MDS1 = mds$points[,1], MDS2 = mds$points[,2], group=metadata$Treatment)
NMDS.mean=aggregate(NMDS[,1:2],list(group=metadata$Treatment),mean)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()
for(g in levels(NMDS$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),
        length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2))))),group=g))
}

  ggplot(data = NMDS, aes(MDS1,MDS2, colour = group,  shape = group, fill = group)) + 
  geom_point(aes(), shape = 21, size = 6, stroke = 1, alpha = 0.7) +
  
  scale_color_manual(values = c("purple","deepskyblue","blue", "black", "orange", "red")) +
    
  scale_fill_manual(values = c("lightslateblue","lightblue","deepskyblue", "gray50", "salmon1", "indianred1")) +
  
  geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size= 1., linetype=4) +
  
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.grid.major=element_line(linetype=3,color="white"))
        
      


  
  



