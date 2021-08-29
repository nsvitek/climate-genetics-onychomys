# environment PCA
# look at loadings
PCA.e.rot<-PCA.environment$rotation

# look at relative amount of variation explained
PCA.e.perc<-round(summary(PCA.environment)$importance[2,]*100,1)

rbind(PCA.e.perc,PCA.e.rot)  %>% write.csv(.,"PCA_environment_loadings.csv")

# source(paste(scriptsdir,"/../../observer-free-morphotype-characterization/anderson.R",sep=""))
# anderson(PCA.environment$sdev) 
#If you want to use only some PCs, first two PCs explain significantly more data than others.
#84% of variation, add PCs 4-5 to get 95% of variation

# phenotype PCA ------------------
PCA.p.perc<-round(summary(PCA.phenotype)$importance[2,]*100,1)
PCA.p.2plot<-cbind(metadata,PCA.phenotype$x)
# create PCA figures ------
#link population metadata to environmental PCs
PCA.e.2plot<-cbind(metadata,PCA.environment$x)

ggplot(data=PCA.e.2plot,aes(x=PC1,y=PC2, color = `Sampling Area`,shape=`Genetic Group`))+
  geom_point(size=2) +
  theme_classic() + theme(legend.position="none") +
  xlab(paste("PC 1 (",PCA.e.perc[1],"%)",sep="")) +
  ylab(paste("PC 2 (",PCA.e.perc[2],"%)",sep=""))
ggsave("PCA_env.pdf", device = cairo_pdf, width = 8, height = 8,units="cm",dpi=600)

ggplot(data=PCA.e.2plot,aes(x=PC1,y=PC2, color = `Sampling Area`,shape=`Genetic Group`))+
  geom_point(size=2) 
ggsave("PCA_legend.pdf", device = cairo_pdf, width = 16, height = 16,units="cm",dpi=300)

# put everything together into a master table for plotting shape --------
ggplot(data=PCA.p.2plot,aes(x=PC1,y=PC2, color = `Sampling Area`,shape=`Genetic Group`))+
  geom_point(size=2) +
  theme_classic() + theme(legend.position="none") +
  xlab(paste("PC 1 (",PCA.p.perc[1],"%)",sep="")) +
  ylab(paste("PC 2 (",PCA.p.perc[2],"%)",sep=""))
ggsave("PCA_phen.pdf", device = cairo_pdf, width = 8, height = 8,units="cm",dpi=600)




#look at a 3D PC plot ----- 
# library(plot3D)
# plot3d(x=master.variables$PC1,y=master.variables$PC2,z=master.variables$PC3,
#        col=c("red","green","blue")[as.factor(master.variables$gengroup)],size=10)

# explore relationships between individual variables, out of curiosity ------
# ggplot(data=master.variables, 
#        aes(x=gengroup,y=season_start_time)) +
#   geom_boxplot(aes(fill=gengroup)) +
#   theme(legend.text = element_text(size = 30))

# # compare shape space of Wyoming Basin to others ---------
# #This is a circle-back section: later on you will find out that differentiation of Wyoming Basin
# #is driving some overall patterns in the dataset. Here, look and see if WB appear to occupy unique shape space
# #unoccupied by specimens from other genetic groups, as you'd predict if WB is phenotypically different
# #(vs. simply being highly clustered in a corner of already-occupied morphospace)
# 
# # Find the convex hull of the points being plotted
# hull.not.1 <- filter(master.variables,gengroup!="wb") %>%
#   slice(chull(PC1, PC2))
# 
# hull.ip.1 <- filter(master.variables,gengroup=="ip") %>%
#   slice(chull(PC1, PC2))
# 
# hull.cp.1 <- filter(master.variables,gengroup=="cp") %>%
#   slice(chull(PC1, PC2))
# 
# hull.wb.1 <- filter(master.variables,gengroup=="wb") %>%
#   slice(chull(PC1, PC2))
# 
# # hull.plot.1<-
# ggplot(data=master.variables, aes(x=PC1,y=PC2)) +
#   # geom_polygon(data = hull.not.1, alpha = 0.5, fill="gray80", linetype = 2) +
#     geom_polygon(data = hull.ip.1, alpha = 0.2, fill="gray80", linetype = 2) +
#     geom_polygon(data = hull.cp.1, alpha = 0.2, fill="gray80", linetype = 3) +
#     geom_polygon(data = hull.wb.1, alpha = 0.5,fill="pink",linetype = 1) + 
#     geom_point(aes(color=population,shape=gengroup)) + 
#     theme_classic() + theme(legend.position="none") +
#     xlab(paste("PC 1 (",PCA.p.perc[1],"%)",sep="")) +
#     ylab(paste("PC 2 (",PCA.p.perc[2],"%)",sep=""))
# 
# 
# hull.not.2 <- filter(master.variables,gengroup!="wb") %>%
#   slice(chull(PC3, PC4))
# 
# hull.ip.2 <- filter(master.variables,gengroup=="ip") %>%
#   slice(chull(PC3, PC4))
# 
# hull.cp.2 <- filter(master.variables,gengroup=="cp") %>%
#   slice(chull(PC3, PC4))
# 
# 
# hull.wb.2 <- filter(master.variables,gengroup=="wb") %>%
#   slice(chull(PC3, PC4))
# 
# # hull.plot.2<-
# ggplot(data=master.variables,aes(x=PC3,y=PC4)) +
#   # geom_polygon(data = hull.not.2, alpha = 0.5, fill="gray80", linetype = 2) +
#   geom_polygon(data = hull.ip.2, alpha = 0.2, fill="gray80", linetype = 2) +
#   geom_polygon(data = hull.cp.2, alpha = 0.2, fill="gray80", linetype = 3) +
#   geom_polygon(data = hull.wb.2, alpha = 0.5,fill="pink",linetype = 1) + 
#   geom_point(aes(color=population,shape=gengroup)) + 
#   theme_classic() + theme(legend.position="none") +
#   xlab(paste("PC 3 (",PCA.p.perc[3],"%)",sep="")) +
#   ylab(paste("PC 4 (",PCA.p.perc[4],"%)",sep=""))
#   
#   
# # # library(gridExtra)
# # grid.arrange(hull.plot.1, hull.plot.2, 
# #              ncol = 1, nrow = 2)
# # ggsave("PCA_phen_hull.pdf", device = cairo_pdf, width = 8, height = 16,units="cm",dpi=600)

