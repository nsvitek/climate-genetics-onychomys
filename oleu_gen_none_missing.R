#Make dataset for PCA, PLS
#handle missing data by replacing with either mean allele frequencies or zero

Y = x
Y@tab<-tab(Y,NA.method="zero")
write.csv(Y,paste(datadir,gen_path_name,'/snps.missing.interpolated.0.csv',sep=""))

Y = x
Y@tab<-tab(Y,NA.method="mean")
write.csv(Y,paste(datadir,gen_path_name,'/snps.missing.interpolated.mean.csv',sep=""))

rm("Y")