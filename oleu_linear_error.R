# functions, dependencies -------
library(dplyr)
library(reshape2)

PME<-function(ANOVA,r){
  # Yezerinac et al. 1992 p. 474 % measurement error
  s.within<-ANOVA$`Mean Sq`[2]
  s.among<-(ANOVA$`Mean Sq`[1]-s.within)/r
  percent_measurement_error<-(s.within/(s.within+s.among))
  return(percent_measurement_error)
}

# read, format data --------
in_df<-dplyr::select(raw_metadata,-c(institutioncode, catalognumber, vial, level, stateprovince, county, 
                                 batch, packing_notes,repacked,is_sequenced, start, end, population, isolated,
                                 crop_lo,crop_hi))

#step 1 in trying to make data "tidy"? 
in_df_long<-melt(in_df,id=c("filename")) 

#unpack variables from label

#break up variable for faceting
in_df_long$tooth<-gsub("(m[1-3]).*","\\1",in_df_long$variable)
in_df_long$orientation<-gsub("m[1-3]_(length|width)_[1-3]","\\1",in_df_long$variable)
in_df_long$replicate<-gsub("m[1-3]_(length|width)_([1-3])","\\2",in_df_long$variable)

# calculate repeatibility ----
#The Roseman & Delezene method for evaluating "within-species repeatibility" (p.233)
#is identical to Yezerinac's equation for percent measurement error, contained in function.

#r is number of repeated measurements per variable. Here, each measurement taken 3 times.
r<-3                               

#make empty matrix to hold percents, rows are tooth positions, columns 
levels.row<-unique(in_df_long$tooth)
levels.col<-unique(in_df_long$orientation)
percents<-matrix(NA,nrow=length(levels.row),
                 ncol=length(levels.col),
                 dimnames=list(levels.row,levels.col))

#fill in each cell in the matrix with a nest of for-loops
for (row in levels.row){
  for (col in levels.col){
    ANOVA<-anova(lm(value~filename, data = filter(in_df_long,
                                                  tooth==row,orientation==col)))
    percents[row,col]<-PME(ANOVA,r)
    
  }
}

#take the average repeatibility for each metric across the tooth row (avg of m1 thru m3)
avg.repeat<-apply(percents,2,function(x) 1-mean(x))

write.csv(percents,paste(datadir,"output/linear_measurement_repeatibility.csv",sep="/"))
#clean workspace
rm(list=c("ANOVA","in_df","in_df_long","col","levels.col","levels.row","r","row"))
