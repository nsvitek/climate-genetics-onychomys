##in silico digest of Peromyscus maniculatus genome to get a sense of numbers for ddRAD-seq project.
##Reference genome downloaded in 5 parts from http://www.ncbi.nlm.nih.gov/Traces/wgs/?val=AYHN01

#install.packages("SimRAD")
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
#source("http://bioconductor.org/biocLite.R")
#biocLite("ShortRead")
library(SimRAD)

##Set location of files
folderinput<-"/scratch/lfs/nvitek/in_silico_digest/genome_p_maniculatus"
folderoutput<-"/scratch/lfs/nvitek/in_silico_digest/"
reportname<-paste(folderoutput,"isd_window_reps",sep="/")

##THE GOAL: Aim for 20x coverage, 6 samples per locality,
##each sample represented by Z loci in a single-end run on 1 lane of an
##Illumina NextSeq500 mid-throughput, producing 110-120M reads. 
coverage<-20
samples_per_locality<-6
reads<-110000000

##Set number of replicates and proportion of genome subsampled per replicate
r<-100
portion<-.10

##Set desired window sizes
windowmin<-c(200,300,400,200,250,300,350,400,450)
windowmax<-c(400,500,600,350,400,450,500,550,600)

##Set restriction enzymes
##EcoRI
cs_5p1<-"G"
cs_3p1<-"AATTC"
##MseI
cs_5p2<-"T"
cs_3p2<-"TAA"


##Analysis starts here
sink(file=paste(reportname,"_report.txt",sep=""),append=FALSE)
Sys.Date()

filename<-list.files(folderinput)
q<-length(filename)
windowreps<-matrix(0,nrow=r,ncol=length(windowmin))
for (a in 1:q){
  for (b in 1:round(r/q)){
    rfsq<-ref.DNAseq(paste(folderinput,filename[a],sep="/"),prop.contigs=(q*portion))
    print(paste("Reference sequence",filename[a],"replicate",b,"of",round(r/q),"loaded."))
    print(paste("Replicate",round(r/q)*a+b-round(r/q),"of",round(r/q)*q))
    digest<-insilico.digest(rfsq,cs_5p1,cs_3p1,cs_5p2,cs_3p2,verbose=TRUE)
    aselect<-adapt.select(digest,type="AB+BA",cs_5p1,cs_3p1,cs_5p2,cs_3p2)
    for (c in 1:length(windowmin)){
      sizeselect<-size.select(aselect,min.size=windowmin[c],max.size=windowmax[c],graph=FALSE,verbose=FALSE)
      windowreps[round(r/q)*a+b-round(r/q),c]<-length(sizeselect@ranges@group)
    }
  }
}

sink()

##write a report for each window frame
sink(file=paste(reportname,"_summary.txt",sep=""),append=FALSE)
print(paste("Summary of",r,"replicates of",portion*100,"% subsampling of genome"))
for (f in 1:length(windowmin)){
  print(paste("**window size:",windowmin[f],"-",windowmax[f],"bp**"))
  avg<-mean(windowreps[,f])
  print(paste("mean # of loci:",round(avg)))
  print(paste("standard deviation:",sd(windowreps[,f])))
  localities<-round(reads/(coverage*samples_per_locality*avg))
  print(paste(localities," localities possible if sampling ",samples_per_locality-1,
     " critters with 1 replicate per locality and ",coverage,"x coverage on one lane achieving ",reads
     ," total reads.",sep=""))
}
colnames(windowreps)<-print(paste(windowmin,"-",windowmax,"bp",sep=""))
write.csv(windowreps,file=paste(reportname,".csv",sep=""),row.names=FALSE,col.names=TRUE)

sink()