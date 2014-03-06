library(methods)
library(lattice)
suppressMessages(library(Matrix))
suppressMessages(library(ClassComparison))
suppressMessages(library(glmnet))
options(warn=-1)

#######################
### arguments input ###
#######################

arg=commandArgs(TRUE)

InputFile=arg[1]          # must be given by user including both 
                          # response variable and siRNA sequences

Seed=as.integer(arg[2])   # Name: "Seed"
                          # specify the seed region to be used in analysis 
                          # an interger out of 1~14.
                          # default is 2, which means 2~7 hexamer is used 

Strand=arg[3]             # Name:"Strand Orientation for Analysis"
                          # specify which strand is used in analysis 
                          # must be given by user, options are "sense", 
                          #"antisense" or "both"

Lambda=as.numeric(arg[4]) # Name:"Lambda Value"
                          # penalty parameter in LASSO regression 
                          # default = 0.001, which can be changed by user
                          # in order to perform differennt LASSO estimation

Library=arg[5]            # Name:"siRNA Library to use"
                          # user has to specify which siRNA library will be 
                          # used in analysis can be "Custom", "Ambion", "New Dharmacon"
                          # or "Old Dharmacon". Default is "Custom"

### new parameter ###

Str=as.numeric(arg[6])    # Name:"Strength of seed-linked effect“
                          # specify cutoff for strenth of seed-linked effet, must be positive
                          # default is 1

Sig=as.numeric(arg[7])    # Name:"Significance (P value)"
                          # specify cutoff for significance (p value)
                          # default is 0.01
                           
#################
### functions ###
#################

# uppercase a sequence and transform T to U ###

TtoU=function(x){
  x=as.character(x)
  x=toupper(x)
  x=chartr("T", "U", x)
  return(x)
}

# reverse a sequence #

StrReverse <- function(x){
  x=as.character(x)
  x=sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
  return(x)
} 

# generate complementary sequence #

DuplexReverse=function(x)
{
  x=as.character(x)
  x=gsub("U","T",x)
  x=gsub("A","U",x)
  x=gsub("T","A",x)
  x=gsub("G","g",x)
  x=gsub("C","G",x)
  x=gsub("g","C",x)
  x=StrReverse(x)
  return(x)
}

# check if input contains AUGC only #

SeqCheck=function(x){
  x=as.character(x)
  x=strsplit(x, "")
  x=unlist(x)
  if(sum(x %in% c("A", "U", "G", "C")) == length(x) ){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

# extract seed sequence from siRNA #

SeedExtract=function(x, start=2, end=7){
  x=as.character(x)
  x=substr(x, start, end)
  return(x)
}

# examine if a pool contain a seed #

SeedExist=function(dfr, s){
  dfr=data.frame(dfr)
  s=as.character(s)
  return(as.numeric(apply(dfr == s, 1, sum)))
}

# seed analysis #

SeedAnalysis=function(Z, data, start=2, end=7){
  
  print("Generating seed.dfr ... ... ")
  data=data.frame(data)
  seed.dfr=data
  for(i in 1:ncol(seed.dfr)){
    seed.dfr[, i]=as.character(SeedExtract(seed.dfr[, i], start=start, end=end))
  }
  colnames(seed.dfr)=paste("seed", 1:ncol(seed.dfr), sep=".")
  print("Generating seed.table ... ... ")
  seed.table=data.frame(table(unlist(seed.dfr)))
  names(seed.table)=c("seed", "family.size")
  
  seed.table$ks.p.value=rep(NA, length(seed.table$seed))
  y=Z
  for(i in 1:length(seed.table$seed)){
    s=SeedExist(seed.dfr, seed.table$seed[i])
    seed.table$ks.p.value[i]=
      ks.test(y[s > 0], y[s== 0])$p.value
    print(paste("seed ", i))
  }
  
  return(list(Z=Z, data=data, seed.table=seed.table, seed.dfr=seed.dfr))
}

# offtarget analysis #

OfftargetAnalysis=function(object, s=0.001){
  
  seed.dfr=object$seed.dfr
  seed.table=object$seed.table
  y=object$Z
  
  print("Generating Matrix X ... ... ")
  s.sub=as.character(seed.table$seed)
  x=matrix(0, length(y), length(s.sub))
  colnames(x)=as.character(s.sub)
  
  for(i in 1:length(s.sub)){
    n=as.character(s.sub[i])
    x[, n]=SeedExist(seed.dfr, n)
    print(paste("Matrix X seed ", i))
  }
  
  print("Fitting ... ... ")
  gc()
  
  fit=glmnet(x, y)
  
  print("offtarget analysis")
  offtarget=data.frame(as.matrix(coef(fit, s=s)))
  offtarget$seed=rownames(offtarget)
  offtarget=offtarget[, c(2, 1)]
  names(offtarget)[2]="coef"
  rownames(offtarget)=NULL
  
  offtarget=merge(offtarget, seed.table, by="seed", all.x=T)
  y.fit=as.numeric(predict(fit, x, s=s))
  y.corrected=y-y.fit
  print("Off-target analysis Done!")
  
  return(list(fit=fit, X=x, offtarget=offtarget, Z.corrected=y.corrected))
  
  
}

# generating data frame of siRNA sequence #

SeqGenerate=function(SeqDat){
  SenseSeq=data.frame(apply(SeqDat, 2, TtoU))
  AntiSeq=data.frame(apply(SenseSeq, 2, DuplexReverse))
  for(i in 1:ncol(SenseSeq)){
    AntiSeq[, i]=DuplexReverse(SenseSeq[, i])
  }
  BothSeq=cbind(SenseSeq, AntiSeq)
  names(BothSeq)=paste("duplex", 1:ncol(BothSeq), sep=".")
  output=list(SenseSeq=SenseSeq, AntiSeq=AntiSeq, BothSeq=BothSeq)
  return(output)
}

# DataAnalysis #

DataAnalysis=function(Z, Seq, SeedStart, Strand, Lambda){
  
  SeedEnd=SeedStart+5
  
  analysis=list()
  
  if(Strand == 'sense'){
    Seq=Seq$SenseSeq
    analysis$seed.analysis=SeedAnalysis(Z=Z, data=Seq, start=SeedStart, end=SeedEnd)
    analysis$offtarget.analysis=OfftargetAnalysis(analysis$seed.analysis, s=Lambda)
  }
  else if(Strand == "antisense"){
    Seq=Seq$AntiSeq
    analysis$seed.analysis=SeedAnalysis(Z=Z, data=Seq, start=SeedStart, end=SeedEnd)
    analysis$offtarget.analysis=OfftargetAnalysis(analysis$seed.analysis, s=Lambda)
  }
  else if(Strand == 'both'){
    Seq=Seq$BothSeq
    analysis$seed.analysis=SeedAnalysis(Z=Z, data=Seq, start=SeedStart, end=SeedEnd)
    analysis$offtarget.analysis=OfftargetAnalysis(analysis$seed.analysis, s=Lambda)
  }
  return(analysis)
}

# new functions #

# FDR calculation # 

fdr.bum=function(p.value, Sig){
  p.sum=attributes(summary(Bum(p.value), tau=Sig))$estimates
  f=p.sum[4]/(p.sum[2]+p.sum[4])
  f=as.numeric(f)
  return(f)
}

# load miRNA database #

load("/home/galaxy/galaxy-dist/tools/offtarget/miRNA.RData")

##################
### data input ###
##################

# data input #

SeedStart=Seed
inputfile=read.csv(InputFile, header=F)
id=inputfile[, 1]
Z=as.numeric(inputfile[, 2]) 

# Z score checking #

if(sum(is.na(Z)) > 0){
  print("Missing Value in Experimental Readouts!")
  quit()
}

if(sum(Z == "") > 0){
  print("Missing Value in Experimental Readouts!")
  quit()
}

#####################
### data analysis ###
#####################

### Custom ###

if(Library == "custom"){
  
  ### check if user selected wrong "Library" parameter ###
  if(ncol(inputfile) < 3){
    print("siRNA sequences not available and please check your Library selection")
    quit()
  }
  ##################################################################
  
  SeqDat=inputfile[, 3:ncol(inputfile)] 
  SeqDat=data.frame(SeqDat)
  
  if(sum(is.na(unlist(SeqDat))) > 0){
    print("Missing Value in siRNA sequences!")
    quit()
  }
  
  if(sum(unlist(SeqDat) == "") > 0){
    print("Missing Value in siRNA sequences!")
    quit()
  }
  
  Seq=SeqGenerate(SeqDat)
  
  SenseSeq=Seq$SenseSeq
  Check=data.frame(apply(SenseSeq, 1, SeqCheck))
  
  if(sum(Check == FALSE) > 0){
    print("Non-Sequence Value in siRNA sequences!")
    quit()
  }
  else{
    analysis=DataAnalysis(Z, Seq, SeedStart, Strand, Lambda)
  }
}

### Ambion ###

if(Library == "ambion"){
  
  id=as.integer(id)
  
  if(sum(is.na(id)) > 0){
    print("Missing Value in Gene ID!")
    quit()
  }
  
  if(sum(id == "") > 0){
    print("Missing Value in Gene ID!")
    quit()
  }
  
  load("/home/galaxy/galaxy-dist/tools/offtarget/Ambion.RData")
  id.ambion=unique(ambion$id)
  
  if(length(id[!(id %in% id.ambion)]) > 0){
    write.csv(c("No Matching siRNA Sequences with Below Gene ID and Removed from Analysis",
                id[!(id %in% id.ambion)]), "unmatched.id.csv")
  }

  Dat=data.frame(id=id, z=Z)
  Dat=merge(Dat, ambion, by="id", sort=F)
  Z=Dat$z
  id=Dat$id
  SeqDat=Dat[, 6:8]
  Seq=SeqGenerate(SeqDat)
  analysis=DataAnalysis(Z, Seq, SeedStart, Strand, Lambda)
  
}

### New Dhar ###

if(Library == "new_dharmacon"){
  
  id=as.integer(id)
  
  if(sum(is.na(id)) > 0){
    print("Missing Value in Gene ID!")
    quit()
  }
  
  if(sum(id == "") > 0){
    print("Missing Value in Gene ID!")
    quit()
  }
  
  load("/home/galaxy/galaxy-dist/tools/offtarget/New_Dhar.RData")
  id.new.dhar=unique(new.dhar$id)
  
  if(length(id[!(id %in% id.new.dhar)]) > 0){
    write.csv(c("No Matching siRNA Sequences with Below Gene ID and Removed from Analysis",
                id[!(id %in% id.new.dhar)]), "unmatched.id.csv")
  }
  
  Dat=data.frame(id=id, z=Z)
  Dat=merge(Dat, new.dhar, by="id", sort=F)
  Z=Dat$z
  id=Dat$id
  SeqDat=Dat[, 7:10]
  Seq=SeqGenerate(SeqDat)
  analysis=DataAnalysis(Z, Seq, SeedStart, Strand, Lambda)
  
}

### Old Dhar ###

if(Library == "old_dharmacon"){
  
  id=as.integer(id)
  
  if(sum(is.na(id)) > 0){
    print("Missing Value in Gene ID!")
    quit()
  }
  
  if(sum(id == "") > 0){
    print("Missing Value in Gene ID!")
    quit()
  }
  
  load("/home/galaxy/galaxy-dist/tools/offtarget/Old_Dhar.RData")
  id.old.dhar=unique(old.dhar$id)
  
  if(length(id[!(id %in% id.old.dhar)]) > 0){
    write.csv(c("No Matching siRNA Sequences with Below Gene ID and Removed from Analysis",
                id[!(id %in% id.old.dhar)]), "unmatched.id.csv")
  }
  
  Dat=data.frame(id=id, z=Z)
  Dat=merge(Dat, old.dhar, by="id", sort=F)
  Z=Dat$z
  id=Dat$id
  SeqDat=Dat[, 7:10]
  Seq=SeqGenerate(SeqDat)
  analysis=DataAnalysis(Z, Seq, SeedStart, Strand, Lambda)
  
}

######################
### results output ###
######################

# code chagned below #

# generating output files #

duplex.dfr=analysis$seed.analysis$data 
names(duplex.dfr)=paste("duplex", 1:length(names(duplex.dfr)), sep=".") 
seed.dfr=analysis$seed.analysis$seed.dfr 
offtarget=analysis$offtarget.analysis$offtarget 
Z.corrected=analysis$offtarget.analysis$Z.corrected
output.corr=cbind(ID=id, Z.Score=Z, Corrected.Z.Score=Z.corrected, 
                  duplex.dfr, seed.dfr)


print("outputing results... ...")
 
if(Library == 'custom'){
  
  print("outputing corrected Z score... ...")
  write.csv(output.corr, "CorrectedZScore.csv", row.names=F)   ### export corrected z score 
  
  print("outputing summary... ...")
  write.csv(offtarget, "Seed_Famliy_Summary.csv", row.names=FALSE)   ### output seed family summary table
  
  print("outputing offtarget seed families... ...")
  
  off.target=offtarget$seed[offtarget$ks.p.value < Sig & (  
    offtarget$coef < -Str | offtarget$coef > Str)] # select off target by custom setup
  
  write.csv(offtarget[offtarget$seed %in% off.target, ], 
            "Off-Target_Seed_Families.csv", row.names=FALSE)   ### export off-target seed families 
  
  print("outputing miRNAs... ...")
  if(nrow(miRNA[miRNA$Seed.sequence %in% off.target, ]) ==  0){
    write.csv("No miRNA detected", "miRNA.csv", row.names=F) 
    fe="No"
  }   ### export miRNAs table sharing common seed
  
  if(nrow(miRNA[miRNA$Seed.sequence %in% off.target, ]) != 0){
    write.csv(miRNA[miRNA$Seed.sequence %in% off.target, ], "miRNA.csv", row.names=F) 
    fe="Yes"
  }   ### export miRNAs talbe sharing common seed 
  
  if(length(off.target) > 0){
    for(i in 1:length(off.target)){
      
      n=off.target[i]
      si=SeedExist(seed.dfr, n)
      si.offtarget=output.corr[si == 1, -3]
      si.offtarget$off.target.seed=n
      si.offtarget=si.offtarget[order(si.offtarget$Z.Score), ]
      print(paste(n, "siRNA pools ... ..."))
      write.csv(si.offtarget, paste("siRNAsPool_OffTargetSeedFamily_", n,  ".csv",  sep=""), row.names=F)
      
    }
  } 
}else{
  
  print("outputing corrected Z score... ...")
  write.csv(output.corr[, 1:3], "CorrectedZScore.csv", row.names=F)
  
  print("outputing summary... ...")
  write.csv(offtarget, "Seed_Famliy_Summary.csv", row.names=FALSE)   ### output seed family summary table
  
  print("outputing offtarget seed families... ...")
  
  off.target=offtarget$seed[offtarget$ks.p.value < Sig & (  
    offtarget$coef < -Str | offtarget$coef > Str)] # select off target by custom setup
  
  write.csv(offtarget[offtarget$seed %in% off.target, ], 
            "Off-Target_Seed_Families.csv", row.names=FALSE)   ### export off-target seed families 
  
  print("outputing miRNAs... ...")
  if(nrow(miRNA[miRNA$Seed.sequence %in% off.target, ]) ==  0){
    write.csv("No miRNA detected", "miRNA.csv", row.names=F) 
    fe="No"
  }   ### export miRNAs table sharing common seed
  
  if(nrow(miRNA[miRNA$Seed.sequence %in% off.target, ]) != 0){
    write.csv(miRNA[miRNA$Seed.sequence %in% off.target, ], "miRNA.csv", row.names=F) 
    fe="Yes"
  }   ### export miRNAs talbe sharing common seed 
  
  if(length(off.target) > 0){
    for(i in 1:length(off.target)){
      
      n=off.target[i]
      si=SeedExist(seed.dfr, n)
      si.offtarget=output.corr[si == 1, -3]
      si.offtarget$off.target.seed=n
      si.offtarget=si.offtarget[order(si.offtarget$Z.Score), ]
      print(paste(n, "siRNA pools ... ..."))
      write.csv(si.offtarget, paste("siRNAsPool_OffTargetSeedFamily_", n,  ".csv",  sep=""), row.names=F)
      
    }
  }
}  


# generating output figures # 

print("outputing figures... ...")

jpeg(file="SeedFamilyOffTarget.jpeg", height=600, width=800) 
par(mar=c(5, 5, 2, 2))
plot(offtarget$coef, -log10(offtarget$ks.p.value), 
     pch=20, col="blue", yaxt="n", cex.axis=1.5, cex.lab=2, 
     xlab="Strength of seed-linked effect", ylab="Significance (P value)")
p.log=ceiling(-log10(min(offtarget$ks.p.value, na.rm=T))/2)
p.range=10^(0:-p.log*2)
axis(2, p.range, 
     at=-log(p.range, 10),
     cex.axis=1.5)
off.target=offtarget$seed[offtarget$coef < -Str & offtarget$ks.p.value < Sig]
nn=length(off.target)
points(offtarget$coef[offtarget$seed %in% off.target], 
      -log10(offtarget$ks.p.value[offtarget$seed %in% off.target]), 
      col="red", pch=20)
off.target=offtarget$seed[offtarget$coef > Str & offtarget$ks.p.value < Sig]
np=length(off.target)
points(offtarget$coef[offtarget$seed %in% off.target], 
       -log10(offtarget$ks.p.value[offtarget$seed %in% off.target]), 
       col="green", pch=20)
abline(v=-Str, col="red", lwd=3)
abline(v=Str, col="green", lwd=3)
abline(h=-log10(Sig), col="blue", lwd=3)
dev.off()

# generating legend #

jpeg(file="SeedFamilyOffTargetLegend.jpeg", width=500, height=100) 
par(mar=c(0, 0, 0, 0))
plot(c(0, 2), c(0.8, 1.2), type="n", bty="n", yaxt="n", xaxt='n')
points(c(0, 0), c(1.1, 0.9), col=c('red', "green"), pch=20, cex=2)
text(1, 1.1, paste('N: ', nn, " (Strength < ", -Str, ", P value < ", Sig, ")", sep=""), cex=2)
text(1, 0.9, paste('N: ', np, " (Strength > ", Str, ",  P value < ", Sig, ")", sep=""), cex=2)
dev.off()

# new output #

# generating summary #

print("Ansysis Summary ... ...")

print("p value ... ...")
p.value=offtarget$ks.p.value[offtarget$coef < -Str | offtarget$coef > Str]

print("FDR ... ...")
f=fdr.bum(p.value, Sig)
print("rounding ... ...")
f=round(f, digits=3)
print(f)
print("FDR done ... ...")


print("off.sum ... ...")
off.sum=matrix(rep("", 30), 15, 2)
off.sum[1,1]="Analysis is successfully done!"
off.sum[3,1]="Parameters set up"
off.sum[4:15, 1]=c("Strand Orientation for Analysis:", "Lambda value:", "Seed:", 
                   "siRNA Library:", "Strength of seed-linked effect:", "P value:", "",
                   "Analysis summary",
                   "Number of negative off-target seed families:", 
                   "Number of positive off-target seed families:", "FDR:", 
                   "Identify miRNA with phenotypical effects:")
off.sum[4:15, 2]=c(Strand, Lambda, Seed, Library, Str, Sig, "", "", nn, np, f, fe)
write.table(off.sum, "Analysis_Summary.csv", row.names=F, col.names=F, sep=",")

print("Analysis done ... ...")

