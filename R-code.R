################################################################################################
############################Collapse tables#####################################################

##Setworking directory
setwd('/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB')

#Read tables
PmA=read.table('almost/Primary-ATAC+CDS.txt',sep="\t", header = F)
rownames(PmA)= as.character(PmA[,1])
PpA=read.table('almost/Primary+ATAC+CDS.txt',sep="\t", header = F)
rownames(PpA)= as.character(PpA[,1])
NmA=read.table('almost/NonPrimary-ATAC+CDS.txt',sep="\t", header = F)
rownames(NmA)= as.character(NmA[,1])
NpA=read.table('almost/NonPrimary+ATAC+CDS.txt',sep="\t", header = F)
rownames(NpA)= as.character(NpA[,1])


OP=ngenes
OP[,9]=(as.numeric(OP[,7]) > 1)
colnames(OP)[9]="InMiddleOperon"

PmA[,ncol(PmA)+1]=OP[rownames(PmA),9]
PpA[,ncol(PpA)+1]=OP[rownames(PpA),9]
NmA[,ncol(NmA)+1]=OP[rownames(NmA),9]
NpA[,ncol(NpA)+1]=OP[rownames(NpA),9]

PmA[,ncol(PmA)+1]=rep("PrimaryNoATAC",nrow(PmA))
PpA[,ncol(PpA)+1]=rep("PrimaryWithATAC",nrow(PpA))
NmA[,ncol(NmA)+1]=rep("NotPrimaryNoATAC",nrow(NmA))
NpA[,ncol(NpA)+1]=rep("NotPrimaryWithATAC",nrow(NpA))

ALL=rbind(PpA,PmA,NpA,NmA)
colnames(ALL)=c('WormBaseID','locus','transcript','GeneBody-CDS-Chr','GeneBody-CDS-0based-Start','GeneBody-CDS-0based-End',
                'GeneBody-Name','GeneBody-Score','GeneBody-Strand','TSSregion-0based-Chr','TSSregion-0based-Start',
                'TSSregion-0based-End','TSSregion-Info','TSSregion-Score','TSSregion-Strand','1stSharedExon-0based-Chr',
                '1stSharedExon-0based-Start','1stSharedExon-0based-End','1stSharedExon-Isoform','1stSharedExon-Score',
                '1stSharedExon-0based-Strand','GeneInsideOperon','Source')

write.table(ALL,'OligoDB_WormBase268-Ensemble95.txt',sep='\t', quote = F, row.names = F)
###############################################################################################


###############################################################################################
#############################Merge Tracks######################################################

##New directory, append,get data
setwd('/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/Rscript')

library("Biostrings")
library("adagio")

##Primers
FwdFile <- readDNAStringSet("/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/OligoDB-primers-Forward.fasta")
seq_name = names(FwdFile)
sequence = paste(FwdFile)
Fwd <- data.frame(seq_name, sequence)

RevFile <- readDNAStringSet("/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/OligoDB-primers-Reverse.fasta")
seq_name = names(RevFile)
sequence = paste(RevFile)
Rev <- data.frame(seq_name, sequence)


##Do for each set, e.g.:
#CEN45ALL Set1Fasta

fastaFile <- readDNAStringSet("/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/normal/perProm/CEN45.fasta")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
id = type = loc = prom = c()
for(name in c(1:length(seq_name))){
  name=seq_name[name]
  tmp=unlist(strsplit(name,"_"))
  id=c(id,tmp[1])
  type=c(type,tmp[2])
  loc=c(loc,unlist(strsplit(tmp[4],"-"))[1])
  prom=c(prom,tmp[length(tmp)])
}
Set <- data.frame(id,type,loc,prom, sequence)

#Order set
tab=data.frame(table(Set[1]))
ktab=tab
con=1
groups=list()
NgRNA=c()
Ngen=c()
while(nrow(ktab)>1){
  knap=knapsack(ktab[,2],rep(1,nrow(ktab)),120)
  groups=append(groups,list(ktab[knap$indices,1]))
  NgRNA=c(NgRNA,knap$capacity)
  Ngen=c(Ngen,knap$profit)
  ktab=ktab[-knap$indices,]
  con = con+1; if(con>1000){print; 'No solution'; exit}
}



###Add primers
##Set boundaries
Tsslim=c(1:22)
ATGlim=c(23:24)
CDSlim=c(25:29)
TAGlim=c(30:31)
Regions=c("TSSRegion","GeneStartbed","1stCDSbed","GeneEndbed")

Wellsnames=list()
Wellseqs=list()
fastanames=c()
fastasequences=c()
keyid=c()

#PerEachgroupWriteSequencesofWell
for(i in 1:length(groups)){
  gcount=1;
  fastanames=c()
  fastasequences=c()
  for(keyid in as.character(groups[[i]])){
    tabtemp=Set[which(as.character(Set[,1])==keyid),]
    for(kind in 1:length(Regions)){
      sgcount=1;
      for(last in which(tabtemp[,3]==Regions[kind])){
        if(kind==1){
        fname=c(paste(">",tabtemp[last,1],"_",tabtemp[last,3],"-",sgcount,"_",tabtemp[last,4],"_",Fwd[gcount,1],"_",Rev[Tsslim[sgcount],1],sep=""));
        fseq=paste(paste(FwdFile[gcount]),tabtemp[last,5],paste(reverseComplement(RevFile[Tsslim[sgcount]])),sep='');
        fastanames=c(fastanames,fname)
        fastasequences=c(fastasequences,fseq)
        sgcount=sgcount+1
        }
        if(kind==2){
          fname=c(paste(">",tabtemp[last,1],"_",tabtemp[last,3],"-",sgcount,"_",tabtemp[last,4],"_",Fwd[gcount,1],"_",Rev[ATGlim[sgcount],1],sep=""));
          fseq=paste(paste(FwdFile[gcount]),tabtemp[last,5],paste(reverseComplement(RevFile[ATGlim[sgcount]])),sep='');
          fastanames=c(fastanames,fname)
          fastasequences=c(fastasequences,fseq)
          sgcount=sgcount+1
          }
        if(kind==3){
          fname=c(paste(">",tabtemp[last,1],"_",tabtemp[last,3],"-",sgcount,"_",tabtemp[last,4],"_",Fwd[gcount,1],"_",Rev[CDSlim[sgcount],1],sep=""));
          fseq=paste(paste(FwdFile[gcount]),tabtemp[last,5],paste(reverseComplement(RevFile[CDSlim[sgcount]])),sep='');
          fastanames=c(fastanames,fname)
          fastasequences=c(fastasequences,fseq)
          sgcount=sgcount+1
          }
        if(kind==4){
          fname=c(paste(">",tabtemp[last,1],"_",tabtemp[last,3],"-",sgcount,"_",tabtemp[last,4],"_",Fwd[gcount,1],"_",Rev[TAGlim[sgcount],1],sep=""));
          fseq=paste(paste(FwdFile[gcount]),tabtemp[last,5],paste(reverseComplement(RevFile[TAGlim[sgcount]])),sep='');
          fastanames=c(fastanames,fname)
          fastasequences=c(fastasequences,fseq)
          sgcount=sgcount+1
          }
      }
    }
    gcount=gcount+1;  
  }
  Wellsnames=append(Wellsnames,list(fastanames))
  Wellseqs=append(Wellseqs,list(fastasequences))
}


##Next sets
fastasets=c("/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/normal/perProm/CEN50-1_crRNA.fasta",
            "/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/normal/perProm/CEN50-1.fasta",
            "/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/normal/perProm/CEN55.fasta")



for(finalfasta in fastasets){

fastaFile <- readDNAStringSet(finalfasta)
seq_name = names(fastaFile)
sequence = paste(fastaFile)
id = type = loc = prom = c()
for(name in seq_name){
  tmp=unlist(strsplit(name,"_"))
  id=c(id,tmp[1])
  type=c(type,tmp[2])
  loc=c(loc,unlist(strsplit(tmp[4],"-"))[1])
  prom=c(prom,tmp[length(tmp)])
}
Set <- data.frame(id,type,loc,prom, sequence)

#Order set
tab=data.frame(table(Set[1]))
ktab=tab
con=1
groups=list()
NgRNA=c()
Ngen=c()
while(nrow(ktab)>1){
  knap=knapsack(ktab[,2],rep(1,nrow(ktab)),120)
  groups=append(groups,list(ktab[knap$indices,1]))
  NgRNA=c(NgRNA,knap$capacity)
  Ngen=c(Ngen,knap$profit)
  ktab=ktab[-knap$indices,]
  con = con+1; if(con>1000){print; 'No solution'; exit}
}


#PerEachgroupWriteSequencesofWell
for(i in 1:length(groups)){
  gcount=1;
  fastanames=c()
  fastasequences=c()
  for(keyid in groups[[i]]){
    tabtemp=Set[which(Set[1]==keyid),]
    for(kind in 1:length(Regions)){
      sgcount=1;
      for(last in which(tabtemp[,3]==Regions[kind])){
        if(kind==1){
          fname=c(paste(">",tabtemp[last,1],"_",tabtemp[last,3],"-",sgcount,"_",tabtemp[last,4],"_",Fwd[gcount,1],"_",Rev[Tsslim[sgcount],1],sep=""));
          fseq=paste(paste(FwdFile[gcount]),tabtemp[last,5],paste(reverseComplement(RevFile[Tsslim[sgcount]])),sep='');
          fastanames=c(fastanames,fname)
          fastasequences=c(fastasequences,fseq)
          sgcount=sgcount+1
        }
        if(kind==2){
          fname=c(paste(">",tabtemp[last,1],"_",tabtemp[last,3],"-",sgcount,"_",tabtemp[last,4],"_",Fwd[gcount,1],"_",Rev[ATGlim[sgcount],1],sep=""));
          fseq=paste(paste(FwdFile[gcount]),tabtemp[last,5],paste(reverseComplement(RevFile[ATGlim[sgcount]])),sep='');
          fastanames=c(fastanames,fname)
          fastasequences=c(fastasequences,fseq)
          sgcount=sgcount+1
        }
        if(kind==3){
          fname=c(paste(">",tabtemp[last,1],"_",tabtemp[last,3],"-",sgcount,"_",tabtemp[last,4],"_",Fwd[gcount,1],"_",Rev[CDSlim[sgcount],1],sep=""));
          fseq=paste(paste(FwdFile[gcount]),tabtemp[last,5],paste(reverseComplement(RevFile[CDSlim[sgcount]])),sep='');
          fastanames=c(fastanames,fname)
          fastasequences=c(fastasequences,fseq)
          sgcount=sgcount+1
        }
        if(kind==4){
          fname=c(paste(">",tabtemp[last,1],"_",tabtemp[last,3],"-",sgcount,"_",tabtemp[last,4],"_",Fwd[gcount,1],"_",Rev[TAGlim[sgcount],1],sep=""));
          fseq=paste(paste(FwdFile[gcount]),tabtemp[last,5],paste(reverseComplement(RevFile[TAGlim[sgcount]])),sep='');
          fastanames=c(fastanames,fname)
          fastasequences=c(fastasequences,fseq)
          sgcount=sgcount+1
        }
      }
    }
    gcount=gcount+1;  
  }
  Wellsnames=append(Wellsnames,list(fastanames))
  Wellseqs=append(Wellseqs,list(fastasequences))
}

}



##Now for other sets
fastasets=c(
  "/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/egl-23.CEN45.fasta",
  "/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/egl-23.CEN50-1_crRNA.fasta",
  "/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/egl-23.CEN50-1.fasta",
  "/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/egl-23.CEN55.fasta",
  "/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/his-72.CEN45.fasta",
  "/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/his-72.CEN50-1_crRNA.fasta",
  "/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/his-72.CEN50-1.fasta",
  "/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/his-72.CEN55.fasta",
  "/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/twk-18.CEN45.fasta",
  "/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/twk-18.CEN50-1_crRNA.fasta",
  "/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/twk-18.CEN50-1.fasta",
  "/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/twk-18.CEN55.fasta"
  )
  



for(finalfasta in fastasets){
  
  fastaFile <- readDNAStringSet(finalfasta)
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  fastanames=c()
  fastasequences=c()
  count=1
##KindOfMatrix20*6
for(i in c(1:6)){
  for(j in c(1:20)){
    if(count > length(seq_name)){next}
    fname=c(paste(">",seq_name[count],"_",count,"_",Fwd[i,1],"_",Rev[j,1],sep=""));
    fseq=paste(paste(FwdFile[i]),sequence[count],paste(reverseComplement(RevFile[j])),sep='');
    fastanames=c(fastanames,fname)
    fastasequences=c(fastasequences,fseq)
    count=count+1
  }
}
  Wellsnames=append(Wellsnames,list(fastanames))
  Wellseqs=append(Wellseqs,list(fastasequences))
}


#Randombinding
fastasets=c(
"/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/extra/Random-binding-set1.fasta",
"/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/extra/Random-binding-set2.fasta"
)


for(finalfasta in fastasets){
  
  fastaFile <- readDNAStringSet(finalfasta)
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  fastanames=c()
  fastasequences=c()
  count=1
  ##KindOfMatrix20*6
  for(i in c(1:6)){
    for(j in c(1:20)){
      if(count > length(seq_name)){next}
      fname=c(paste(">",seq_name[count],"_",count,"_",Fwd[i,1],"_",Rev[j,1],sep=""));
      fseq=paste(paste(FwdFile[i]),sequence[count],paste(reverseComplement(RevFile[j])),sep='');
      fastanames=c(fastanames,fname)
      fastasequences=c(fastasequences,fseq)
      count=count+1
    }
  }
  Wellsnames=append(Wellsnames,list(fastanames))
  Wellseqs=append(Wellseqs,list(fastasequences))
}



###Repairs
fastasets=c(
"/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/HomsArmsCassete/Repairs-1.fasta",
"/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/HomsArmsCassete/Repairs-2.fasta",
"/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/HomsArmsCassete/Repairs-3.fasta"
)

for(finalfasta in fastasets){
  
  fastaFile <- readDNAStringSet(finalfasta)
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  fastanames=c()
  fastasequences=c()
  count=1
  ##KindOfMatrix20*6
  for(i in c(1:6)){
    for(j in c(1:20)){
      if(count > length(seq_name)){next}
      fname=c(paste(">",seq_name[count],"_",count,"_",Fwd[i,1],"_",Rev[j,1],sep=""));
      fseq=paste(paste(FwdFile[i]),sequence[count],paste(reverseComplement(RevFile[j])),sep='');
      fastanames=c(fastanames,fname)
      fastasequences=c(fastasequences,fseq)
      count=count+1
    }
  }
  Wellsnames=append(Wellsnames,list(fastanames))
  Wellseqs=append(Wellseqs,list(fastasequences))
}

##Add random no binding in spots left with space
Well2fill=which(unlist(lapply(Wellsnames,length)) != 120)
#6,6,6
#3,7,37
#47
##Fill with no binding
#/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/extra/Random_setGRNA_NObinding.CEN45.fasta
#/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/extra/Random_setGRNA_NObinding.CEN50-1_crRNA.fasta
#/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/extra/Random_setGRNA_NObinding.CEN50-1.fasta
#/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/extra/Random_setGRNA_NObinding.CEN55.fasta
Nbin45 <- readDNAStringSet("/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/extra/Random_setGRNA_NObinding.CEN45.fasta")
Nbin50cr <- readDNAStringSet("/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/extra/Random_setGRNA_NObinding.CEN50-1_crRNA.fasta")
Nbin55 <- readDNAStringSet("/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/extra/Random_setGRNA_NObinding.CEN55.fasta")
Nbin50 <- readDNAStringSet("/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/extra/Random_setGRNA_NObinding.CEN50-1.fasta")

One23=function(x){
  tmp=x
fastanames=c()
fastasequences=c()
for(set in 1:3){
  fname=c(paste(names(tmp)[set],"_",Fwd[6,1],"_",Rev[set,1],sep=""));
  fseq=paste(paste(FwdFile[6]),paste(tmp)[set],paste(reverseComplement(RevFile[set])),sep='');
  fastanames=c(fastanames,fname)
  fastasequences=c(fastasequences,fseq)
}
return(data.frame(fastanames,fastasequences))
}

Four210=function(x){
  tmp=x
  con=1
  fastanames=c()
  fastasequences=c()
for(set in 4:10){
  fname=c(paste(names(tmp)[set],"_",Fwd[6,1],"_",Rev[con,1],sep=""));
  fseq=paste(paste(FwdFile[6]),paste(tmp)[set],paste(reverseComplement(RevFile[con])),sep='');
  fastanames=c(fastanames,fname)
  fastasequences=c(fastasequences,fseq)
  con=con+1
}
  return(data.frame(fastanames,fastasequences))
}

Eleven247=function(x){
  tmp=x
  con=1
  fastanames=c()
  fastasequences=c()
  for(set in 11:41){
    fname=c(paste(names(tmp)[set],"_",Fwd[6,1],"_",Rev[con,1],sep=""));
    fseq=paste(paste(FwdFile[6]),paste(tmp)[set],paste(reverseComplement(RevFile[con])),sep='');
    fastanames=c(fastanames,fname)
    fastasequences=c(fastasequences,fseq)
    con=con+1
  }
  con=1
  for(set in 42:47){
    fname=c(paste(names(tmp)[set],"_",Fwd[7,1],"_",Rev[con,1],sep=""));
    fseq=paste(paste(FwdFile[7]),paste(tmp)[set],paste(reverseComplement(RevFile[con])),sep='');
    fastanames=c(fastanames,fname)
    fastasequences=c(fastasequences,fseq)
    con=con+1
  }
  return(data.frame(fastanames,fastasequences))
}
#40 41 42 CEN45
coso=One23(Nbin45)
Wellsnames[[40]]=append(Wellsnames[[40]],as.character(coso[,1]))
Wellseqs[[40]]=append(Wellseqs[[40]],as.character(coso[,2]))
coso=Four210(Nbin45)
Wellsnames[[41]]=append(Wellsnames[[41]],as.character(coso[,1]))
Wellseqs[[41]]=append(Wellseqs[[41]],as.character(coso[,2]))
coso=Eleven247(Nbin45)
Wellsnames[[42]]=append(Wellsnames[[42]],as.character(coso[,1]))
Wellseqs[[42]]=append(Wellseqs[[42]],as.character(coso[,2]))

#82 83 84 CEN50cr
coso=One23(Nbin50cr)
Wellsnames[[82]]=append(Wellsnames[[82]],as.character(coso[,1]))
Wellseqs[[82]]=append(Wellseqs[[82]],as.character(coso[,2]))
coso=Four210(Nbin50cr)
Wellsnames[[83]]=append(Wellsnames[[83]],as.character(coso[,1]))
Wellseqs[[83]]=append(Wellseqs[[83]],as.character(coso[,2]))
coso=Eleven247(Nbin50cr)
Wellsnames[[84]]=append(Wellsnames[[84]],as.character(coso[,1]))
Wellseqs[[84]]=append(Wellseqs[[84]],as.character(coso[,2]))

#124 125 126 CEN50
coso=One23(Nbin50)
Wellsnames[[124]]=append(Wellsnames[[124]],as.character(coso[,1]))
Wellseqs[[124]]=append(Wellseqs[[124]],as.character(coso[,2]))
coso=Four210(Nbin50)
Wellsnames[[125]]=append(Wellsnames[[125]],as.character(coso[,1]))
Wellseqs[[125]]=append(Wellseqs[[125]],as.character(coso[,2]))
coso=Eleven247(Nbin50)
Wellsnames[[126]]=append(Wellsnames[[126]],as.character(coso[,1]))
Wellseqs[[126]]=append(Wellseqs[[126]],as.character(coso[,2]))

#166 167 168 CEN55
coso=One23(Nbin55)
Wellsnames[[166]]=append(Wellsnames[[166]],as.character(coso[,1]))
Wellseqs[[166]]=append(Wellseqs[[166]],as.character(coso[,2]))
coso=Four210(Nbin55)
Wellsnames[[167]]=append(Wellsnames[[167]],as.character(coso[,1]))
Wellseqs[[167]]=append(Wellseqs[[167]],as.character(coso[,2]))
coso=Eleven247(Nbin55)
Wellsnames[[168]]=append(Wellsnames[[168]],as.character(coso[,1]))
Wellseqs[[168]]=append(Wellseqs[[168]],as.character(coso[,2]))

##LastFIN
##Next sets
fastasets=c("/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/FIN/SetChromFIN.CEN50-1.fasta")



for(finalfasta in fastasets){
  
  fastaFile <- readDNAStringSet(finalfasta)
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  id = type = loc = prom = c()
  for(name in seq_name){
    tmp=unlist(strsplit(name,"_"))
    id=c(id,tmp[1])
    type=c(type,tmp[2])
    loc=c(loc,unlist(strsplit(tmp[4],"-"))[1])
    prom=c(prom,tmp[length(tmp)])
  }
  Set <- data.frame(id,type,loc,prom, sequence)
  
  #Order set
  tab=data.frame(table(Set[1]))
  ktab=tab
  con=1
  groups=list()
  NgRNA=c()
  Ngen=c()
  while(nrow(ktab)>1){
    knap=knapsack(ktab[,2],rep(1,nrow(ktab)),120)
    groups=append(groups,list(ktab[knap$indices,1]))
    NgRNA=c(NgRNA,knap$capacity)
    Ngen=c(Ngen,knap$profit)
    ktab=ktab[-knap$indices,]
    con = con+1; if(con>1000){print; 'No solution'; exit}
  }
  
  
  #PerEachgroupWriteSequencesofWell
  for(i in 1:length(groups)){
    gcount=1;
    fastanames=c()
    fastasequences=c()
    for(keyid in groups[[i]]){
      tabtemp=Set[which(Set[1]==keyid),]
      for(kind in 1:length(Regions)){
        sgcount=1;
        for(last in which(tabtemp[,3]==Regions[kind])){
          if(kind==1){
            fname=c(paste(">",tabtemp[last,1],"_",tabtemp[last,3],"-",sgcount,"_",tabtemp[last,4],"_",Fwd[gcount,1],"_",Rev[Tsslim[sgcount],1],sep=""));
            fseq=paste(paste(FwdFile[gcount]),tabtemp[last,5],paste(reverseComplement(RevFile[Tsslim[sgcount]])),sep='');
            fastanames=c(fastanames,fname)
            fastasequences=c(fastasequences,fseq)
            sgcount=sgcount+1
          }
          if(kind==2){
            fname=c(paste(">",tabtemp[last,1],"_",tabtemp[last,3],"-",sgcount,"_",tabtemp[last,4],"_",Fwd[gcount,1],"_",Rev[ATGlim[sgcount],1],sep=""));
            fseq=paste(paste(FwdFile[gcount]),tabtemp[last,5],paste(reverseComplement(RevFile[ATGlim[sgcount]])),sep='');
            fastanames=c(fastanames,fname)
            fastasequences=c(fastasequences,fseq)
            sgcount=sgcount+1
          }
          if(kind==3){
            fname=c(paste(">",tabtemp[last,1],"_",tabtemp[last,3],"-",sgcount,"_",tabtemp[last,4],"_",Fwd[gcount,1],"_",Rev[CDSlim[sgcount],1],sep=""));
            fseq=paste(paste(FwdFile[gcount]),tabtemp[last,5],paste(reverseComplement(RevFile[CDSlim[sgcount]])),sep='');
            fastanames=c(fastanames,fname)
            fastasequences=c(fastasequences,fseq)
            sgcount=sgcount+1
          }
          if(kind==4){
            fname=c(paste(">",tabtemp[last,1],"_",tabtemp[last,3],"-",sgcount,"_",tabtemp[last,4],"_",Fwd[gcount,1],"_",Rev[TAGlim[sgcount],1],sep=""));
            fseq=paste(paste(FwdFile[gcount]),tabtemp[last,5],paste(reverseComplement(RevFile[TAGlim[sgcount]])),sep='');
            fastanames=c(fastanames,fname)
            fastasequences=c(fastasequences,fseq)
            sgcount=sgcount+1
          }
        }
      }
      gcount=gcount+1;  
    }
    Wellsnames=append(Wellsnames,list(fastanames))
    Wellseqs=append(Wellseqs,list(fastasequences))
  }
  
}
###Take Well 193 and put top 40 into well 185 (what a coincidence :D)
##Better make a subset fasta


fastaFile <- readDNAStringSet('/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/FIN/SetChromFIN.tmp')
seq_name = names(fastaFile)
sequence = paste(fastaFile)
id = type = loc = prom = c()
for(name in seq_name){
  tmp=unlist(strsplit(name,"_"))
  id=c(id,tmp[1])
  type=c(type,tmp[2])
  loc=c(loc,unlist(strsplit(tmp[4],"-"))[1])
  prom=c(prom,tmp[length(tmp)])
}
Set <- data.frame(id,type,loc,prom, sequence)

#Order set
tab=data.frame(table(Set[1]))
ktab=tab
con=1
groups=list()
NgRNA=c()
Ngen=c()
while(nrow(ktab)>1){
  knap=knapsack(ktab[,2],rep(1,nrow(ktab)),120)
  groups=append(groups,list(ktab[knap$indices,1]))
  NgRNA=c(NgRNA,knap$capacity)
  Ngen=c(Ngen,knap$profit)
  ktab=ktab[-knap$indices,]
  con = con+1; if(con>1000){print; 'No solution'; exit}
}


#PerEachgroupWriteSequencesofWell
for(i in 1:length(groups)){
  gcount=5;
  fastanames=c()
  fastasequences=c()
  for(keyid in groups[[i]]){
    tabtemp=Set[which(Set[1]==keyid),]
    for(kind in 1:length(Regions)){
      sgcount=1;
      for(last in which(tabtemp[,3]==Regions[kind])){
        if(kind==1){
          fname=c(paste(">",tabtemp[last,1],"_",tabtemp[last,3],"-",sgcount,"_",tabtemp[last,4],"_",Fwd[gcount,1],"_",Rev[Tsslim[sgcount],1],sep=""));
          fseq=paste(paste(FwdFile[gcount]),tabtemp[last,5],paste(reverseComplement(RevFile[Tsslim[sgcount]])),sep='');
          fastanames=c(fastanames,fname)
          fastasequences=c(fastasequences,fseq)
          sgcount=sgcount+1
        }
        if(kind==2){
          fname=c(paste(">",tabtemp[last,1],"_",tabtemp[last,3],"-",sgcount,"_",tabtemp[last,4],"_",Fwd[gcount,1],"_",Rev[ATGlim[sgcount],1],sep=""));
          fseq=paste(paste(FwdFile[gcount]),tabtemp[last,5],paste(reverseComplement(RevFile[ATGlim[sgcount]])),sep='');
          fastanames=c(fastanames,fname)
          fastasequences=c(fastasequences,fseq)
          sgcount=sgcount+1
        }
        if(kind==3){
          fname=c(paste(">",tabtemp[last,1],"_",tabtemp[last,3],"-",sgcount,"_",tabtemp[last,4],"_",Fwd[gcount,1],"_",Rev[CDSlim[sgcount],1],sep=""));
          fseq=paste(paste(FwdFile[gcount]),tabtemp[last,5],paste(reverseComplement(RevFile[CDSlim[sgcount]])),sep='');
          fastanames=c(fastanames,fname)
          fastasequences=c(fastasequences,fseq)
          sgcount=sgcount+1
        }
        if(kind==4){
          fname=c(paste(">",tabtemp[last,1],"_",tabtemp[last,3],"-",sgcount,"_",tabtemp[last,4],"_",Fwd[gcount,1],"_",Rev[TAGlim[sgcount],1],sep=""));
          fseq=paste(paste(FwdFile[gcount]),tabtemp[last,5],paste(reverseComplement(RevFile[TAGlim[sgcount]])),sep='');
          fastanames=c(fastanames,fname)
          fastasequences=c(fastasequences,fseq)
          sgcount=sgcount+1
        }
      }
    }
    gcount=gcount+1;  
  }
}

Wellsnames[[185]]=append(Wellsnames[[185]],fastanames)
Wellseqs[[185]]=append(Wellseqs[[185]],fastasequences)
Wellsnames[[193]]=c()
Wellseqs[[193]]=c()



##New experiments at the end
##Remove two last spots
Wellsnames[[192]]=c()
Wellseqs[[192]]=c()
Wellsnames[[191]]=c()
Wellseqs[[191]]=c()


#Randombinding
fastasets=c(
"/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/lastones/unc-119.CEN50-1.fasta",
"/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/lastones/CRISPorguides_FP.CEN50-1.fasta"
)


for(finalfasta in fastasets){
  
  fastaFile <- readDNAStringSet(finalfasta)
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  fastanames=c()
  fastasequences=c()
  count=1
  ##KindOfMatrix30*4
  for(i in c(1:4)){
    for(j in c(1:30)){
      if(count > length(seq_name)){next}
      fname=c(paste(">",seq_name[count],"_",count,"_",Fwd[i,1],"_",Rev[j,1],sep=""));
      fseq=paste(paste(FwdFile[i]),sequence[count],paste(reverseComplement(RevFile[j])),sep='');
      fastanames=c(fastanames,fname)
      fastasequences=c(fastasequences,fseq)
      count=count+1
    }
  }
  Wellsnames=append(Wellsnames,list(fastanames))
  Wellseqs=append(Wellseqs,list(fastasequences))
}

##FInal dpy



dpy45 <- readDNAStringSet("/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/lastones/dpy/dpy-10.CEN45.fasta")
dpy50c <- readDNAStringSet("/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/lastones/dpy/dpy-10.CEN50-1_crRNA.fasta")
dpy50 <- readDNAStringSet("/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/lastones/dpy/dpy-10.CEN50-1.fasta")
dpy55 <- readDNAStringSet("/home/velazqam/Documents/Projects/Main-sgRNA-discovery/datasets/PipeLine2GenerateDB/GenerateSets/AppendData/add/lastones/dpy/dpy-10.CEN55.fasta")

Onedpy=function(x){
  tmp=x
  fastanames=c()
  fastasequences=c()
  for(set in 1){
    fname=c(paste(names(tmp)[set],"_",Fwd[7,1],"_",Rev[set,1],sep=""));
    fseq=paste(paste(FwdFile[6]),paste(tmp)[set],paste(reverseComplement(RevFile[set])),sep='');
    fastanames=c(fastanames,fname)
    fastasequences=c(fastasequences,fseq)
  }
  return(data.frame(fastanames,fastasequences))
}

coso=Onedpy(dpy45)
Wellsnames[[169]]=append(Wellsnames[[169]],as.character(coso[,1]))
Wellseqs[[169]]=append(Wellseqs[[169]],as.character(coso[,2]))
Wellsnames[[173]]=append(Wellsnames[[173]],as.character(coso[,1]))
Wellseqs[[173]]=append(Wellseqs[[173]],as.character(coso[,2]))
Wellsnames[[177]]=append(Wellsnames[[177]],as.character(coso[,1]))
Wellseqs[[177]]=append(Wellseqs[[177]],as.character(coso[,2]))

coso=Onedpy(dpy50c)
Wellsnames[[170]]=append(Wellsnames[[170]],as.character(coso[,1]))
Wellseqs[[170]]=append(Wellseqs[[170]],as.character(coso[,2]))
Wellsnames[[174]]=append(Wellsnames[[174]],as.character(coso[,1]))
Wellseqs[[174]]=append(Wellseqs[[174]],as.character(coso[,2]))
Wellsnames[[178]]=append(Wellsnames[[178]],as.character(coso[,1]))
Wellseqs[[178]]=append(Wellseqs[[178]],as.character(coso[,2]))

coso=Onedpy(dpy50)
Wellsnames[[171]]=append(Wellsnames[[171]],as.character(coso[,1]))
Wellseqs[[171]]=append(Wellseqs[[171]],as.character(coso[,2]))
Wellsnames[[175]]=append(Wellsnames[[175]],as.character(coso[,1]))
Wellseqs[[175]]=append(Wellseqs[[175]],as.character(coso[,2]))
Wellsnames[[179]]=append(Wellsnames[[179]],as.character(coso[,1]))
Wellseqs[[179]]=append(Wellseqs[[179]],as.character(coso[,2]))

coso=Onedpy(dpy55)
Wellsnames[[172]]=append(Wellsnames[[172]],as.character(coso[,1]))
Wellseqs[[172]]=append(Wellseqs[[172]],as.character(coso[,2]))
Wellsnames[[176]]=append(Wellsnames[[176]],as.character(coso[,1]))
Wellseqs[[176]]=append(Wellseqs[[176]],as.character(coso[,2]))
Wellsnames[[180]]=append(Wellsnames[[180]],as.character(coso[,1]))
Wellseqs[[180]]=append(Wellseqs[[180]],as.character(coso[,2]))




#write first files with headers to replace appends
##Write 'for' for final output
write(paste("pool_num","oligo_num","oligo_seq","oligo_len",	"oligo_ext_name",sep=","), file = 'OligoDBv1List.csv')
for(well in 1:length(Wellseqs)){
  for(seqid in 1:length(unlist(Wellsnames[[well]]))){
    write(paste(well,seqid,unlist(Wellseqs[[well]])[seqid],nchar(unlist(Wellseqs[[well]])[seqid]),unlist(Wellsnames[[well]])[seqid],sep=","),
          file = 'OligoDBv1List.csv',append = TRUE,sep='')
  }
}




