###################################
##  simulate the sequencing procedure and genotype likelihood estimation
##  x is the true genotype
##  d is the average depth for all the individuals or it can be a vector of average depths, one for each individuals
##  e is an error rate in the sequencing procedure.
##  norm is an optional normalization of the likelihood. This will no affect the results but can help to avoid underflow problems.
##  liusy99@mail.sysu.edu.cn March 26th, 2024
########################################################

options(scipen=999)


getGenoProb<-function(x){
  c(
    x[1]^2,     ##AA
    2*x[1]*x[2],##AC
    2*x[1]*x[3],##AG
    2*x[1]*x[4],##AT
    x[2]^2     ,##CC
    2*x[2]*x[3],##CG
    2*x[2]*x[4],##CT
    x[3]^2,     ##GG
    2*x[3]*x[4],##GT
    x[4]^2      ##TT
    )
}


getpileup<-function(tgeno,d=0.06,errMat,norm=FALSE){
  bases=c("A","C","G","T")
  n<-length(tgeno) # the true genotype interger vector of all individuals
  dep<-rpois(n,d) # the individual depth
  bq=5 #ASCII 5 corresponding to 20 in Sanger 33 Phred-scale system
  cells=rep("*",n*3)
  for(i in 1:n){
    cells[3*i-2]=dep[i]
    if(dep[i]!=0){
        multinom=rmultinom(1,dep[i],errMat[tgeno[i],])
        b=rep("NA",4)
        for(k in 1:4){ # 1:length(bases) or 1:length(multinom)
            if(multinom[k,1]!=0){
                b[k]=paste(rep(bases[k],multinom[k,1]),collapse="")
            }
        }
        cells[3*i-1]=paste(b[b!="NA"],collapse="")
        cells[3*i-0]=paste(rep(bq,nchar(cells[3*i-1])),collapse="")
    }
  }
  return(cells)
  #SY assuming the genotype is AA: rmultinom(1,4,c(0.99,0.01/3,0.01/3,0.01/3))
}
      

getBaseProb<-function(e){
  HO<-(1-e)
  HE<-0.5*(1-2/3*e)
  EE<-e/3
  errMat<- rbind(
       #A   C   G   T #SY
        #assuming the genoytpe is AA and error rate is 0.01, c(0.09,0.01/3,0.01/3, 0.01/3)
       c(HO,EE,EE,EE), ##AA 
       c(HE,HE,EE,EE), ##AC
       c(HE,EE,HE,EE), ##AG
       c(HE,EE,EE,HE), ##AT
       c(EE,HO,EE,EE), ##CC
       c(EE,HE,HE,EE), ##CG
       c(EE,HE,EE,HE), ##CT
       c(EE,EE,HO,EE), ##GG
       c(EE,EE,HE,HE), ##GT
       c(EE,EE,EE,HO) ##TT
       )
  errMat
}


run<-function(af,d,emat,n){
  #the probability of the 10 genotype assuming HWE
  # AA AC AG AT CC CG CT GG GT TT
  Gf<-getGenoProb(af[1:4])
  ##simulate the true genotypes
  TrueGeno<-sample(1:10,n,prob=Gf,replace=TRUE)  
  ##generate pileup based on the True Genotype and the error rate 
  line=getpileup(TrueGeno,d=depth,errMat=emat)
   
   return(line)
  
}


af_table<-function(interval=0.00001){
    af1=rep(1,100)
    af_mono=cbind(af1,0,0,0,1)
    af=seq(0.00001,0.5,interval) #the minor allele frequency will always be given to the non-A alleles
    af_bi=cbind(1-af,af,0,0,2)
    af_tri=cbind(1-af,af/2,af/2,0,3)
    af_tetra=cbind(1-af,af/3,af/3,af/3,4)
    all_af=rbind(af_mono,af_bi,af_tri,af_tetra) ## the true frequency of A C G T
    all_loci_af=cbind(1:nrow(all_af),all_af)
    all_af_sim=all_loci_af[loci_start:loci_end,]
#    write.table(all_af_sim,aftable,col.names=F,row.names=F,quote=F)
    return(all_af_sim)
}

##################### arguments
args=commandArgs(T)
#number of individuals e.g. 1000
ind<-as.numeric(args[1])
#depth e.g. 1X
depth <- as.numeric(args[2])
#error rate
error <- as.numeric(args[3])
loci_start<-as.numeric(args[4])
loci_end<-as.numeric(args[5])
outfile<-args[6]
#aftable<-args[7]
outputdir=dirname(outfile)

##################### af_tables
all_af_sim=af_table(interval=0.00001)

#################### error matrix
#the probabilities of the four bases given a genotype
errMat<-getBaseProb(error) 

####################
loci=length(all_af_sim[,1])
result=matrix("*",nrow=loci,ncol=3*(ind+1),byrow=T)
for (l in 1:loci){
    start=loci_start+l-1
    header=c("chrN",start,"A")
    body=as.vector(run(all_af_sim[l,2:5],d=depth,emat=errMat,n=ind))
    result[l,]=c(header,body)
}

#header=matrix(rep(c("chrN",1,"A"),loci),c(loci,3),byrow=T)
#result=cbind(header,t(apply(all_af_sim[,1:4],1,run,d=depth,emat=errMat,n=ind)))
write.table(result,outfile,col.names=F,row.names=F,quote=F,sep="\t")

