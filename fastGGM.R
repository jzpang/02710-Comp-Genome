###data process

#########t-stat
myttest<-function(x,y,s0){
	t=NULL
for (indexiii in 1:length(x[,1])){

    if(indexiii%%1000==0){cat (indexiii,"\n")}
#splitgrp=split(as.matrix(x[indexiii,]),grp,drop=false)
x1=as.numeric(x[indexiii,])
x2=as.numeric(y[indexiii,])
n1=length(x1)
n2=length(x2)

sx1x2=sqrt(((n1-1)*var(x1)+(n2-1)*var(x2))/(n1+n2-2))

t1=(mean(x1)-mean(x2))/(sx1x2*sqrt((1/n1)+(1/n2))+s0)

t=c(t,t1)

}

return(t)

}

####
pValue<-function(t.stat){
    pval=NULL
    for(i in 1:length(t.stat)){
    if(t.stat[i]>0){
       pval=c(pval,2*pt(t.stat[i],386,lower.tail=FALSE))
    }else{
      pval=c(pval,2*pt(t.stat[i],386))
    }
    }
    return (pval)
}

###perumutation
##random the group for 10 times and get the t-stat for each set
groupStat<-function(x,x.grp,s0){
    groupP=NULL
    grp=x.grp
    set.seed(1)
    for(i in 1:10){
        grp=x.grp[sample(1:length(x.grp))]
        t.stat=myttest(x,grp,s0)
        groupP=cbind(groupP,t.stat)
    }
    return (groupP)
}

findingFdr<-function(data,data.grp){
	original.Tstat=myttest(data,data.grp,0)
	groupT.stat=groupStat(data,data.grp,0)
   fdr=NULL
	for(i in 1:length(original.Tstat[])){
		numExceed=length(groupT.stat[abs(groupT.stat)>=abs(original.Tstat[i])])
		fdr=c(fdr,numExceed/10/length(original.Tstat[abs(original.Tstat)>=abs(original.Tstat[i])]))

	}
	return (fdr)
}


data.all=read.table('GSE63061_normalized.txt',sep='\t',stringsAsFactors=F)
data=data.all[-1,1];

info=read.table('patient_info-1.txt',stringsAsFactors=F)

data.name=NULL
for(i in data.all[1,-1]){
name=gsub('X','',i)
line=grep(name,info[,1])
group=0
if(info[line,2]=="yes"){
	group=2
}else{
	group=1
}
data.name=c(data.name,group)
}

t.statistic=myttest(data,data.name,0)
pvalue=pValue(t.statistic)

##BH ajustment 
adjustPvalue=p.adjust(pvalue,'BH')

index=which(adjustPvalue<0.05)

data.GGM=t(data[index,])
data.GGM=as.matrix(data.GGM)
data.GGM=apply(data.GGM,1,function(x) as.numeric(x))


###Bonferroni ajustment
adjust_Bon=p.adjust(pvalue,'bonferroni')

index_Bon=which(adjustPvalue<0.05)

Bon.GGM=t(data[index_Bon,])
Bon.GGM=as.matrix(data.GGM)
Bon.GGM=apply(data.GGM,2,function(x) as.numeric(x))


###permutation test fdr
fdr_value=findingFdr(data,data.name)
index_permu=which(fdr_value<0.05)

###running GGM
library(FastGGM)

outlist1 <- FastGGM(data.GGM)

partialCor=outlist1$partialCor

