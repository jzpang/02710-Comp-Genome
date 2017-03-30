####Perform a significant test

myperm<-function(data,data.grp){
	mat=matrix(1:10,2,5)
	layout(mat)
	layout.show(10)
	set.seed(1)

	for (i in 1:10){
		datar=t(apply(data,1,sample))
		T=apply(datar,1,myttest,grp=data.grp)
		p=P(T)
		hist(p,col="lightblue",labels=T,xlab="p-value")
	}
}


####perform a t-test
myttest <- function(x, grp,s0){
	c=sort(unique(grp))
	c1=grep(c[1], grp)
	c2=grep(c[2], grp)
	x1=x[c1]
	x2=x[c2]
	s1=var(x1)
	s2=var(x2)
	m1=mean(x1)
	m2=mean(x2)
	l1=length(x1)
	l2=length(x2)
	s=sqrt(((l1-1)*s1+(l2-1)*s2)/(l1+l2-2))
	t=(m1-m2)/(s*sqrt(1/l1+1/l2)+s0)
	return (t)
}

###Calculate adjusted p-value
P<-function(T){
	P=c()
for (i in T){
	if (i>0){
		j=(1-pt(i,df=45))*2
	}
	else{
		j=2*pt(i,df=45)
	}
	P=c(P,c(j))
}
return (P)
}


####Calculate emperical FDR
epfdr<-function(data, data.grp){}
	M=c()
	set.seed(1)
	for (i in 1:10){
		mygrp=data.grp[sample(1:length(data.grp))]
		T=apply(data,1,myttest,grp=mygrp)
		M=cbind(M,T)
	}

	T=apply(data,1,myttest,grp=data.grp)

	epfdr=c()
	for (i in T){
		n0=length(T[abs(T)>=abs(i)])
		np=(length(M[abs(M)>=abs(i)]))/10
		c=np/n0
		epfdr=c(epfdr,c(c))
	}
	return (epfdr)
}