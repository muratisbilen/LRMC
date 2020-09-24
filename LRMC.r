library(survival)
setwd("Input Files")
inputFiles=try(shell("dir /B ",intern=T,wait=T))
setwd("../")

for(files in inputFiles){
	output = paste(substring(files, 1, nchar(files)-4), " LRMC Graphs.pdf", sep="")
	output2 = paste(substring(files, 1, nchar(files)-4), " LRMC Statistics.txt", sep="")
	
	setwd("Input Files")
	d = as.matrix(read.delim(file=files,sep="\t",header=T))
	setwd("../")
	
	surv = as.numeric(d[nrow(d)-1,-c(1,ncol(d))])
	stat = as.numeric(d[nrow(d),-c(1,ncol(d))])
	gs = as.vector(as.matrix(d[1:(nrow(d)-2),ncol(d)]))
	ps = as.vector(as.matrix(d[1:(nrow(d)-2),1]))
	
	mat = c()

	for(i in 1:(nrow(d)-2)){
		mat = cbind(mat,as.numeric(as.matrix(d[i,-c(1,ncol(d))])))
	}

	colnames(mat) = ps
	mat2 = data.frame(mat)
	a = apply(mat2,1,sum)
	mat3 = mat2[which(!is.na(a)),]

	setwd("Output Files")
	pdf(output)
	setwd("../")
	thrs = c()
	genes = c()
	probesets = c()
	pvals = c()
	hrs = c()
	
	for(j in 1:ncol(mat3)){
		k = mat3[,j]
		sk = sort(unique(k[which(!is.na(k))]))
		sk2 = sk[2:(length(sk)-2)]
		p = c()
		p2 = c()
		color=c()
		
		for(i in 2:(length(sk)-2)){
			g = k
			thrs = c(thrs,sk[i])
			genes = c(genes,gs[j])
			probesets = c(probesets,ps[j])
			pos1 = which(k<=sk[i])
			pos2 = which(k>sk[i])
			g[pos1] = 0
			g[pos2] = 1

			lr = survdiff(Surv(surv,stat)~g)
			cp = coxph(Surv(surv,stat)~g)
			KM = survfit(Surv(surv,stat)~g)
			plr	= pchisq(lr$chisq, 1, lower.tail=FALSE)
			hr = exp(cp$coef)
			hrs = c(hrs,hr)
			p = c(p,-log10(plr))
			p2 = c(p2,plr)
			pvals = c(pvals,plr)
			
			if(!is.na(plr) && !is.na(hr)){
				if(plr<1 && hr>1){
					color = c(color,"red")
				}else if(plr<1 && hr<1){
					color = c(color,"blue")
				}else{
					color = c(color,"grey")
				}
			}else{
				color = c(color,NA)
			}
		}
		
		possel = which(!is.na(p))
		sk3 = sk2[possel]
		pp = p[possel]
		color2 = color[possel]
		
		quan = quantile(k)
		plot(c(quan[2],quan[2]),c(-0.5,max(pp,-log10(0.05))+0.5),ylim=c(0,max(pp,-log10(0.05))),xlim=c(min(sk3),max(sk3)),type="l",lty=2,col="darkgrey",lwd=2,main=paste(gs[j],ps[j],"LRMC Graph",sep="\n"), ylab="-log10(p)",xlab="Expression Threshold",cex.lab=1.5,cex.axis=1.5,cex=1.3)
		lines(c(quan[3],quan[3]),c(-0.5,max(pp,-log10(0.05))+0.5),lty=2,col="darkgrey",lwd=2)
		lines(c(quan[4],quan[4]),c(-0.5,max(pp,-log10(0.05))+0.5),lty=2,col="darkgrey",lwd=2)
		lines(sk3,pp,type="b",pch=16,ylim=c(0,max(p,-log10(0.05))))
		points(sk3,pp,col=color2,pch=16,cex=1.3)
		abline(h=-log10(0.05),lty=2,lwd=2)
	}
	res = c("Gene Symbol","Probeset","HR","Threshold","Log-Rank p value")
	res = rbind(res,cbind(genes,probesets,hrs,thrs,pvals))
	setwd("Output Files")
	write(file=output2,t(res),ncol=ncol(res),sep="\t")
	dev.off()
	setwd("../")
}
