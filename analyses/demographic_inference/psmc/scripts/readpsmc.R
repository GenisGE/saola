
# read pscm output psmc.result() adapted from function in https://datadryad.org/stash/dataset/doi:10.5061/dryad.0618v

##-------Rescale the ith iteration result of PSMC, and make ready for plotting
# file: result file from PSMC
# i.iteration: the ith iteration
# mu: mutation rate
# s: bin size
# g: years per generation

psmc.result<-function(file,i.iteration=25,mu=1e-8,s=100,g=1)
{
	X<-scan(file=file,what="",sep="\n",quiet=TRUE)
	
	START<-grep("^RD",X)
	END<-grep("^//",X)
	
	X<-X[START[i.iteration+1]:END[i.iteration+1]]
	
	TR<-grep("^TR",X,value=TRUE)
	RS<-grep("^RS",X,value=TRUE)
	
	write(TR,"temp.psmc.result")
	theta0<-as.numeric(read.table("temp.psmc.result")[1,2])
	N0<-theta0/4/mu/s # this is the same as theta0/(4*mu)/s, don't panic next time
	
	write(RS,"temp.psmc.result")
	a<-read.table("temp.psmc.result")
	Generation<-as.numeric(2*N0*a[,3]) # a[,3] is t_k
	Ne<-as.numeric(N0*a[,4]) #a[,4] is lambda_k
	
	file.remove("temp.psmc.result")
	
	n.points<-length(Ne)
	YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
		Generation[n.points])*g
	Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
		Ne[n.points])
	
	data.frame(YearsAgo,Ne)
}
