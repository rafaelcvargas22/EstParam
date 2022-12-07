#if(!(require(edgeR))) install.packages('edgeR')
#library('edgeR')

input <- function(inputfile) {
	  X <<- data.matrix(read.table(file=inputfile, sep="\t", header=F))
}

run <- function() {
	type=1                            # DEFAULTED TO 1 AS 2 IS NOT IMPLEMENTED
	if (!class(X)%in%c("ExpressionSet","matrix"))
		stop(class(X))
    
    	if(!is.matrix(X))  X=exprs(X)
    		if (!(type==1|type==2)) stop("type must be either 1 or 2")

		if(type==1){# use simple estimation with minimum assumption
			# Code from function getDisp1, seed defaulted to 2015
	          	seed=2015
	    		seqDepth=colSums(X)
    			k=seqDepth/median(seqDepth)
    			X2=sweep(X, 2, k, FUN="/")##pseudo Counts
	            	
  			m=rowMeans(X2)
			v=rowVars(X2)
			phi.g0 = phi.g = (v-m)/m^2
			## only keep those with good coverage
			Good=m>30 & rowMeans(X>0)>0.8
			phi.g0=phi.g0[Good]
			phi.g0=replace(phi.g0,is.na(phi.g0),.001)
			phi.g0=replace(phi.g0,phi.g0<.001,0.001)

			## for those with unobserved dispersion, sample from phi.g0
			ii=(phi.g<=0.001) | is.na(phi.g)
			set.seed(seed)
			phi.g[ii]=sample(phi.g0, sum(ii), replace=TRUE)
							            
			res <<- list(seqDepth=seqDepth,lmeans=log(m),lOD=log(phi.g))
		}
      
      #if(type==2){ #use edgeR          # TO BE IMPLEMENTED
        #require(edgeR)
        #res <<- getDisp2(X)
      #}
}

rowVars=function(x) {
	n0 <- ncol(x)
      	EX <- rowMeans(x, na.rm=TRUE)
	vv <- rowSums(sweep(x,1,EX)^2)/(n0-1)
	vv
}

output <- function(outputfile){
	fileConn <- file(outputfile);
  	writeLines(paste(res), fileConn)
  	close(fileConn);
}
