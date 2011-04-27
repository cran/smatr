
sma <- function(formula, data, subset, na.action, log='',
	 method=c("SMA","MA"), type=c("elevation","shift"), alpha=0.05, 
	 slope.test=NA, elev.test=NA,
	 ...)
{
	method <- match.arg(method)
	type <- match.arg(type)

	# Model frame (code borrowed from 'lm')
	mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())

	# Log-transform 
	log <- tolower(log)
	if(log == "x")mf[,2] <- log10(mf[,2])
	if(log == "y")mf[,1] <- log10(mf[,1])
	if(log == "xy"){
		mf[,2] <- log10(mf[,2])
		mf[,1] <- log10(mf[,1])
	}
	if(!(log %in% c("","x","y","xy")))
		warning("Log transformation ignored! Use one of '', 'x','y' or 'xy'")
	
	#REMOVE INCOMPLETE DATA (DAN)
	mf <- mf[(!is.na(mf[,2])&!is.na(mf[,1])),]
	
	#n <- nrow(mf)  #not needed (DAN)

	# Prepare testing by group.
	mt <- attr(mf, "terms")
	vn <- names(mf)  # variable names
	tn <- attr(mt, "term.labels")   # term names
	
	# Check if there is group testing, and check that formula is compatible.
	if(length(tn) == 1)grouptest <- "none"
	if(length(tn) == 2){
		formulacheck <- all(tn == c(vn[2], vn[3]))
		if(formulacheck){
			if(type == "elevation")grouptest <- "elevcom"
			if(type == "shift")grouptest <- "shiftcom"
		} else grouptest <- "malformed"
	}
	if(length(tn) == 3){
		formulacheck <- all(tn == c(vn[2], vn[3], paste(vn[2],":",vn[3],sep="")))
		if(formulacheck)grouptest <- "slopecom" else grouptest <- "malformed"
	}
	if(length(tn) > 3 || grouptest == "malformed"){
		warning("Formula not supported by sma() and/or ma(), and is ignored.")
		grouptest <- "none"
	}
	
	# Check for intercept. Also note that it is not allowed to drop the intercept for some group tests,
	# but this is not yet tested here!
	if(attr(mt, "intercept") == 0)intercept <- FALSE
	if(attr(mt, "intercept") == 1)intercept <- TRUE
	
	#CHANGES HERE  - REMOVED CLAULCATIONS FORM IF LOOPS, SO AS TO NOT DUPLICATE CODE (DAN)
	
	#Determine grouping 
	
	if(grouptest %in% c("elevcom","shiftcom","slopecom")){
		ngroups <- nlevels(mf[,3])
		grps<-mf[,3]
		lv <- levels(grps)
		commoncoef <- line.cis(mf[,1], mf[,2], intercept=intercept, method=method, alpha=alpha, ...)
		
		#run group tests
		if(grouptest == "elevcom"){
			if(!intercept)stop("Cannot perform elevation test without fitted intercept.")
			grouptestresult <- elev.com(mf[,1], mf[,2], mf[,3], alpha=alpha, method=method, group.names=lv)
		}
		if(grouptest == "slopecom"){
			grouptestresult <- slope.com(mf[,1], mf[,2], mf[,3], alpha=alpha, intercept=intercept, method=method, group.names=lv)
		}
		if(grouptest == "shiftcom"){
			grouptestresult <- shift.com(mf[,1], mf[,2], mf[,3], intercept=intercept, method=method, group.names=lv)
		}	
	 }
	else{ # single group
		ngroups<-1	
		grps<-as.factor(rep("all", length(mf[,1])))
		lv <- levels(grps)
		commoncoef <- NA
		grouptestresult <- ""
	}
	 
	#Calculate stuff for each groupGet the sma coefficients 
	coeff <- list(); n<- list(); r2<- list(); pval <- list(); from <- list(); to<-list(); slopetest <-list(); elevtest <-list();
		
	for(i in 1:ngroups){
		X <- mf[grps == lv[i],2]
		Y <- mf[grps == lv[i],1]
		
		#groupsize
		n[[i]] <- length(X); 
		
		#sma coefficients 
		coeff[[i]] <- line.cis(Y,X,intercept=intercept, method=method, alpha=alpha, ...)   
		B=coeff[[i]][2,1]; 
				
		if(!intercept){
			coeff[[i]] <- coeff[[i]][2,]
			a= 0;
		} else {
			a= coeff[[i]][1,1]		
		}
		
		#correlation	---> what about when no intercept, need to change?
		r2[[i]]<- cor(X, Y)^2;    
        pval[[i]] <- cor.test(X, Y, method = "pearson")$p.value;
      
      	# Test slope against some specified value
     	if(!is.na(slope.test)){
			slopetest[[i]] <- slope.test(Y,X,  test.value=slope.test, method=method, alpha=alpha, intercept=intercept)
		} else {
			slopetest[[i]] <-slope.test(Y,X,  test.value=NA, method=method, alpha=alpha, intercept=intercept)
		}
	
		# Test elevation against some specified value
		if(!is.na(elev.test)){
			if(!intercept)stop("Cannot perform elevation test without fitted intercept.")
				elevtest[[i]] <- elev.test( Y,X, test.value=elev.test, method=method, alpha=alpha)
		} else {
				elevtest[[i]] <- elev.test( Y,X, test.value=NA, method=method, alpha=alpha)
		}
      	
      	#determine range of fitted values (as X value)
        if(method=="SMA"){
  	    	from[[i]] <- (min(Y+B*X) - a)/(2.0*B) 
  	    	to[[i]] <-(max(Y+B*X)-a)/(2.0*B)
  	    } else if (method =="MA"){
  	    	from[[i]] <- (min(X+B*Y) - B*a)/(1+B)
  	    	to[[i]] <-(max(X+B*Y) - B*a)/(1+B)
  	    }
  	   	if(log %in% c("x","xy"))
		{
			from[[i]] = 10^from[[i]]
			to[[i]] = 10^to[[i]]
		}    
	}
	#apply names to new variables
	names(coeff) <- lv; names(n) <- lv; names(r2) <- lv; names(pval) <- lv; names(from) <- lv; names(to) <- lv;
	
	l <- list()
	l$coef <- coeff
	l$commoncoef <- commoncoef
	l$alpha <- alpha
	l$method <- method
	l$intercept <- intercept
	l$call <- match.call()
	l$data <- mf
	l$log <- log
	l$variables <- names(mf)
	l$origvariables <- all.vars(match.call()$formula)
	l$groups <-lv;
	l$gt <- grouptest
	l$gtr <- grouptestresult
	
	l$slopetest <- slopetest
	l$elevtest <- elevtest
	l$n <- n;
	l$r2 <- r2;
	l$pval <- pval;
	l$from <- from;
	l$to <-to;
	
	l$summary <-NULL;
#	for(i in 1:ngroups){
#		l$summary <-rbind(l$summary, data.frame(group=l$groups[i], n= l$n[i] , r2 <-l$r2[i], pval <-l$pval[i], Slope = l$coef[[i]][2,1], Slope_lowCI = l$coef[[i]][2,2], Slope_highCI = l$coef[[i]][2,3],  Int = l$coef[[i]][1,1], Int_lowCI = l$coef[[i]][1,2], Int_highCI = l$coef[[i]][1,3], Slope_test = l$slopetest[[i]]$test.value, Slope_test_p= l$slopetest[[i]]$p, Elev_test = l$elevtest[[i]]$test.value, Elev_test_p= l$elevtest[[i]]$p))
#	}
		
	class(l) <- "sma"
	
    return(l)
}




