## ---------------------------
## Script name: isoTracker
## Purpose of script: Finds groups of features that could correspond to the same isotopic envelope and that are statistically siginicant dysregulated between labeled and non-labeled groups.
##
## Authors: Xavier Domingo-Almenara, Ph.D (author); Markus M. Rinschen, M.D. (contributer)
##
## Date Created: 2019-10-09
##
## Copyright (c) Xavier Domingo-Almenara, 2020
## Email: xavier.domingoa@eurecat.org
##
## ---------------------------
##
## Notes: this is script was designed to work with specifically formatted data from the following paper:
##   Rinschen M.M., Palygin O., Golosova D., Domingo-Almenara X., et al. Lysine’s metabolic activity conveys kidney protection in hypertension. (2020) Submitted
##
## ---------------------------

library(igraph)
library(openxlsx)
library(magrittr)
	

getPPM <- function(mass.error, mass){
	(mass.error * 10^6)/round(mass)
}


isotopeDec <- function(massVal, MZ.vector, rtVal, RT.vector, maxRTdiff, isotope.dis, ppm.error, nPositions)
{	
	mzVect <- massVal + cumsum(rep(isotope.dis, nPositions))
	rtWinInd <- which(abs(rtVal - RT.vector) < maxRTdiff)
	mzVectorInd <- sapply(mzVect, function(j) which.min(abs(MZ.vector[rtWinInd]-j)))
	ppmErrors <- mapply(function(y,z){getPPM(abs(z-y),z)}, MZ.vector[rtWinInd][mzVectorInd], mzVect)
	outV <- rep(0, length(mzVect))
	outV[which(ppmErrors<ppm.error)] <- rtWinInd[mzVectorInd][which(ppmErrors<ppm.error)]
	
	outV	
}


getLabelSignificance <- function(parentInt, isoInt, labelInds)
{	
	isoInt[is.na(isoInt)] <- 0
	parentInt[is.na(parentInt)] <- 0
	
	difUnlab <- parentInt[-labelInds] - isoInt[-labelInds] 
	difLab <- parentInt[labelInds] - isoInt[labelInds] 

	pVal <- try(2*t.test(difUnlab, difLab, alternative='greater')$p.value, silent=T)
	if(class(pVal)=='try-error') pVal <- 1
	if(pVal>1) pVal <- 1

	pVal	
}

lf <- list.files(pattern='xlsx')

for(fln in lf){
		
	feature.table <- openxlsx::read.xlsx(fln)
	isoNumber <- 23
	
	experimentClass <- (strsplit(fln, "_")[[1]] %>% strsplit(split='\\.'))[[2]][1]
	
	gGroups <- data.frame(matrix(c(rep('con',3),rep('1w',3),rep('2w',3),rep('3w',3), 27:ncol(feature.table)), ncol=2), stringsAsFactors=FALSE)
	colnames(gGroups) <- c('class','index')
	gGroups$index <- as.numeric(as.vector(gGroups$index))
	
		mzVector <- as.numeric(as.vector(feature.table[,'mzmed']))
		mzOrd <- order(mzVector)
		mzVectorInd <- rep(1, length(mzVector))
		mzVector <- mzVector[mzOrd]
		featureIdx <- feature.table[mzOrd,'featureidx']
		featureName <- feature.table[mzOrd,'name']
		rtVector <- 60*as.numeric(as.vector(feature.table[mzOrd,'rtmed']))
		isoList <- list()
		for(i in 1:length(mzVector)){
			if(is.na(mzVectorInd[i])) next
			y <- isotopeDec(mzVector[i], mzVector, rtVector[i], rtVector, maxRTdiff=2, ppm.error=10, isotope.dis=1.003355, nPositions=isoNumber-1)
			if(all(y==0)) next
			mzVectorInd[y] <- NA
			isoList <- c(isoList, 1)
			isoList[[length(isoList)]] <- c(i, y)
		}
		
		isoTab <- matrix('', nrow=nrow(feature.table), ncol=(4+(4*isoNumber)))
		cTagP <- c(rep('con',isoNumber),rep('1w',isoNumber),rep('2w',isoNumber),rep('3w',isoNumber))
		iTagP <- c(0:(isoNumber-1), 0:(isoNumber-1), 0:(isoNumber-1), 0:(isoNumber-1))
		cNames <- sapply(1:length(cTagP), function(x) paste(cTagP[x], '+', iTagP[x], sep=''))
		
		colnames(isoTab) <- c('id','isotope', 'nIsotopes', 'pval', cNames)
		
		imgDirName <- paste0('IMG_', experimentClass, '/')
		unlink(imgDirName, recursive=T)
		dir.create(imgDirName)		
		nonLabClass <- unique(gGroups$class)[1]	
		kId <- 1
		for(i in 1:length(isoList)){
			isoPosInds <- which(isoList[[i]]!=0)
			xInd <- isoList[[i]][isoPosInds]
			yCon <- feature.table[mzOrd[xInd],]
						
			yCon[yCon==0] <- 1
			dCon <- sapply(unique(gGroups$class), function(k) {
				oCon <- yCon[,gGroups$index[which(gGroups$class==k)]]
				oConD <- as.vector(unlist(apply(oCon, 2, function(x) dist(x/x[1]))))
				oConD
			})
	
			pVals <- apply(combn(1:ncol(dCon),2), 2, function(x) t.test(dCon[,x[1]], dCon[,x[2]])$p.value)
			pVals[is.na(pVals)] <- 1
			if(!any(pVals<0.05)) next
			gCon <- sapply(unique(gGroups$class), function(k) rowSums(yCon[,gGroups$index[which(gGroups$class==k)]]))
	
			if(length(xInd)==2 & gCon[2,1]>gCon[1,1]) next
	
			gMaxs <- apply(gCon,2,max)
			isoAbund <- round(100*sweep(gCon, 2, gMaxs, '/'),0)
			minPVal <- round(min(pVals, na.rm=T),5)
			isoPos <- sapply(isoPosInds-1, function(x) paste('M+', x, sep=''))
			
			isoTab[xInd, 1:4] <- cbind(rep(kId,length(xInd)), isoPos, rep(length(xInd),length(xInd)), rep(minPVal,length(xInd)))
			isoMatPos <- c(4+isoPosInds, 4+unlist(sapply(cumsum(rep(isoNumber,3)), function(x) x + isoPosInds, simplify=FALSE)))
			isoTab[xInd[1],isoMatPos] <- as.vector(isoAbund)	
								
			pName <- paste(featureName[xInd], collapse=', ')			
			pdf(file=paste(imgDirName, featureIdx[xInd] , '.pdf', sep=''))
			plot(mzVector[xInd], gCon[,1]/max(gCon[,1]), ylim=c(-1, 1), type='h', lwd=2, xlim=c(min(mzVector[xInd])-5, max(mzVector[xInd]) + 5), main=pName, xlab='m/z', ylab='Rel. Intensity (%)')	
			lines(mzVector[xInd], -1*gCon[,4]/max(gCon[,4]), type='h', col='#FFD5D5', lwd=2)
			lines(mzVector[xInd], -1*gCon[,3]/max(gCon[,3]), type='h', col='#FF8282', lwd=2)
			lines(mzVector[xInd], -1*gCon[,2]/max(gCon[,2]), type='h', col='#FF0000', lwd=2)
			dev.off()
			kId <- kId + 1
		}
		
		fileOutput <- paste0(strsplit(fln, split='\\.')[[1]][1], '_IsoTrkr.xlsx')
		featTabOutput <- cbind(feature.table[mzOrd,], isoTab)
		write.xlsx(featTabOutput, file=fileOutput)

}
	