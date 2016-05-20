#
# Code obtained from [2] under the GPL2 license.  Original code was modified for [2] from [3] by authors of [2].
#
# [2] Smith, C. A., Want, E. J., O'Maille, G., Abagyan, R., & Siuzdak, G. (2006). XCMS: processing mass spectrometry data for metabolite profiling using nonlinear peak alignment, matching, and identification. Analytical chemistry, 78(3), 779-787. \link[dest=10.1021/ac051437y]{10.1021/ac051437y}
# [3] Du, P., Kibbe, W. A., & Lin, S. M. (2006). Improved peak detection in mass spectrum by incorporating continuous wavelet transform-based pattern matching. Bioinformatics, 22(17), 2059-2065. \link[dest=10.1093/bioinformatics/btl355]{10.1093/bioinformatics/btl355}
#

getLocalMaximumCWT = function(wCoefs, minWinSize=5, amp.Th=0)
  {        ## from package MassSpecWavelet
    localMax <- NULL
    scales <- as.numeric(colnames(wCoefs))
    
    for (i in 1:length(scales)) {
      scale.i <- scales[i]
      winSize.i <- scale.i * 2 + 1
      if (winSize.i < minWinSize) {
        winSize.i <- minWinSize
      }
      temp <- localMaximum(wCoefs[,i], winSize.i)
      localMax <- cbind(localMax, temp)
    }
    ## Set the values less than peak threshold as 0
    localMax[wCoefs < amp.Th] <- 0
    colnames(localMax) <- colnames(wCoefs)
    rownames(localMax) <- rownames(wCoefs)
    localMax
  }


localMaximum = function (x, winSize = 5)
  {   ## from package MassSpecWavelet
    len <- length(x)
    rNum <- ceiling(len/winSize)
    
    ## Transform the vector as a matrix with column length equals winSize
    ##		and find the maximum position at each row.
    y <- matrix(c(x, rep(x[len], rNum * winSize - len)), nrow=winSize)
    y.maxInd <- apply(y, 2, which.max)
    ## Only keep the maximum value larger than the boundary values
    selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > x[winSize]))
    
    ## keep the result
    localMax <- rep(0, len)
    localMax[(selInd-1) * winSize + y.maxInd[selInd]] <- 1
    
    ## Shift the vector with winSize/2 and do the same operation
    shift <- floor(winSize/2)
    rNum <- ceiling((len + shift)/winSize)
    y <- matrix(c(rep(x[1], shift), x, rep(x[len], rNum * winSize - len - shift)), nrow=winSize)
    y.maxInd <- apply(y, 2, which.max)
    ## Only keep the maximum value larger than the boundary values
    selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > x[winSize]))
    localMax[(selInd-1) * winSize + y.maxInd[selInd] - shift] <- 1
    
    ## Check whether there is some local maxima have in between distance less than winSize
    maxInd <- which(localMax > 0)
    selInd <- which(diff(maxInd) < winSize)
    if (length(selInd) > 0) {
      selMaxInd1 <- maxInd[selInd]
      selMaxInd2 <- maxInd[selInd + 1]
      temp <- x[selMaxInd1] - x[selMaxInd2]
      localMax[selMaxInd1[temp <= 0]] <- 0
      localMax[selMaxInd2[temp > 0]] <- 0
    }
    
    localMax
  }

getRidge = function(localMax, iInit=ncol(localMax), step=-1, iFinal=1, minWinSize=3, gapTh=3, skip=NULL)
  {  ## modified from package MassSpecWavelet
    
    scales <- as.numeric(colnames(localMax))
    if (is.null(scales))  scales <- 1:ncol(localMax)
    
    maxInd_curr <- which(localMax[, iInit] > 0)
    nMz <- nrow(localMax)
    
    if (is.null(skip))	{
      skip <- iInit + 1
    }
    
    ## Identify all the peak pathes from the coarse level to detail levels (high column to low column)
    ## Only consider the shortest path
    if ( ncol(localMax) > 1 ) colInd <- seq(iInit+step, iFinal, step)
    else colInd <- 1
    ridgeList <- as.list(maxInd_curr)
    names(ridgeList) <- maxInd_curr
    peakStatus <- as.list(rep(0, length(maxInd_curr)))
    names(peakStatus) <- maxInd_curr
    
    ## orphanRidgeList keep the ridges disconnected at certain scale level
    ## Changed by Pan Du 05/11/06
    orphanRidgeList <- NULL
    orphanRidgeName <- NULL
    nLevel <- length(colInd)
    
    for (j in 1:nLevel) {
      col.j <- colInd[j]
      scale.j <- scales[col.j]
      
      if (colInd[j] == skip) {
        oldname <- names(ridgeList)
        ridgeList <- lapply(ridgeList, function(x) c(x, x[length(x)]))
        ##peakStatus <- lapply(peakStatus, function(x) c(x, x[length(x)]))
        names(ridgeList) <- oldname
        ##names(peakStatus) <- oldname
        next
      }
      
      if (length(maxInd_curr) == 0) {
        maxInd_curr <- which(localMax[, col.j] > 0)
        next
      }
      
      ## The slide window size is proportional to the CWT scale
      ## winSize.j <- scale.j / 2 + 1
      winSize.j <- floor(scale.j/2)
      if (winSize.j < minWinSize) {
        winSize.j <- minWinSize
      }
      
      selPeak.j <- NULL
      remove.j <- NULL
      for (k in 1:length(maxInd_curr)) {
        ind.k <- maxInd_curr[k]
        start.k <- ifelse(ind.k-winSize.j < 1, 1, ind.k-winSize.j)
        end.k <- ifelse(ind.k+winSize.j > nMz, nMz, ind.k+winSize.j)
        ind.curr <- which(localMax[start.k:end.k, col.j] > 0) + start.k - 1
        ##ind.curr <- which(localMax[, col.j] > 0)
        if (length(ind.curr) == 0) {
          status.k <- peakStatus[[as.character(ind.k)]]
          ## bug  work-around
          if (is.null(status.k)) status.k <- gapTh +1
          ##
          if (status.k > gapTh & scale.j >= 2) {
            temp <- ridgeList[[as.character(ind.k)]]
            orphanRidgeList <- c(orphanRidgeList, list(temp[1:(length(temp)-status.k)]))
            orphanRidgeName <- c(orphanRidgeName, paste(col.j + status.k + 1, ind.k, sep='_'))
            remove.j <- c(remove.j, as.character(ind.k))
            next
          } else {
            ind.curr <- ind.k
            peakStatus[[as.character(ind.k)]] <- status.k + 1
          }
        } else {
          peakStatus[[as.character(ind.k)]] <- 0
          if (length(ind.curr) >= 2)  ind.curr <- ind.curr[which.min(abs(ind.curr - ind.k))]
        }
        ridgeList[[as.character(ind.k)]] <- c(ridgeList[[as.character(ind.k)]], ind.curr)
        selPeak.j <- c(selPeak.j, ind.curr)
      }
      ## Remove the disconnected lines from the currrent list
      if (length(remove.j) > 0) {
        removeInd <- which(names(ridgeList) %in% remove.j)
        ridgeList <- ridgeList[-removeInd]
        peakStatus <- peakStatus[-removeInd]
      }
      
      ## Check for duplicated selected peaks and only keep the one with the longest path.
      dupPeak.j <- unique(selPeak.j[duplicated(selPeak.j)])
      if (length(dupPeak.j) > 0) {
        removeInd <- NULL
        for (dupPeak.jk in dupPeak.j) {
          selInd <- which(selPeak.j == dupPeak.jk)
          selLen <- sapply(ridgeList[selInd], length)
          removeInd.jk <- which.max(selLen)
          removeInd <- c(removeInd, selInd[-removeInd.jk])
          orphanRidgeList <- c(orphanRidgeList, ridgeList[removeInd.jk])
          orphanRidgeName <- c(orphanRidgeName, paste(col.j, selPeak.j[removeInd.jk], sep='_'))
        }
        selPeak.j <- selPeak.j[-removeInd]
        ridgeList <- ridgeList[-removeInd]
        peakStatus <- peakStatus[-removeInd]
      }
      
      ## Update the names of the ridgeList as the new selected peaks
      ##if (scale.j >= 2) {
      if (length(ridgeList) > 0) names(ridgeList) <- selPeak.j
      if (length(peakStatus) > 0) names(peakStatus) <- selPeak.j
      ##}
      
      ## If the level is larger than 3, expand the peak list by including other unselected peaks at that level
      if (scale.j >= 2) {
        maxInd_next <- which(localMax[, col.j] > 0)
        unSelPeak.j <- maxInd_next[!(maxInd_next %in% selPeak.j)]
        newPeak.j <- as.list(unSelPeak.j)
        names(newPeak.j) <- unSelPeak.j
        ## Update ridgeList
        ridgeList <- c(ridgeList, newPeak.j)
        maxInd_curr <- c(selPeak.j, unSelPeak.j)
        ## Update peakStatus
        newPeakStatus <- as.list(rep(0, length(newPeak.j)))
        names(newPeakStatus) <- newPeak.j
        peakStatus <- c(peakStatus, newPeakStatus)
      } else {
        maxInd_curr <- selPeak.j
      }
    }
    
    ## Attach the peak level at the beginning of the ridge names
    if (length(ridgeList) > 0) names(ridgeList) <- paste(1, names(ridgeList), sep='_')
    if (length(orphanRidgeList) > 0) names(orphanRidgeList) <- orphanRidgeName
    ## Combine ridgeList and orphanRidgeList
    ridgeList <- c(ridgeList, orphanRidgeList)
    if (length(ridgeList) == 0) return(NULL)
    
    ## Reverse the order as from the low level to high level.
    ridgeList <- lapply(ridgeList, rev)
    ## order the ridgeList in increasing order
    ##ord <- order(selPeak.j)
    ##ridgeList <- ridgeList[ord]
    
    ## Remove possible duplicated ridges
    ridgeList <- ridgeList[!duplicated(names(ridgeList))]
    
    attr(ridgeList, 'class') <- 'ridgeList'
    attr(ridgeList, 'scales') <- scales
    return(ridgeList)
  }
