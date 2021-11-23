load("XPN.AFX.rda")
library(fastcluster)
library(skmeans)
library(varbvs)
library(phonTools)
library(calibrate)
dyn.load("HARMONY_MLE.so")
is.loaded("HARMONY_MLE_C_C")
harmony.afx.static <- function(x,
                               K=10,
                               L=4,
                               p1.names=0,
                               p2.names=0,
                               gene.cluster="skmeans",
                               assay.cluster="hclust",
                               corr="pearson",
                               iterations=1,
                               skip.match=FALSE,
                               is.assays.identical=FALSE,
                               output.iterations=FALSE) {
  harmony.result <- harmony.static(
    x,
    XPN.AFX.RDA,
    K,
    L,
    p1.names,
    p2.names,
    gene.cluster,
    assay.cluster,
    corr,
    iterations,
    skip.match,
    is.assays.identical,
    output.iterations)
  
  x.vector <- as.vector(as.matrix(harmony.result$x))
  y.vector <- as.vector(as.matrix(XPN.AFX.RDA))
  x.mean <- mean(x.vector)
  y.mean <- mean(y.vector)
  x.sd <- sd(x.vector)
  y.sd <- sd(y.vector)
  
  return (as.matrix((harmony.result$x - x.mean) / x.sd * y.sd + y.mean))
}

harmony.static = function(platform1.data,
                          platform2.data,
                          K=10,
                          L=4,
                          p1.names=0,
                          p2.names=0,
                          gene.cluster="kmeans",
                          assay.cluster="kmeans",
                          corr="pearson",
                          iterations=30,
                          skip.match=FALSE,
                          is.assays.identical=FALSE,
                          output.iterations=FALSE){
  #If K or L is not a single value, it is taken as a list of possible values.
  
  if (is.assays.identical) {
    cat("Switching to assays identical mode\n")
  }
  
  if (is.assays.identical && ncol(platform1.data) != ncol(platform2.data)) {
    stop("Number of columns must be equal for is.assays.identical = TRUE")
  }
  
  #Match names
  input = processplatforms(list(x=platform1.data,y=platform2.data),namesvec = c(p1.names, p2.names), skip.match=skip.match)
  x <- input$x
  y <- input$y
  input <- NA
  
  #Remove the medians
  x_demed = x - rowMedians(as.matrix(x))
  y_demed = y - rowMedians(as.matrix(y))
  
  #Get the dimensions
  nx=dim(x)[2]
  ny=dim(y)[2]
  mx=dim(x)[1]
  my=dim(y)[1]
  
  #Create the combined dataframe for clustering purposes
  combined <- cbind(x_demed,y_demed)
  
  #Detect K and L if necessary
  K = detect_K(combined,K,corr)
  L = detect_L(x_demed,y_demed,L,corr)
  
  xout = 0
  yout = y
  
  iterations.data <- list()
  for(iter in 1:iterations){
    cat("HARMONY iteration",iter,"out of",iterations,"\n")
    
    #Do the assay clustering
    if (is.character(assay.cluster))
    {
      if (is.assays.identical) {
        assayclusters = harmonyassaycluster(.5*(x_demed + y_demed), y_demed, L, assay.cluster, corr)[1:nx]
      } else {
        assayclusters = harmonyassaycluster(x_demed, y_demed, L, assay.cluster, corr)
      }
    }
    else
    {
      assayclusters = assay.cluster
    }
    
    if (is.assays.identical) {
      assayclusters = c(assayclusters, assayclusters)
      names(assayclusters) <- c(colnames(x), colnames(y))
    }
    
    #Do the gene clustering
    if (is.character(gene.cluster))
    {
      geneclusters = harmonygenecluster(combined, K, gene.cluster, corr)
    }
    else
    {
      geneclusters = gene.cluster
    }
    
    #Estimate the HARMONY parameters
    xparams = harmonyparamsp(x,geneclusters,assayclusters[1:nx])
    yparams = harmonyparamsp(y,geneclusters,assayclusters[(nx+1):(nx+ny)])
    
    #Calcuate the weighted averages
    nor = 1/(nx + ny)
    Aav = try((1/(xparams$nblock+yparams$nblock))*(xparams$nblock*xparams$A + yparams$nblock*yparams$A))
    bav = nor*(nx*xparams$b + ny*yparams$b)
    cav = nor*(nx*xparams$c + ny*yparams$c)
    sigmaav = sqrt(nor*(nx*xparams$s2 + ny*yparams$s2))
    sigmax = sqrt(xparams$s2)
    #sigmay = sqrt(yparams$s2)
    
    
    #Calculate the expanded A
    expAx = harmonyclusterexpand(xparams$A,geneclusters,assayclusters[1:nx])
    #expAy = harmonyclusterexpand(yparams$A,geneclusters,assayclusters[(nx+1):(nx+ny)])
    
    #Calculate the residuals
    epsilonx = as.matrix(x) - (as.vector(xparams$b) * expAx + kronecker(ones(1,nx),xparams$c))
    #epsilony = as.matrix(y) - (as.vector(yparams$b) * expAy + kronecker(ones(1,ny),yparams$c))
    
    #Calculate the expanded average A
    expAavx = harmonyclusterexpand(Aav,geneclusters,assayclusters[1:nx])
    #expAavy = harmonyclusterexpand(Aav,geneclusters,assayclusters[(nx+1):(nx+ny)])
    
    #Calculate the output values 
    xcurout = ((as.vector(bav) * expAavx) + kronecker(ones(1,nx),cav) + as.vector(sigmaav/sigmax) * epsilonx)
    ycurout = 0 #((as.vector(bav) * expAavy) + kronecker(ones(1,ny),cav) + as.vector(sigmaav/sigmay) * epsilony)
    
    xout = xout + (1/iterations) * xcurout
    yout = yout
    #+ (1/iterations) * ycurout
    
    if (output.iterations) {
      iterations.data[[iter]] <- list(
        assay.cluster = assayclusters,
        gene.cluster = geneclusters,
        xparams = xparams,
        yparams = yparams,
        x = xcurout,
        y = ycurout
      )
    }
    
  }#end of the enclosing for loop
  
  #Put the rownames back in and convert to data frames
  xout = as.data.frame(xout,row.names=rownames(x))
  yout = as.data.frame(yout,row.names=rownames(y))
  
  #All done!
  if (!output.iterations)
    return(list(x=xout, y=yout))
  
  return(list(x=xout, y=yout, iterations = iterations.data))
}

processplatforms = function(datalist, namesvec=NULL, skip.match=FALSE){
  #Convert data from various formats to the proper format for use 
  #with all the crossnorm normalization functions
  
  for(i in 1:length(datalist)){
    if(is.matrix(datalist[[i]])){
      datalist[[i]] <- as.data.frame(datalist[[i]])
    }
  }
  
  if (is.null(namesvec)){
    namesvec <- numeric(length(datalist))
    for (i in 1:length(datalist)){
      namesvec[i] <- 0
    }
  }
  
  #Put the row names in their places
  for (i in 1:length(namesvec)){
    if(namesvec[i] != 0){
      rownames(datalist[[i]]) = datalist[[i]][,namesvec[i]]
      datalist[[i]] = datalist[[i]][,-1*namesvec[i]]
    }	
  }
  
  if(!skip.match){
    #Create the common genes list
    commongenes <- rownames(datalist[[1]])
    for (i in 2:length(datalist)){
      commongenes <- intersect(commongenes,rownames(datalist[[i]]))
    }
    
    #Put it all together
    for (i in 1:length(datalist)){
      l <- which(rownames(datalist[[i]])%in%commongenes)
      datalist[[i]] <- datalist[[i]][l,]
    }
  }
  return(datalist)
}

rowMedians = function(aMatrix, na.rm = FALSE){
  #Also works on data.frames
  m = dim(aMatrix)[1]
  ret = numeric(m)
  for (i in 1:m){
    ret[i] = mean(median(aMatrix[i,],na.rm=na.rm),na.rm=TRUE)
  }
  return(ret)
}

detect_K = function(x,K,corr="pearson"){
  
  #Determine K if a range was given
  if(length(K)>1 & length(K) != dim(x)[1]){
    xdiss = 1 - cor(t(x),method=corr)
    pamklist=pamk(xdiss,krange=K,diss=TRUE, keep.diss=FALSE, keep.data=FALSE)
    K = pamklist$nc
    genepamobj = pamklist$pamobject
  }else{
    genepamobj = NULL	
  }
  return(K)
}

detect_L = function(x,y,L,corr="pearson"){
  
  #Determine L if a range has been given
  if(length(L)>1 & length(L) != (dim(x)[2]+dim(y)[2])){
    xdiss = 1 - cor(as.matrix(x), method = corr)
    ydiss = 1 - cor(as.matrix(y), method = corr)
    Lx=pamk(xdiss,krange=L,diss=TRUE, keep.diss=FALSE, keep.data=FALSE)$nc
    Ly = pamk(ydiss,krange=L,diss=TRUE, keep.diss=FALSE, keep.data=FALSE)$nc
    L = max(Lx,Ly)
    cat(as.character(L), "assay clusters detected...\n")
  }
  return(L)
}

harmonyassaycluster = function(x,y,L,cluster="kmeans",corr="pearson"){
  
  #The numbers of assays
  nx = ifelse(is.null(x), 0, dim(x)[2])
  ny = ifelse(is.null(y), 0, dim(y)[2])
  
  #The number of genes
  m = ifelse(is.null(x), dim(y)[1], dim(x)[1])
  
  if (L == 1)
    return (rep(1, nx + ny))
  
  if(is.numeric(cluster) | is.factor(cluster)){
    return(as.numeric(cluster))
  }
  
  if ((nx == 0 || ny == 0)
      && cluster != "kmeans"
      && cluster != "skmeans"
      && cluster != "hclust"
      && cluster != "random")
    stop(paste0("Assay clusterization type \"", cluster, "\" is not supported"))
  
  #Compute the two correlation matrices if necessary
  if (cluster == 'pam' | (length(L)>1)){
    if (nx > 0)
      xdiss = 1 - cor(as.matrix(x), method = corr)
    if (ny > 0)
      ydiss = 1 - cor(as.matrix(y), method = corr)
  }
  
  #Determine L if a range has been given
  if(length(L)>1){
    Lx=pamk(xdiss,krange=L,diss=TRUE, keep.diss=FALSE, keep.data=FALSE)$nc
    Ly = pamk(ydiss,krange=L,diss=TRUE, keep.diss=FALSE, keep.data=FALSE)$nc
    L = max(Lx,Ly)
    cat(as.character(L), "assay clusters detected...\n")
  }
  
  #Do the clustering
  if (cluster == "classic"){
    
    done = FALSE
    while (!done){
      xyclust <- ckmeans(t(cbind(x,y)),centers=L, iter.max=20, distance=corr)
      #print(xyclust$cluster)
      xclust <- xyclust$cluster[1:nx]
      yclust <- xyclust$cluster[(nx+1):(nx+ny)]
      #print(xclust)
      #			print(yclust)
      #			print(nlevels(as.factor(xclust)))
      #			print(nlevels(as.factor(yclust)))
      done <- nlevels(as.factor(xclust)) == L & nlevels(as.factor(yclust)) == L
      if(!done){
        cat("Not all platforms contained all assay clusters.  Trying again (only \"classic\" mode has this problem)...\n")
      }
    }
    return(xyclust$cluster)
  }
  else if (cluster == "pam"){
    xclust = pam(xdiss, k=L, diss=TRUE, keep.diss=FALSE, keep.data=FALSE,cluster.only=TRUE)
    yclust = pam(ydiss, k=L, diss=TRUE, keep.diss=FALSE, keep.data=FALSE,cluster.only=TRUE)
  }else if (cluster == "kmeans"){
    
    if (nx > 0) {
      #This loop ensures there are no empty clusters.
      done = FALSE
      while (!done){
        xclust = ckmeans(t(x), centers=L, iter.max = 20, distance = corr)
        
        done = !xclust$empty
      }
      xclust = xclust$cluster
    }
    
    if (ny > 0) {
      done = FALSE
      while (!done){
        yclust = ckmeans(t(y), centers=L, iter.max = 20, distance = corr)
        
        done = !yclust$empty
      }
      yclust = yclust$cluster	
    }
  } else if (cluster == "flexclust"){
    
    distCor1 = function (x, centers) 
    {
      z <- matrix(0, nrow(x), ncol = nrow(centers))
      for (k in 1:nrow(centers)) {
        z[, k] <- 1 - cor(t(x), centers[k, ], method=corr)
      }
      z
    }
    #This loop ensures there are no empty clusters.
    done = FALSE
    while (!done){
      xclust = kcca(t(x), k=L, family=kccaFamily(dist=distCor1,cent="centMean"), simple=TRUE)
      yclust = kcca(t(y), k=L, family=kccaFamily(dist=distCor1,cent="centMean"), simple=TRUE)
      
      done = (dim(xclust@clusinfo)[1] == L) &(dim(yclust@clusinfo)[1] == L) 
    }
    xclust = xclust@cluster
    yclust = yclust@cluster	
  }else if(cluster == "random"){
    if (nx > 0)
      xclust = sample(c(1:L,sample(1:L,nx-L,replace = TRUE)),nx,replace=FALSE)
    
    if (ny > 0)
      yclust = sample(c(1:L,sample(1:L,ny-L,replace = TRUE)),ny,replace=FALSE)
  } else if (cluster == 'hclust') {
    sd.zeros = (nx > 0 && any(apply(x, 2, sd) < .Machine$double.eps ^ 0.5)) ||
      (ny > 0 && any(apply(y, 2, sd) < .Machine$double.eps ^ 0.5))
    if (sd.zeros) {
      cat("Identical values in columns are detected (after substracting row medians).\n")
      cat("Switching to assays random clustering.\n")
      return (harmonyassaycluster(x, y, L, "random", corr))
    }
    
    if (nx > 0)
      xclust = cutree(hclust.vector(t(x), method="ward"), k = L)
    
    if (ny > 0)
      yclust = cutree(hclust.vector(t(y), method="ward"), k = L)
  }
  
  if (nx > 0 && ny > 0) {
    #Compute cluster averages
    xave = matrix(NA,m,L)
    yave = matrix(NA,m,L)
    for (i in 1:L){
      xinds = xclust==i
      yinds = yclust==i
      xave[,i] = rowMeans(as.matrix(x[,xinds],nrow=m,ncol=sum(as.numeric(xinds))))
      yave[,i] = rowMeans(as.matrix(y[,yinds],nrow=m,ncol=sum(as.numeric(yinds))))
    }
    
    #Compute the cluster correlation matrix
    clustercor = cor(xave,yave, method = corr)
    clustercor[is.na(clustercor)] <- 1
    
    #Map clusters
    xtoymap = matrix(NA,L,1)
    for (i in 1:L){
      highest = which.max(clustercor)
      xind = highest%%L
      
      yind = highest%/%L + 1
      
      if (xind==0){
        xind = L
        yind = yind-1
      }
      xtoymap[xind]=yind
      clustercor[xind,] = -2
      clustercor[,yind] = -2
    }
    
    #Change the x clusters using the map
    newxclust = xclust
    for (i in 1:L){
      newxclust[xclust==i] = xtoymap[i]
    }
    xclust = newxclust
  }
  
  rr <- c()
  if (nx > 0)
    rr <- c(rr, xclust)
  
  if (ny > 0)
    rr <- c(rr, yclust)
  
  #Return the clustering in a combined vector
  return (rr)
}

harmonygenecluster = function(x,K,cluster="pam",corr="pearson"){
  if (K == 1)
    return (rep(1, nrow(x)))
  
  #Dimensions
  m = dim(x)[1]
  n = dim(x)[2]
  
  #Compute the dissimilarity matrix
  if (cluster=="pam"|(length(K)>1)){
    xdiss = 1 - cor(t(x),method=corr)
  }
  
  #Determine K if a range was given
  if(length(K)>1 ){
    pamklist=pamk(xdiss,krange=K,diss=TRUE, keep.diss=FALSE, keep.data=FALSE)
    K = pamklist$nc
    genepamobj = pamklist$pamobject
  }else{
    genepamobj = NULL	
  }
  
  
  #Do the gene clustering
  if(cluster=="pam"){
    if (!is.null(genepamobj)){
      geneclusters = genepamobj$clustering
    }else{
      geneclusters = pam(xdiss, k=K, diss=TRUE, keep.diss=FALSE, keep.data=FALSE, cluster.only=TRUE)
    }
  }else if(cluster=="kmeans"){
    #This loop ensures there are no empty clusters.
    done = FALSE
    while (!done){
      geneclusters = ckmeans(x,centers=K,iter.max=1000,distance=corr)
      done = !geneclusters$empty			
    }
    geneclusters = geneclusters$cluster
  }else if(cluster=="flexclust"){
    
    
    distCor1 = function (x, centers) #change
    {
      z <- matrix(0, nrow(x), ncol = nrow(centers))
      for (k in 1:nrow(centers)) {
        z[, k] <- 1 - cor(t(x), centers[k, ], method=corr)
      }
      z
    }
    #This loop ensures there are no empty clusters.
    done = FALSE
    while (!done){
      geneclusters = kcca(x, k=K, family=kccaFamily(dist=distCor1,cent="centMean"), simple=TRUE) #change
      #geneclusters = ckmeans(x,centers=K,iter.max=1000,distance=corr)
      done = (dim(geneclusters@clusinfo)[1]==K)#change		
    }
    geneclusters = geneclusters@cluster
  }else if (cluster == "random"){
    geneclusters = sample(c(1:K,sample(1:K,m-K,replace = TRUE)),m,replace=FALSE)
  }else if (cluster == "hclust"){
    geneclusters = cutree(hclust.vector(x, method="ward"), k = K)
  }else if (cluster == "skmeans"){
    zero.genes = names(which(apply(x, 1, sd) < .Machine$double.eps ^ 0.5))
    non.zero.genes = setdiff(rownames(x), zero.genes)
    non.zero.clusters = skmeans(as.matrix(x[non.zero.genes,]), K, method = "genetic", control=list(start="S"))$cluster
    
    geneclusters = rep(0, nrow(x))
    names(geneclusters) = rownames(x)
    geneclusters[zero.genes] = 1
    geneclusters[non.zero.genes] = non.zero.clusters
  }else if (cluster == "kmeanspp"){
    geneclusters = kmeanspp(x, K)$cluster
  }else if (cluster == "hclustmem"){
    stop("Non supported clusterization method")
    
    # find min and max
    x.center = apply(x, 2, mean)
    x.dist.to.center = apply(x - x.center, 1, function(x) sum(x ^ 2))
    x.centroid.1 = x[which.min(x.dist.to.center),]
    x.centroid.2 = x[which.max(x.dist.to.center),]
    x.centroids = as.matrix(rbind(x.centroid.1, x.centroid.2))
    
    distCor1 = function (x, centers) #change
    {
      z <- matrix(0, nrow(x), ncol = nrow(centers))
      for (k in 1:nrow(centers)) {
        z[, k] <- 1 - cor(t(x), centers[k, ], method=corr)
      }
      z
    }
    
    geneclusters = kcca(x, k=x.centroids, family=kccaFamily(dist=distCor1,cent="centMean"), simple=TRUE) #change
    geneclusters = geneclusters@cluster
    
    first.genes = names(which(geneclusters == 1))
    second.genes = names(which(geneclusters == 2))
    
    first.K = round(K * length(first.genes) / m)
    second.K = K - first.K
    
    hgeneclusters.1 = cutree(hclust.vector(x[first.genes,], method="ward"), k = first.K)
    hgeneclusters.2 = cutree(hclust.vector(x[second.genes,], method="ward"), k = second.K)
    
    geneclusters = rep(0, nrow(x))
    names(geneclusters) = rownames(x)
    geneclusters[first.genes] = hgeneclusters.1
    geneclusters[second.genes] = hgeneclusters.2 + first.K
  } else{
    cat("Unknown gene clustering method:",cluster,"\n")
  }
  
  return(geneclusters)
  
}

harmonyparamsp = function(x, geneclusters, assayclusters, method){
  #x should be a dataframe or hashdataframe	
  
  #Get K and L
  K=max(geneclusters)
  L=max(assayclusters)
  
  #Get number of assays and genes
  numassays=length(x)
  numgenes=length(rownames(x))
  
  #Number of assays in each cluster
  nj = matrix(table(as.factor(assayclusters)),1,L)
  
  #Set up the output variables
  A = matrix(nrow=K,ncol=L)
  nblock = matrix(nrow=K,ncol=L)
  b = matrix(nrow=numgenes, ncol=1)
  c = matrix(nrow=numgenes, ncol=1)
  sigma = matrix(nrow=numgenes, ncol=1)
  mu = matrix(nrow=numgenes,ncol=L)
  
  
  #For each gene cluster, estimate the HARMONY parameters
  for (a in 1:K){
    
    #Get the logical indexing vector for gene cluster a 
    geneinds = geneclusters==a
    numgenesa = sum(as.numeric(geneinds))
    
    #Get the data for this gene cluster
    xa = as.matrix(x[geneinds,])
    
    #Get the number of genes in this cluster
    Ga = dim(xa)[1]
    S = dim(xa)[2]
    
    #For each assay cluster and each gene, get the (relatively) unconstrained
    #MLE for mu
    mua = matrix(nrow=Ga,ncol=L)
    for (bb in 1:L){
      assayinds = assayclusters==bb
      xab = matrix(xa[,assayinds],numgenesa,sum(as.numeric(assayinds)))
      mua[,bb]=t(t(rowMeans(xab)))
    }
    
    
    #Calculate sigma2
    expmua = harmonyclusterexpand(mua, colclusters = assayclusters)
    sigma2 = xa - expmua
    sigma2 = sigma2*sigma2
    sigma2 = harmonyclustercollapse(sigma2,colclusters = assayclusters)
    ## HARMONY_MLE
    solna = HARMONY_MLE_C(mua,sigma2,nj)
    
    
    ba = solna$b
    Aa = solna$A
    ca = solna$c
    s2a = solna$s2
    
    #Write the values for this gene group into the larger arrays
    sigma[geneinds,] = s2a
    mu[geneinds,] = mua
    A[a,] = Aa
    nblock[a,] = nj
    b[geneinds,] = ba
    c[geneinds,] = ca
  }
  
  return(list(A=A,nblock=nblock,b=b,c=c,s2=sigma,mu=mu))
  
}

harmonyclusterexpand = function(aMatrix, rowclusters = 1:(dim(aMatrix)[1]), colclusters = 1:(dim(aMatrix)[2]), noise = 0){
  #This is basically a fancy repmat.  rowclusters and colclusters 
  #as produced by pam.  aMatrix must have the same number of rows
  #and columns as there are row and column clusters.
  
  l = dim(aMatrix)[1]
  w = dim(aMatrix)[2]
  
  #Get the dimensions of the output
  m = length(rowclusters)
  n = length(colclusters)
  
  #Get the number of row and column clusters
  nrowclust = max(rowclusters)
  ncolclust = max(colclusters)
  
  #initialize the vertically expanded matrix
  vert = matrix(0,m,dim(aMatrix)[2])
  
  #Do the vertical expansion
  for (i in 1:m){
    vert[i,]=aMatrix[rowclusters[i],] + (noise*randn(1,w))
  }
  
  #initialize the fully expanded matrix
  output = matrix(NA,m,n)
  
  #Do the horizontal expansion
  for (j in 1:n){
    output[,j] = vert[,colclusters[j]] + (noise*randn(m,1))
  }
  
  #That is it!
  return(output)
}

harmonyclustercollapse = function(aMatrix, rowclusters = 1:(dim(aMatrix)[1]), colclusters = 1:(dim(aMatrix)[2])){
  #This undoes the work of harmonyclusterexpand
  
  l = dim(aMatrix)[1]
  w = dim(aMatrix)[2]
  
  #collapse the rows
  rowidx = split(1:l,rowclusters)
  lc = length(rowidx)
  collapse1 = matrix(NA,lc,w)
  for ( i in 1:lc){
    collapse1[i,] =  colMeans(matrix(aMatrix[rowidx[[i]],],length(rowidx[[i]]),w))
  }
  
  #collapse the columns
  colidx = split(1:w,colclusters)
  wc = length(colidx)
  collapse2 = matrix(NA,lc,wc)
  for ( i in 1:wc){
    collapse2[,i] =  rowMeans(matrix(collapse1[,colidx[[i]]],lc,length(colidx[[i]])))
  }
  
  return(collapse2)
  
}

HARMONY_MLE_C = function(xbar,sigma2,nj){
  
  
  n = sum(nj)
  I = dim(xbar)[1]
  J = dim(xbar)[2]
  
  A = zeros(1,J)
  b = ones(I,1)
  c = zeros(I,1)
  s2 = ones(I,1)
  
  cout = .C("HARMONY_MLE_C_C",as.double(xbar),A=as.double(A),b=as.double(b),c=as.double(c),s2=as.double(s2),as.double(sigma2),as.integer(I),as.integer(J),as.integer(n),as.integer(nj))
  
  A = matrix(cout$A,1,J)
  b = matrix(cout$b,I,1)
  c = matrix(cout$c,I,1)
  s2 = matrix(cout$s2,I,1)
  
  return(list(A=A,b=b,c=c,s2=s2))
}