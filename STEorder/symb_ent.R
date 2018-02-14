symbol <- function(x) {
  # x is a short vector
  s <- c()
  for (i in 1:(length(x)-1)) {
    s <- c(s, as.character(as.numeric((x[i] > x[i+1]))))
  }
  paste(s,collapse="")
}

symbolTransform <- function(x, ord) {
  # takes a numeric vector
  # and returns a sequence of symbols
  # ord is the markov order
  symbList <- c()
  for (i in 1:(length(x)-ord)) {
    symbList <- c(symbList, symbol(x[i:(i+ord)]))    
  }
  symbList
}

###################
# Symbolic Transfer Entropy
# following paper: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.100.158101

buildJoint <- function(xsym, ysym, n, l) {
  jointCounts <- table(xsym[1:(n-l)], ysym[1:(n-l)], xsym[l:n])
  jointProb <- jointCounts/sum(jointCounts)
  jointProb
}

safediv <- function(a,b) {
  # divide a/b
  # if b != 0
  # a can be a vector
  if (b == 0) {
    return(rep(0, length(a)))
  } else {
    a/b
  }
}

buildCond1 <- function(xsym, ysym, n, l){
  a <- xsym[1:(n-l)]
  b <- ysym[1:(n-l)]
  c <- xsym[l:n]
  countxy <- table(a,b)  # given x_i and y_i
  countxyz <- table(a,b,c)  # prob of x_(i+l)
  for (si in unique(a)){
    for (sj in unique(b)) {
      countxyz[si,sj,] <- safediv(countxyz[si,sj,], countxy[si,sj])
    }
  }
  countxyz
}

buildCond2 <- function(xsym, n, l){
  # want P(x_l | x)
  a <- xsym[1:(n-l)]
  b <- xsym[l:n]
  countxx <- table(a,b)
  countx  <- rowSums(countxx) 
  for (si in unique(a)) {
    countxx[si,] <- safediv(countxx[si,], countx[si])
  }
  countxx
}


ste <- function(y,x,ord,l) {
  # x,y are time series
  # ord is the markov order or embedding dimension
  # here we assume the direction is y -> x (y influencing x)
  # l is the time delay (for y -> time delay l -> x)
  # predicting x out into the future using current x and y.
  xsym <- symbolTransform(x,ord)
  ysym <- symbolTransform(y,ord)
  n <- length(xsym)
  jointProb <- buildJoint(xsym,ysym,n,l)  
  condProbTop <- buildCond1(xsym,ysym,n,l) 
  condProbBot <- buildCond2(xsym,n,l)
  a <- xsym[1:(n-l)]
  b <- ysym[1:(n-l)]
  c <- xsym[l:n]
  te <- 0
  for (i in 1:(n-l+1)) {
    #     sum of  joint(a,b,c) * log2  P(c | a,b) /  P(c | a)
    te <- te + jointProb[a[i],b[i],c[i]] * log2( safediv(condProbTop[a[i],b[i],c[i]], condProbBot[a[i],c[i]]) )
  }
  te
}
