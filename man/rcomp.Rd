\name{rcomp}
\alias{rcomp}
\title{rcomp}
\description{
Returns a random sample of size n of a user specified multisection composite distribution
}
\usage{
rcomp(nn,dists,par,borders,par.pos,buffer)
}
\arguments{
\item{nn}{Desired random sample size}
\item{dists}{A vector of strings stating the desired partial distributions in order}
\item{par}{A list of parameters, a vector of parameters for each partial distribution, with the first two being the interval limits and the second argument being he weights to be used}
\item{borders}{Optional: If the distribution has to have continuous and differentiable catentation points, the user specifies a list for each of area following the first, containing a range for a parameter of the following partial distribution to lie within.}
\item{par.pos}{Optional: If 'borders' is non empty, e.g a smooth function is desired, here the vector of parameter positions that need to be changed is specified. Default value is 1, meaning the first parameter for each partial distribution is amended}
\item{buffer}{Optional: A two dimensional vector, containing the values for upper and lower buffer from the respective catenation points during optimization}
}
\value{An object of length  \code{nn}, giving the random numbers}
\author{Martin Wiegand, Saralees Nadarajah}
\examples{par<-list()

distvec<-c("lnorm","gamma")

par[[1]]<-c(0,1,Inf)

par[[2]]<-c(1)

par[[3]]<-c(0,1)

par[[4]]<-c(1,1)


n<-1000

# non-continuous case

y1<-rcomp(n,distvec,par)

# continuous case

y2<-rcomp(n,distvec,par,borders=list(c(0.00001,10)),buffer=c(10e-5,0))


par(mfrow=c(1,2),oma=rep(0,4))

hist(y1,nclass=10,xlab="x",ylab="Frequency",main="")

hist(y2,nclass=10,xlab="x",ylab="Frequency",main="")}
\keyword{Random number generation}
