\name{dcomp}
\alias{dcomp}
\title{dcomp}
\description{
Returns a density function of a user specified multisection composite distribution
}
\usage{
dcomp(xx,dists,par,borders,par.pos,buffer)
}
\arguments{
\item{xx}{Evaluation locations}
\item{dists}{A vector of strings stating the desired partial distributions in order}
\item{par}{A list of parameters, a vector of parameters for each partial distribution, with the first two being the interval limits and the second argument being the weights to be used}
\item{borders}{Optional: If the distribution has to have continuous and differentiable catentation points, the user specifies a list for each of area following the first, containing a range for a parameter of the following partial distribution to lie within.}
\item{par.pos}{Optional: If 'borders' is non empty, e.g a smooth function is desired, here the vector of parameter positions that need to be changed is specified. Default value is 1, meaning the first parameter for each partial distribution is amended}
\item{buffer}{Optional: A two dimensional vector, containing the values for upper and lower buffer from the respective catenation points during optimization}
}
\value{An object of the same length as \code{xx}, giving the density values}
\author{Martin Wiegand, Saralees Nadarajah}
\examples{par<-list()
distvec<-c("lnorm","gamma")
par[[1]]<-c(0,1,Inf)
par[[2]]<-c(1)
par[[3]]<-c(0,1)
par[[4]]<-c(1,1)

x<-seq(0,3,0.01)
# non-continuous case
y1<-dcomp(x,distvec,par)
# continuous case
y2<-dcomp(x,distvec,par,borders=list(c(0.00001,10)),buffer=c(10e-5,0))
 
par(mfrow=c(1,2),oma=rep(0,4))
xrange<-range(x)
yrange<-range(y1,y2)
plot(x,y1,type="l",xlab="x",ylab="Density function",xlim=xrange,ylim=yrange)
abline(v=1)
plot(x,y2,type="l",xlab="x",ylab="Density function",xlim=xrange,ylim=yrange)
abline(v=1)}
\keyword{Density function}
