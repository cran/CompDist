\name{par.fit}
\alias{par.fit}
\title{par.fit}
\description{
Returns the parameters fitted to a random sample along with a number of error measures, such as the log likelihood, AIC, BIC, AICc, CAIC and HQC.
}
\usage{
par.fit(data,dists,par,borders,par.pos,optit,buffer,cont)
}
\arguments{
\item{data}{Data set to be fitted to the distribution}
\item{dists}{A vector of strings stating the desired partial distributions in order}
\item{par}{A list of parameters, a vector of parameters for each partial distribution, with the first two being the interval limits and the second argument being he weights to be used}
\item{borders}{Optional: If the distribution has to have continuous and differentiable catentation points, the user specifies a list for each of area following the first, containing a range for a parameter of the following partial distribution to lie within.}
\item{par.pos}{Optional: If 'borders' is non empty, e.g a smooth function is desired, here the vector of parameter positions that need to be changed is specified. Default value is 1, meaning the first parameter for each partial distribution is amended}
\item{optit}{Number of iteration loops over the parameter optimisation}
\item{buffer}{Optional: A two dimensional vector, containing the values for upper and lower buffer from the respective catenation points during optimization}
\item{cont}{Logical value for smooth catenation points. Default FALSE.}
}
\value{Gives parameter estimates and values of the log likelihood, AIC, BIC, AICc, CAIC and HQC.}
\author{Martin Wiegand, Saralees Nadarajah}
\examples{# Generate random data


par<-list()

distvec<-c("lnorm","gamma")

par[[1]]<-c(0,1,Inf)

par[[2]]<-c(1)

par[[3]]<-c(0,1)

par[[4]]<-c(1,1)


n<-1000

# non-continuous case

r1<-rcomp(n,distvec,par)

# continuous case

r2<-rcomp(n,distvec,par,borders=list(c(0.00001,10)),buffer=c(10e-5,0))


# Initial Guess



par<-list()

distvec<-c("lnorm","gamma")

par[[1]]<-c(0,1,Inf)

par[[2]]<-c(1)

par[[3]]<-c(0,0.5)

par[[4]]<-c(0.5,1)



# Fitting



# non-continuous case

estimate1<-par.fit(r1,distvec,par,optit=1)

# continuous case

estimate2<-par.fit(r2,distvec,par,borders=list(c(0.00001,10)),optit=1,buffer=c(10e-5,0),cont=TRUE)



x<-seq(0,30,0.01)

# non-continuous case

y1<-dcomp(x,distvec,estimate1$Parameter)

# continuous case

y2<-dcomp(x,distvec,estimate2$Parameter,borders=list(c(0.00001,10)),buffer=c(10e-5,0))



par(mfrow=c(1,2),oma=rep(0,4))

hist(r1,probability=TRUE,breaks=40,main="",xlab="Data",ylab="Fitted density")

lines(x,y1,col="red")

hist(r2,probability=TRUE,breaks=40,main="",xlab="Data",ylab="Fitted density")

lines(x,y2,col="red")



estimate1

estimate2}
\keyword{Data fitting}
