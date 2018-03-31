\name{cytofData}
\docType{data}
\alias{cytofData}
\alias{cytofData-package}
\alias{bodenmiller2012bcrxl}
\title{SummarizedExperiment with marker abundances from Bodenmiller et al (2012)}
\description{
  A SummarizedExperiment object with marker abundances for 35 markers, prepared 
  according to the steps outlined in \code{inst/scripts/make-data.R}
}
\usage{
  bodenmiller2012bcrxl(metadata=FALSE)
}
\arguments{
  \item{metadata}{
    \code{logical} value indicating whether metadata only should be returned
    or if the resource should be loaded. Default behavior(metadata=FALSE) 
    loads the data.
  }
}
\examples{
  bodenmiller2012bcrxl()
}
\value{A SummarizedExperiment object}
\format{SummarizedExperiment}
\source{FCS files}
\references{
  Bodenmiller et al
}
\keyword{datasets}
