% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/profile_gene.R
\name{profile_gene}
\alias{profile_gene}
\title{Plot Gene Expression Profiles in a Data Frame}
\usage{
profile_gene(
  data,
  scale = "none",
  x.title = "Sample/conditions",
  y.title = "Value",
  text.size = 15,
  text.angle = 45
)
}
\arguments{
\item{data}{A data frame, where rows are genes and columns are features/conditions.}

\item{scale}{The way to to scale the data. If \code{none} (default), no scaling. If \code{row}, the data is scalaed independently. If \code{all}, all the data is scaled as a whole.}

\item{x.title, y.title}{X-axis title and Y-axis title respectively.}

\item{text.size}{The size of axis title and text.}

\item{text.angle}{The angle of axis text.}
}
\value{
An image of ggplot.
}
\description{
Plot Gene Expression Profiles in a Data Frame
}
\examples{
# See examples in the function "spatial_enrich".
}
\references{
H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
\cr Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/.
}
\seealso{
\code{spatial_enrich}
}
\author{
Jianhai Zhang \email{jzhan067@ucr.edu; zhang.jianhai@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@ucr.edu}
}
