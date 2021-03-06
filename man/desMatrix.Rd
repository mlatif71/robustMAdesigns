% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/desMatrix.R
\name{desMatrix}
\alias{desMatrix}
\title{Generating a Design Matrix}
\usage{
desMatrix(layout, design = "DS")
}
\arguments{
\item{layout}{experimetnal layout of interest, e.g. "1x3" is for one-factor experiment with three treatments and "2x2" is for two-factor experiments with both the factors with two levels}

\item{design}{design from the layout of interest, e.g. dye-swarp ("DS"), circular loop ("CL"), etc.}
}
\value{
A data frame of which first two columnscorrespond to the dye effects cy3 and cy5, respectively and the remaining columns correspond to the number of treatments to be compared.
}
\description{
Generates a design matrix corresponding to an experimental layout
and a design
}
\examples{

# design matrix for dye-swap design from 1x3 layout
    desMatrix(layout="1x3", design="DS")

# design matrix for circular loop design from 2x2 layout
    desMatrix(layout="2x2", design="CL")

}
