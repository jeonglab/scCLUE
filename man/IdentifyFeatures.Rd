% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scCLUE.R
\name{IdentifyFeatures}
\alias{IdentifyFeatures}
\title{Identify feature genes (highly variable genes)}
\usage{
IdentifyFeatures(inData = data, pFeat = 0.2)
}
\arguments{
\item{inData}{The gene expression matrix in log scale. Each row corresponds to the gene and column corresponds to the cell}

\item{pFeat}{Percentage of potential features to be selected}
}
\description{
Identify feature genes (highly variable genes)
}
