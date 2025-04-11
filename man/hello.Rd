\name{LF_VCR}
\alias{LF_VCR}
\alias{LF_VCR.ae}
\title{LF-VCR: Latent Factor Extraction with Sparse Group Lasso and Autoencoder Feature Extraction}
\usage{
LF_VCR(X, number.K, covariate, matrix, categorical, y)
LF_VCR.ae(X, number.K, covariate, matrix, categorical, y)
}
\description{
\strong{LF_VCR} extracts latent factors from the data and fits them using a sparse group lasso approach. \cr\cr
\strong{LF_VCR.ae} uses an autoencoder to perform feature extraction prior to subsequent analysis.
}
\arguments{
  \item{X}{A data matrix.}
  \item{number.K}{The number of latent factors to extract.}
  \item{covariate}{The covariate to be adjusted.}
  \item{matrix}{A character string indicating whether thresholding is applied to the correlation or covariance matrix.}
  \item{categorical}{A logical flag denoting whether the outcome \code{y} is continuous (\code{FALSE}) or categorical (\code{TRUE}).}
  \item{y}{The outcome variable.}
}
\details{
The function \code{LF_VCR} implements a two-stage procedure: first, it extracts latent factors from the input data, and then it fits these factors using a sparse group lasso model. In contrast, \code{LF_VCR.ae} begins with an autoencoder-based feature extraction step before continuing with the modeling process, providing a useful alternative when deep feature representation is beneficial.
}
\examples{
# Example usage for LF_VCR:
result1 <- LF_VCR(X = data_matrix,
                  number.K = 3,
                  covariate = cov,
                  matrix = "covariance",
                  categorical = FALSE,
                  y = outcome)

# Example usage for LF_VCR.ae:
result2 <- LF_VCR.ae(X = data_matrix,
                     number.K = 3,
                     covariate = cov,
                     matrix = "covariance",
                     categorical = FALSE,
                     y = outcome)
}
