\name{sampleData}
\alias{sampleData}
\alias{sampledata}
\title{Sample data from the posterior (or prior) distribution.}
\description{
  The purpose of the \code{sampleData()} function is to simulate new data
  from a model that has already been estimated. This can faciliate
  posterior predictive checks, as well as prior predictive checks (setting
  prisamp = TRUE during model estimation).}
\usage{
sampleData(object, nrep = NULL, conditional = FALSE, type = "response",
           simplify = FALSE, ...)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{blavaan}}.}
\item{nrep}{How many datasets to generate? If not supplied, defaults to
  the total number of posterior samples.}
\item{conditional}{Logical indicating whether to sample from the
  distribution that is marginal over latent variables (\code{FALSE};
  default) or from the distribution that conditions on latent variables
  (\code{TRUE}). For \code{TRUE}, you must set \code{save.lvs = TRUE}
  during model estimation.}
\item{type}{The type of data desired (only relevant to ordinal
  data). The \code{type = "response"} option generates ordinal data. The
\code{type = "link"} option generates continuous variables underlying
ordinal data (which would be cut by thresholds to yield ordinal data).}
\item{simplify}{For single-group models, should the list structure be
  simplified? This makes each dataset a single list entry, instead of a
  list within a list (which reflects group 1 of dataset 1). Defaults to \code{FALSE}.}
\item{...}{Other arguments, which for now is only \code{parallel}. Parallelization via \code{future_lapply()} is available by setting \code{parallel = TRUE}.}
}
\details{
This is a convenience function to generate data for posterior or prior
predictive checking. The underlying code is also used to generate data
for posterior predictive p-value computation.
}
\seealso{
This function overlaps with \code{blavPredict()}. The
\code{blavPredict()} function is more focused on generating pieces of
data conditioned on other pieces of observed data (i.e., latent
variables conditioned on observed variables; missing variables
conditioned on observed variables). In contrast, the \code{sampleData()}
function is more focused on generating new data given the sampled model parameters.
}
\examples{
\dontrun{
data(HolzingerSwineford1939, package = "lavaan")

## fit model
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- bcfa(HS.model, data = HolzingerSwineford1939)

## 1 dataset generated from the posterior
out <- sampleData(fit, nrep = 1)

## nested lists: 1 list entry per nrep.
## then, within a rep, 1 list entry per group
## so our dataset is here:
dim(out[[1]][[1]])

## 1 posterior dataset per posterior sample:
out <- sampleData(fit)

## obtain the data on x1 across reps and summarize:
x1dat <- sapply(out, function(x) x[[1]][,1])
summary( as.numeric(x1dat) )
}
}
