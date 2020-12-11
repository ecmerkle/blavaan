## predictions from blavaan object; similar to lavPredict, but lavPredict is never called
## overload standard R function `predict'
setMethod("predict", "blavaan",
function(object, newdata = NULL) {
    blavPredict(object = object, newdata = newdata, type = "lv")
})

blavPredict <- function(blavobject, newdata = NULL, type = "lv") {

  stopifnot(inherits(blavobject, "blavaan"))
  blavmodel <- blavobject@Model
  blavpartable <- blavobject@ParTable
  blavdata <- blavobject@Data
  standata <- blavobject@external$mcmcdata
  
  type <- tolower(type)
  if(type %in% c("latent", "lv", "factor", "factor.score", "factorscore"))
      type <- "lv"
  if(type %in% c("ov","yhat"))
      type <- "yhat"
  if(type %in% c("ypred", "ydist"))
      type <- "ypred"
  if(type %in% c("ymis", "ovmis"))
      type <- "ymis"
  
  
  stantarget <- lavInspect(blavobject, "options")$target == "stan"

  if(!is.null(newdata)) stop("blavaan ERROR: posterior predictions for newdata are not currently supported")
  
  ## lv: posterior dist of lvs (use blavInspect functionality); mcmc list
  ## lvmeans: use blavInspect functionality; matrix
  ## yhat: posterior expected value of ovs conditioned on lv samples; mcmc list
  ## ypred: posterior predictive distribution of ovs conditioned on lv samples; mcmc list
  ## ymis: posterior predictive distribution of missing values conditioned on observed values; matrix
  if(type == "lv") {
    out <- blavInspect(blavobject, 'lvs')
  } else if(type == "lvmeans") {
    out <- blavInspect(blavobject, 'lvmeans')
  } else if(type %in% c("yhat", "ypred", "ymis")) {
    if(!stantarget) stop(paste0("blavaan ERROR: '", type, "' is only supported for target='stan'"))

    if(type %in% c("yhat", "ypred")) {
      ## TODO; ypred is yhat plus noise
      out <- NULL
    }

    if(type == "ymis") {
      out <- samp_data(blavobject@external$mcmcout, blavmodel, blavpartable, standata, blavdata)
    }
  } else {
    stop("blavaan ERROR: unknown type supplied; use lv lvmeans yhat ypred ymis")
  }
  
  out
}
