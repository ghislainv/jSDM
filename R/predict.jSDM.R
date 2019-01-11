# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ecology.ghislainv.fr
# license         :GPLv3
# ==============================================================================

predict.jSDM <- function(object, newdata=NULL, type="mean", probs=c(0.025,0.975), ...) {

    ##= Check
    if (!(type %in% c("mean","quantile","posterior"))) {stop("type must be \"mean\", \"quantile\" or \"posterior\"")}
    if (sum(probs<0)>0 | sum(probs>1)>0) {stop("probs must be a vector of probabilities between (0,1)")}
    
    ##= Link function depending on family
    if (object$family=="binomial") {
        inv.link <- function (x) {inv_logit(x)}
    }
    if (object$family=="poisson") {
        inv.link <- function (x) {exp(x)}
    }
    
    ##= Model specifications
    model.spec <- object$model.spec
    
    ##= Matrix for predictions
    if (is.null(newdata)) {
        newdata <- model.spec$data
    }
    suitability <- model.spec$suitability
    mf.pred <- model.frame(formula=suitability,data=newdata)
    X.pred <- model.matrix(attr(mf.pred,"terms"),data=mf.pred)
    npred <- nrow(X.pred)      
    ##= Model parameters
    np <- ncol(X.pred)
    ngibbs <- model.spec$mcmc+model.spec$burnin
    nthin <- model.spec$thin
    nburn <- model.spec$burnin
    nsamp <- model.spec$mcmc/nthin
    ##= Matrix of MCMC for beta parameters
    beta.mat <- as.matrix(object$mcmc)[,c(1:np)]
    ##= Posterior mean
    if (type=="mean") {
        term <- rep(0, npred)
        ##= Loop on samples
        for (i in 1:nsamp) {
            beta <- beta.mat[i,]
            link.term <- X.pred %*% beta
            term <- term + inv.link(link.term)
        }
        term.pred <- as.numeric(term/nsamp)
    }
    ##= Full posterior
    if (type %in% c("quantile","posterior")) {
        term <- matrix(NA,nrow=nsamp,ncol=npred)
        ##= Loop on samples
        for (i in 1:nsamp) {
            beta <- beta.mat[i,]
            link.term <- X.pred %*% beta
            term[i,] <- inv.link(link.term)
        }
        if (type=="quantile") {
            term.mean <- apply(term,2,mean)
            term.quant <- apply(term,2,quantile,probs)
            term.pred <- as.data.frame(t(rbind(term.mean,term.quant)))
            names(term.pred)[1] <- c("mean")
        }
        if (type=="posterior") {
            term.pred <- coda::mcmc(term,start=nburn+1,end=ngibbs,thin=nthin)
        }
    }

    ## Output
    return(term.pred)
}

# End