#' Generate Counterfactual Couples
#'
#' @description  This function will sample a set of alternate partners from a data.frame for a set of actual
#'               unions in order to generate a set of counterfactual unions that can be used in a conditional
#'               logit model to predict the characteristics that are likely to lead to a union match.
#' @param n the number of counterfactual unions to create for each real union.
#' @param actual a data.frame identifying identifying the actual couples. This dataset should have the variable identified
#'              by \code{geo} and variables for husbands and wives ids that end in "h" and "w" respectively. It should also
#'              contain all the same variables ending in "h" and "w" that should be kept in the final analysis.
#' @param men a data.frame object identifying all the potential "male" partners. This data.frame must include the
#'           variable identified by \code{geo} above and it must also contain an id variable whose name starts with
#'           the string given by \code{id} and end with "h". Any person-specific variable that should be kept in the
#'           results should end with an "h" (e.g. ageh, raceh). These variables should correspond exactly to variables
#'           in \code{actual} and should correspond to variables in \code{women} that end in "w".
#' @param women a data.frame object identifying all the potential "female" partners. Its format should be the
#'             same as \code{men}, except variables should end in "w" rather than "h".
#' @param geo a character string giving the name of the variable that identifies the clusters that alternate
#'            partners should be sampled from. This variable must be named the same way in \code{men},
#'            \code{women}, and \code{actual}.
#' @param id a character string giving part of the variable name to identify partners in \code{actual}, \code{men}, and
#'           \code{women}. In each of the datasets, the actual variable name should be appended with "h" and "w" for
#'           husbands and wives respectively.
#' @param weight A character string that identifies a weight variable for sampling of potential partners. If left
#'              \code{NULL}, partners will be sampled with equal probability. This weight variable must exist in both the
#'              \code{men} and \code{women} data.
#' @param keep a vector of character strings identifying additional variables in \code{actual} that should be kept in the
#'            results.
#' @param verbose if true, the function will report its progress.
#' @details This program uses the \pkg{wrswoR} package to do fast weighted sampling without replacement. Large sets of
#'          alternate spouses and actual unions can be sampled relatively quickly, but processing time will begin to
#'          increase exponentially with very large datasets within each cluster.
#'
#' For each actual union, the program randomly determines one of the two spouses to sample with 50/50 odds. This is to
#' ensure that all characteristics of a single spouse are fixed within the fixed effects conditional logit model.
#'
#' Currently the format of the datasets must be followed exactly for the function to work correctly.
#' @return The output of this program is a data.frame of actual and counterfactual unions. It will keep all
#' variables in the three datasets that end in an "h" or "w" as well as:
#' \item{geo}{description the cluster identifier used in the function.}
#' \item{group}{a unique identifier based on the id of the spouse from the actual union who had partners
#'              sampled for them. This should be used as the fixed effect in fixed effects models.}
#' \item{choice}{a boolean variables that is \code{TRUE} if this is an actual union and \code{FALSE} if this is a
#' counterfactual union. This variable should be used as the dependent variable in a fixed effect conditional logit model.}
#' @examples
#' #generate three counterfactual couples for each real couples
#' #in example ACS data
#' market <- generateCouples(3,acs.couples,
#'                           acs.malealters,acs.femalealters,
#'                           "state",weight="perwt",keep="hhwt")
#'
#' #check that there is one real marriage and three counterfactual
#' #marriages for each case
#' summary(tapply(market$choice,market$group,sum))
#' summary(tapply(!market$choice,market$group,sum))
#'
#' \dontrun{
#' #load survival function and run clogit command to estimate how age
#' #differences and racial exogamy affect the log-odds of union formation
#' require(survival)
#' model <- clogit(choice~I(ageh-agew)+I((ageh-agew)^2)+I(raceh!=racew)
#'                        +strata(group), data=market)
#' summary(model)
#' }
generateCouples <- function(n, actual, men, women, geo, id="id",
                            weight=NULL, keep=NULL, verbose=TRUE) {

  if(!require(wrswoR)) {
    stop("package wrswoR must be installed.")
  }

  if(verbose) {
    cat("\tSet up....")
  }

  #remove cases that are missing values on clustering variable
  actual <- actual[!is.na(actual[,geo]),]
  men <- men[!is.na(men[,geo]),]
  women <- women[!is.na(women[,geo]),]

  #remove men and women that come from clusters that are not present in actual
  markets <- unique(actual[,geo])
  men <- men[men[,geo] %in% markets,]
  women <- women[women[,geo] %in% markets,]

  #check if sample size exceeds the number of alternates in market area
  small.areas <- names(which(table(men[,geo])<n | table(women[,geo])<n))
  if(any(small.areas %in% markets)) {
    warning(paste("The following areas had less than the required sample for men and women and the entire list was drawn:",
                  paste(small.areas[small.areas %in% markets],collapse=", "),collapse=" "))
  }

  #decide whether to sample husband or wife
  actual$sampled <- sample(c("h","w"),nrow(actual),replace=TRUE)

  if(verbose) {
    cat("Done\n")
  }

  if(verbose) {
    cat("\tSampling from markets...")
  }

  #split the datasets up by the clustering variable
  markets.actual <- split(actual, actual[,geo], drop=TRUE)
  markets.men <- split(men, men[,geo], drop=TRUE)
  markets.women <- split(women, women[,geo], drop=TRUE)

  if(length(markets.actual)!=length(markets.men)) {
    stop("The number of clusters for alternate male partners does not equal the number of clusters for actual unions")
  }
  if(length(markets.actual)!=length(markets.women)) {
    stop("The number of clusters for alternate female partners does not equal the number of clusters for actual unions")
  }
  if(any(names(markets.actual)!=names(markets.men))) {
    stop("The clusters for actual unions do not match the clusters for alternate male partners")
  }
  if(any(names(markets.actual)!=names(markets.women))) {
    stop("The clusters for actual unions do not match the clusters for alternate female partners")
  }

  #sample the counterfactual unions. The use of the mapply and do.call dramatically speeds up processing
  #time compared to for-looping through each cluster
  fakes <- do.call("rbind",
                   mapply(function(market.actual, market.men, market.women) {
                                rbind(organizeColumns(samplePartners(subset(market.actual,sampled=="h"),
                                                                     market.women,n,"w", weight, id),
                                                      geo, keep),
                                      organizeColumns(samplePartners(subset(market.actual,sampled=="w"),
                                                                     market.men,n,"h", weight, id),
                                                      geo, keep))
                          },
                          markets.actual, markets.men, markets.women,SIMPLIFY=FALSE))
  row.names(fakes) <- NULL

  if(verbose) {
    cat("Done.\n")
  }

  if(verbose) {
    cat("\tCombining data....")
  }

  actual$choice <- TRUE
  actual$group <- ifelse(actual$sampled=="h", actual$idh, actual$idw)
  partners <- rbind(organizeColumns(actual,geo,keep), fakes)
  partners$group <- as.factor(partners$group)
  partners <- partners[order(partners$group),]

  if(verbose) {
    cat("Done\n")
  }

  return(partners)
}


#' Sample partners
#'
#'@description This function will draw a set of counterfactual partners for a given set of real unions.
#' It is primarily intended for internal use by the \code{generateCouples} functiom.
#' @param actual a data.frame of actual unions identical in structure to that for \code{generateCouples}.
#' @param eligibles a data.frame of eligible partners to be sampled. This should
#'             have the same structure as that for \code{generateCouples}.
#' @param n The number of alternate partners to sample for each real spouse.
#' @param partner a character sting of  either "h" or "w" indicating which type of partner should be sampled.
#' @param  weight optional name of weights variable in eligibles that
#'          allows for different probabilities of sampling eligible partners. If left out
#'          then sampling will be done with equal probability
#' @param id a character string indicating the first part of the variable name for the id to identify indidividuals
#'        in actual and eligibles. Identical to the argument in generateCouples.
#' @details This function uses the \pkg{wrswoR} package for fast weighted sampling without replacement. In addition,
#'      for speed, the initial matching is done without concern for potentially sampling the correct partner from
#'      the set of alternate partners. This is then checked and removed in a secondary stage for the rare cases where
#'      such duplication occurs.
#' @return a data.frame object similar to \code{actual} but containing only counterfactual unions.
samplePartners <- function(actual, eligibles, n, partner, weight=NULL, id="id") {

  if(nrow(actual)==0) {
    return(NULL)
  }

  #set up ego and partner string ids
  ego <- "h"
  if(partner=="h") {
    ego <- "w"
  }
  if(partner!="h" & partner!="w") {
    warning("partner argument must be h or w")
    break
  }

  id.ego <- paste(id,ego,sep="")
  id.partner  <- paste(id,partner,sep="")

  #do the sampling
  if(is.null((weight))) {
    eligibles$weight <- rep(1,nrow(eligibles))
  } else {
    eligibles$weight <- eligibles[,weight]
  }

  #if n is greater than the number of eligibles, then reduce n, but produce warning
  if(n>nrow(eligibles)) {
    n <- nrow(eligibles)
  }

  idx.fakes <- as.vector(replicate(nrow(actual),
                                   sample_int_expj(nrow(eligibles),n,eligibles$weight)))

  #combine the actual spouse with the fakes
  fakes <- cbind(actual[rep(1:nrow(actual),each=n),!grepl(paste(partner,"$",sep=""),colnames(actual))],
                 eligibles[idx.fakes,grep(paste(partner,"$",sep=""),colnames(eligibles))],
                 row.names=NULL)

  #address any cases where we might have accidently drawn the real match
  dupe.idx <- which(actual[rep(1:nrow(actual),each=n),id.partner]==fakes[,id.partner])

  toRemove <- NULL

  for(j in dupe.idx) {
    #resample, while explicitly removing all already assigned spouses from pool (including real one)
    alreadyassigned <- fakes[which(fakes[,id.ego]==fakes[j,id.ego]),id.partner]
    validEligibles <- eligibles[which(!(eligibles[,id.partner] %in% alreadyassigned)),]
    if(nrow(validEligibles)==0) {
      toRemove <- c(toRemove,j)
    } else {
      #newPartner <- eligibles[sample(which(!(eligibles[,id.partner] %in% alreadyassigned)),1),]
      newPartner <- validEligibles[sample_int_expj(nrow(validEligibles),1,validEligibles$weight),]
      fakes[j,grep(paste(partner,"$",sep=""),colnames(fakes))] <- newPartner[grep(paste(partner,"$",sep=""),
                                                                                  colnames(eligibles))]
    }
  }

  #remove any duplicates where we couldn't resample
  if(length(toRemove)>0) {
    fake <- fakes[-toRemove]
  }

  fakes$group <- fakes[,id.ego]
  fakes$choice <- FALSE

  return(fakes)
}


#' Organize columns for final union data
#' @description This is an internal function used by \code{generateCouples} to organize the final
#'            dataset returned from the function.
#' @param couples A data frame of the actual and counterfactual couples.
#' @param geo A character string giving the name of the clustering variable.
#' @param keep A vector of character strings giving additional variables to keep that don't end in
#'             "h" or "w".
organizeColumns <- function(couples, geo, keep=NULL) {
  return(couples[c(which(colnames(couples) %in% c(geo,"group","choice","sampled",keep)),
                   grep("h$",colnames(couples)),
                   grep("w$",colnames(couples)))])
}




#' Run fixed effects conditional logit model on multiple datasets and pool results
#'
#' @description Given a list of multiple datasets of real and counterfactual unionts and a model formula, this function
#'              will calculate a discrete choice model for each dataset and pool the results, taking account of both
#'              within and between variance in calculating the standard error.
#' @param formula an object of class \code{\link[stats]{formula}} specifying the \code{\link[survival]{clogit}} model to be
#'                performed on each dataset.
#' @param datasets a list of datasets where each dataset is produced from the \code{\link{generateCouples}} function.
#' @param method A character string indicating the estimation method to use in the \code{\link[survival]{clogit}} model.
#' @param parallel A boolean indicating whether to use the parallel package to increase the speed of estimation via multiple cores.
#' @details Because the dataset of real and counterfactual unions is created by sampling among all the possible alternate partners,
#'          the results of models will vary as a result of this sampling process. Therefore, it may be useful to generate multiple
#'          datasets and pool model results across these datasets in a manner identical to multiple imputation, where the standard
#'          errors of estimates are adjusted for the variance in coefficient estimates across datasets.
#'
#'          This function is a convenience function that will perform this pooling and produce properly adjusted results.
#'          The reported coefficients from the model are given by taking the mean across all datasets. The reported variance \eqn{V}
#'          for each parameter is given by:
#'          \deqn{V=W+(1+1/m)B}
#'          Where \eqn{m} is the number of datasets, \eqn{W} is the within variance, estimated by the square of the mean standard
#'          error across datasets, and \eqn{B} is the between variance estimated by the variance of coefficient estimates across
#'          datasets.
#'
#'          Models are estimated using the \code{\link[survival]{clogit}} function from the \pkg{survival} package. This package
#'          must be installed.
#' @return a list containing the following objects:
#' \item{coefficients}{a data.frame object with the following elements:}
#' \itemize{
#'     \item{b.pool: }{The average coefficient across datasets.}
#'     \item{se.pool: }{standard error that combined within and between variance}
#'     \item{z.pool: }{z-statistic from dividing b by se}
#'     \item{pvalue.pool: }{p-value for the hypothesis test that the coefficient is zero in the population}
#'     \item{within.var: }{The square of the mean standard error across datasets}
#'     \item{between.var: }{the variance of the coefficient across datasets}
#'     }
#' \item{deviance}{A vector of deviances for each model.}
#' \item{bic}{A vector of BIC statistic for each dataset relative to the null model.}
#' @examples
#' markets <- replicate(5, generateCouples(3,acs.couples,
#'                            acs.malealters,acs.femalealters,
#'                            "state",weight="perwt",verbose=FALSE),
#'                      simplify=FALSE)
#'
#' poolChoiceModel(choice~ageh+I(ageh^2)+I(ageh-agew)+I((ageh-agew)^2)+strata(group),
#'                 markets)
poolChoiceModel <- function(formula, datasets, method="exact", parallel=FALSE) {
  if(!require(survival)) {
    stop("survival package must be installed.")
  }
  if(parallel & require(parallel)) {
    models <- mclapply(datasets, function(dataset) {
      clogit(formula, data=dataset, method=method)
    })
  } else {
    models <- lapply(datasets, function(dataset) {
      clogit(formula, data=dataset, method=method)
    })
  }

  m <- length(models)
  b <- sapply(models, coef)
  se <- sapply(models, function(model) {summary(model)$coef[,3]})
  deviance <- sapply(models, function(model) {-2*model$loglik[2]})
  bic <- sapply(models, function(model) {diff(-2*model$loglik)+length(model$coef)*log(sum(model$y[,2]))})

  b.pool <- apply(b,1,mean)
  between.var <- apply(b,1,var)
  within.var <- apply(se^2,1,mean)
  se.pool <- sqrt(within.var+between.var+between.var/m)
  z.pool <- b.pool/se.pool
  pvalue.pool <- (1-pnorm(abs(z.pool)))*2

  return(list(coefficients=data.frame(b.pool,se.pool,z.pool,pvalue.pool,within.var,between.var),
              deviance=deviance, bic=bic))
}


#this will take a husband and wife categorical variable and generate a set of
#exogamy categories, where endogamy is the reference. It can also make these
#terms symmetric

#' Create exogamy terms from partner's characteristics
#'
#' @description This function will calculate a variable that combines two identical factor
#' variables for each partner, replacing all cases where partners match to a reference category
#' of endogamous unions.
#' @param varh A factor vector giving the categories for the "husband."
#' @param varw A factor vector giving the categories for the "wife."
#' @param symmetric If set to true, then the ordering of the categories by partner will be ignored.
#' @details This function operates much like the \code{\link[base]{interaction}} function to combine partner
#'  characteristics. The function assumes that the two variables being combined are identical categorical variables
#'  for each partner (e.g. husband's and wife's education, husband's and wife's race). It differs from a basic interaction
#'  in two ways. First, any cases where husband and wife have the same category will be recorded as the "Endog" category and
#'  will be set to the reference. Second, users can specify that terms should be symmetric such that couples are put in the
#'  same exogamous category without respect to whether it was the husband or wife who belonged to the specific categories.
#'  For example, if symmetric, then the category for white husbands and black wives would be the same as for black husbands
#'  and white wives.
#'
#'
#' @return a factor vector containing the combined categories, with "Endog" set as the reference.
#' @examples
#'  market <- generateCouples(3,acs.couples,
#'                           acs.malealters,acs.femalealters,
#'                           "state",weight="perwt")
#'  market$rexog <- createExogamyTerms(market$raceh, market$racew)
#'  market$rexog.sym <- createExogamyTerms(market$raceh, market$racew, symmetric=TRUE)
#'  levels(market$rexog)
#'  levels(market$rexog.sym)
#'
#'  \dontrun{
#'  require(survival)
#'  clogit(choice~rexog+strata(group), data=market)
#'  clogit(choice~rexog.sym+strata(group), data=market)
#'  }
createExogamyTerms <- function(varh, varw, symmetric=FALSE) {
  ncat <- nlevels(varh)
  if(nlevels(varw) != ncat) {
    stop("partner variables have different number of categories")
  }
  if(any(!(levels(varw) %in% levels(varh)))) {
    stop("partner variables have different sets of categories")
  }
  exog <- interaction(varh, varw, sep=".", lex.order=TRUE)
  #replace all couples with same characterstics with a single endogamous ("Endog") level
  l <- levels(exog)
  l[seq(from=1, by=(ncat+1),to=length(l))] <- "Endog"
  if(symmetric) {
    #ignore the ordering of husbands and wives
    for(i in 1:length(l)) {
      if(l[i]=="Endog") {
        next
      }
      x <- strsplit(l[i],"\\.")
      reverse <- paste(x[[1]][2],x[[1]][1],sep=".")
      if(which(l==reverse)<i) {
        l[i] <- reverse
      }
    }
  }
  levels(exog) <- l
  exog <- relevel(exog, ref="Endog")
  return(exog)
}
