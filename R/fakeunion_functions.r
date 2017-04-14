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
#' markets <- generateCouples(3,acs.couples,
#'                            acs.malealters,acs.femalealters,
#'                            "state",weight="perwt",keep="hhwt")
#'
#' #check that there is one real marriage and three counterfactual
#' #marriages for each case
#' summary(tapply(markets$choice,markets$group,sum))
#' summary(tapply(!markets$choice,markets$group,sum))
#'
#' \dontrun{
#' #load survival function and run clogit command to estimate how age
#' #differences and racial exogamy affect the log-odds of union formation
#' require(survival)
#' model <- clogit(choice~I(ageh-agew)+I((ageh-agew)^2)+I(raceh!=racew)
#'                        +strata(group), data=markets)
#' summary(model)
#' }
generateCouples <- function(n, actual, men, women, geo, id="id",
                            weight=NULL, keep=NULL, verbose=TRUE) {

  if(!require(wrswoR)) {
    cat("package wrswoR not found")
    return()
  }

  if(verbose) {
    cat("\tSet up....")
  }

  #get all unique markets
  markets <- unique(na.omit(actual[,geo]))

  #initialize dataset of fake marriages
  fakes <- NULL

  #set up groups for actual marriages
  actual$group <- NA

  if(verbose) {
    cat("Done\n")
  }

  if(verbose) {
    cat("\tSampling from markets:")
  }

  #loop through markets and create fakes
  for(i in markets) {

    if(verbose) {
      cat("\n\t\t")
      cat(i)
    }

    #get the index for actual marriages
    idx <- which(actual[,geo]==i)

    #split the actual unions so that roughly half get wives randomly selected
    #and the other half get husbands randomly selected
    tmp <- sample(1:2,length(idx),replace=TRUE)
    idx.h <- idx[tmp==1]
    idx.w <- idx[tmp==2]
    rm(tmp)

    if(length(idx.h)>0) {
      actual$group[idx.h] <- actual[idx.h,paste(id,"h",sep="")]
      fakes.men <- samplePartners(actual[idx.h,],women[which(women[,geo]==i),],n,"w", weight, id)
      fakes <- rbind(fakes, organizeColumns(fakes.men,geo,keep))
    }

    if(length(idx.w)>0) {
      actual$group[idx.w] <- actual[idx.w,paste(id,"w",sep="")]
      fakes.women <- samplePartners(actual[idx.w,],men[which(men[,geo]==i),],n,"h", weight, id)
      fakes <- rbind(fakes, organizeColumns(fakes.women,geo,keep))
    }
  }

  if(verbose) {
    cat("\n\tSampling complete.")
  }

  if(verbose) {
    cat("\n\tCombining data....")
  }

  real <- actual[!is.na(actual[,geo]),]
  real$choice <- TRUE

  #combine fakes with reals
  partners <- rbind(organizeColumns(real,geo,keep), fakes)
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
    weight <- rep(1,nrow(eligibles))
  } else {
    weight <- eligibles[,weight]
  }

  idx.fakes <- as.vector(replicate(nrow(actual),
                                   sample_int_expj(nrow(eligibles),n,weight)))

  #combine the actual spouse with the fakes
  fakes <- cbind(actual[rep(1:nrow(actual),each=n),!grepl(paste(partner,"$",sep=""),colnames(actual))],
                 eligibles[idx.fakes,grep(paste(partner,"$",sep=""),colnames(eligibles))],
                 row.names=NULL)

  #address any cases where we might have accidently drawn the real match
  dupe.idx <- which(actual[rep(1:nrow(actual),each=n),id.partner]==fakes[,id.partner])

  for(j in dupe.idx) {
    #resample, while explicitly removing all already assigned spouses from pool (including real one)
    alreadyassigned <- fakes[which(fakes[,id.ego]==fakes[j,id.ego]),id.partner]
    newPartner <- eligibles[sample(which(!(eligibles[,id.partner] %in% alreadyassigned)),1),]
    fakes[j,grep(paste(partner,"$",sep=""),colnames(fakes))] <- newPartner[grep(paste(partner,"$",sep=""),colnames(eligibles))]
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
  return(couples[c(which(colnames(couples) %in% c(geo,"group","choice",keep)),
                   grep("h$",colnames(couples)),
                   grep("w$",colnames(couples)))])
}
