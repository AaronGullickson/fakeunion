#' Married Couples from ACS 2011-2015
#'
#' A dataset containing information on a sample of married couples from the American Community Survey 2011-2015
#' data for the states of New York and California. Data were downloaded from IPUMS.
#'
#' @format A data frame containing data on 1000 married couples from the 2011-2015 ACS data drawn randomly from
#'        the states of New York and California. The data frame contains the following variables for each couple:
#'        \describe{
#'          \item{state}{State of residence. This can be used to select appropriate cluster to sample alternate spouses}
#'          \item{hhwt}{Household weight for the given couple}
#'          \item{idw}{unique id of wife}
#'          \item{agew}{age of wife}
#'          \item{racew}{race of wife. Either "W" (white) or "B" (black)}
#'          \item{idh}{unique id of husband}
#'          \item{ageh}{age of husband}
#'          \item{raceh}{race of husband. Either "W" (white) or "B" (black)}
#'        }
#' @source \url{https://usa.ipums.org/usa/}
"acs.couples"

#' Alternate male partners from ACS 2011-2015
#'
#' A dataset containing information on a sample of unmarried and recently married men from the American Community
#' urvey 2011-2015 data for the states of New York and California. Data were downloaded from IPUMS.
#'
#' @format A data frame containing data on 1000 unmarried and recently married men from the 2011-2015 ACS data drawn
#'         randomly from the states of New York and California. The data frame contains the following variables:
#'        \describe{
#'          \item{state}{State of residence. This can be used to select appropriate cluster to sample alternate spouses}
#'          \item{perwt}{Person weight for the given respondent}
#'          \item{idh}{unique id of husband}
#'          \item{ageh}{age of husband}
#'          \item{raceh}{race of husband. Either "W" (white) or "B" (black)}
#'        }
#' @source \url{https://usa.ipums.org/usa/}
"acs.malealters"

#' Alternate female partners from ACS 2011-2015
#'
#' A dataset containing information on a sample of unmarried and recently married women from the American Community
#' urvey 2011-2015 data for the states of New York and California. Data were downloaded from IPUMS.
#'
#' @format A data frame containing data on 1000 unmarried and recently married women from the 2011-2015 ACS data drawn
#'         randomly from the states of New York and California. The data frame contains the following variables:
#'        \describe{
#'          \item{state}{State of residence. This can be used to select appropriate cluster to sample alternate spouses}
#'          \item{perwt}{Person weight for the given respondent}
#'          \item{idw}{unique id of wife}
#'          \item{agew}{age of wife}
#'          \item{racew}{race of wife. Either "W" (white) or "B" (black)}
#'        }
#' @source \url{https://usa.ipums.org/usa/}
"acs.femalealters"

