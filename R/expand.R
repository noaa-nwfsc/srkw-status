#' Function expand expands a data frame from animals as observations to animal-years as obserations
#'
#' @param filename A file location (csv file) with births and deaths of each animal
#' @param start_year Starting year for the expanded data, starts at 1976
#' @param current_year End year for the expanded data, defaults to current calendar year
#' @param fem_age_mat Female age at maturity, 10 years of age by default
#' @param fem_age_senesc Female age at reproductive senescence, defaults to 43
#'
#' @import dplyr
#' @importFrom utils read.csv
#'
#' @export
#' @examples
#' \dontrun{
#' library(kwdemog)
#' data(orca)
#' expanded_data = expand(orca)
#' # or use a filename
#' expanded_data = expand("use_this_file.csv", start_year = 1979, current_year = 2018,
#' fem_age_mat = 10, fem_age_senesc = 43)
#' }
expand <- function(filename, start_year = 1976, current_year = NULL,
  fem_age_mat = 10, fem_age_senesc = 43) {
  if(is.null(current_year)) current_year = as.numeric(substr(Sys.Date(),1,4))
  directory = find.package("kwdemog")
  ages2stages = read.csv(paste0(directory, "/extdata/ages2stages.csv"), stringsAsFactors = FALSE)
  if(class(filename) == "data.frame") orca = filename
  if(class(filename) == "character") orca = read.csv(filename, stringsAsFactors = FALSE)
  orca$population = ifelse(orca$pod %in% c("J001","K001","L001"), "SRKW", "NRKW")
  # create expanded data frame of all years and animals
  alldat = expand.grid("year" = start_year:current_year,
    "animal" = unique(orca$animal), stringsAsFactors = FALSE)

  alldat = left_join(alldat, orca)
  alldat$age = alldat$year - alldat$birth
  alldat$age[which(alldat$age < 0)] = NA
  alldat$age[which(alldat$age > alldat$death)] = NA

  # add alive column for live / died. NAs before an animal was born, or after they died
  alldat$alive = NA
  alldat$alive[which(alldat$year == alldat$death)] = 0
  alldat$alive[which(alldat$year >= alldat$birth & alldat$year < alldat$death)] = 1
  alldat$alive[which(alldat$year >= alldat$birth & is.na(alldat$death))] = 1

  # add birth information. This is somewhat complicated because females
  # can't give birth in the years before / after they give birth (18 mo gestation)
  alldat$gave_birth = NA
  fems = unique(alldat$animal[which(alldat$alive == 1 & alldat$age %in% seq(fem_age_mat,fem_age_senesc) & alldat$sexF1M2==1)])
  for(i in 1:length(fems)) {
    alldat$gave_birth[which(alldat$animal == fems[i] & alldat$age %in% seq(fem_age_mat,fem_age_senesc) & alldat$sexF1M2==1)] = 0
    birth_years = alldat$year[which(alldat$mom == fems[i] & alldat$age == 0)]
    alldat$gave_birth[which(alldat$animal==fems[i] & alldat$year %in%birth_years)] = 1
    alldat$gave_birth[which(alldat$animal==fems[i] & alldat$year %in%(birth_years-1))] = NA
    alldat$gave_birth[which(alldat$animal==fems[i] & alldat$year %in%(birth_years+1))] = NA
  }

  # Fill in NAs for births in 1979 because of births in 1978
  alldat$gave_birth[alldat$animal%in%c("J010","K003","L011","L032","A11","D07","G12") &
      alldat$year=="1979" & alldat$gave_birth != 1] = NA

  # Delete some NRKW individuals that were not seen in certain periods
  # I01 matriline not seen 1994-1996
  alldat$gave_birth[alldat$pod=="I01" & alldat$year%in%c(seq(1994,1996)) & alldat$gave_birth != 1] = NA
  alldat$gave_birth[alldat$pod=="I18" & alldat$year%in%c(seq(1994,1996)) & alldat$gave_birth != 1] = NA
  #X$Birth[X$matriline=="H01" & X$time1%in%c(seq(1994,1999)) & X$Birth != 1] = NA
  #X$Birth[X$matriline=="H03" & X$time1%in%c(seq(1994,1999)) & X$Birth != 1] = NA
  alldat$gave_birth[alldat$pod=="H01" & alldat$year%in%c(seq(1994,1999)) & alldat$gave_birth != 1] = NA

  alldat$alive[alldat$animal=="I70" & alldat$year > 1993] = NA
  alldat$alive[alldat$animal=="I59" & alldat$year > 1993] = NA
  alldat$alive[alldat$animal=="I60" & alldat$year > 1989] = NA
  alldat$alive[alldat$animal=="I66" & alldat$year > 1991] = NA
  alldat$alive[alldat$animal=="I10" & alldat$year > 1975] = NA
  alldat$alive[alldat$animal=="I84" & alldat$year > 2002] = NA
  alldat$alive[alldat$animal=="R10" & alldat$year > 1975] = NA
  alldat$alive[alldat$animal=="R11" & alldat$year > 1975] = NA
  alldat$alive[alldat$animal=="I17" & alldat$year%in%c(seq(1992,1996))] = NA
  alldat$alive[alldat$animal=="I26" & alldat$year%in%c(seq(1992,1996))] = NA
  alldat$alive[alldat$animal=="R07" & alldat$year%in%c(seq(1976,1980))] = NA
  alldat$alive[alldat$animal=="R05" & alldat$year%in%c(seq(1976,1981))] = NA

  return(alldat)
}
