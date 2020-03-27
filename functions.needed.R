# 16.03.2017 :  - changed use of do.trunc in rtriang.perc and triangulize. Now it should either
#                 be NA for no truncation, or it should have two values. The minimum of these two
#                 values will be used as lower truncation mark, and the maximum as upper truncation
#                 mark.
#               - introduced possibility of having negative TC
#               - introduced argument "unit" in triangulize, to enable skipping distribution-making
#                 of TCs = 1, but not preventing having inputs = 1
# 
# History of modifications:   2017.05.02: corrected mistake in triangulize.TC
#                             2017.05.03: changed message from triangulize.TC
#                                         added rtrapez.perc function

# =================================================================================================

# introduce symmetric triangular distribution with percentage (default being 50%)
rtriang.perc <- function(mmode, perc = 0.5, N = SIM, do.trunc = NA)
{
  require(mc2d)
  
  stopifnot(perc <= 1, perc >= 0)
  
  if(is.na(mmode)){
    return(rep(NA,N))
  }
  
  mmin <- mmode - abs(mmode*perc)
  mmax <- mmode + abs(mmode*perc)
  
  if(!is.na(do.trunc[1])){
    # if wished, truncate the distribution under 0 and above 1
    return(rtrunc(distr = "rtriang", n = N,
                  linf = min(do.trunc), lsup = max(do.trunc),
                  min = mmin, mode = mmode, max = mmax))
  } else {
    return(rtriang(n = N, min = mmin, mode = mmode, max = mmax))
  }
  
  return(out)
}

# introduce trapezoidal distribution with percentage (default being 50%)
rtrapez.perc <- function(data, perc = 0.5, N = SIM, do.trunc = NA)
{
  require(trapezoid)
  
  stopifnot(length(data) == 2, perc <= 1, perc > 0)
  if(any(is.na(data))){
    return(rep(NA,N))
  }
  
  data <- sort(data)
  
  mmin <- data[1] - abs(data[1]*perc)
  mmax <- data[2] + abs(data[2]*perc)
  
  if(!is.na(do.trunc[1])){
    # if wished, truncate the distribution under 0 and above 1
    return(rtrunc(distr = "rtrapezoid", n = N,
                  linf = min(do.trunc), lsup = max(do.trunc),
                  min = mmin, mode1 = data[1], mode2 = data[2], max = mmax))
  } else {
    return(rtrapezoid(n = N, min = mmin, mode1 = data[1], mode2 = data[2], max = mmax))
  }
  
  return(out)
}


# =================================================================================================

# introduce function that creates a triangular distribution from a modal value and an uncertainty
# for the input vector
triangulize.input <- function(modalvalues, uncertainty, N, do.trunc = NA){
  
  # 1. Test the input
  if(!all(is.numeric(modalvalues)) & !all(is.numeric(uncertainty))){
    stop("modalvalues and uncertainty have to be numeric.")
  }
  
  if((is.vector(modalvalues) & !is.vector(uncertainty)) |
     (!is.vector(modalvalues) & is.vector(uncertainty))){
    stop("modalvalues and uncertainty are not both vectors, cannot proceed to calculation.")
  }
  
  if(!(all(N == floor(N)) & length(N) == 1)){
    stop("N needs to be a single integer number.")
  }
  
  # 2. create a triangular distribution for every vector element of modalvalues
  
  # create empty matrix
  Distr <- matrix(NA, length(modalvalues), N)
  if(!is.null(names(modalvalues))){
    rownames(Distr) <- names(modalvalues)
  }
  
  # fill matrix with distributions
  for(i in 1:length(modalvalues)){
    if(is.na(modalvalues[i])){ # if the modal value is NA, everything is set to NA
      next
    } else if(modalvalues[i] == 0){ # if no input (mass), no distribution is created
      Distr[i,] <- 0
    } else {
      # else create a triangular distribution around the modal value with the uncertainty given
      Distr[i,] <- rtriang.perc(mmode = modalvalues[i], perc = uncertainty[i], N = N,
                                do.trunc = do.trunc)
    }
  }
  
  message(paste(format(Sys.time(), "%H:%M:%S"),"Input triangulization complete."))
  message(paste("-------- There are", length(which(modalvalues != 0)), "sources."))
  
  return(Distr)
}

# =================================================================================================

# introduce function that creates a triangular distribution from a modal value and an uncertainty
triangulize.TC <- function(modalvalues, uncertainty, N, do.trunc = NA){
  
  # 1. Test the input
  if(!all(is.numeric(modalvalues)) & !all(is.numeric(uncertainty))){
    stop("modalvalues and uncertainty have to be numeric.")
  }
  
  if((is.matrix(modalvalues) & !is.matrix(uncertainty)) |
     (!is.matrix(modalvalues) & is.matrix(uncertainty))){
    stop("modalvalues and uncertainty are not both matrices, cannot proceed to calculation.")
  }
  
  if(!(all(N == floor(N)) & length(N) == 1)){
    stop("N needs to be a single integer number.")
  }
  
  # 2. create a triangular distribution for every matrix/vector element of modalvalues
    
  # create empty named list with one element per compartment. Every list element will contain the outflows
  # to the different compartments.
  Distr <- sapply(colnames(modalvalues), function(x) NULL)
  
  # fill list with distributions
  for(i in 1:length(Distr)){
    # if there is no outflow from the compartment considered, don't create any data
    if(length(which(modalvalues[,i] != 0)) == 0){
      next
    } else {
      # find out what destination compartments there are as index or as name
      dest.index <- which(modalvalues[,i] != 0 | is.na(modalvalues[,i]))
      destinations <- rownames(modalvalues)[dest.index]
      # create a named list, with one element per compartment where it flows to
      Distr[[i]] <- sapply(destinations, function(x) NULL)
      
      if(length(destinations) == 1){
        # if there is only one outflow (only one value that is non-zero or NA), set all data to 1
        Distr[[i]][[1]] <- rep(1, N)
        cat(paste0("-------- TC from '", names(Distr)[[i]], "' to '", destinations, "' set to 1. \n"))
      } else {
        # else, create triangular distribution
        for(j in 1:length(destinations)){
          Distr[[i]][[j]] <- rtriang.perc(mmode = modalvalues[dest.index[j],i],
                                          perc = uncertainty[dest.index[j],i],
                                          N = N, 
                                          do.trunc = do.trunc)
        }
      }
      
      rm(destinations, dest.index)
    }
  }
  
  message(paste(format(Sys.time(), "%H:%M:%S"),"TC triangulization complete."))
  message(paste("-------- There are", length(which(sapply(Distr, is.null))), "sinks."))
  message(paste("-------- There are", length(which(is.na(modalvalues))), "undefined flows."))
  
  return(Distr)
}

# =================================================================================================

# introduce function that normalizes outflows
normalize <- function(Distr){
  
  # 1. Test the input
  if(!is.list(Distr)){
    stop("Distr should be a list of lists.")
  }
  
  # 2. Normalization step: make sure that all the flows going out of a compartment sum up to 1
  Distr.Norm <- Distr
  
  for(j in 1:length(Distr)){
    
    # if there are no outflows, leave as it is
    if(is.null(Distr[[j]])){
      next
      
    } else {
      # sum the flows that leave a compartment
      thesum <- apply(sapply(Distr[[j]], cbind), 1, sum)
      
      for(i in 1:length(Distr[[j]])){
        Distr.Norm[[j]][[i]] <- Distr[[j]][[i]]/thesum
      }
    }
  }
 
  message(paste(format(Sys.time(), "%H:%M:%S"),"Normalization complete.")) 
  
  return(Distr.Norm)
}

# =================================================================================================

# Monte-Carlo simulation
solve.MC <- function(TC.Distr, inp.Distr, N){
  
  # prepare output
  Mass <- matrix(NA, length(TC.Distr), N, dimnames = list(names(TC.Distr), NULL))
  
  # Monte-Carlo iterations
  for(k in 1:N){
    
    # create empty TC matrix
    themat <- matrix(0, length(TC.Distr), length(TC.Distr), dimnames = list(names(TC.Distr), names(TC.Distr)))
    # prepare k-th matrix for equation solving
    for(j in 1:length(TC.Distr)){
      if(is.null(TC.Distr[[j]])){
        next
      } else {
        for(i in 1:length(TC.Distr[[j]])){
          themat[names(TC.Distr[[j]])[i],j] <- TC.Distr[[j]][[i]][k]
        }
      }
    }
    
    # transform the matrix
    themat <- -themat
    diag(themat) <- 1
    
    Mass[,k] <- solve(themat, inp.Distr[,k])
  }
  
  message(paste(format(Sys.time(), "%H:%M:%S"),"Simulation complete."))
  
  return(Mass)
}


