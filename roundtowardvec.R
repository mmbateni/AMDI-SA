roundtowardvec <- function(X, roundvec = NULL, type = "round") {
  if (missing(X)) {
    cat("Usage: roundtowardvec(X, [roundvec], [type])\n")
    return(NULL)
  } else if (is.null(X)) {
    return(NULL)
  }
  
  type <- tolower(type)
  if (is.null(roundvec)) {
    if (type == "round") {
      newnums <- round(X)
    } else if (type == "away") {
      newnums <- sign(X) * ceiling(abs(X))
    } else if (type == "fix") {
      newnums <- trunc(X)
    } else if (type == "floor") {
      newnums <- floor(X)
    } else if (type == "ceil") {
      newnums <- ceiling(X)
    } else {
      stop("Round type not recognized. Options are:\n'round' - round to nearest value\n'floor' - round toward -Inf\n'ceil'  - round toward Inf\n'fix'   - round toward 0\n'away'  - round away from 0")
    }
  } else {
    roundvec <- unique(roundvec)
    newnums <- X
    if (type == "round") {
      roundvec <- sort(roundvec)
      for (i in rev(seq_along(X))) {
        if (!(X[i] %in% roundvec)) {
          diffs <- abs(roundvec - X[i])
          if (X[i] >= 0) {
            ind <- length(diffs) - which.min(diffs)
          } else {
            ind <- which.min(diffs)
          }
          newnums[i] <- roundvec[ind]
        }
      }
    } else if (type == "fix") {
      roundvec <- sort(c(roundvec, 0))
      for (i in rev(seq_along(X))) {
        if (!(X[i] %in% roundvec)) {
          if (X[i] > 0) {
            if (X[i] > min(roundvec)) {
              newnums[i] <- roundvec[max(which(X[i] > roundvec))]
            } else {
              newnums[i] <- 0
            }
          } else if (X[i] < 0) {
            if (X[i] < max(roundvec)) {
              newnums[i] <- roundvec[min(which(X[i] < roundvec))]
            } else {
              newnums[i] <- 0
            }
          }
        }
      }
    } else if (type == "ceil") {
      roundvec <- sort(roundvec)
      for (i in rev(seq_along(X))) {
        if (!is.na(X[i]) && !(X[i] %in% roundvec)) {
          if (X[i] < max(roundvec)) {
            newnums[i] <- roundvec[min(which(X[i] < roundvec))]
          } else {
            newnums[i] <- Inf
          }
        }
      }
    } else if (type == "floor") {
      roundvec <- sort(roundvec)
      for (i in rev(seq_along(X))) {
        if (!is.na(X[i]) && !(X[i] %in% roundvec)) {
          if (X[i] > min(roundvec)) {
            newnums[i] <- roundvec[max(which(X[i] > roundvec))]
          } else {
            newnums[i] <- -Inf
          }
        }
      }
    } else if (type == "away") {
      roundvec <- sort(roundvec)
      for (i in rev(seq_along(X))) {
        if (!(X[i] %in% roundvec)) {
          if (X[i] > 0) {
            if (X[i] < max(roundvec)) {
              newnums[i] <- roundvec[min(which(X[i] < roundvec))]
            } else {
              newnums[i] <- Inf
            }
          } else if (X[i] < 0) {
            if (X[i] > min(roundvec)) {
              newnums[i] <- roundvec[max(which(X[i] > roundvec))]
            } else {
              newnums[i] <- -Inf
            }
          } else if (X[i] == 0) {
            diffs <- abs(roundvec - X[i])
            ind <- length(diffs) - which.min(diffs)
            newnums[i] <- roundvec[ind]
          }
        }
      }
    } else {
      stop("Round type not recognized. Options are:\n'round' - round to nearest value\n'floor' - round toward -Inf\n'ceil'  - round toward Inf\n'fix'   - round toward 0\n'away'  - round away from 0")
    }
  }
  return(newnums)
}
