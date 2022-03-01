rarefy_putback <- function(x, rnumber)
  {
    if (!identical(all.equal(x, round(x)), TRUE)) 
      stop("function is meaningful only for integers (counts)")
    x <- as.matrix(x)
    if (ncol(x) == 1)
      x <- t(x)
    if (length(rnumber) > 1 && length(rnumber) != nrow(x))
      stop(gettextf(
        "length of 'sample' and number of rows of 'x' do not match"))
    rnumber <- rep(rnumber, length=nrow(x))
    colnames(x) <- colnames(x, do.NULL = FALSE)
    nm <- colnames(x)
    for (i in 1:nrow(x)) {
      row <- sample(rep(nm, times=x[i,]), rnumber[i],replace=T)
      row <- table(row)
      ind <- names(row)
      x[i,] <- 0
      x[i,ind] <- row
    }
    x
  }