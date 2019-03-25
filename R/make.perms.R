#' Making permutated index
#' 
#' 
#' @param dat1 an element of \code{ydata1}
#' @param dat2 an element of \code{ydata2}
#' @param nperm a scalar, number of permutation
#' @param .index see \code{grid} in \code{DVDtest}
#' @param adj see \code{permadj} in \code{DVDtest}
#' @return a matrix, permuted indices
#' @author Philip Reiss, Meng Xu
#' @seealso \code{DVDtest}
#' @keywords internal
make.perms <-
function(dat1, dat2, nperm, .index, adj) {
  n1 <- length(unique(dat1$.obs))
  n2 <- length(unique(dat2$.obs))
  matt <- matrix(NA,nperm, n1)
  for (i in 1:nperm) {
    bothdat <- rbind(dat1, dat2)
    check <- 0
    if (adj == T){
      while (check == 0) {
        which1 <- sample(n1+n2, n1)
        perm1 <- bothdat[bothdat$.obs %in% unique(bothdat$.obs)[which1],]
        perm2 <- bothdat[!bothdat$.obs %in% unique(bothdat$.obs)[which1],]
        if (min(perm1$.index) <= min(.index) & max(perm1$.index) >= max(.index) 
            & min(perm2$.index) <= min(.index) & max(perm2$.index) >= max(.index)) check <- 1
      }
    }
    else which1 <- sample(n1+n2, n1)
    matt[i,] <- which1
  }
  matt
}
