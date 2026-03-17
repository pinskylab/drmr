##' @title Generates an adjacency matrix
##' @param x an \code{sf} object representing the patches.
##' @return An adjacency \code{matrix}.
##' @author lcgodoy
##' @export
gen_adj <- function(x) {
  if (!inherits(x, "sfc")) {
    message("Assuming patches are 1d and ordered in such a way that the first and last patches have only one neighbor.")
    aux <-
      rbind(c(0, 0), c(1, 0),
            c(1, 1), c(0, 0)) |>
      list() |>
      sf::st_polygon() |>
      sf::st_sfc()
    out <- sf::st_make_grid(aux,
                            n = c(1, length(x)),
                            what = "polygons")
  } else
    out <- x
  spdep::poly2nb(out) |>
    spdep::nb2mat(style = "B") |>
    as.matrix()
}

##' @title Get nodes for ICAR spatial random effects
##' @param adj adjacency \code{matrix}
##' @return A \code{list} containing the number of edges and the neighbors.
##' @author lcgodoy
get_nodes <- function(adj) {
  ladj <- apply(adj > 0, MARGIN = 1,
                which)
  nodes <-
    Map(\(row, col) {
      matrix(c(rep(as.integer(row),
                   length(col)),
               col), ncol = 2)
    }, names(ladj), ladj)
  nodes <- do.call(rbind, nodes)
  return(list("N_edges" = NROW(nodes),
              "neighbors" = t(nodes)))
}

##' @title Scaling factor for ICAR
##'
##' @description Using results from Rue and Held 2005 and Morris et al. 2019.
##' 
##' @param adj adjacency \code{matrix}
##' @return A \code{scalar} representing the scaling factor.
##' @author lcgodoy
get_scaling <- function(adj) {
  Q <- matrix(0, nrow = nrow(adj), ncol = ncol(adj))
  diag(Q) <- apply(adj, 1, sum)
  Q <- Q - adj
  jitter <- max(diag(Q)) * sqrt(.Machine$double.eps)
  Q <- Q +
    diag(jitter,
         ncol = NCOL(Q), nrow = NROW(Q))  
  Q_inv <- tryCatch(
      solve(Q),
      error = function(e) ginv(Q)
  )
  A <- matrix(1, nrow = 1, ncol = NCOL(Q))
  QA <- tcrossprod(Q_inv, A)
  Q_inv_const <- Q_inv - tcrossprod(QA, QA) / as.numeric(A %*% QA)
  vars <- diag(Q_inv_const)
  exp(mean(log(vars[vars > 0])))
}
