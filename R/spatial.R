##' @title Generates an adjacency matrix
##' @param x the patches: an \code{sf} or \code{sfc} object with one polygon
##'   per patch. For 1-D patches, \code{x} can also be a vector with one
##'   element per patch, in their order along the line.
##' @return An adjacency \code{matrix}. A patch with no neighbors, such as an
##'   island, gets a row and a column of zeros.
##' @author lcgodoy
##' @export
gen_adj <- function(x) {
  if (inherits(x, "sf"))
    x <- sf::st_geometry(x)
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
    spdep::nb2mat(style = "B", zero.policy = TRUE) |>
    as.matrix()
}

##' @title Get nodes for ICAR spatial random effects
##' @description Two sites are neighbors if either entry of \code{adj} for
##'   the pair is positive. The diagonal is ignored.
##' @param adj adjacency \code{matrix}
##' @return A \code{list} containing the number of edges and the neighbors.
##'   Each pair of neighbors appears once, so the ICAR precision is the graph
##'   Laplacian, the matrix that [get_scaling()] uses.
##' @author lcgodoy
get_nodes <- function(adj) {
  adj <- as.matrix(adj)
  stopifnot(NROW(adj) == NCOL(adj))
  a <- (adj > 0) | t(adj > 0)
  ## one column per pair of neighbors (i < j), ordered by i
  nodes <- which(a & upper.tri(a), arr.ind = TRUE)
  nodes <- nodes[order(nodes[, 1], nodes[, 2]), , drop = FALSE]
  return(list("N_edges" = NROW(nodes),
              "neighbors" = unname(t(nodes))))
}

##' @title Connected components of a graph
##' @param a a symmetric \code{logical} adjacency \code{matrix}.
##' @return An \code{integer} vector with the component of each site.
##' @author lcgodoy
get_components <- function(a) {
  n <- NROW(a)
  comp <- integer(n)
  k <- 0L
  for (s in seq_len(n)) {
    if (comp[s] == 0L) {
      k <- k + 1L
      comp[s] <- k
      queue <- s
      while (length(queue) > 0) {
        nb <- which(a[queue[1], ] & comp == 0L)
        comp[nb] <- k
        queue <- c(queue[-1], nb)
      }
    }
  }
  return(comp)
}

##' @title Data for the ICAR spatial random effects
##' @description Follows Freni-Sterrantino et al. (2018) for graphs that are
##'   not connected. Each connected component with two or more sites gets its
##'   own sum-to-zero constraint and its own scaling factor. The scaling
##'   makes the geometric mean of the marginal variances of the scaled effect
##'   equal to 1 in each component. The marginal variance of a site depends
##'   on its position in the graph. Each site with no neighbors gets an
##'   independent standard normal effect. So \code{sigma_s} is a typical
##'   marginal standard deviation of the spatial effect.
##' @param adj adjacency \code{matrix}. Two sites are neighbors if either
##'   entry for the pair is positive. Weights and the diagonal are ignored.
##' @return A \code{list} with the ICAR inputs for the Stan models:
##'   \code{N_edges}, \code{neighbors}, \code{scaling} (one value per site:
##'   the scaling factor of its component, or 1 for a site with no
##'   neighbors), \code{n_comp_s} (number of components with two or more
##'   sites), \code{comp_size_s} (their sizes), and \code{comp_order_s} (the
##'   sites, component by component, with the sites without neighbors last).
##' @references Freni-Sterrantino, A., Ventrucci, M., and Rue, H. (2018). A
##'   note on intrinsic conditional autoregressive models for disconnected
##'   graphs. \emph{Spatial and Spatio-temporal Epidemiology}, 26, 25-34.
##' @author lcgodoy
icar_data <- function(adj) {
  adj <- as.matrix(adj)
  stopifnot(NROW(adj) == NCOL(adj))
  a <- (adj > 0) | t(adj > 0)
  diag(a) <- FALSE
  n <- NROW(a)
  comp <- get_components(a)
  sizes <- tabulate(comp)
  multi <- which(sizes > 1)
  scaling <- rep(1, n)
  for (k in multi) {
    idx <- which(comp == k)
    scaling[idx] <- get_scaling(1 * a[idx, idx, drop = FALSE])
  }
  comp_order <- c(unlist(lapply(multi, function(k) which(comp == k))),
                  which(sizes[comp] == 1))
  nodes <- get_nodes(a)
  return(list(N_edges = array(as.integer(nodes$N_edges), dim = 1),
              neighbors = nodes$neighbors,
              scaling = array(scaling, dim = n),
              n_comp_s = length(multi),
              comp_size_s = array(as.integer(sizes[multi]),
                                  dim = length(multi)),
              comp_order_s = array(as.integer(comp_order), dim = n)))
}

##' @title Placeholder ICAR data for models without ICAR effects
##' @return A \code{list} with the same elements as [icar_data()], with the
##'   sizes that the Stan models expect when \code{sp_re} is off.
##' @author lcgodoy
icar_data_none <- function() {
  return(list(N_edges = integer(0),
              neighbors = matrix(1L, ncol = 1, nrow = 1),
              scaling = numeric(0),
              n_comp_s = 0L,
              comp_size_s = integer(0),
              comp_order_s = integer(0)))
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

##' @title Check an adjacency matrix
##' @param adj_mat adjacency \code{matrix}.
##' @param n_sites number of sites.
##' @return \code{adj_mat}, invisibly. The function stops if \code{adj_mat}
##'   does not have one row and one column per site, or if it has values that
##'   are missing, infinite, or negative.
##' @author lcgodoy
check_adj_mat <- function(adj_mat, n_sites) {
  a <- as.matrix(adj_mat)
  if (NROW(a) != n_sites || NCOL(a) != n_sites) {
    stop("'adj_mat' must have one row and one column per site.",
         call. = FALSE)
  }
  if (any(!is.finite(a)) || any(a < 0)) {
    stop("'adj_mat' must contain finite, non-negative values only.",
         call. = FALSE)
  }
  return(invisible(adj_mat))
}
