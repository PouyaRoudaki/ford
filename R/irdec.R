#' Title
#'
#' @param Y
#' @param X
#' @param dist.type.X
#' @param na.rm
#'
#' @returns
#' @export
#'
#' @examples
irdec <- function(Y, X, dist.type.X = "continuous", na.rm = TRUE){

  # if inputs are not in proper matrix format change if possible
  # otherwise send error
  if(!is.vector(Y)) {
    Y = as.vector(Y)
  }
  if(!is.matrix(X)) {
    X = as.matrix(X)
  }

  if((length(Y) != nrow(X))) stop("Number of rows of Y and X should be equal.")

  if (na.rm == TRUE) {
      # cases without na
      wo_na = complete.cases(Y,X)
      # NAs are removed here:
      X = as.matrix(X[wo_na,])
      Y = Y[wo_na]
  }

  n = length(Y)

  if (n < 3) {
    stop("Not enough points to find 3 nearest neighbors which is required because of j != i and j != N(i) condition.")
  }

  if(dist.type.X == "continuous"){ # for continuous data repetition is very unlikely so we do not need randomized neighbor extraction

    # Extract the matrix of first three neighbors(including itself) for all observations. Note that by convention of RANN::nn2() for any x_i the first neighbor is itself.
    nn_X <- RANN::nn2(X, query = X, k = 3)

    # Extract the vector of "first" neighbors as N(i) for all data
    nn1_index_X = nn_X$nn.idx[, 2]

    # Extract the vector of "second" neighbors as N'(i) for all data where N'(i) is useful when N(i) = j and since we are looking for N^{(-j)}(i) we can switch to this one.
    nn2_index_X = nn_X$nn.idx[, 3]

  } else if(dist.type.X == "discrete"){ # for discrete data repetition occurs with high probability we need randomized neighbor extraction

    # Extract the matrix of first two neighbors(excluding itself) for all observations.
    rnn_X <- randomized_nn(X)$two_neighbors

    # Extract the vector of "first" neighbors as N(i) for all data
    nn1_index_X <- rnn_X[, 1]

    # Extract the vector of "second" neighbors as N'(i) for all data where N'(i) is useful when N(i) = j and since we are looking for N^{(-j)}(i) we can switch to this one.
    nn2_index_X <- rnn_X[, 2]

  } else{
    stop("Please choose your dist.type.X from 'continuous' or 'discrete'.")
  }

  # Find Y's ranks
  R_Y <- rank(Y, ties.method = "max")

  # vector of F_{n,j}(Y_j)
  G_Y <- (R_Y - 1) / (n - 1)

  # vector of denominators:  F_{n,j}(Y_j) * (1 - F_{n,j}(Y_j))
  D <- G_Y * (1 - G_Y)

  # initials vector of numerator:  \sum_{j, R_j != 1 or n} \sum_{i!=j and N(i) != j} \bone{ Y_{j} \in K_{i}} where K{i} = [min(Y_i, Y_{N(i)}),max(Y_i, Y_{N(i)})]
  C_Y <- rep(0, n)

  # Find the maximum and minimum which is required for K(i)
  Rmatrix <- rbind(R_Y, R_Y[nn1_index_X])
  Rmax <- apply(Rmatrix, 2, max)
  Rmin <- apply(Rmatrix, 2, min)

  # Find the vector of numerator:  \sum_{i!=j and N(i) != j} \bone{ Y_{j} \in K_{i}} where K{i} = [min(Y_i, Y_{N(i)}),max(Y_i, Y_{N(i)})]
  for (j in 1:n) {

    # First find:  \sum_{i != j} \bone{ Y_{j} \in K_{i}}
    # The minus 1 is for the case that i == j. In this case \bone{ Y_{j} \in K_{i}} = 1 so we need to reduce the sum by 1.
    C_Y[j] <- length(which((Rmax >= R_Y[j]) & (Rmin <= R_Y[j]))) - 1

    # Now for each i we should subtract cases that N(i) = j. In these cases (possibly several) \bone{ Y_{j} \in K_{i}} = 1.
    # So we subtract all of them and check what is the next neighbor since the formula is written based on N^{(-j)}(i) which is the closest neighbor of i without i and j.
    idx <- which(nn_index_X == j)


    if(length(idx) > 0) {

      # we subtract all N(i) = j
      C_Y[j] <- C_Y[j] - length(idx)
      for (i in idx) {
        # we add 1 if Y_{j} \in [min(Y_i, Y_{N'(i)}),max(Y_i, Y_{N'(i)})]
        C_Y[j] <- C_Y[j] + as.numeric((R_Y[j] <= max(R_Y[i], R_Y[nn2_index_X[i]])) &
                                        (R_Y[j] >= min(R_Y[i], R_Y[nn2_index_X[i]])))
      }
    }
  }

  # remove bad indices: first rank and last rank of Y.
  bad_id <- which(R_Y == 1 | R_Y == n)

  # Find irdec: 1 - 1/(2 * (n - 1) * (n - 2)) \sum_{j, R_j != 1 or n} 1/ F_{n,j}(Y_j) * (1 - F_{n,j}(Y_j)) * \sum_{i!=j and N(i) != j} \bone{ Y_{j} \in K_{i}} where K{i} = [min(Y_i, Y_{N(i)}),max(Y_i, Y_{N(i)})]
  nu_n = 1 - sum(C_Y[-bad_id] / t1[-bad_id]) / (2 * (n - 1) * (n - 2))

  # Return the irdec
  return(nu_n)
}
