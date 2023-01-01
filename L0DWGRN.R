

#' Title
#'
#' @param X 
#' @param Y 
#' @param penalty 
#' @param nFolds 
#' @param seed 
#' @param maxSuppSize 
#' @param nGamma 
#' @param gammaMin 
#' @param gammaMax 
#'
#' @return
#' @export
#'
#' @examples
LO_fit <- function(X, Y,
                   penalty = penalty,
                   nFolds = 10,
                   seed = 1,
                   maxSuppSize = maxSuppSize,
                   nGamma = 5,
                   gammaMin = 0.0001, gammaMax = 10) {
  tryCatch(
    {
      fit <- L0Learn.cvfit(X, Y,
        penalty = penalty,
        maxSuppSize = maxSuppSize,
        nFolds = 10,
        seed = 1,
        nGamma = 5,
        gammaMin = 0.0001,
        gammaMax = 10
      )
      fit_inf <- print(fit)
      optimalGammaIndex <- which(unlist(lapply(fit$cvMeans, min)) == min(unlist(lapply(fit$cvMeans, min))))
      gamma <- fit$fit$gamma[optimalGammaIndex]
      lambda_list <- fit_inf[which(fit_inf$gamma == gamma), ]
      if (is.null(maxSuppSize)) {
        lambda <- min(lambda_list$lambda)
      } else {
        if (maxSuppSize %in% lambda_list$maxSuppSize) {
          lambda <- lambda_list$maxSuppSize[which(lambda_list$maxSuppSize == maxSuppSize)]
        } else {
          lambda <- min(lambda_list$lambda)
        }
      }
      temp <- coef(fit, lambda = lambda, gamma = gamma)
    },
    error = function(e) {
      fit <- L0Learn.fit(X, Y,
        penalty = penalty,
        maxSuppSize = maxSuppSize,
        nGamma = 5,
        gammaMin = 0.0001,
        gammaMax = 10
      )
      fit_inf <- print(fit)
      fit_inf <- fit_inf[order(fit_inf$suppSize, decreasing = TRUE), ]
      lambda <- fit_inf$lambda[1]
      gamma <- fit_inf$gamma[1]
      temp <- coef(fit,
        lambda = lambda,
        gamma = gamma
      )
    }
  )
}

#' Title
#'
#' @param matrix 
#' @param penalty 
#' @param regulators 
#' @param targets 
#' @param maxSuppSize 
#'
#' @return
#' @export
#'
#' @examples
L0DWGRN <- function(matrix,
                  penalty = NULL,
                  regulators = NULL,
                  targets = NULL,
                  maxSuppSize = NULL) {
  library(L0Learn)
  matrix <- as.data.frame(t(matrix))
  weightdf <- c()
  if (is.null(penalty)) {
    penalty <- "L0"
  }
  if (is.null(maxSuppSize)) {
    maxSuppSize <- dim(matrix)[2]
  }
  if (is.null(targets)) {
    targets <- colnames(matrix)
  }
  if (!is.null(regulators)) {
    matrix_reg <- matrix[, regulators]
    for (i in 1:length(targets)) {
      if (targets[i] %in% regulators) {
        X <- as.matrix(matrix_reg[, -which(colnames(matrix_reg) == targets[i])])
      } else {
        X <- as.matrix(matrix_reg)
      }
      Y <- matrix[, targets[i]]
      temp <- LO_fit(X, Y,
        penalty = penalty,
        nFolds = 10,
        seed = 1,
        maxSuppSize = maxSuppSize,
        nGamma = 5,
        gammaMin = 0.0001,
        gammaMax = 10
      )
      temp <- as.vector(temp)
      wghts <- temp[-1]
      wghts <- abs(wghts)
      # wghts <- wghts / max(wghts)
      if (F) {
        wghts <- wghts / max(wghts)
        indices <- sort.list(wghts, decreasing = TRUE)
        zeros <- which(wghts <= 0.8)
        # wghts[1:length(wghts)] <- 1
        wghts[zeros] <- 0
      }
      if (length(wghts) != dim(X)[2]) {
        weightd <- data.frame(regulatoryGene = colnames(X), targetGene = targets[i], weight = 0)
        # weightd <- data.frame(regulatoryGene = targets[i], targetGene = colnames(X), weight = 0)
      } else {
        weightd <- data.frame(regulatoryGene = colnames(X), targetGene = targets[i], weight = wghts)
        # weightd <- data.frame(regulatoryGene = targets[i], targetGene = colnames(X), weight = wghts)
      }
      # weightd$weight <- weightd$weight / max(weightd$weight)
      weightdf <- rbind.data.frame(weightdf, weightd)
      if (i == length(regulators)) {
        weightdf <- weightdf[order(weightdf$weight, decreasing = TRUE), ]
      }
    }
  } else {
    regulators <- colnames(matrix)
    for (i in 1:length(regulators)) {
      X <- as.matrix(matrix[, -which(colnames(matrix) == regulators[i])])
      Y <- matrix[, regulators[i]]
      temp <- LO_fit(X, Y,
        penalty = penalty,
        nFolds = 10, seed = 1,
        maxSuppSize = maxSuppSize,
        nGamma = 5,
        gammaMin = 0.0001, gammaMax = 10
      )
      temp <- as.vector(temp)
      wghts <- temp[-1]
      wghts <- abs(wghts)
      # wghts <- wghts / max(wghts)
      if (F) {
        wghts <- wghts / max(wghts)
        indices <- sort.list(wghts, decreasing = TRUE)
        zeros <- which(wghts <= 0.8)
        # wghts[1:length(wghts)] <- 1
        wghts[zeros] <- 0
      }
      if (length(wghts) != dim(X)[2]) {
        weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulators[i], weight = 0)
        # weightd <- data.frame(regulatoryGene = regulators[i], targetGene = colnames(X), weight = 0)
      } else {
        weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulators[i], weight = wghts)
        # weightd <- data.frame(regulatoryGene = regulators[i], targetGene = colnames(X), weight = wghts)
      }
      # weightd$weight <- weightd$weight / max(weightd$weight)
      weightdf <- rbind.data.frame(weightdf, weightd)
      if (i == length(regulators)) {
        weightdf <- weightdf[order(weightdf$weight, decreasing = TRUE), ]
      }
    }
  }
  return(weightdf)
}

#' Plot function --------------------------------------------------
#' Title
#'
#' @param data 
#' @param plotType 
#'
#' @return
#' @export
#'
#' @examples
L0Plot <- function(data, plotType = NULL) {
  if (is.null(plotType)) {
    plotType <- boxplot
  }
  if (plotType == boxplot) {
    p <- ggplot(
      data,
      aes(x = Method, y = AUPRC)
    ) +
      # geom_violin(aes(fill = Method),
      #     trim = FALSE
      # ) +
      geom_boxplot(aes(fill = Method),
        width = 0.8
      ) +
      stat_compare_means(
        method = "wilcox.test",
        label = "p.signif",
        comparisons = my_comparisons,
        bracket.size = 0.6,
        sizen = 4,
        color = "#6699cc"
      ) +
      scale_fill_manual(values = mycol) +
      scale_x_discrete(labels = methods) +
      labs(x = "Methods", y = "AUPRC") +
      theme(legend.position = "bottom") +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = 10
        )
      )
  }
  p
}
