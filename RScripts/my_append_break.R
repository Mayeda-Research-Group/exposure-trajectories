my_append_break <- function(data, brk, warp.model = warp.model, id = NULL, 
                            typ = "pred"){
  k <- length(brk)
  app <- data[data$first, ]
  if (!is.null(id)) {
    idx <- app$HHIDPN %in% id
    app <- app[idx, ]
  }
  nap <- nrow(app)
  app$first <- FALSE
  app$typ <- typ
  app$occ <- NA
  app <- app[rep.int(seq_len(nap), length(brk)), ]
  app$age <- rep(brk, each = nap)
  app$age2 <- predict(warp.model, newdata = app)
  X <- splines::bs(app$age, knots = brk, Boundary.knots = c(brk[1], 
                                                            brk[k] + 1e-04), degree = 1)
  X <- X[, -(k + 1)]
  app[, paste0("x", seq_len(ncol(X)))] <- X
  app[, c("hgt.z", "wgt.z", "bmi.z")] <- NA
  app <- rbind(data, app)
  data <- app[order(app$id, app$age), ]
  return(data)
}