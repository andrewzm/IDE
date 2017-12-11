IDE <- function(f, data, dt, process_basis = NULL, kernel_basis = NULL, grid_size = 41, forecast = 0, hindcast = 0) {

  if(!is(f,"formula"))
    stop("f needs to be of class formula")

  if(!is(data, "ST"))
    stop("data needs to be of class ST")

  if(!is(dt, "difftime"))
    stop("dt needs to be of class difftime")

  if(is.null(process_basis)) {
        process_basis <- auto_basis(manifold = plane(),
                                    data = data,
                                    regular = 1,
                                    nres = 2)
    }

    G_const <- new("Basis",
                 manifold = plane(),
                 fn = list(function(s) rep(1, nrow(s))),
                 pars = list(),
                 df = data.frame(),
                 n = 1)

    if(is.null(kernel_basis)) {
        kernel_basis <- list(G_const, G_const, G_const, G_const)
    }

    ## Set functions
    set <- function(k = NULL, sigma2_eps = NULL, sigma2_eta = NULL) {
        if(!is.null(sigma2_eps)) {
            sigma2_eps <<- sigma2_eps
            Q_eps <<- 1/(sigma2_eps) * Diagonal(m)
        }

        if(!is.null(sigma2_eta)) {
            sigma2_eta <<- sigma2_eta
            Q_eta <<- 1/(sigma2_eta) * Diagonal(r)
        }

        if(!is.null(k)) {
            if(!(length(k) == length(kernel_basis)))
                stop("Need to supply as many parameters as kernel basis functions")
            k <<- k
            M <<- Mfun(kernel_basis, k)
        }

        if(exists("M") & exists("Q_eta"))
            if(!is.null(sigma2_eta) | !is.null(k)) {
                Q <<- construct_Q(Q_eta, M, T)
                if(max(abs(eigen(M)$values)) < 2) { # safegaurd against blow-up
                    Qpost <<- crossprod(chol(Q_eps) %*% PHI_obs) + Q
                    Qpostchol <<- FRK:::cholPermute(Qpost)
                    update_betahat()
                    update_alpha()
                }
            }

    }

    get <- function(obj) {
        switch(obj, "sigma2_eps" = sigma2_eps,
               "sigma2_eta" = sigma2_eta,
               "Q_eps" = Q_eps,
               "Q_eta" = Q_eta,
               "alphahat" = alphahat,
               "alphavar" = alphavar,
               "betahat" = betahat,
               "coordnames" = coordnames,
               "data" = data,
               "PHI_obs" = PHI_obs,
               "plausible_ranges" = plausible_ranges,
               "process_basis" = process_basis,
               "kernel_basis" = kernel_basis,
               "Qpost" = Qpost,
               "X_obs" = X_obs,
               "Q" = Q,
               "M" = M,
               "k" = k,
               "f" = f,
               "r" = r,
               "s" = s,
               "m" = m,
               "T" = T,
               "Z" = Z)
    }

    update_alpha <- function() {
        alphahat <<-  FRK:::cholsolve(Q = Qpost,
                                      y = t(PHI_obs) %*% Q_eps %*% (Z - X_obs %*% betahat),
                                      perm = TRUE,
                                      cholQp = Qpostchol$Qpermchol,
                                      P = Qpostchol$P)
    }

    #update_beta <- function(X_obs, PHI_obs, Q_eps, Qpost_cholsolve, Z) {
    update_betahat <- function() {
        Qpost_cholsolve <- function(y) {
            FRK:::cholsolve(Q = Qpost,
                            y = y,
                            perm = TRUE,
                            cholQp = Qpostchol$Qpermchol,
                            P = Qpostchol$P)
        }

        tPHIQepsZ <- t(PHI_obs) %*% Q_eps %*% Z
        tXQeps <- (t(X_obs) %*% Q_eps)
        tXQepsPHI <- tXQeps %*% PHI_obs
        Part1 <- solve(tXQeps %*% X_obs - tXQepsPHI %*% Qpost_cholsolve(t(tXQepsPHI)))
        Part2 <- tXQeps %*% Z - tXQepsPHI %*% Qpost_cholsolve(y = tPHIQepsZ)
        betahat <<- Part1 %*% Part2
    }


    ## Log likelihood
    negloglik <- function() {
        if(max(abs(eigen(M)$values)) < 1 &
           all(sapply(3:length(k), function(i) all(abs(k[[i]]) < axis_ranges[i-2])))) {
            Ztilde <- Z - X_obs %*% betahat

            Qchol <- FRK:::cholPermute(Q)
            -((0.5*logdet(Qchol$Qpermchol) +
                   0.5*logdet(chol(Q_eps)) -
                   0.5*logdet(Qpostchol$Qpermchol) -
                   0.5*t(Ztilde) %*% (Q_eps %*% Ztilde) +
                   0.5*t(Ztilde) %*% Q_eps %*% PHI_obs %*% alphahat) %>% as.numeric())
        } else {
            1e10
        }
    }

    # Cast to STIDF
    if(is(data, "STFDF"))
        data <- as(data, "STIDF")

    # Remove NAs
    data <- subset(data, !is.na(data@data[,all.vars(f)[1]]))

    r <- nbasis(process_basis)
    m <- length(data)
    bbox <- data@sp@bbox
    coordnames <- colnames(coordinates(data))
    s <- construct_grid(bbox, grid_size,
                        coordnames = coordnames)
    axis_ranges <- apply(bbox, 1, diff)
    time_points <- seq(min(time(data)) - dt*hindcast,
                       max(time(data)) + dt*forecast,
                       by = dt)

      #sort(unique(time(data)))
    if(!all(time(data) %in% time_points))
      stop("Data time points need to be equidistant on chose time interval")
    T <- length(time_points)

    Z <- data[["z"]]
    PHI_obs_list <- lapply(1:T, function(i) {
        if(any(time(data) %in% time_points[i])) {
          eval_basis(process_basis, coordinates(data[,time_points[i]]))
        } else { Zeromat(0,r)}})
    PHI_obs <- do.call("bdiag", PHI_obs_list)
    X_obs <-  model.matrix(f, data)
    betahat <- solve(crossprod(X_obs,X_obs)) %*% t(X_obs) %*% Z # OLS


    Ztilde <- Z - X_obs %*% betahat
    set(sigma2_eps = var(Ztilde)/2,
        sigma2_eta = var(Ztilde)/2)

    plausible_ranges <- data.frame(k1 = 150 / (s$area*grid_size^2) * c(0.01,10),
                                   k2 = max((axis_ranges / 2)^2)*c(0.001,1),
                                   k3 = max(axis_ranges)*c(-0.5,0.5),
                                   k4 = max(axis_ranges)*c(-0.5,0.5),
                                   sigma2_eps =  c(var(Ztilde)) * c(0.01, 2),
                                   sigma2_eta =  c(var(Ztilde)/100, var(Ztilde)*2))


    Mfun <- construct_M(process_basis, s)
    kinit <- c(150 / (s$area*grid_size^2),
               0.002 * (s$area*grid_size^2),
               0, 0)
    k <- lapply(1:length(kernel_basis),
                function(i) rep(kinit[i], nbasis(kernel_basis[[i]])))
    # 150, 0.002, -0.1, 0.1
    set(k = k)

    update_alpha()        # Initial guess
    update_betahat()      # GLS

    IDEobj <- list(set = set,
                   get = get,
                   update_alpha = update_alpha,
                   update_betahat = update_betahat,
                   negloglik = negloglik)
    class(IDEobj) <- "IDE"
    IDEobj
}

fit.IDE <- function(object, method = "DEoptim") {

    ## Optimise log likelihood
    optimfun <- function(theta, IDEmodel) {
        ki <- theta[1:4]
        sigma2_eps <- exp(theta[5])
        sigma2_eta <- exp(theta[6])
        ki[1] <- ki[1]*1000
        ki[2] <- exp(ki[2]*10)
        IDEmodel$set(k =  lapply(1:4,function(i) ki[i]),
                     sigma2_eps = sigma2_eps,
                     sigma2_eta = sigma2_eta)
        IDEmodel$negloglik()
    }

    P <- IDEmodel$get("plausible_ranges")


    if(method == "DEoptim") {

        O <- DEoptim(fn = optimfun,
                     lower = c(P$k1[1]/1000,log(P$k2[1])/10, P$k3[1], P$k4[1],
                               log(P$sigma2_eps[1]), log(P$sigma2_eta[1])),
                     upper = c(P$k1[2]/1000,log(P$k2[2])/10, P$k3[2], P$k4[2],
                               log(P$sigma2_eps[2]), log(P$sigma2_eta[2])),
                     control = list(parallelType = 1,
                                    packages = c("Matrix","FRK", "sp", "dplyr"),
                                    parVar = c("IDEmodel","construct_kernel",
                                               "rep.col", "construct_Q","Zeromat",
                                               "logdet")),
                     IDEmodel = object)
        theta <- O$optim$bestmem
    } else {
        stop("Only DEoptim implemented for now")
    }

    ki <- theta[1:4]
    sigma2_eps <- exp(theta[5])
    sigma2_eta <- exp(theta[6])
    ki[1] <- ki[1]*1000
    ki[2] <- exp(ki[2]*10)
    object$set(k =  lapply(1:4,function(i) ki[i]),
                 sigma2_eps = sigma2_eps,
                 sigma2_eta = sigma2_eta)

    list(optim_results = O,
         IDEmodel = object)
}

predict.IDE <- function(object, newdata = NULL) {
    alphahat <- object$get("alphahat")
    betahat <- object$get("betahat")
    coordnames <- object$get("coordnames")
    process_basis <- object$get("process_basis")
    data <- object$get("data")
    Qpost <- object$get("Qpost")
    f <- object$get("f")

    if(is.null(newdata)) {
        s <- object$get("s")
        T <- object$get("T")
        r <- object$get("r")
        PHI_pred_1 <- eval_basis(process_basis, s$s_grid_mat)
        PHI_pred <- do.call("bdiag", lapply(1:T, function(x) PHI_pred_1))
        newdata <- s$s_grid_df[, coordnames] %>%
            expand.grid.df(data.frame(t = 0:(T-1)))

    } else {
        dates <- unique(time(data))
        PHI_pred_1 <- list()
        for(i in seq_along(dates)) {
            newdata_1 <- subset(newdata, time(newdata) == dates[i])
            if(length(newdata_1) == 0){
                PHI_pred_1[[i]] <- Zeromat(0,r)
            } else {
                PHI_pred_1[[i]] <- eval_basis(process_basis, newdata_1)
            }
        }
        PHI_pred <- do.call("bdiag", PHI_pred_1)
    }
    if(is(newdata, "ST")) {
        newdata[[all.vars(f)[1]]] <- 0            # dummy data
    } else {
        newdata[all.vars(f)[1]] <- 0            # dummy data
    }
    X_pred <- model.matrix(f, newdata)

    newdata$Ypred <- (X_pred %*% betahat + PHI_pred %*% alphahat) %>% as.numeric()
    Qpost_dense <- densify(Qpost,t(PHI_pred) %*% PHI_pred)
    Qpostchol <- FRK:::cholPermute(Qpost_dense)
    Ssparseinv <- FRK:::Takahashi_Davis(Q = Qpost_dense,
                                        cholQp = Qpostchol$Qpermchol,
                                        P = Qpostchol$P)
    newdata$Ypredse <- rowSums(PHI_pred * (PHI_pred %*% Ssparseinv)) %>%
        as.numeric() %>% sqrt()
    newdata
}

show_kernel <- function(IDEmodel) {
    kernel_basis <- IDEmodel$get("kernel_basis")
    k <- IDEmodel$get("k")
    s <- IDEmodel$get("s")
    ndim <- dimensions(kernel_basis[[1]])
    K <- construct_kernel(kernel_basis, k)
    centred_grid <- scale(s$s_grid_mat, scale = FALSE)
    Kmat <- K(matrix(rep(0,ndim),1,ndim),
              centred_grid)
    s$s_grid_df$m <- as.numeric(Kmat)
    if(ndim == 2) {
        ggplot(s$s_grid_df) + geom_tile(aes(centred_grid[,1],centred_grid[,2],fill=m)) +
            theme_bw() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
            scale_fill_continuous(low = "white", high = "black") +
            xlab("x1 - s1") + ylab("x2 - s2")
    } else {
        stop("Plotting only implemented for 2D spatial fields")
    }
}

simIDE <- function(T = 9, nobs = 100) {

    set.seed(1)
    zlocs <- data.frame(s1 = runif(nobs),
                        s2 = runif(nobs))

    ## Spatial decomposition
    Y_basis <- auto_basis(manifold = plane(),
                          data = SpatialPoints(zlocs),
                          regular = 1,
                          nres = 2)
    r <- nbasis(Y_basis)

    ## Kernel decomposition
    G_const <- new("Basis",
                   manifold = plane(),
                   fn = list(function(s) rep(1, nrow(s))),
                   pars = list(),
                   df = data.frame(),
                   n = 1)

    ## Regression coeffocients
    beta <- c(0.2,0.2,0.2)

    ## Other parameters
    nobs <- 100
    sigma2_eta <- 0.01^2
    sigma2_eps <- 0.01^2
    Sigma_eta <- sigma2_eta * Diagonal(r)
    Sigma_eps <- sigma2_eps * Diagonal(nobs * T)
    Q_eta <- Sigma_eta %>% solve()
    Q_eps <- Sigma_eps %>% solve()

    ## Simulate process
    bbox <- matrix(c(0,0,1,1),2,2)
    s <- construct_grid(bbox, 41)
    Mfun <- construct_M(Y_basis, s)
    K_basis <- list(G_const, G_const, G_const, G_const)
    M <- Mfun(K_basis, list(150, 0.002, -0.1, 0.1))
    alpha <- matrix(0,nbasis(Y_basis),T)
    alpha[65,1] <- 1
    PHI <- eval_basis(Y_basis, s$s_grid_mat)
    s$s_grid_df$Y0 <- (PHI %*% alpha[,1]) %>% as.numeric()
    for(i in 1:(T-1)) {
        alpha[,i+1] <- (M %*% alpha[,i]) %>% as.numeric() + sqrt(sigma2_eta)*rnorm(nbasis(Y_basis))
        s$s_grid_df[paste0("Y",i)] <- (PHI %*% alpha[,i+1]) %>% as.numeric()
    }
    s_long <- gather(s$s_grid_df, time, val, -s1, -s2)
    g_truth <- ggplot(s_long) + geom_tile(aes(s1,s2,fill=val)) + facet_wrap(~time) +
        scale_fill_distiller(palette="Spectral", limits = c(-0.1,1.1)) +
        coord_fixed(xlim=c(0,1), ylim = c(0,1))

    ## Observe process
    PHI_obs_1 <- eval_basis(Y_basis, zlocs[,1:2] %>% as.matrix())
    PHI_obs <- do.call("bdiag", lapply(1:T, function(x) PHI_obs_1))
    X_obs <-  cbind(1, do.call("rbind", lapply(1:T, function(x) zlocs))) %>% as.matrix()
    Z <- X_obs %*% beta + PHI_obs %*% c(alpha) + sqrt(sigma2_eps) * rnorm(nrow(PHI_obs))
    z_df <- data.frame(expand.grid.df(zlocs, data.frame(t = 0:(T-1))))
    z_df$z <- Z %>% as.numeric()
    g_obs <- ggplot(z_df) + geom_point(aes(s1, s2, colour = z)) +
        facet_wrap(~t) +
        scale_colour_distiller(palette = "Spectral") +
        coord_fixed(xlim=c(0,1), ylim = c(0,1))

    z_STIDF <- STIDF(sp = SpatialPoints(z_df[,c("s1","s2")]),
                     time = as.Date(z_df$t, origin = "2017-12-01"),
                     data = z_df)

    list(s_df = s_long,
         z_df = z_df,
         z_STIDF = z_STIDF,
         g_truth = g_truth,
         g_obs = g_obs)
}

construct_grid <- function(bbox, ngrid, coordnames = NULL) {
    ndim <- nrow(bbox)
    if(length(ngrid) == 1)
        ngrid <- rep(ngrid, ndim)
    if(!(length(ngrid) == ndim) | !is.numeric(ngrid))
        stop("ngrid needs to be a numeric (which will be rounded if not an integer)
             with length one or equal to the number of columns in bbox")
    ngrid <- round(ngrid)

    s <- lapply(1:nrow(bbox), function(i)
        seq(bbox[i,1], bbox[i,2], length.out = ngrid[i]))
    s_grid <- do.call("expand.grid", s)
    if(is.null(coordnames)) {
        names(s_grid) <- paste0("s",1:ndim)
    } else {
        names(s_grid) <- coordnames
    }
    list(s_grid_df = s_grid,
         s_grid_mat = s_grid %>% as.matrix(),
         ds = (bbox[,2] - bbox[,1])/ngrid,
         area = prod((bbox[,2] - bbox[,1])/ngrid))
}

construct_M <- function(Y_basis, s) {
    PHI <- eval_basis(Y_basis, s$s_grid_mat)
    GRAM <- crossprod(PHI)*s$area
    GRAM_inv <- solve(GRAM)
    ndim <- dimensions(Y_basis)
    function(K_basis, ki) {
        K <- construct_kernel(K_basis, ki)
        Kmat <- K(s$s_grid_mat, s$s_grid_mat)
        M <- GRAM_inv %*% crossprod(t(Kmat) %*% PHI, PHI)*s$area^2
    }
}

construct_Q <- function (Q_eta, M, T)
{
    n <- nrow(Q_eta)
    QM <- Q_eta %*% M
    MtQM <- crossprod(chol(Q_eta) %*% M) + Q_eta
    MtQ <- t(QM)
    for (i in 0:(T - 3)) {
        if (i == 0) {
            Q <- cBind(-QM, MtQM, -MtQ, Zeromat(n, ((T - 3) - i) *
                                                    n))
        }
        else if (i == (T - 3)) {
            Q <- rBind(Q, cBind(Zeromat(n, n * i), -QM, MtQM,
                                -MtQ))
        }
        else {
            Q <- rBind(Q, cBind(Zeromat(n, n * i), -QM, MtQM,
                                -MtQ, Zeromat(n, ((T - 3) - i) * n)))
        }
    }
    Q <- rBind(cBind(Q_eta  + MtQM, -MtQ, Zeromat(n, (T - 2) * n)),
               Q,
               cBind(Zeromat(n, n *(T - 2)), -QM, Q_eta))
    return(Q)
}


Zeromat <- function (ni, nj = NULL)
{
    if (is.null(nj))
        nj <- ni
    return(as(sparseMatrix(i = {
    }, j = {
    }, dims = c(ni, nj)), "dgCMatrix"))
}

rep.col<-function(x,n){
    l <- lapply(1:n, function(i) x)
    y <- do.call("c", l)
    matrix(y, ncol = n, byrow = FALSE)
}

construct_kernel <- function(Basis, ki) {
    if(!is.list(Basis)) stop("Basis needs to be of class list")
    if(!all(sapply(Basis, function(x) is(x,"Basis"))))
        stop("All Basis functions need to be of class Basis")
    ndim <- dimensions(Basis[[1]])

    function(s, r) {
        if(1) { ## Actual basis
            theta_s <- list()
            for(i in 1:(2 + ndim)) {
                theta_s[[i]] <- (eval_basis(Basis[[i]], s) %*% ki[[i]]) %>%
                    as.numeric() %>%
                    rep.col(nrow(r))
            }
            theta_s_1 <- lapply(theta_s, function(x) x[,1])
            D <- fields::rdist(s + do.call("cbind", theta_s_1[3:(2 + ndim)]), r)
            theta_s[[1]] * exp(-D^2/theta_s[[2]])
        } else {
            D <- fields::rdist(t(t(s) + c(ki[[3]], ki[[4]])), r)
            ki[[1]] * exp(-D^2/ki[[2]])
        }
    }
}

logdet <- function (L)
{
    diagL <- diag(L)
    return(2 * sum(log(diagL)))
}

#' @title Densify with explicit zeroes
#'
#' @description This function takes two sparse matrices and returns the first matrix padded with explicit zeros so that it is at least dense (probably denser) than the second matrix. This function only works with matrices of class Matrix
#' #'
#' @param A object of class Matrix
#' @param B object of class Matrix
#' @return object of class Matrix
#' @export
#' @examples
#' require(Matrix)
#' Q1 <- sparseMatrix(i=c(1,2,2),j=c(1,1,2),x=c(0.1,0.2,1))
#' Q2 <- sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.3,0.2,1))
#' Q1dens <- densify(Q1,Q2)
#' Q1
#' Q1dens
densify <- function(A,B) {
    ## Makes A at least as dense as B
    As <- symb(A)
    Bs <- symb(B)
    delta <- as(As - Bs,"dgTMatrix")
    idx <- which(delta@x == -1)
    addon <- sparseMatrix(delta@i+1,delta@j+1,x=0)
    A <- A + addon
    A
}

#' @title Return the symbolic representation of a Matrix
#'
#' @description This function takes an object of class Matrix and returns the same Matrix with all elements replaced with 1
#' #'
#' @param A object of class Matrix
#' @return object of class Matrix
#' @export
#' @examples
#' require(Matrix)
#' Q <- sparseMatrix(i=c(1,2,2),j=c(1,1,2),x=c(0.1,0.2,1))
#' Qsymb <- symb(Q)
#' Qsymb
symb <- function(A) {
    A@x <- rep(1,length(A@x))
    A
}


expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
