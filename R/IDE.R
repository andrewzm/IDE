# IDE: An R Software package for implementing integro-difference
# equation models
# Copyright (c) 2018 Andrew Zammit-Mangion
# Author: Andrew Zammit-Mangion, azm (at) uow.edu.au
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#'
#' Integro-difference equation
#'
#' The Integro-Difference Equation model is a linear, dynamical model used to model phenomena that evolve in space and in time. At the heart of the model is the kernel, which dictates how the process evolves from one time point to the next. Both process and parameter reduction are used to facilitate computation, and spatially-varying kernels are allowed. Data used to estimate the parameters are assumed to be readings of the process corrupted by Gaussian measurement error. Parameters are fitted by maximum likelihood,  and estimation is carried out using an evolution algorithm.
#' @name IDE-package
#' @docType package
#' @import methods
#' @import ggplot2
#' @import Matrix
#' @import sp
#' @import spacetime
#' @import parallel
#' @import dplyr
#' @import sparseinv
#' @importFrom FRK distR auto_basis eval_basis nbasis plane
#' @importFrom DEoptim DEoptim
#' @importFrom tidyr gather spread
#' @importFrom stats .getXlevels coefficients dist kmeans lm median model.extract model.frame model.matrix na.fail optim runif sd terms var time rnorm formula coef predict
#' @importFrom utils data
NULL
