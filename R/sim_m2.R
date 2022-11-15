#' sim_m2
#'
#' @param pars
#'
#' @return
#' @export
#' @import matrixcalc
#' @import  MASS
#' @import tidyverse
#' @examples
sim_m2=function(pars=list(I=2,J=3,w_a=w_a,sigma_a=1,u_mean=0,sigma_mean=1,u_env=0,sigma_env=1, u_gen=0,sigma_env=1,sigma_ge=2,u_ge=0,sigma_ge=1,u_epsilon=0,sigma_epsilon=1)){
  # setting the dimensions
  I=pars$I
  J=pars$J
  env_labels=sprintf("environment_%s",seq(1:J))
  gen_label=sprintf("genotype_%s",seq(1:I))
  df_m2=expand.grid(environment_label=env_labels,genotype_label=gen_label)
  ##creating incidence matrix

  library(modelr)
  inc=model_matrix(~.-1, data=df_m2,  contrasts.arg =  lapply(df_m2[, sapply(df_m2, is.factor), drop = FALSE],
                                                              contrasts, contrasts = FALSE))
  ## getting the incidence matrix of environments Z_E.
  z_e=inc %>%
    dplyr::select(contains("environment"))
  ## getting the incidence matrix of genotypes
  z_a=inc %>%
    dplyr::select(contains("genotype"))
  ##
  #generate matrix Ga which will be used to generate the genomic effects and interactions effects to
  G_a=tcrossprod(pars$w_a/ncol(pars$w_a))
  sigma_a=pars$sigma_a
  rownames(G_a)=sprintf("genotype_%s",seq(1:I))
  colnames(G_a)=sprintf("genotype_%s",seq(1:I))
  ## obtaining the interactions
  # we start by obtaining the varcovarstructure
  g_a_part=as.matrix(z_a)%*%G_a%*%t(as.matrix(z_a))
  z_e_part=(as.matrix(z_e)%*%t(as.matrix(z_e)))*pars$sigma_ge
  # make the hadarmar product

  var_covar_ge=hadamard.prod(g_a_part, z_e_part)

  train_test=c("train","test")
  # we need to produce a training and a test dataset.
  # first we create a list to store the train and
  # test dataframes.
  df_train_test=vector("list",2)
  # now we loop through the simulations to produce a # training and a testing dataset
  for (i in 1:length(df_train_test)){
    #generate the interaction effects.
    library(MASS)
    ge=mvrnorm(n=1,mu = rep(pars$u_ge,I*J),
               Sigma = var_covar_ge)
    #generate the genetic effets
    g_i <- mvrnorm(n = 1,
                   mu = rep(pars$u_gen,I),
                   Sigma = G_a*sigma_a)
    #### now we simulate the environmental effect
    e_j=rnorm(J,mean=pars$u_env,sd=pars$sigma_env)
    ### now we simulate the errors
    epsilon_ij=rnorm(n=I*J,mean=pars$u_epsilon, sd=pars$sigma_epsilon)

    ## we need now to generate the interaction terms which are
    ##$M N\left(0,\left[Z_{a} G_{a} Z_{a}^{\prime}\right] \odot\left[Z_{E} ## Z_{E}^{\prime}\right] \sigma_{a e}^{2}\right)$

    df_m2$gen_effect=rep(g_i,each=J)
    df_m2$env_effect=rep(e_j,I)
    df_m2$epsilon=epsilon_ij
    df_m2$mu=rep(rnorm(1),I*J)
    df_m2$ge=ge
    df_m2$y=(df_m2$gen_effect+df_m2$env_effect+df_m2$epsilon+df_m2$mu+df_m2$ge)
    df_train_test[[i]]<-df_m2
  }

  names(df_train_test)<-c("train","test")
  return(list(pars=pars,
              K_G=G_a,
              test_df=df_train_test$test,
              train_df=df_train_test$train))

}
