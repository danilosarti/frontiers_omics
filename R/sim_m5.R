#' sim_m5
#'
#' Simulates from model 5 presented in
#'
#' @param pars
#'
#' @return
#' @export
#'
#' @examples
sim_m5=function(pars=list(I=2,J=3,K=5,w_a=w_a,sigma_a=1,u_mean=0,sigma_mean=1,u_env=0,sigma_env=1, u_gen=0,sigma_env=1,u_epsilon=0,sigma_epsilon=1, x_env=x_env,sigma_w=1,u_w=0,sigma_ge=1,u_ge=0,sigma_aw=1)){
  # setting the dimensions
  I=pars$I
  J=pars$J
  env_labels=sprintf("environment_%s",seq(1:J))
  gen_label=sprintf("genotype_%s",seq(1:I))
  df_m5=expand.grid(environment_label=env_labels,genotype_label=gen_label)
  ##creating incidence matrix


  inc=model_matrix(~.-1, data=df_m5,  contrasts.arg =  lapply(df_m5[, sapply(df_m5, is.factor), drop = FALSE],
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
  sigma_ag=pars$sigma_a
  ## obtaining the interactions
  # we start by obtaining the varcovarstructure
  g_a_part=as.matrix(z_a)%*%G_a%*%t(as.matrix(z_a))
  z_e_part=(as.matrix(z_e)%*%t(as.matrix(z_e)))*pars$sigma_aw
  # make the hadarmar product

  var_covar_gw=hadamard.prod(g_a_part, z_e_part)
  ## creating the matrix ke of env covariates:
  #k_e=(as.matrix(pars$x_env)%*%t(as.matrix(pars$x_env)))/(sum(diag((as.matrix(pars$x_env)))/nrow(pars$x_env)))
  k_e=tcrossprod(pars$x_env)/ncol(pars$x_env)
  train_test=c("train","test")
  # we need to produce a training and a test dataset.
  # first we create a list to store the train and
  # test dataframes.
  df_train_test=vector("list",2)
  # now we loop through the simulations to produce a # training and a testing dataset
  for (i in 1:length(df_train_test)){
    #generate the interaction effects.
    #library(MASS)
    gw=mvrnorm(n=1,mu = rep(pars$u_ge,I*J),
               Sigma = var_covar_gw)
    #generate the genetic effets
    g_i <- mvrnorm(n = 1,
                   mu = rep(pars$u_gen,I),
                   Sigma = G_a*pars$sigma_a)
    #### now we simulate the environmental effect
    e_j=rnorm(J,mean=pars$u_env,sd=pars$sigma_env)
    ### now we simulate the errors
    epsilon_ij=rnorm(n=I*J,mean=pars$u_epsilon, sd=pars$sigma_epsilon)
    ## generating w
    w_k <- mvrnorm(n = 1,
                   mu = rep(pars$u_w,pars$J),
                   Sigma = k_e*pars$sigma_w)
    ## we need now to generate the interaction terms which are
    ##$M N\left(0,\left[Z_{a} G_{a} Z_{a}^{\prime}\right] \odot\left[Z_{E} ## Z_{E}^{\prime}\right] \sigma_{a e}^{2}\right)$

    df_m5$gen_effect=rep(g_i,each=J)
    df_m5$env_effect=rep(e_j,I)
    df_m5$epsilon=epsilon_ij
    df_m5$mu=rep(rnorm(1),I*J)
    df_m5$w_k=w_k
    df_m5$gw=gw
    df_m5$y=(df_m5$gen_effect+df_m5$env_effect+df_m5$epsilon+df_m5$mu+df_m5$w_k+df_m5$gw)
    df_train_test[[i]]<-df_m5
  }

  names(df_train_test)<-c("train","test")
  return(list(pars=pars,
              test_df=df_train_test$test,
              train_df=df_train_test$train))

}
