#' sim_m1
#'
#' A function to simulate from Model one proposed in
#'
#' @param pars a list containing the parameter required to simulate from model 1
#'
#' @return
#' @export
#'
#' @examples

sim_m1<-function(pars=list(w_a=w_a,u_mean=0,sigma_mean=1, sigma_ag=1,u_env=0,sigma_env=1, u_epsilon=0, sigma_epsilon=1,J=3,I=2,u_gen=1)){
  #declaring the dimensions
  I=pars$I
  J=pars$J
  #generate matrix Ga which will be used to
  G_a=tcrossprod(pars$w_a/ncol(pars$w_a))
  dim(G_a)
  rownames(G_a)=sprintf("genotype_%s",seq(1:I))
  colnames(G_a)=sprintf("genotype_%s",seq(1:I))
  sigma_ag=pars$sigma_ag
  library(MASS)
  train_test=c("train","test")
  # we need to produce a training and a test dataset.
  # first we create a list to store the train and
  # test dataframes.
  df_train_test=vector("list",2)
  # now we loop through the simulations to produce a # training and a testing dataset
  for (i in 1:length(df_train_test)){
    #generate the genetic effets
    g_i <- mvrnorm(n = 1,
                   mu = rep(pars$u_gen,I),
                   Sigma = G_a*sigma_ag)
    #### now we simulate the environmental effect
    e_j=rnorm(J,mean=pars$u_env,sd=pars$sigma_env)
    ### now we simulate the errors
    epsilon_ij=rnorm(n=I*J,mean=pars$u_epsilon, sd=pars$sigma_epsilon)
    env_labels=sprintf("environment_%s",seq(1:J))
    gen_label=sprintf("genotype_%s",seq(1:I))
    df_m1=expand.grid( environment_label=env_labels,genotype_label=gen_label)
    df_m1$gen_effect=rep(g_i,each=J)
    df_m1$env_effect=rep(e_j,I)
    df_m1$epsilon=epsilon_ij
    df_m1$mu=rep(rnorm(1,mean=pars$u_mean,sd=pars$sigma_mean),I*J)
    df_m1$y=(df_m1$gen_effect+df_m1$env_effect+df_m1$epsilon+df_m1$mu)
    df_train_test[[i]]<-df_m1
  }
  names(df_train_test)<-c("train","test")
  return(list(pars=pars,
              K_G=G_a,
              test_df=df_train_test$test,
              train_df=df_train_test$train))
}

