#' sim_blues
#'
#'A function to simulate blues in a block design context.
#'
#' @param pars a list containing the number of genotypes I, environments J, Z, marginal mean of treatmens
#' @import stats
#' @return
#' @export
#'
#' @examples
sim_blue_model=function(pars = list(I = 4,J=3,Z=3,R=2,u_treat=0,sigma_treat=1,u_gen=0,
                                    sigma_gen=1,u_block=0,sigma_block=1,u_pb=0,sigma_pb=1, u_gp=0, sigma_gp=1,u_mean=0,sigma_mean=1,u_epsilon=0, sigma_epsilon=1
)){
  # defining the dimensions
  I=pars$I
  Z=pars$Z
  R=pars$R
  J=pars$J
  train_test=c("train","test")
  # we need to produce a training and a test dataset.
  # first we create a list to store the train and
  # test dataframes.
  df_train_test=vector("list",2)
  # now we loop through the simulations to produce a # training and a testing dataset
  for (i in 1:length(df_train_test)){
    #generating effects of treatment
    treatment=rnorm(n=Z,mean=pars$u_treat,sd=pars$sigma_treat)
    treat_label <- sprintf("treatment_%s",seq(1:Z))
    #effects of genotype
    gen=rnorm(n=I,mean=pars$u_gen,sd=pars$sigma_gen)
    gen_label=sprintf("genotype_%s",seq(1:I))
    #effects of block
    block=rnorm(n=R,mean=pars$u_block,sd=pars$sigma_block)
    block_label=sprintf("block_%s",seq(1:R))
    # nested effect between treatment and blocks
    pb=c()
    for(z in 1:R){
      pbzr=rnorm(n=Z,mean=pars$u_pb,sd=pars$sigma_pb)
      pb=c(pb,pbzr)
    }
    # organize the nestings:
    df_pb=expand.grid(treat_lable=treat_label,block_lable=block_label)
    df_pb$pb=pb

    # nested effect between genotype and treatment.
    gp=c()
    for(w in 1:Z){
      gpiz=rnorm(n=I,mean=pars$u_gp, sd=pars$sigma_gp)
      gp=c(gp,gpiz)
    }
    #arranging the nestings
    df_gp=expand.grid(gen_lable=gen_label,treat_lable=treat_label)
    df_gp$gp=gp

    #overall mean \mu
    mu=rnorm(n=1,mean=pars$u_mean, sd=pars$sigma_mean)
    #error
    epsilon=rnorm(n=I*Z*R, mean=pars$u_epsilon,sd=pars$sigma_epsilon)
    df=expand.grid( treatment_label=treat_label,genotype_label=gen_label,block_label=block_label)

    #we must now populate the data set with the effects.
    df$treatment_effect=rep(treatment,I*R)
    df$gen_effect=rep(rep(gen,each=Z),R)
    df$block_effect=rep(block,each=I*Z)
    df$epsilon=epsilon
    df$mu=mu

    ##we need to insert the nesting into df.
    df$pb=df_pb$pb
    df$gp=df_gp$gp
    ## generating y's
    df$y=df$treatment_effect+df$gen_effect+df$block_effect+df$mu+df$epsilon+df$pb+df$gp
    df_train_test[[i]]<-df
  }
  names(df_train_test)<-c("train","test")
  return(list(pars=pars,
              test_df=df_train_test$test,
              train_df=df_train_test$train))
}
