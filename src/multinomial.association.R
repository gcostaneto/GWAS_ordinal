multinomial.association <- function(home.path = NULL,  # home directory. If null, getwd()
                                    output.path = NULL,  # output directory. If null, getwd()
                                    pheno.trait = NULL, # data.frame for taxa x traits
                                    geno.marker = NULL, # data.frame for taxa x markers
                                    covariate = NULL, # data.frame for taxa x covariates
                                    #   partial = 'full_model' # partial = c('full_model', 'envo_model','geno_model')
                                    n.core = NULL, #  numeric value denoting the numbe of cores. 
                                    output_name = NULL, # charcter denoting the name of the model under test
                                    parallel = TRUE, # if is TRUE, and n.core is null, then use detectCores() - 1
                                    verbose = FALSE) #,save.meta=TRUE)
{
  message(paste0('Start at................', format(Sys.time(), "%a %b %d %X %Y")))
  message(paste0('1. Data Check ..........', format(Sys.time(), "%a %b %d %X %Y")))
  
  if(is.null(home.path)) home.path = getwd()
  if(is.null(   output.path))    output.path = getwd()
  if(is.null( output_name)) output_name <-'multinomial'
  
  source('https://raw.githubusercontent.com/gcostaneto/envFeatures/main/src/MakeMyLogRegEq.R')
  
  if (!requireNamespace("ordinal", quietly = TRUE)) {
    utils::install.packages("ordinal")
    require("ordinal")
  }
  
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    utils::install.packages("doParallel")
    require("doParallel")
  }
  
  if (!requireNamespace("parallel", quietly = TRUE)) {
    utils::install.packages("parallel")
    require("parallel")
  }
  
  if (!requireNamespace("foreach", quietly = TRUE)) {
    utils::install.packages("foreach")
    require("foreach")
  }
  
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    utils::install.packages("reshape2")
    require("reshape2")
  }
  
  if(is.null(pheno.trait)) stop('Missing the Phenotypic Traits')
  if(is.null(geno.marker)) stop('Missing the Genomic Features (markers)')
  
  SNP_Q_model=Q_model=SNP_model = '1+SNP'
  
  if(!is.null(covariate))
  {
    SNP_Q_model = c(SNP_Q_model,colnames(covariate))
    Q_model = colnames(covariate)
  }
  
  if(is.null(covariate)) covariate = data.frame(matrix(1,nrow = nrow(geno.marker),ncol = 1))
  
  #'-------------------------------------------------------------------------------------------------
  # Building the Models     #####
  #'-------------------------------------------------------------------------------------------------
  
  
  
  
  form_null = MakeMyLogRegEq(Y =  'trait',X_tested =    '1' ) 
  form_full = MakeMyLogRegEq(Y =  'trait',X_tested =    SNP_Q_model)
  form_geno = MakeMyLogRegEq(Y =  'trait',X_tested =    SNP_model ) 
  form_cova = MakeMyLogRegEq(Y =  'trait',X_tested =    Q_model ) 
  
  
  n.geno.features  = ncol(geno.marker )
  
  
  message(paste0('2. makeCluster..........', format(Sys.time(), "%a %b %d %X %Y")))
  
  if(isTRUE(parallel))
  {
    if(is.null( n.core))
    {
      if(isTRUE(parallel))
      {
        cl <- parallel::makeCluster(parallel::detectCores() - 1)
        doParallel::registerDoParallel(cl)
      }
    }
    if(!is.null(n.core))
    {
      cl <-  parallel::makeCluster(n.core)
      doParallel::registerDoParallel(cl)
    }
    message(paste0('   threads = ',n.core))
  }
  if(isFALSE(parallel))
  {
    message(paste0('   parallel = FALSE',n.core))
  }
  message(paste0('Done!...................', format(Sys.time(), "%a %b %d %X %Y")))
  
  #'-------------------------------------------------------------------------------------------------
  # Running the models     #####
  #'-------------------------------------------------------------------------------------------------
  
  start = format(Sys.time(), "%a %b %d %X %Y")
  message(paste0('3. Running MLR .........',start))
  
  output_GWAS =
    try(
      foreach::foreach(geno.feature.j = 1:n.geno.features , 
                       .packages = c('plyr','reshape2','data.table','tidyverse','ordinal'),.combine = "rbind") %dopar%
        {
          
          #'---------------------------------------------------------
          # Model Step 1: organizing the models and inputs
          #'---------------------------------------------------------
          
          df_model = 
            cbind(
              data.frame(
                SNP = geno.marker[,geno.feature.j] ,
                trait = as.factor(pheno.trait[,1])),
              covariate)
          
          df_model = data.frame(   df_model)
          
          # running ordinal multinomial log regression
          
          null_model      <- try(ordinal::clm(formula = form_null,data = df_model)) # M0: null model (~1 model)
          full_model      <- try(ordinal::clm(formula = form_full,data = df_model)) # M1: markers + covariates (~1+SNP+Q model)
          geno_model      <- try(ordinal::clm(formula = form_geno,data = df_model)) # M2 : markers (~1 + SNP model)
          cova_model      <- try(ordinal::clm(formula = form_cova,data = df_model)) # M3 : markers (~1 + covariates)
          
          if(!is.null(full_model))
          {
            
            #             table_PA =     table(df_model$gFeature)
            #'---------------------------------------------------------
            # Model Step 2: LRT  test and McFadden's pseudo R-sq
            #'---------------------------------------------------------
            
            H1 = as.data.frame(anova(   null_model,  geno_model)) # testing SNP effect with no covariates
            H2 = as.data.frame(anova(   cova_model,  full_model)) # testing SNP effect over the effect of covariates
            
            
            # general model statistics
            output.stat = data.frame(GenoFeature    = colnames(geno.marker)[geno.feature.j],
                                     
                                     # Model AIC
                                     AIC_null      =   NA,
                                     AIC_SNP       =   NA,
                                     AIC_SNP_Q     =   NA,
                                     
                                     # Model p-value
                                     chisqp_H1 = NA,
                                     chisqp_H2 = NA,
                                     
                                     # Model p-value
                                     log_chisqp_H1 = NA,
                                     log_chisqp_H2 = NA,
                                     
                                     # McFadden's pseudo R-sq
                                     McF.pR2_H1   = NA,
                                     McF.pR2_H2   = NA)
            
            
            # Model AIC
            output.stat$AIC_null   =   H1$AIC[1]
            output.stat$AIC_SNP    =   H1$AIC[2]
            output.stat$AIC_SNP_Q  =   H2$AIC[2]
            
            # Model p-value
            
            output.stat$chisqp_H1 = H1$`Pr(>Chisq)`[2]
            output.stat$chisqp_H2 = H2$`Pr(>Chisq)`[2]
            
            output.stat$log_chisqp_H1 = -log10(H1$`Pr(>Chisq)`[2])
            output.stat$log_chisqp_H2 = -log10(H2$`Pr(>Chisq)`[2])
            
            # McFadden's pseudo R-sq
            output.stat$ McF.pR2_H1   = try(100*round(1 - geno_model$logLik/null_model$logLik,4))
            output.stat$ McF.pR2_H2   = try(100*round(1 - full_model$logLik/cova_model$logLik,4))
            
            #'---------------------------------------------------------
            #' Model Step 3: coef and partial p-values for each predictor
            #'---------------------------------------------------------
            
            output.partial       = as.data.frame(summary(full_model)$coefficients)
            output.partial       = output.partial[which(rownames( output.partial) %in% 'SNP'),]
            
            output.partial <-  data.frame(Feature = rownames( output.partial), output.partial[,c(1:2,4)] )
            names( output.partial)[2:4] = c('effect','se_effect','p_value')
            output.partial$log10 = -log10( output.partial$p_value)
            
            
            # maximum log10p
            nulllog =  any(output.partial$log10[which( output.partial$log10 == Inf)])
            if(isTRUE( nulllog ))  
            {
              max_value = max(output.partial$log10[-which( output.partial$log10 == Inf)],na.rm=T)
              output.partial$log10[which( output.partial$log10 == Inf)] =   max_value 
            }
            
            output.partial <- 
              output.partial %>% 
              reshape2::melt() %>%
              mutate(GenoFeature    = colnames(geno.marker)[geno.feature.j]) %>% 
              mutate(variable=paste0('partial_',variable)) %>% 
              reshape2::dcast(GenoFeature~variable+Feature,value.var = 'value')
            
            #'---------------------------------------------------------
            #' Model Step 4: Exporting ;)
            #'---------------------------------------------------------
            
            if(isTRUE(save))
            {
              data.table::fwrite(       output.partial, paste0(output.path,'/model_partial',output_name ,".csv"),append = T)
              data.table::fwrite(       output.stat, paste0(output.path,'/model_statistics',output_name ,".csv"),append = T)
            }
            
            
            output <-
              output.stat %>% 
              merge(output.partial ,by='GenoFeature')
          }
          
          
          
          
          return(output)
        }
    )
  
  if(isTRUE(parallel)) {parallel::stopCluster(cl) }
  
  end = Sys.time()
  message(paste0('Done!...................', format(Sys.time(), "%a %b %d %X %Y")))
  message(paste0('End at..................', format(Sys.time(), "%a %b %d %X %Y")))
  
  if(isTRUE(save))
  {
    message(paste0('saving files at......',output.path,'\n'))
    data.table::fwrite(   output_GWAS, paste0(output.path,'/Full_Results_',output_name ,".csv"))
  }
  
  return(   output_GWAS)
  
}