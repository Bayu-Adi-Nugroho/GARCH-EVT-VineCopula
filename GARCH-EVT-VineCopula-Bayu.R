
#My interpretation to construct GARCH - EVT - Vine Copula
#Please cite my paper if you use these codes:
#Mappadang, A., Nugroho, B. A., Lestari, S. D., Elizabeth, & Lestari, T. K. (2024). 
#Measuring value-at-risk and expected shortfall of newer cryptocurrencies: new insights. 
#Cogent Business & Management, 11(1). 
#https://doi.org/10.1080/23311975.2024.2416096

#The first step is to select the appropriate GARCH (thanks to Dr. David Gabauer):

library(rugarch)
library(xts)
# Load or prepare your time series data (e.g., 'returns_data')
# Example: 

data("EuStockMarkets")
index_data = tsbox::ts_xts(EuStockMarkets)

DAX <- index_data$DAX
FTSE <- index_data$FTSE
SMI <- index_data$SMI
all_index <- cbind(DAX,FTSE,SMI)

# Create stationary dataset
Z1 = diff(log(all_index[, c("DAX", "FTSE", "SMI")]))
Z1 = na.omit(Z1)

# Step 1 : Selecting the appropriate ARMA 

get_pq = function(data,max_p=5,max_q=5){
  final.aic = 0
  final.order = c(0,0,0)
  for(p in 0:max_p) for(q in 0:max_q){
    if(p == 0 && q == 0){
      next
    }
    arimaFit = tryCatch(arima(data,order = c(p,0,q)),
                        error = function(err) FALSE,
                        warning = function(err) FALSE)
    
    if(!is.logical(arimaFit)){
      
      current.aic <- AIC(arimaFit)
      if(current.aic < final.aic){
        
        final.aic <- current.aic
        final.order <- c(p,0,q)
        final.arima <- arima(data,order = final.order)
      }
    }else{
      next
    }
  }
  return(final.order)
}

optimal_lag= apply(Z1,MARGIN = 2,FUN = get_pq)
## function for converting matrics to list
mat2list = function(data){
  return(as.list(as.data.frame(Z1,stringAsFactors = FALSE)))
}

# Step 2 : Selecting the best-fitting GARCH model using the ARMA from step 1

source("GARCH selection.R")

GARCHselection = function(x, distributions=c("norm","snorm","std","sstd","ged","sged"), models=c("sGARCH","eGARCH","gjrGARCH","iGARCH"), prob=0.05, conf.level=0.90, lag=20, ar=0, ma=0) {
  message("A dynamic version of the optimal univariate GARCH selection procedure is implemented according to:\n Antonakakis, N., Chatziantoniou, I., & Gabauer, D. (2021). The impact of Euro through time: Exchange rate dynamics under different regimes. International Journal of Finance & Economics, 26(1), 1375-1408.")
  if (!is(x, "zoo")) {
    stop("Data needs to be of type 'zoo'")
  }
  GARCH_IC = matrix(Inf, nrow=length(distributions), ncol=length(models))
  colnames(GARCH_IC) = models
  rownames(GARCH_IC) = distributions
  spec_list = list()
  table_list = list()
  for (i in 1:length(models)) {
    spec_list[[i]] = list()
    table_list[[i]] = list()
    for (j in 1:length(distributions)) {
      spec_list[[i]][[j]] = list()
    }
    names(spec_list[[i]]) = distributions
  }
  names(spec_list) = names(table_list) = models
  
  for (j in 1:length(models)) {
    message(paste0("-",models[j]))
    for (i in 1:length(distributions)) {
      message(paste0("--",distributions[i]))
      {
        ugarch.spec = rugarch::ugarchspec(mean.model=list(armaOrder=c(ar,ma)),
                                          variance.model=list(model=models[j], garchOrder=c(1,1)),
                                          distribution.model=distributions[i])
      }
      ugarch.fit = rugarch::ugarchfit(ugarch.spec, data=x, solver="hybrid", solver.list=list(outer.iter=10, inner.iter=1000, eval.se=FALSE, tol=1e-12))
      if (ugarch.fit@fit$convergence==0) {
        fit = GARCHtests(ugarch.fit, prob=prob, conf.level=conf.level, lag=lag)
        GARCH_IC[i,j] = fit$InformationCriterion
        spec_list[[models[j]]][[distributions[i]]] = ugarch.spec
        table_list[[j]][[distributions[i]]] = fit
      }
    }
  }
  GARCH_selection = which(GARCH_IC==min(GARCH_IC),arr.ind=TRUE)
  best_ugarch = spec_list[[GARCH_selection[2]]][[GARCH_selection[1]]]
  best_table = table_list[[GARCH_selection[2]]][[GARCH_selection[1]]]
  return = list(best_ugarch=best_ugarch, best_table=best_table, GARCH_IC=GARCH_IC, spec_list=spec_list, table_list=table_list)
}

#for example DAX
Best_GARCH= GARCHselection(Z1$DAX, ar=5, ma=4) #See optimal_lag for ar and ma
Best_GARCH[["best_ugarch"]]


# After getting the appropriate ARMA and GARCH, we compute the residuals
# Assume that we have armaOrder is c(5,4), eGARCH, and std 
spec <- ugarchspec(
  mean.model = list(armaOrder = c(5, 4), include.mean = TRUE), 
  #the armaOrder is from STEP TWO, lets assume we have c(0,0)
  variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
  ##the model = "gjrGARCH" is from STEP ONE, lets assume we have gjrGARCH and std
  distribution.model = "std" # Use Student's t-distribution
)

# Fit the model
garch_fit1 <- ugarchfit(spec = spec, data = Z1$DAX)
garch_fit2 <- ugarchfit(spec = spec, data = Z1$FTSE)
garch_fit3 <- ugarchfit(spec = spec, data = Z1$SMI)

#Note that my paper used ugarchroll since it focuses on out-of-sample data.

# Extract standardized residuals (zt)
# These should be approximately i.i.d., with mean 0 and variance 1
standardized_residuals1 <- residuals(garch_fit1, standardize = TRUE)
standardized_residuals2 <- residuals(garch_fit2, standardize = TRUE)
standardized_residuals3 <- residuals(garch_fit3, standardize = TRUE)

#The third step is to fit EVT to the tails of the residuals:

#To transform the results of a Generalized Pareto Distribution (gpd.fit) 
#to a copula in R, you must convert the GPD-fitted tail exceedances 
#into uniform [0, 1] marginal distributions (or "pseudo-observations") and then fit a copula model to these transformed values. 

#This process is typically done in three stages: 
#fitting the tails, transforming to uniform, and fitting the copula. 

#fitting the tails
library(evd)
# Set thresholds (e.g., 90th percentile)
u1 <- quantile(standardized_residuals1, 0.9)
u2 <- quantile(standardized_residuals2, 0.9)
u3 <- quantile(standardized_residuals3, 0.9)

# Fit GPD
fit1 <- fpot(standardized_residuals1, threshold = u1)
fit2 <- fpot(standardized_residuals2, threshold = u2)
fit3 <- fpot(standardized_residuals3, threshold = u3)

#transforming to uniform
# Function to transform data using GPD fit
transform_to_unif <- function(data, fit, threshold) {
  # Empirical distribution function for data below threshold
  below_thresh <- data[data <= threshold]
  emp_cdf <- ecdf(below_thresh)
  
  # Calculate uniform values
  n <- length(data)
  u <- numeric(n)
  
  # Proportion of data below threshold
  p_below <- mean(data <= threshold)
  
  # Transform
  u[data <= threshold] <- emp_cdf(data[data <= threshold]) * p_below
  
  # Use GPD to transform upper tail
  # pgpd is from evd package
  u[data > threshold] <- p_below + (1 - p_below) * 
    pgpd(data[data > threshold], loc = threshold,
         scale = fit$estimate["scale"],
         shape = fit$estimate["shape"])
  return(u)
}

u1_trans <- transform_to_unif(standardized_residuals1, fit1, u1)
u2_trans <- transform_to_unif(standardized_residuals2, fit2, u2)
u3_trans <- transform_to_unif(standardized_residuals3, fit3, u3)
u_data <- cbind(u1_trans, u2_trans, u3_trans)

library(VineCopula)
library(copula)

Pobs_Matrix <- as.matrix(u_data)

Copula_Fit <- RVineStructureSelect(Pobs_Matrix, progress = TRUE)

summary(Copula_Fit)

#the best fit is C-Vine Copula

# selecting the best family

count_par = function(input_mat){
  num = 0
  for(i in 1:3){ #----------------------------> 2 if there are two assets
    for(j in 1:3){ #----------------------------> 2 if there are two assets
      if(input_mat[i,j] != 0){
        num = num + 1
      }  
    }
  }
  return(num)
}

cores = detectCores()
par(mfrow=c(1,1)) 

# all family
cores = detectCores()
par(mfrow=c(1,1)) 
gaussian = RVineStructureSelect(data = u_data,
                                familyset = 1,
                                progress = F,
                                type = 'CVine',
                                treecrit = 'tau',
                                cores = cores)
LLH_gaus = gaussian$logLik
AIC_gaus = gaussian$AIC
BIC_gaus = gaussian$BIC
par_count_gaus = count_par(gaussian$par)
res1 = c(par_count_gaus,LLH_gaus,AIC_gaus,BIC_gaus)



tstudent = RVineStructureSelect(data = u_data,
                                familyset = 2,
                                progress = F,
                                type = 'CVine',
                                treecrit = 'tau',
                                core = cores)
LLH_t = tstudent$logLik
AIC_t = tstudent$AIC
BIC_t = tstudent$BIC
par_count_t = count_par(tstudent$par)   #+ count_par(all_tstudent$par2)  ?
res2 = c(par_count_t,LLH_t,AIC_t,BIC_t)

clayton = RVineStructureSelect(data = u_data,
                               familyset = 3,
                               progress = F,
                               type = 'CVine',
                               treecrit = 'tau',
                               core = cores)
LLH_clayton = clayton$logLik
AIC_clayton = clayton$AIC
BIC_clayton = clayton$BIC
par_count_clayton = count_par(clayton$par)   #+ count_par(all_tstudent$par2)  ?
res3 = c(par_count_clayton,LLH_clayton,AIC_clayton,BIC_clayton)


gumbel = RVineStructureSelect(data = u_data,
                              familyset = 4,
                              progress = F,
                              type = 'CVine',
                              treecrit = 'tau',
                              core = cores)
LLH_gumbel = gumbel$logLik
AIC_gumbel = gumbel$AIC
BIC_gumbel = gumbel$BIC
par_count_gumbel = count_par(gumbel$par)   #+ count_par(all_tstudent$par2)  ?
res4 = c(par_count_gumbel,LLH_gumbel,AIC_gumbel,BIC_gumbel)


frank = RVineStructureSelect(data = u_data,
                             familyset = 5,
                             progress = F,
                             type = 'CVine',
                             treecrit = 'tau',
                             core = cores)
LLH_frank = frank$logLik
AIC_frank = frank$AIC
BIC_frank = frank$BIC
par_count_frank = count_par(frank$par)   #+ count_par(all_tstudent$par2)  ?
res5 = c(par_count_frank,LLH_frank,AIC_frank,BIC_frank)


joe = RVineStructureSelect(data = u_data,
                           familyset = 6,
                           progress = F,
                           type = 'CVine',
                           treecrit = 'tau',
                           core = cores)
LLH_joe = joe$logLik
AIC_joe = joe$AIC
BIC_joe = joe$BIC
par_count_joe = count_par(joe$par)   #+ count_par(all_tstudent$par2)  ?
res6 = c(par_count_joe,LLH_joe,AIC_joe,BIC_joe)

res_all = rbind(res1,res2,res3,res4,res5,res6)
rownames(res_all) = c('Gaussian','t-student','clayton','gumbel','frank','joe')
colnames(res_all) = c('Number of assets','LogLik','AIC','BIC')

# Joe copula is the best family

# next is to simulate the copula then transform it to returns
set.seed(2024)
PseudoScenarios <- data.frame(RVineSim(N = 10000, RVM = joe))

library(gamlss)
library(gamlss.dist)
library(gamlss.add)
fit1=fitDist(Z1[,"DAX"], k = 3, type = "realline", trace = FALSE, try.gamlss = TRUE)#####Identify the suitable distribution

SimDAX = qLO(PseudoScenarios[,1],mu=fit1$mu, sigma=fit1$sigma)####Simulated DAX

