## strat_score.R

## Assume input datafile with following format:
##   Column 1:     Disease Indicator
##   Column 2:     Genotype at test locus
##   Column 3-12:  First 10 principal components

phenofile <- "input.dat"
pheno <- matrix(scan(phenofile), ncol=12, byrow=T)

## dis: vector of disease outcomes (1=case, 0=control)
## g: test-locus genotype (0,1,2)
dis  <- pheno[,1]
g    <- pheno[,2]
pc   <- pheno[,-c(1:2)]
Ntot <- length(dis)


## form stratification score based on 10 principal components:
sscore_pc <- glm(dis~pc[,1:10], family=binomial())

## assign subjects to 5 strata based on predicted odds:
sslp_pc    <- sscore_pc$linear.predictors
qsslp_pc   <- quantile(sslp_pc, probs=c(0.2,0.4,0.6,0.8), names=T)
stratid_pc <- array(dim=Ntot)

for(j in 1:Ntot){
    if(sslp_pc[j] <= qsslp_pc[1]){
		stratid_pc[j] <- 1
    }

    if ((sslp_pc[j] > qsslp_pc[1]) && (sslp_pc[j] <= qsslp_pc[2])){
		stratid_pc[j] <- 2
    }

    if ((sslp_pc[j] > qsslp_pc[2]) && (sslp_pc[j] <= qsslp_pc[3])){
		stratid_pc[j] <- 3
    }

    if ((sslp_pc[j] > qsslp_pc[3]) && (sslp_pc[j] <= qsslp_pc[4])){
		stratid_pc[j] <- 4
    }

    if (sslp_pc[j] > qsslp_pc[4]){
		stratid_pc[j] <- 5
    }
}


## construct indicator functions for each stratum (stratum 5 is baseline category):
nstrat_ss <- 5
stind_pc <- matrix(nrow=Ntot, ncol=(nstrat_ss - 1))
for(k in 1:(nstrat_ss - 1)){
	stind_pc[,k] <- ifelse(stratid_pc==k,1,0)
}


## construct stratified logistic-regression model testing for association
## between disease and test genotype adjusting for the stratification-score
## indicators:
adjust_ss_pc <- summary.glm(glm(dis ~ g + stind_pc[, 1:(nstrat_ss - 1)], family=binomial()))
p_ss_pc <- adjust_ss_pc$coefficients[2,4]
