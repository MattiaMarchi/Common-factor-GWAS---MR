# R code for the main model in the paper #
# The impact of poverty on mental illness: Emerging evidence of a causal relationship #
# R code by Mattia Marchi (mattiamarchimd@gmail.com) 
# May 18, 2023

####----------------------------------------------------------------------------------------------------------------------------###
###-------------------------------------------1. Load the required packages-----------------------------------------------------###
####----------------------------------------------------------------------------------------------------------------------------###
library(stringr); library(TwoSampleMR); library(MendelianRandomization); library(tidyverse); library(data.table); library(cause)
library(ggplot2); require(GenomicSEM); library(ggpubr); library(ieugwasr); library(MRPRESSO); library(qqman)

####----------------------------------------------------------------------------------------------------------------------------###
###-----------------------------------------2. Latent poverty-factor estimation-------------------------------------------------###
####----------------------------------------------------------------------------------------------------------------------------###
#This codes will produce multivariable GWAS sumstats for latent phenotype called POVERTY made of:
#Household Income (HI) + Occupational Income (OI) + Social Deprivation (SD)
#Ref to the method: https://github.com/GenomicSEM/GenomicSEM/wiki/4.-Common-Factor-GWAS; to the paper: 10.1038/s41562-019-0566-x

#Steps:
#1. Munge the sumstats
#2. Run multivariable LDSC
#3. Prepare the sumstats for GWAS
#4. Combine the sumstats and LDSC output and run the common factor GWAS

#-----Munge
hi <- "HI_forGenomicSEM.txt" #Change the file name as needed #N=379598; sample.prev=NA; pop.prev=NA (because its a continuous trait - if it were binary, look at wikihow on github for calculating Neff)
sd <- "SD_forGenomicSEM.txt" #N=440350; sample.prev=NA; pop.prev=NA
oi <- "OI_forGenomicSEM.txt" #N=282963; sample.prev=NA; pop.prev=NA
#create vector of the summary statistics files
files <- c(hi,sd,oi)
#define the reference file being used to allign alleles across summary stats
#here we are using hapmap3
hm3 <- "eur_w_ld_chr/w_hm3.snplist"
#name the traits 
trait.names <- c("HI", "SD", "OI")
#list the sample sizes. All have SNP-specific sum of effective sample sizes (Neff), SD is continuous so it's NA
N = c(379598, 440350, 282963)
#define the imputation quality filter
info.filter = 0.9
#define the MAF filter
maf.filter = 0.01
#run munge
munge(files = files, hm3 = hm3, trait.names = trait.names, N = N, info.filter = info.filter, maf.filter = maf.filter)

#-----Run LDSC
traits <- c("HI.sumstats.gz", "SD.sumstats.gz", "OI.sumstats.gz")
sample.prev <- c(NA, NA, NA)
population.prev <- c(NA, NA, NA)
ld <- "eur_w_ld_chr/"
wld <- "eur_w_ld_chr/"
hi_sd_oi_LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
hi_sd_oi_LDSCoutput$V #is the sampling covariance matrix in the format expected by lavaan
hi_sd_oi_LDSCoutput$S #is the covariance matrix (on the liability scale for case/control designs)
hi_sd_oi_LDSCoutput$m #number of SNPs used to construct the LD score
#Standard Errors
k <- nrow(hi_sd_oi_LDSCoutput$S)
SE <- matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <- sqrt(diag(hi_sd_oi_LDSCoutput$V))
#-----Prepare the sumstats for GWAS
files <- c(hi,sd,oi)
ref <- "reference.1000G.maf.0.005.txt"
trait.names <- c("HI", "SD", "OI")
se.logit <- c(F,F,F)
OLS <- c(T,T,T)
info.filter <- 0.6
maf.filter <- 0.01
poverty_sumstats <- sumstats(files=files, ref=ref, trait.names=trait.names, se.logit=se.logit, OLS=OLS, linprob=NULL,
                             N=NULL, betas=NULL, info.filter=info.filter, maf.filter=maf.filter, keep.indel=FALSE,
                             parallel=FALSE,cores=NULL)
#-----Run common factor GWAS
#run the multivariate GWAS using parallel processing - Note: running that on standalone PCs may exceed computation power, you can avoind that by running serial processing, but computation time largely increase.
#Run this code on HPC
poverty_factor <- commonfactorGWAS(covstruc = hi_sd_oi_LDSCoutput, SNPs = poverty_sumstats,
                                   estimation = "DWLS", cores = NULL, toler = FALSE, SNPSE = FALSE,
                                   parallel = TRUE, Output = NULL, GC="standard", MPI=FALSE)
#Manhattan plot
mn <- manhattan(poverty_factor2, chr="CHR", bp="BP", snp="SNP", p="Pval_Estimate", genomewideline = 8,
                logp = T) #, highlight = snpsOfInterest)
mn

####----------------------------------------------------------------------------------------------------------------------------###
###-------------------------------------------3. Mendelian Randomization--------------------------------------------------------###
####----------------------------------------------------------------------------------------------------------------------------###
####Be sure to convert all to Betas - see example script below
outcome  <- fread("your_outcome_GWAS_sumstat.txt")
outcome$BETA <- log(outcome$OR)
outcome$Z <- qnorm(1 - outcome$P/2) 
outcome$SE <- abs(outcome$BETA/outcome$Z)
write.table(outcome, file = "your_outcome_GWAS_sumstat.txt", row.names = F, col.names = T, quote = F)
###Below the example script for one Twosample MR analysis of "your_exposure" against "your_outcome", replace them as needed
#Import exposure data
exp <- "your_exposure_GWAS_sumstat.txt"
exposure_data <- read_exposure_data(exp,
                                    clump = FALSE,
                                    sep = " ",
                                    snp_col = "SNP",
                                    beta_col = "beta",
                                    se_col = "se",
                                    effect_allele_col = "effect_allele",
                                    other_allele_col = "other_allele",
                                    pval_col = "pval",
                                    ncase_col = "ncase",
                                    ncontrol_col = "ncontrol",
                                    samplesize_col = "samplesize",
                                    min_pval = 1e-200,
                                    log_pval = F)
exposure_data_small <- subset(exposure_data, exposure_data$pval.exposure<5e-8)
exposure_data_small <- clump_data(exposure_data_small, clump_kb = 10000, clump_r2 = 0.001,
                                  clump_p1 = 5e-8, clump_p2 = 1, pop = "EUR")
#Import outcome data
out <- "your_outcome_GWAS_sumstat.txt"
outcome_data <- read_outcome_data(out,
                                  snps = NULL,
                                  sep = "\t",
                                  snp_col = "ID",
                                  beta_col = "beta",
                                  se_col = "se",
                                  effect_allele_col = "effect_allele",
                                  other_allele_col = "other_allele",
                                  pval_col = "pval",
                                  ncase_col = "ncase",
                                  ncontrol_col = "ncontrol",
                                  samplesize_col = "samplesize",
                                  min_pval = 1e-200,
                                  log_pval = FALSE)
#Harmonization
dat <- harmonise_data(exposure_data_small, outcome_data, action = 2)
#---------standard MR
mr_results <- mr(dat, method_list = c('mr_ivw_fe','mr_weighted_median','mr_egger_regression'))
mr_results
generate_odds_ratios(mr_results)
#egger intercept 
egger <- mr_pleiotropy_test(dat)
egger
#egger I2
object <- subset(dat, select = c(beta.exposure, se.exposure, beta.outcome, se.outcome))
colnames(object) <- c('bx','bxse','by','byse')
object_in <- mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure, by = dat$beta.outcome, byse = dat$se.outcome)
res <- mr_egger(object_in)
res
Fstat <- mr_divw(object_in)
mr_heterogeneity(dat)
#Sensitivity analyses
res_single <- mr_singlesnp(dat)
res_loo <- mr_leaveoneout(dat)
#Plots
p1 <- mr_scatter_plot(mr_results, dat) # mrplot
p1 <- mr_forest_plot(res_single) #forest plot
p1 <- mr_funnel_plot(res_single) #funnel plot
p1 <- mr_leaveoneout_plot(res_loo) # leave one out plot
#Direction of causality
steiger <- directionality_test(dat)
steiger
#MR-presso
presso <- run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)
presso
#wrap results altogether
a <- mr_wrapper(dat)
a
#Eventually write a report
mr_report(dat)

####----------------------------------------------------------------------------------------------------------------------------###
###---------------------------4. Causal Analysis Using Summary Effect Estimates (CAUSE)-----------------------------------------###
####----------------------------------------------------------------------------------------------------------------------------###
#CAUSE consist in 4 steps:
#1 - Format the data for use with CAUSE
#2 - Calculate nuisance parameters
#3 - LD pruning
#4 - Fit CAUSE and look at the results
###Below the example script for one CAUSE analysis of "your_exposure" against "your_outcome", replace them as needed
#-----Format the data for use with CAUSE
#Import exposure data
your_exposure_file <- "your_exposure_GWAS_sumstat.txt"
exp <- fread(your_exposure_file)
#Import outcome data
your_outcome_file <- "your_outcome_GWAS_sumstat.txt"
out <- fread(your_outcome_file)
#Merge the two in a single dataset
exp_out <- gwas_merge(exp, out, snp_name_cols = c("SNP", "SNP"), 
                      beta_hat_cols = c("beta", "beta"), 
                      se_cols = c("se", "se"), 
                      A1_cols = c("effect_allele", "effect_allele"), 
                      A2_cols = c("other_allele", "other_allele"), 
                      pval_cols = c("pval", "pval"))

#-----Calculate nuisance parameters
set.seed(100)
varlist <- with(exp_out, sample(snp, size = 1000000, replace = FALSE))
params <- est_cause_params(exp_out, varlist)
#-----LD pruning
r2_thresh  <-  0.01
pval_thresh <-  1e-3
exp_out_clump <- ld_clump_local(exp_out, clump_kb = 10000,  clump_r2 = r2_thresh, clump_p = pval_thresh,
                                bfile = "/~/REF/EUR",
                                plink_bin = "/~/plinkbinr/bin/plink_Darwin")
top_vars <- exp_out_clump$rsid
#-----Fit CAUSE and look at the results
res <- cause(X = exp_out, variants = top_vars, param_ests = params)
#look at results
res$elpd
#It returns "delta_elpd": Estimated difference in elpd. If delta_elpd is negative, model 2 has a better fit
plot(res$sharing)
plot(res$causal)
summary(res, ci_size=0.95)
plot(res)
plot(res, type = "data")
res$causal
res$sharing

####----------------------------------------------------------------------------------------------------------------------------###
###-----------------------------------------------5. Multivariable MR (MVMR)----------------------------------------------------###
####----------------------------------------------------------------------------------------------------------------------------###
###Below the example script for one MVMR analysis of "your_exposure1" and "your_exposure2" against "your_outcome", replace them as needed
#Import data exposure1
exp1 <- "your_exposure1_GWAS_sumstat.txt"
exposure1 <- read_exposure_data(exp1,
                                clump = FALSE,
                                sep = " ",
                                snp_col = "SNP",
                                beta_col = "beta",
                                se_col = "se",
                                effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele",
                                pval_col = "pval",
                                ncase_col = "ncase",
                                ncontrol_col = "ncontrol",
                                samplesize_col = "samplesize",
                                min_pval = 1e-200,
                                log_pval = F)
exposure1_small <- subset(exposure1, exposure1$pval.exposure<5e-8)
exposure1_small <- clump_data(exposure1_small, clump_kb = 10000, clump_r2 = 0.001,
                              clump_p1 = 5e-8, clump_p2 = 1, pop = "EUR")
#Import data exposure2
exp2 <- "your_exposure2_GWAS_sumstat.txt"
exposure2 <- read_exposure_data(exp2,
                                clump = FALSE,
                                sep = " ",
                                snp_col = "SNP",
                                beta_col = "beta",
                                se_col = "se",
                                effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele",
                                pval_col = "pval",
                                ncase_col = "ncase",
                                ncontrol_col = "ncontrol",
                                samplesize_col = "samplesize",
                                min_pval = 1e-200,
                                log_pval = F)
exposure2_small <- subset(exposure2, exposure2$pval.exposure<5e-8)
exposure2_small <- clump_data(exposure2_small, clump_kb = 10000, clump_r2 = 0.001,
                              clump_p1 = 5e-8, clump_p2 = 1, pop = "EUR")
#Import outcome data
out <- "your_outcome_GWAS_sumstat.txt"
outcome_data <- read_outcome_data(out,
                                  snps = NULL,
                                  sep = "\t",
                                  snp_col = "ID",
                                  beta_col = "beta",
                                  se_col = "se",
                                  effect_allele_col = "effect_allele",
                                  other_allele_col = "other_allele",
                                  pval_col = "pval",
                                  ncase_col = "ncase",
                                  ncontrol_col = "ncontrol",
                                  samplesize_col = "samplesize",
                                  min_pval = 1e-200,
                                  log_pval = FALSE)
#Format the data for MVMR
#Harmonize each exposure with outcome data
out_exp1 <- harmonise_data(exposure1_small, outcome_data, action = 2)
out_exp2 <- harmonise_data(exposure2_small, outcome_data, action = 2)
#make a unique dataset with all the SNPs in either exposure and in the outcome GWAS
mvmr <- rbind(out_exp1, out_exp2)
#Add beta and se columns for each exposure relative to the SNPs
e1 <- exposure1[exposure1$SNP %in% mvmr$SNP, ]
str(e1)
e1 <- e1[,c(1,4,5)]
colnames(e1) <- c("SNP", "beta.exp1", "se.exp1")
mvmr <- merge(mvmr, e1, by="SNP")

e2 <- exposure2[exposure2$SNP %in% mvmr$SNP, ]
str(e2)
e2 <- e2[,c(1,4,5)]
colnames(e2) <- c("SNP", "beta.exp2", "se.exp2")
mvmr <- merge(mvmr, e2, by="SNP")
#Run MVMR
MRMVInputObject <- mr_mvinput(bx = cbind(mvmr$beta.exp1, mvmr$beta.exp2),
                              bxse = cbind(mvmr$se.exp1, mvmr$se.exp2),
                              by = mvmr$beta.outcome,
                              byse = mvmr$se.outcome,
                              exposure = c("exp1","exp2"),
                              outcome = "outcome")
MRMVObject <- mr_mvivw(MRMVInputObject,
                       model = "default",
                       correl = FALSE,
                       distribution = "normal",
                       alpha = 0.05)
MRMVObject