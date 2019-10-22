library(tidyverse)
chr <- 22
p <- 9537
n <- 2000

tmp.dir <- "/scratch/PI/candes/ukbiobank_tmp/fastphase/"
K.list <- c(10,20,30,40,50,60,70)
    
Results <- tibble()
for( K in K.list ){
    filename <- paste(tmp.dir, "phased_K",K,"/ukb_hap_chr",chr,"_finallikelihoods", sep="")
    dat <- read_delim(filename, skip=1, delim="\t", col_names=c("Seed", "Likelihood"))
    dat$K <- as.integer(K)
    Results <- rbind(Results,dat)
}

Results$Parameters <- p*(Results$K+1)+Results$K
Results$AIC <- 2*Results$Parameters - 2*Results$Likelihood
Results$BIC <- log(n) * Results$Parameters - 2*Results$Likelihood
    
Results %>% ggplot(aes(x=K, y=Likelihood)) + geom_point() + theme_bw()
Results %>% ggplot(aes(x=K, y=AIC)) + geom_point() + theme_bw()
