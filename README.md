-as.data.frame(data)

OR.meta <- data$OR.meta

data <- within (data, Beta= log(OR.meta))

library(dplyr)
data <- data %>% mutate(Beta = log(OR.meta))

names(data)[names(data)== "OA"] <- "other_allele"
names(data)[names(data)== "EA"] <- "effect_allele"

OR.meta <- data$OR.meta
data <- within (data, Beta= log(OR.meta))
library(dplyr)
data <- data %>% mutate(Beta = log(OR.meta))
as.in.data(alcohol_outcome_data)
as.data.frame(alcohol_outcome_data)

read_exposure_data <- function(anxiety_data, clump=FALSE, snp_col="SNPID", beta_col="Effect", se_col="StdErr", eaf_col="Freq1", effect_allele_col="Allele1", other_allele_col="Allele2", pval_col="P.value", samplesize_col="TotalN", min_pval=5e-8, log_pval=FALSE, chr_col="chr", pos_col="BP")
{
  exposure_dat <- data.table::fread(anxiety_data, header=TRUE, sep=sep)
  exposure_dat <- format_data(
    as.data.frame(exposure_dat),
    type="exposure",
    snps=NULL,
    snp_col=SNPID,
    beta_col=Effect,
    se_col=StdErr,
    eaf_col=Freq1,
    effect_allele_col=Allele1,
    other_allele_col=Allele2,
    pval_col=P.value,
    samplesize_col=TotalN,
    min_pval=min_pval,
    log_pval=log_pval,
    chr_col=chr,
    pos_col=BP
  )
  exposure_dat$data_source.exposure <- "textfile"
  if(clump)
  {
    exposure_dat <- clump_data(exposure_dat)
  }
  return(exposure_dat)
}
head(read_exposure_data)






anxiety.exposure <- read_exposure_data(
  filename = anxiety_data,
  snp_col = "SNPID",
  sep="/t",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "p.value",
  samplesize_col = "TotalN"
)
my_data <- read.csv("anxiety_data.txt")


anxiety.exposure <- read_exposure_data(
  filename = my_data
 
)
exposure_dat <- extract_instruments(c('ukb-b-18336'))
exposure_dat <- extract_instruments(c('finn-b-F5_ANXIETY'))
exposure_dat <- extract_instruments(c('ukb-b-17243'))


anxiety_exp_dat <- read_exposure_data(
  filename = my_data,
  snp_col = "SNPID",
  sep=",",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "p.value",
  samplesize_col = "TotalN"
)


exposure_dat <- data.table::fread(my_data, header=TRUE, sep=sep)
exposure_dat <- format_data(
  as.data.frame(my_data),
  type="exposure",
  snps=NULL,
  snp_col=SNPID,
  beta_col=Effect,
  se_col=StdErr,
  eaf_col=Freq1,
  effect_allele_col=Allele1,
  other_allele_col=Allele2,
  pval_col=P.value,
  samplesize_col=TotalN,
  min_pval=min_pval,
  log_pval=log_pval,
  chr_col=chr,
  pos_col=BP
)
my_data <- system.file("anxiety_data.csv", package = "TwoSampleMR")
my_data <- system.file("extdata/anxiety_data.csv", package = "TwoSampleMR")
anxiety_data <- system.file("extdata/anxiety_data.csv", package = "TwoSampleMR")
head(anxiety_data)

anxiety_exp_dat <- read_exposure_data(
  filename = my_data,
  snp_col = "SNPID",
  sep="\t",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "P.value",
  samplesize_col = "TotalN"
)
head(my_data)

anxiety_exp_dat$exposure <- "Anxiety"
my_data_exp_dat$exposure <- "Anxiety"

as.data.frame(my_data)
anxiety2_file <- system.file("extdata/my_data.csv", package = "TwoSampleMR")

anxiety2_file <- read_exposure_data(
  filename = my_data,
  snp_col = "SNPID",
  sep="\t",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "P.value",
  samplesize_col = "TotalN"
)

my_data2_file <- read_exposure_data(
  filename = my_data,
  snp_col = "SNPID",
  sep="\t",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "P.value",
  samplesize_col = "TotalN"
)
exposure_dat <- extract_instruments(c('ukb-b-18336'))
anxiety_exp_dat <- format_data(anxiety_data, type="exposure")

my_data <- data.frame(
  SNP = c("SNPID"),
  beta = c("Effect"),
  se = c("StdErr"),
  effect_allele = c("Allele1")
)

write.csv2(exposure_dat, exposure_dat.csv)

outcome_dat <- read_outcome_data(
  snps = exposure_dat,
  filename = "outcome_anxiety_data.xlsx",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "P.val",
  samplesize_col = "samplesize"
)
as.data.frame(outcome_anxiety_data)
as.data.frame(my_data)

exposure_dat <- data.table::fread(my_data, header=TRUE, sep = "\t")
n = read.csv('my_data', sep ="\\|", row.names = NULL ,check.names = FALSE, header = TRUE, stringsAsFactors = FALSE )
head(my_data)

exposure_dat

write.csv(my_data, "my_data.csv")
