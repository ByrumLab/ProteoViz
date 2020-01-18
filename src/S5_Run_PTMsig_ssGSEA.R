source('src/ssGSEA2.0.R')
rel_path <- "data/PTMsig"
gct_file <- paste0(rel_path, "/phospho_PTMsig_input.gct")
dir.create(paste0(rel_path, "/output"))
output_location <- paste0(rel_path, "/output/output")


## ##########################################################
##  define parameters below:
## ##########################################################

## sGSEA / PSEA parameters
sample.norm.type    = "rank"              ## "rank", "log", "log.rank", "none" 
weight              = 0.75                ## value between 0 (no weighting) and 1 (actual data counts)
statistic           = "area.under.RES"    ## "Kolmogorov-Smirnov"
output.score.type   = "NES"               ## 'ES' or 'NES'
nperm               = 1e3                 ## No. of permutations
min.overlap         = 3                  ## minimal overlap between gene set and data
correl.type         = "z.score"           ## 'rank', 'z.score', 'symm.rank'
par                 = T                   ## use 'doParallel' package?
spare.cores         = 1                   ## No. of cores to leave idle
export.signat.gct   = F                   ## if TRUE gene set GCT files will be exported 
extended.output     = T                   ## if TRUE the GCT files will contain stats on gene set overlaps etc.   

sig_db <- "databases/ptm.sig.db.all.flanking.human.v1.8.1.gmt"

## import signature database
signat.all <- unlist(lapply(sig_db, readLines))
signat.all <- strsplit(signat.all, '\t')
names(signat.all) <- sapply(signat.all, function(x)x[1])
signat.all <- lapply(signat.all, function(x) x[-c(1,2)])

names(gct_file) <- paste(  sub('\\.gct$', '', sub('.*/','', gct_file)), 'ssGSEA', sep='_' )
input.ds <- gct_file


## ########################################
## ssGSEA

i <- 1
gsea.res <-
  ssGSEA2(
    gct_file,
    gene.set.databases = sig_db,
    sample.norm.type = sample.norm.type,
    weight = weight,
    statistic = statistic,
    output.score.type = output.score.type,
    nperm  = nperm,
    min.overlap  = min.overlap,
    correl.type = correl.type,
    output.prefix = output_location,
    par = par,
    spare.cores = spare.cores,
    param.file = F,
    export.signat.gct = export.signat.gct,
    extended.output = extended.output
  )


