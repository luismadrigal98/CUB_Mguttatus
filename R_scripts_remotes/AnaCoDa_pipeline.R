#!/usr/bin/env Rscript
#
# ******************************************************************************
# Runner of the MCMC for the AnaCoDa-based analysis
# 
# @author Luis Javier Madrigal-Roca
# 
# @date 11/20/2025
# ______________________________________________________________________________

library(argparse)
library(AnaCoDa)
library(stringr)

# ******************************************************************************
# 1) Parsing arguments ----
# ______________________________________________________________________________

parser <- ArgumentParser()

parser$add_argument("-i",
                    "--input",
                    help = "FASTA file with CDSs (transcript-based fasta file)",
                    type = "character",
                    default = "./")
parser$add_argument("-o",
                    "--output",
                    help = "Directory for writing out the results. Will automatically generate lowest-level directory if not already generated.",
                    type = "character",
                    default = ".")
parser$add_argument("-d",
                    "--div",
                    help = "Number of steps to diverge from starting values. Will be applied at beginning of each run, with the exception of the last.",
                    type = "integer",
                    default = 0)
parser$add_argument("-s",
                    "--samp",
                    help = "Number of samples (how far the MCMC goes...)",
                    type = "integer",
                    default = 5000)
parser$add_argument("-a",
                    "--adapt",
                    help = "Adaptive Width By Samples, i.e. will adapt every i samples",
                    type = "integer",
                    default = 50)
parser$add_argument("-t",
                    "--thin",
                    help = "Thinning value. Total number of iterations will be samples * thinning",
                    type = "integer",
                    default = 20)
parser$add_argument("-p",
                    "--percentage_to_keep",
                    help = "Percentage of samples used to estimate posterior for parameters. Default: last 50%",
                    default = 0.5)
parser$add_argument("-n",
                    "--threads",
                    help = "Number of threads to use for MCMC",
                    type = "integer",
                    default = 1)
parser$add_argument("--obs_phi",
                    type = "character",
                    default = NULL)
parser$add_argument("--dEta",
                    help = "Initial dEta values. Assumes csv format with columns AA,Codon,DEta. First line should be a header.",
                    type = "character")
parser$add_argument("--dM",
                    help = "Initial dM values. Assumes csv format with columns AA,Codon,DM. First line should be a header.",
                    type = "character")
parser$add_argument("--sphi_initial_values",
                    help = "Initial values for the parameter sphi. Notice that if input vector is shorter than number of mixtures the vector is going to be recycled.",
                    default = 1)
parser$add_argument("--phi",
                    help = "Initial Phi values (expression). Assumes csv format with Gene IDs in first column and Phi values in second column.", 
                    type = "character")
parser$add_argument("--est_phi",
                    help = "Use this flag to indicate estimation of Phi. Otherwise, Phi will not be estimated.",
                    action = "store_true")
parser$add_argument("--est_csp",
                    help = "Use this flag to indicate estimation of CSP. Otherwise, CSP will not be estimated.",
                    action = "store_true")
parser$add_argument("--est_hyp",
                    help = "Use this flag to indicate estimation of Hyperparameters. Otherwise, Hyperparameters will not be estimated.",
                    action = "store_true")
parser$add_argument("--est_mix",
                    help = "Use this flag to indicate estimation of mixture assignment. Otherwise, assignment will not be estimated.",
                    action = "store_true")
parser$add_argument("--max_num_runs",
                    help = "Max number of runs to perform.",
                    type = "integer",
                    default  =  10)
parser$add_argument("--fix_dEta",
                    help = "Use this flag to fix dEta at starting value.",
                    action = "store_true")
parser$add_argument("--fix_dM",
                    help = "Use this flag to fix dM at starting value.",
                    action = "store_true")
parser$add_argument("--number_of_mixtures",
                    type = "integer",
                    default = 1)
parser$add_argument("--mix_def",
                    type = "character",
                    default = "allUnique")
parser$add_argument("--mixture_assignment",
                    type = "character",
                    default = NULL)
parser$add_argument("--codon_table",
                    type = "integer",
                    default = 1)
parser$add_argument("--restart_file",
                    type = "character",
                    default = NULL)
parser$add_argument("--init_sphi",
                    type = "character",
                    default = NULL)

args <- parser$parse_args()

message("Setting up AnaCoDa run with the following arguments ...")
message("========================================================")
message("          AnaCoDa MCMC Configuration Summary            ")
message("========================================================")

for (name in names(args)) {
  val <- args[[name]]
  # Handle NULLs so they print nicely
  if (is.null(val)) {
    val_str <- "NULL"
  } else {
    val_str <- as.character(val)
  }
  # Print aligned key-value pairs
  message(sprintf("%-25s : %s", name, val_str))
}
message("========================================================")

input <- args$input
directory <- args$output
div <- args$div
samples <- args$samp
adaptiveWidth <- args$adapt
thinning <- args$thin
percentage.to.keep <- args$percentage_to_keep
num.threads <- args$threads
dEta.file <- args$dEta
dM.file <- args$dM
sphi.init <- args$sphi_initial_values
phi.files <- args$phi
est.csp <- args$est_csp
est.phi <- args$est_phi
est.hyp <- args$est_hyp
est.mix <- args$est_mix
max.num.runs <- args$max_num_runs
fix.dEta <- args$fix_dEta
fix.dM <- args$fix_dM
number.of.mixtures <- args$number_of_mixtures
mix.def <- args$mix_def
mix.assign <- args$mixture_assignment
restart.file <- args$restart_file
codon.table <- args$codon_table
init.sphi <- args$init_sphi

# Expression
obs.phi <- args$phi 
with.phi <- !is.null(obs.phi)

# ******************************************************************************
# 2) Auxiliary functions ----
# ______________________________________________________________________________

calcDeltaMFromIntronsLocal <- function(...) {
  return(0)
}

## Outputs CSP estimates
createParameterOutput <- function(parameter, dir_name, numMixtures, 
                                  samples, mixture.labels, 
                                  samples.percent.keep=1, 
                                  relative.to.optimal.codon=F, 
                                  report.original.ref=T)
{
  for (i in 1:numMixtures)
  {
    getCSPEstimates(parameter, paste(dir_name, "Parameter_est", 
                                     mixture.labels[i], sep="/"), i, 
                    samples*samples.percent.keep, 
                    relative.to.optimal.codon=relative.to.optimal.codon, 
                    report.original.ref = report.original.ref)
  }
}

## Outputs traces for CSPs and plots expected frequencies as function of log10(\phi) 
createTracePlots <- function(trace, 
                             model,
                             genome,
                             numMixtures,
                             samples,
                             mixture.labels,
                             samples.percent.keep = 1)
{
  for (i in 1:numMixtures)
  {
    plot(trace, what = "Mutation", mixture = i)
    plot(trace, what = "Selection", mixture = i)
    
    plot(model, genome, samples = samples*samples.percent.keep, 
         mixture = i,main = mixture.labels[i])
  }
  plot(trace, what="AcceptanceRatio")
}

# ******************************************************************************
# 3) Setting up the MCMC ----
# ______________________________________________________________________________

if (with.phi && !is.null(obs.phi))
{
  genome <- initializeGenomeObject(file = input,
                                   match.expression.by.id = TRUE,
                                   observed.expression.file = obs.phi)
} else
  {
  genome <- initializeGenomeObject(file = input)
  }

size <- length(genome)
index <- c(1:size)

# 3.1) Mixture ----

if (!is.null(mix.assign))
{
  tmp <- read.csv(mix.assign,sep="\t",header=T,stringsAsFactors=F)
  geneAssignment <- tmp[,2]
  numMixtures <- length(unique(tmp[,2]))
  mixture.labels <- as.character(sort(unique(tmp[,2])))
} else if (is.null(mix.assign) && number.of.mixtures > 1)
{
  warning("Number of mixtures greater than 1 but no assignment was provided. \n
  Gene assignments are going to be estimated")
  est.mix = TRUE # Notice that this overrides the user input.
} else
{
  geneAssignment <- rep(1,size)
  mixture.labels <- paste0("Cluster_", geneAssignment)
  numMixtures <- 1
}

# 3.2) Phi and sphi ----

init_phi <- NULL
sphi_init <- rep(sphi.init, numMixtures)

if (!is.null(phi.files))
{
  segment_exp <- read.table(file = phi.files,
                            sep = ",",
                            header = TRUE)
  
  init_phi <- c(init_phi,
                segment_exp[,"Mean"])
  
  sphi_init <- rep(sd(log(init_phi)),
                   numMixtures)
  
  if(length(genome) !=  length(init_phi))
  {
    stop("length(genomeObj) !=  length(init_phi), but it should.")
  } else{
    message("Initial Phi values successfully files loaded.");
  }
}

if (!is.null(obs.phi))
{
  obs.phi <- read.csv(obs.phi,
                      header=T,
                      row.names=1)
  n.obs.phi <- ncol(obs.phi)
  
  s_eps <- rep(0.1,
               n.obs.phi)
} else {
  s_eps <- 0.1
}

mutation.prior.mean <- 0

# ******************************************************************************
# 4) Layout for running the MCMC iteratively ----
# ______________________________________________________________________________

if (!dir.exists(directory)) dir.create(directory)

# Check if input FASTA file exists
if (!file.exists(input)) {
  stop(paste("Input FASTA file not found:", input))
}

# Helper function for safe directory creation
safe_dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(path)) {
      stop(paste("Failed to create directory:", path))
    }
  }
}

safe_dir_create(directory)

done <- FALSE
done.adapt <- FALSE
run_number <- 1
param.conv <- TRUE

while((!done) && (run_number <= max.num.runs))
{
  if (adaptiveWidth == 0)
  {
    percent.to.keep <- 1
    adaptiveWidth <- 20
    div_run <- 0
  } else {
    percent.to.keep <- 0.5
    div_run <- div
  }
  if (is.null(restart.file))
  {
    parameter <- initializeParameterObject(genome,sphi_init,numMixtures, geneAssignment,init.sepsilon = s_eps,split.serine = TRUE, mixture.definition = mix.def, initial.expression.values = init_phi,init.w.obs.phi=with.phi,mutation.prior.mean=mutation.prior.mean)
    if (length(dM.file) > 0)
    {
      if (!file.exists(dM.file)) {
        stop(paste("dM file not found:", dM.file))
      }
      parameter$initMutationCategories(dM.file,1,fix.dM)
    } 
    if (length(dEta.file) > 0)
    {
      if (!file.exists(dEta.file)) {
        stop(paste("dEta file not found:", dEta.file))
      }
      parameter$initSelectionCategories(dEta.file,1,fix.dEta)
    }
    
  } else {
    previous <- stringr::str_extract(pattern="run_[0-9]+",string=restart.file)
    run_number <- as.numeric(stringr::str_extract(pattern="[0-9]+",string=previous)) + 1
    parameter<-initializeParameterObject(init.with.restart.file = restart.file,model="ROC")
  }
  if (!is.null(init.sphi))
  {
    tmp <- as.numeric(readLines(init.sphi))
    parameter$setStdDevSynthesisRate(tmp,0)
  }
  steps.to.adapt <- (samples*thinning)*(1-percent.to.keep)
  dir_name <- paste0(directory,"/run_",run_number)
    safe_dir_create(dir_name)
    safe_dir_create(paste(dir_name,"Graphs",sep="/"))
    safe_dir_create(paste(dir_name,"Restart_files",sep="/"))
    safe_dir_create(paste(dir_name,"Parameter_est",sep="/"))
    safe_dir_create(paste(dir_name,"R_objects",sep="/"))
  
  mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, 
                               adaptive.width = adaptiveWidth,
                               est.expression = est.phi, 
                               est.csp = est.csp, 
                               est.hyper = est.hyp,
                               est.mix = est.mix)
  
  mcmc$setStepsToAdapt(steps.to.adapt)
  
  model <- initializeModelObject(parameter, 
                                 "ROC", 
                                 with.phi,
                                 fix.observation.noise=F)
  setRestartSettings(mcmc, 
                     paste(dir_name,"Restart_files/rstartFile.rst",sep="/"), 
                     adaptiveWidth, F)
  
  
  sys.runtime <- system.time(
    runMCMC(mcmc, genome, model, num.threads, div = div_run)
  )
  
  # Output runtime for mcmc
  sys.runtime <- data.frame(Value = names(sys.runtime),
                            Time = as.vector(sys.runtime))
  write.table(sys.runtime,file = paste(dir_name,"mcmc_runtime.csv",sep="/"),
              sep=",", col.names = T, row.names = T, quote = F)
  
  # Creates R objects, which can be later loaded for re-analzying already completed runs
  writeParameterObject(parameter, paste(dir_name,"R_objects/parameter.Rda",
                                        sep="/"))
  writeMCMCObject(mcmc, 
                  file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))
  
  # Output CSP file
  createParameterOutput(parameter = parameter, dir_name = dir_name, numMixtures = numMixtures, samples = samples, mixture.labels = mixture.labels, samples.percent.keep = percent.to.keep, relative.to.optimal.codon = F, report.original.ref = T)
  
  # Output phi file
  expressionValues <- getExpressionEstimates(parameter,c(1:size),
                                             samples*percent.to.keep,
                                             genome = genome)
  write.table(expressionValues,file=paste(dir_name,
                                          "Parameter_est/gene_expression.txt",
                                          sep="/"), sep=",", col.names = T, 
              quote = F, row.names = F)
  
  # Plots different aspects of trace
  trace <- parameter$getTraceObject()
  pdf(paste(dir_name,"Graphs/mcmc_traces.pdf",sep="/"))
  plot(mcmc,what = "LogPosterior")
  plot(mcmc,what="LogLikelihood")
  plot(trace, what = "ExpectedPhi")
  if (est.hyp)
  {
    plot(trace,what="Sphi")
  }
  if (with.phi)
  {
    plot(trace,what="Sepsilon")
  }
  if (est.csp)
  {
    ## Calculate auto-correlation and convergence of CSP traces
    param.conv <- TRUE
    if (!fix.dEta)
    {
      acfCSP(parameter,
             csp = "Selection",
             numMixtures = numMixtures,
             samples = samples * percent.to.keep)
      for (i in 1:numMixtures)
      {
        param.diag<-convergence.test(trace,
                                     samples = samples * percent.to.keep,
                                     thin = thinning,
                                     what = "Selection",
                                     mixture = i,
                                     frac1 = 0.25,frac2 = 0.5)
        z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
        if (length(z.scores) > 5)
        {
          param.conv <- FALSE
        }
        write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_eta_",i,".txt"),ncolumns = 1)
      }
    }
    if (!fix.dM)
    {
      acfCSP(parameter,csp = "Mutation",
             numMixtures = numMixtures,
             samples = samples * percent.to.keep)
      for (i in 1:numMixtures)
      {
        param.diag<-convergence.test(trace,samples=samples*percent.to.keep,
                                     thin = thinning,what="Mutation",
                                     mixture = i,frac1 = 0.25, frac2 = 0.5)
        z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
        if (length(z.scores) > 5)
        {
          param.conv <- FALSE
        }
        write(param.diag$z,
              paste0(dir_name,
                     "/Parameter_est/convergence_delta_M_",i,".txt"),
              ncolumns = 1)
      }
    }
  }
  dev.off()
  
  pdf(paste(dir_name,"Graphs/CSP_traces_CUB_plot.pdf",sep="/"), width = 11, height = 12)
  createTracePlots(trace=trace,model=model,genome=genome,numMixtures=numMixtures,samples=samples,samples.percent.keep = percent.to.keep,mixture.labels = mixture.labels)
  dev.off()
  
  diag <- convergence.test(mcmc,
                           samples = samples * percent.to.keep,
                           thin = thinning, 
                           frac1 = 0.1,
                           frac2 = 0.5)
  z <- abs(diag$z)
  
  ## Can end if overall log(posterior) and CSP parameters have converged
  #done <- (z < 1.96) && param.conv
  done <- F
  rm(parameter)
  rm(trace)
  rm(model)
  run_number <- run_number + 1
}