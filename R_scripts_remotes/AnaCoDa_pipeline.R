#!/usr/bin/env Rscript
#
# ******************************************************************************
# Runner of the MCMC for the AnaCoDa-based analysis
# 
# @author Luis Javier Madrigal-Roca
# 
# @date 11/20/2025
#
# IMPORTANT NOTE ON ANACODA BUG WORKAROUND:
# There is a known bug in AnaCoDa where using with.phi=TRUE in the model
# (i.e., including observed expression in the likelihood) causes a segmentation
# fault during MCMC execution. This occurs regardless of whether dM is fixed.
#
# WORKAROUND IMPLEMENTED:
# When expression data is provided (--phi) AND mutation parameters are fixed
# (--fix_dM), this script:
#   1. Initializes phi values FROM observed expression (init.w.obs.phi = TRUE)
#   2. Does NOT estimate phi (est.expression = FALSE) - phi stays at empirical values
#   3. Sets with.phi = FALSE in the model to avoid the segfault
#
# This means phi is FIXED at empirical expression values, and only selection
# coefficients (dEta) are estimated conditional on those fixed phi values.
# This is a valid analysis approach when you trust your expression measurements.
#
# See GitHub issue: https://github.com/clandere/AnaCoDa/issues/XXX
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
parser$add_argument("--dEta",
                    help = "Initial dEta values. Assumes csv format with columns AA,Codon,DEta. First line should be a header.",
                    type = "character")
parser$add_argument("--dM",
                    help = "Initial dM values. Assumes csv format with columns AA,Codon,DM. First line should be a header.",
                    type = "character")
parser$add_argument("--sphi_init",
                    help = "Initial values for the parameter sphi. Notice that if input vector is shorter than number of mixtures the vector is going to be recycled.",
                    default = 1)
parser$add_argument("--sepsilon_init",
                    help = "Initial values for the sepsilon hyperparameter. Notice that AnaCoDa expects one value per data source in the phi csv. Example: '0.5, 0.1', always as a string",
                    default = NULL)
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
                    help = "sphi value per mixture as a string. If multiple values, separate them using `,`",
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
sphi.init <- args$sphi_init
sepsilon.init <- args$sepsilon_init
phi.file <- args$phi
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

# Expression
obs.phi <- args$phi 
with.phi <- !is.null(obs.phi)

# ******************************************************************************
# FIXED PHI WORKAROUND (modeling choice, not a bug workaround)
# When using fixed dM with expression data, we treat observed expression as the
# "true" phi values and estimate only selection coefficients conditional on them:
#   - Initialize phi FROM observed expression
#   - Do NOT estimate phi (keep it fixed at empirical values)
#   - Set with.phi=FALSE in the model
#
# NOTE ON ANACODA BUG (now fixed in AnaCoDA_fixed/AnaCoDa):
# The original segfault (address 0x29) when withPhi=TRUE was caused by a missing
# argument in Trace::initializeROCTrace() — it passed estimateSynthesisRate (bool)
# in the position of numObservedPhiSets, initializing SynthesisOffset traces with
# size 1 instead of the actual number of phi groupings. This has been fixed in
# AnaCoDA_fixed/AnaCoDa/src/Trace.cpp. Once the fixed package is installed, the
# phi-only model (withPhi=TRUE, est.phi=TRUE) will work correctly.
# ******************************************************************************
use.fixed.phi.workaround <- with.phi && fix.dM

if (use.fixed.phi.workaround) {
  message("")
  message("============================================================")
  message("  FIXED-PHI MODE: Fixed dM with expression data detected    ")
  message("============================================================")
  message("Phi will be FIXED at empirical expression values (modeling")
  message("choice). Only CSP and hyperparameters will be estimated.")
  message("============================================================")
  message("")

  # Override est.phi to FALSE - phi will not be estimated
  est.phi <- FALSE
}

# ******************************************************************************
# 2) Auxiliary functions ----
# ______________________________________________________________________________

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

# --- 3.1 Genome Initialization ---
# CRITICAL: Validate input files before passing to C++ initializer

if (!file.exists(input)) {
  stop(paste("Input FASTA file not found:", input))
}

if (with.phi && !is.null(obs.phi)) {
  if (!file.exists(obs.phi)) {
    stop(paste("Expression file not found:", obs.phi))
  }

  # Read phi file using fread (auto-detects TSV or CSV)
  phi_raw <- tryCatch(
    data.table::fread(obs.phi, header = TRUE, data.table = FALSE),
    error = function(e) stop(paste("Failed to read expression file:", e$message))
  )
  if (ncol(phi_raw) < 2) {
    stop("Expression file must have at least 2 columns (GeneID + expression values)")
  }

  # Filter organellar genes (gene IDs with .O followed by digits).
  # These genes have codon usage and expression levels incompatible with the
  # nuclear ROC model: e.g., O006900 has phi up to 6792 TPM, which causes
  # exp(phi * selection) overflow -> NaN in calculateLogCodonProbabilityVector.
  organellar_mask <- grepl("\\.O[0-9]", phi_raw[[1]])
  if (any(organellar_mask)) {
    n_org <- sum(organellar_mask)
    org_ids <- phi_raw[[1]][organellar_mask]
    message(sprintf("Filtering %d organellar gene(s) from phi file: %s",
                    n_org, paste(org_ids, collapse = ", ")))
    phi_raw <- phi_raw[!organellar_mask, , drop = FALSE]
  }

  message(paste("Expression file has", ncol(phi_raw) - 1, "expression column(s),",
                nrow(phi_raw), "genes after filtering"))

  # Normalize phi values so their mean ≈ 1 (consistent with AnaCoDa's LogNormal prior,
  # which has mean = exp(-sigma^2/2) ≈ 1).  Without normalization, CPM values up to
  # ~37 000 cause numerical instability: extreme phi enters calculateLogCodonProbability-
  # Vector and calculateLogLikelihoodRatioForHyperParameters in ways the overflow guards
  # do not fully cover (log(phi) → -Inf when phi ≈ 0 → NaN via -Inf - (-Inf)).
  phi_cols <- 2:ncol(phi_raw)
  all_phi_vals <- unlist(phi_raw[, phi_cols])
  phi_norm_mean <- mean(all_phi_vals[all_phi_vals > 0])
  phi_raw[, phi_cols] <- phi_raw[, phi_cols] / phi_norm_mean
  message(sprintf("Normalized phi by mean of non-zero values (%.4f); max phi after norm: %.2f",
                  phi_norm_mean, max(phi_raw[, phi_cols])))

  # Write as CSV — AnaCoDa's C++ readObservedPhiValues requires comma-delimited input
  obs.phi.filtered <- tempfile(fileext = ".csv")
  on.exit(unlink(obs.phi.filtered), add = TRUE)
  data.table::fwrite(phi_raw, obs.phi.filtered)

  # Initialize genome with filtered expression data
  genome <- initializeGenomeObject(file = input,
                                   match.expression.by.id = TRUE,
                                   observed.expression.file = obs.phi.filtered)
} else {
  genome <- initializeGenomeObject(file = input)
}

# Remove organellar genes from the genome object (gene IDs matching .O followed by digits).
# These genes may be present in the FASTA and are incompatible with the nuclear ROC model:
# their extreme expression values (e.g., O006900 up to 6792 TPM) cause exp(phi*selection)
# overflow -> NaN in calculateLogCodonProbabilityVector during MCMC.
gene_names_all <- getNames(genome)
organellar_in_genome <- grepl("\\.O[0-9]", gene_names_all)
if (any(organellar_in_genome)) {
  n_removed <- sum(organellar_in_genome)
  removed_ids <- gene_names_all[organellar_in_genome]
  message(sprintf("Removing %d organellar gene(s) from genome: %s",
                  n_removed, paste(removed_ids, collapse = ", ")))
  nuclear_indices <- which(!organellar_in_genome)  # 1-indexed for AnaCoDa
  genome <- genome$getGenomeForGeneIndices(nuclear_indices, FALSE)
}

size <- length(genome)
message(paste("Genome loaded with", size, "genes"))

if (size == 0) {
  stop("No genes loaded from FASTA file. Check file format.")
}

index <- c(1:size)

# --- 3.2 Mixture Setup ---

if (!is.null(mix.assign)) {
  tmp <- read.csv(mix.assign, sep="\t", header=T, stringsAsFactors=F)
  geneAssignment <- tmp[,2]
  numMixtures <- length(unique(tmp[,2]))
  mixture.labels <- as.character(sort(unique(tmp[,2])))
} else if (is.null(mix.assign) && number.of.mixtures > 1) {
  warning("Number of mixtures > 1 but no assignment provided. Estimating assignment.")
  geneAssignment <- rep(1, size)
  numMixtures <- number.of.mixtures
  mixture.labels <- paste0("Cluster_", 1:numMixtures)
  est.mix <- TRUE 
} else {
  geneAssignment <- rep(1, size)
  mixture.labels <- c("Cluster_1")
  numMixtures <- 1
}

message(paste("Number of mixtures:", numMixtures))
message(paste("Mixture labels:", paste(mixture.labels, collapse=", ")))

# --- 3.3 Phi (Initial Values) Setup ---

init_phi <- NULL

# --- 3.4 Hyperparameters (sphi & sepsilon) ---

# SPHI Handling
# Parse the comma-separated string from arguments
if (!is.null(sphi.init)) {
  sphi_vals <- as.numeric(unlist(strsplit(as.character(sphi.init), ",")))
  
  # If 1 mixture, ensure scalar. If multiple, ensure vector length matches.
  if (numMixtures == 1) {
    sphi_input <- sphi_vals[1] 
  } else {
    sphi_input <- rep(sphi_vals, length.out = numMixtures)
  }
} else {
  sphi_input <- rep(1, numMixtures) # Default fallback
}

# SEPSILON Handling
if (!is.null(sepsilon.init)) {
  # Use user input vector
  s_eps <- as.numeric(unlist(strsplit(as.character(sepsilon.init), ",")))
} else if (!is.null(obs.phi)) {
  # Auto-detect number of columns if no sepsilon provided
  # Read just the header to count columns efficiently
  temp_df <- data.table::fread(obs.phi, nrows = 1, header = TRUE, data.table = FALSE)
  n_sources <- ncol(temp_df) - 1 # Minus GeneID column
  s_eps <- rep(0.1, n_sources)
} else {
  s_eps <- 0.1 # Fallback
}

mutation.prior.mean <- 0

# ******************************************************************************
# 4) Layout for running the MCMC iteratively ----
# ______________________________________________________________________________

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
    # WORKAROUND: When using fixed dM with expression data, initialize phi
    # FROM observed expression values. Since est.phi is set to FALSE above,
    # these values will remain fixed throughout the MCMC.
    use.obs.phi.for.init <- use.fixed.phi.workaround
    
    parameter <- initializeParameterObject(genome = genome, 
                                           sphi = sphi_input,
                                           num.mixtures = numMixtures, 
                                           gene.assignment = geneAssignment,
                                           init.sepsilon = s_eps,
                                           split.serine = TRUE, 
                                           mixture.definition = mix.def, 
                                           initial.expression.values = NULL,
                                           init.w.obs.phi = use.obs.phi.for.init,
                                           mutation.prior.mean = mutation.prior.mean)
    
    # Initialize dM (mutation) categories if file provided
    if (!is.null(dM.file) && nchar(dM.file) > 0)
    {
      if (!file.exists(dM.file)) {
        stop(paste("dM file not found:", dM.file))
      }
      message(paste("Initializing mutation categories from:", dM.file))
      parameter$initMutationCategories(c(dM.file), 1, fix.dM)
    } 
    
    # Initialize dEta (selection) categories if file provided
    if (!is.null(dEta.file) && nchar(dEta.file) > 0)
    {
      if (!file.exists(dEta.file)) {
        stop(paste("dEta file not found:", dEta.file))
      }
      message(paste("Initializing selection categories from:", dEta.file))
      parameter$initSelectionCategories(c(dEta.file), 1, fix.dEta)
    }
    
  } else {
    previous <- stringr::str_extract(pattern="run_[0-9]+",string=restart.file)
    run_number <- as.numeric(stringr::str_extract(pattern="[0-9]+",string=previous)) + 1
    parameter<-initializeParameterObject(init.with.restart.file = restart.file,model="ROC")
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
  
  # WORKAROUND: Set with.phi=FALSE in model when using the fixed phi workaround
  # to avoid AnaCoDa segfault. In this mode:
  #   - Phi is initialized from and fixed at empirical expression values
  #   - Phi is NOT estimated (est.expression = FALSE)
  #   - The likelihood does not include observation error model
  #   - Only selection coefficients (dEta) are estimated
  use.phi.in.model <- with.phi && !use.fixed.phi.workaround
  
  model <- initializeModelObject(parameter, 
                                 "ROC", 
                                 use.phi.in.model,
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
  # Only plot Sepsilon trace if we're actually using phi in the model
  # (not when using the fixed phi workaround)
  if (with.phi && !use.fixed.phi.workaround)
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