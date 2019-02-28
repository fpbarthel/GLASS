suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(plyr))
suppressPackageStartupMessages(require(dplyr))

#' Arguments for seqz #####
outdir        <- snakemake@params[["outdir"]]
prefix        <- snakemake@params[["prefix"]]
seqz_gz_file  <- snakemake@input[[1]]
break_method  <- snakemake@params[["break_method"]]
kmin          <- snakemake@params[["kmin"]]
mc.cores      <- snakemake@threads

cat("\n\n#######################################################\n\n")
cat("...... prefix= ",prefix,"\n")
cat("...... outdir= ",outdir,"\n")
cat("...... seqz_gz_file= ",seqz_gz_file,"\n")
cat("...... kmin= ",kmin,"\n")
cat("...... break_method= ",break_method,"\n")
cat("...... mc.cores= ",mc.cores,"\n")
cat("#######################################################\n\n")

#setwd(outdir)

library(sequenza)
is.error <- function(x) inherits(x, "try-error")

#' ###################################################
#' ### GC-normalization
#' ###################################################

#' Read the seqz file ----------
seqz.data   <- read.seqz(seqz_gz_file)

gc.stats    <- gc.sample.stats(seqz_gz_file)

gc.stats    <- gc.norm(x = seqz.data$depth.ratio, gc = seqz.data$GC.percent)

gc.vect     <- setNames(gc.stats$raw.mean, gc.stats$gc.values)

seqz.data$adjusted.ratio <- seqz.data$depth.ratio / gc.vect[as.character(seqz.data$GC.percent)]
saveRDS(seqz.data,file=file.path(outdir,paste(prefix,'_sequenza_gcnorm.rds',sep="")))

pdf(file.path(outdir,paste(prefix,'_sequenza_gcnorm.pdf',sep="")),width=11,height=6)
par(mfrow = c(1,2), cex = 1, las = 1, bty = 'l')
matplot(gc.stats$gc.values, gc.stats$raw,type = 'b', col = 1, pch = c(1, 19, 1), lty = c(2, 1, 2), xlab = 'GC content (%)', ylab = 'Uncorrected depth ratio')
legend('topright', legend = colnames(gc.stats$raw), pch = c(1, 19, 1))
hist2(seqz.data$depth.ratio, seqz.data$adjusted.ratio,breaks = prettyLog, key = vkey, panel.first = abline(0, 1, lty = 2),xlab = 'Uncorrected depth ratio', ylab = 'GC-adjusted depth ratio')
dev.off()

    
#' ###################################################
#' ### Extract the information from the seqz file.
#' ###################################################
cat('\n....run sequenza.extract with a kmin',kmin,'break_method=',break_method,'\n')
z <-try(sequenza_extract_data <- sequenza.extract(seqz_gz_file,kmin = kmin, breaks.method=break_method))
if (is.error(z)) {
  stop('error: sequenza.extract')
} else {          
    names(sequenza_extract_data)
    saveRDS(sequenza_extract_data,
        file=file.path(outdir,paste(prefix,'_sequenza_extract_data.rds',sep="")))
}

#' ###################################################
#' ### Inference of cellularity and ploidy
#' ###################################################
z <-try(cp_table_data <- sequenza.fit(sequenza_extract_data,mc.cores=mc.cores))
if (is.error(z)) {
    stop("error in sequenza.fit")
} else {
    saveRDS(cp_table_data,
            file=file.path(outdir,paste(prefix,'_sequenza_fit_data.rds',sep="")))
}

#' ###################################################
#' ### Results of model fitting
#' ###################################################
#The last part of the workflow is to apply the estimated parameters. There is an
#all-in-one function that plots and saves the results, giving control on file names
#and output directory:

# sequenza.results--------
## cint, cellularity, ploidy are calculated in the same way in sequenza.results function,
## will be used for plotting cp plot.
cint        <- get.ci(cp_table_data)
cellularity <- cint$max.cellularity
ploidy      <- cint$max.ploidy
z <-try(sequenza.results(sequenza.extract = sequenza_extract_data, 
                        cp.table = cp_table_data,
                        sample.id = prefix, 
                        out.dir=outdir,
                        cellularity = cellularity, 
                        ploidy = ploidy))
if (is.error(z)) {
  stop('ERROR: sequenza.results') 
}

# Write purity/ploidy------------
outf  <-file.path(outdir,paste(prefix,'_cellularity.ploidy.txt',sep=""))
write.table(data.frame(cellularity,ploidy),file=outf, sep="\t",col.names=TRUE,row.names=FALSE)


  
#' ### Get results from alternative solution ########################
.runSequenza_altSol_v1 <-function(outdir, prefix) {
  
  require(sequenza)
  require(plyr)
  require(dplyr)
  
  seqz_extract_data_rds 	<-file.path(outdir,paste0(prefix,'_sequenza_extract_data.rds'))
  alt_solu_file           <-file.path(outdir,paste0(prefix,'_alternative_solutions.txt'))
  cp_fn 									<-file.path(outdir,paste0(prefix,'_cellularity.ploidy.txt'))
  seqz_fit_fn 						<-file.path(outdir,paste0(prefix,'_sequenza_fit_data.rds'))
  
  stopifnot(file.exists(seqz_extract_data_rds))
  stopifnot(file.exists(alt_solu_file))
  stopifnot(file.exists(cp_fn))
  
  cat("\n\n..........................................\n\n")
  cat("... Getting alternative solution. prefix",prefix,"\toutdir",outdir,"\n")
  cat("..........................................\n\n")
  
  #setwd(outdir)
  
  ###################################################
  
  #' Get seqz extract data ########
  seqz_ext_data <-readRDS(seqz_extract_data_rds)
  
  #' Get the best solution for cellularity and ploidy ########
  cp_sol <-read.delim(cp_fn,sep="\t",header=TRUE,as.is=TRUE)
  head(cp_sol)
  cp_sol <-mutate(cp_sol,
                  cellularity=as.numeric(as.matrix(cellularity)),
                  ploidy=as.numeric(as.matrix(ploidy)))
  
  #' Get alternative solution for cellularity and ploidy ########
  alt_sol <-read.delim(alt_solu_file,sep="\t",header=TRUE,as.is=TRUE)
  head(alt_sol)
  alt_sol <-mutate(alt_sol,
                   cellularity=as.numeric(as.matrix(cellularity)),
                   ploidy=as.numeric(as.matrix(ploidy)))
  #stopifnot(nrow(alt_sol)>=2)
  
  #' 20180720: 1st line in alt sol may or may not be identical to cp_sol ########
  #stopifnot(identical(alt_sol$cellularity[1],cp_sol$cellularity))
  #stopifnot(identical(alt_sol$ploidy[1],cp_sol$ploidy))
  
  for (iix in 1:nrow(alt_sol)) {
    
    altsol_outdir <-file.path(outdir,'alt_sols',paste0('alt_sol_',iix))
    
    if (!file.exists(dirname(altsol_outdir))) {
      dir.create(dirname(altsol_outdir))
    }
    
    #if (file.exists(altsol_outdir)) unlink(altsol_outdir)
    if (!file.exists(altsol_outdir)) dir.create(altsol_outdir)
    
    write.table(alt_sol[iix,],
                file=file.path(altsol_outdir,"cellularity.ploidy.txt"),
                sep="\t",
                col.names=TRUE,
                row.names=FALSE,
                quote=FALSE)
    
    ## run sequenza.results with the alternative solution ---------
    # function (sequenza.extract, cp.table = NULL, sample.id, out.dir = getwd(),
    #           cellularity = NULL, ploidy = NULL, female = TRUE, CNt.max = 20,
    #           ratio.priority = FALSE, XY = c(X = "X", Y = "Y"), chromosome.list = 1:24)
    z <-try(sequenza.results(sequenza.extract = seqz_ext_data, 
                             sample.id = "alt_sol", 
                             out.dir=altsol_outdir,
                             cellularity=alt_sol[iix,"cellularity"],
                             ploidy=alt_sol[iix,"ploidy"]))
    if (is.error(z)) {
      stop(sprintf('failed to run sequenza.results with %d th alternative solution, cellularity=%f, ploidy=%f',iix, alt_sol[iix,"cellularity"],alt_sol[iix,"ploidy"])) 
    }
    
  }
  
  number_alt_sol <-nrow(alt_sol)
  return(number_alt_sol)
  
}

n_altsol <-.runSequenza_altSol_v1(outdir, prefix)

cat('\nDONE',prefix,'\n')
print(sessionInfo())

