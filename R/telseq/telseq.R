library(data.table)
library(dplyr)
library(tidyr)
library(parallel)

## Constants/parameters
K = 7
G = 332720800
c = 46000

resultsdir  = "results/telseq"
ptrn        = "telseq.txt$"

##################################################################################################################################################################

message("Merging telseq output")

# list of telseq data files
tsfiles = list.files(resultsdir, pattern=ptrn, recursive=T, full.names=T)

tsdat = mclapply(tsfiles, function(fn) {
  aliquot_barcode = unlist(strsplit(basename(fn),"\\."))[1]
  
  f = tryCatch(data.table::fread(fn), error=function(e) e)
  
  if(inherits(f, "error"))
    message(fn)
  
  f = data.frame(aliquot_barcode,f)
  
  total_reads = as.numeric(as.character(f$Total))
  total_reads_wm = round(sum(total_reads,na.rm=T))
  
  mapped_reads = as.numeric(as.character(f$Mapped))
  mapped_reads_wm = round(sum(mapped_reads,na.rm=T))
  
  duplicate_reads = as.numeric(as.character(f$Duplicates))
  duplicate_reads_wm = round(sum(duplicate_reads,na.rm=T))
  
  tel = apply(f[,na.omit(match(paste('TEL',K:99,sep=''),colnames(f)))],1,sum,na.rm=T)
  tel_wm = weighted.mean(tel, total_reads, na.rm=T)
  
  gc = apply(f[,match(paste('GC',4:5,sep=''),colnames(f))],1,sum)
  gc_wm = weighted.mean(gc, total_reads, na.rm=T)
  
  len = (tel/gc)*(G/c)
  len_wm = weighted.mean(len, total_reads, na.rm=T)
  
  out = data.frame(aliquot_barcode,
                   total_reads=total_reads_wm,
                   mapped_reads=mapped_reads_wm,
                   duplicate_reads=duplicate_reads_wm,
                   tel=tel_wm,
                   K,
                   G,
                   c,
                   gc=gc_wm,
                   length=len_wm)
  
  return(out)
})

tsdat = rbindlist(tsdat) %>% as.data.frame()

con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")
DBI::dbWriteTable(con, DBI::Id(schema="analysis",table="telseq"), tsdat)
