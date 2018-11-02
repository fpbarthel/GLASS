## Add aligned files

library(odbc)
library(DBI)
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

bamfiles = list.files("/fastscratch/verhaak-lab/GLASS-C7/results/align/bqsr/", pattern = "bam$", full.names = T)
md5files = list.files("/fastscratch/verhaak-lab/GLASS-C7/results/align/bqsr/", pattern = "md5$", full.names = T)
filesizes = sapply(bamfiles, function(f) file.info(f)$size)
filemd5s  = sapply(md5files, function(f) readLines(f, warn=F))

files_add = data.frame(aliquot_barcode = gsub(".realn.mdup.bqsr.bam", "", basename(bamfiles)), file_name = basename(bamfiles), file_size = unname(filesizes), file_md5sum= unname(filemd5s), file_format = "aligned BAM", stringsAsFactors = F)
                       

dbWriteTable(con, Id(schema="analysis",table="files"), files_add, append=T)

tmp = dbReadTable(con, Id(schema="biospecimen",table="aliquots"))
write.csv(tmp, file = "aliquots.csv")
