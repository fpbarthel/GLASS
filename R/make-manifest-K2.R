
cases_file      = "data/manifest/K2/cases"
samples_file    = "data/manifest/K2/samples"
aliquots_file   = "data/manifest/K2/aliquots"
readgroups_file = "data/manifest/K2/readgroups"
files_file      = "data/manifest/K2/files"
pairs_file      = "data/manifest/K2/pairs"

json_ext = "json"
text_ext = "tsv"

rginfo = read.delim("/projects/verhaak-lab/metastatic_gliosarcoma/Sample_Identifiers.txt", as.is=T)
fqfiles = list.files("/projects/verhaak-lab/metastatic_gliosarcoma/N06057_MK_MJ1801151/raw_data", pattern = "[1-2].fq.gz$", full.names = T)

## Create random IDs for these samples
set.seed(27072018)
uuids = stringi::stri_rand_strings(5, 6, "[A-Z0-9]")
## "WQTP4O" "NYD1BJ" "PA9D9P" "SW2KKR" "LY5QYK"

files = data.frame(aliquot_id = sprintf("GLSS-K2-0001-%s-%s", c("TP", "R1", "R2", "TM", "NB"), uuid = uuids),
                   file_path = paste(fqfiles[seq(1,9,2)], fqfiles[seq(2,10,2)], sep=","),
                   file_name = paste(basename(fqfiles)[seq(1,9,2)], basename(fqfiles)[seq(2,10,2)], sep=","),
                   file_uuid = sprintf("file-%s", uuids),
                   file_size = NA,
                   file_md5sum = NA,
                   file_format = "FQ",
                   stringsAsFactors = F)

readgroups = data.frame(file_uuid = sprintf("file-%s", uuids),
                        aliquot_id = sprintf("GLSS-K2-0001-%s-%s", c("TP", "R1", "R2", "TM", "NB"), uuid = uuids),
                        RGID = sprintf("%s.%s", substr(rginfo$Flowcell.ID[seq(1,9,2)],1,4), rginfo$Flowcell.Lane[seq(1,9,2)]),
                        RGPL = "ILLUMINA",
                        RGPU = sprintf("%s.%s", rginfo$Flowcell.ID[seq(1,9,2)], rginfo$Flowcell.Lane[seq(1,9,2)]),
                        RGLB = uuids,
                        RGPI = 0,
                        RGDT = NA,
                        RGSM = sprintf("GLSS-K2-0001-%s", c("TP", "R1", "R2", "TM", "NB")),
                        RGCN = "GENEWIZ",
                        stringsAsFactors = F)


### aliquots
aliquots = data.frame(sample_id = sprintf("GLSS-K2-0001-%s", c("TP", "R1", "R2", "TM", "NB")),
                      aliquot_uuid = uuids,
                      aliquot_id = sprintf("GLSS-K2-0001-%s-%s", c("TP", "R1", "R2", "TM", "NB"), uuid = uuids),
                      portion = 1,
                      analyte_type = "DNA",
                      analysis_type = "WGS",
                      stringsAsFactors = F)

## Samples
samples = data.frame(case_id = "GLSS-K2-0001",
                     sample_id = sprintf("GLSS-K2-0001-%s", c("TP", "R1", "R2", "TM", "NB")),
                     legacy_sample_id = gsub(".fq.gz", "", basename(fqfiles)[seq(1,9,2)]),
                     sample_type = c("TP", "R1", "R2", "TM", "NB"),
                     stringsAsFactors = F)

### Cases
cases = data.frame(case_id = "GLSS-K2-0001",
                   project_id = "K2",
                   stringsAsFactors = F)

## Pairs
pairs = data.frame(case_id = "GLSS-K2-0001",
                   pair_id = aliquots$aliquot_id[1:4],
                   tumor_aliquot_id = aliquots$aliquot_id[1:4],
                   normal_aliquot_id = "GLSS-K2-0001-NB-LY5QYK")

## Make manifest
