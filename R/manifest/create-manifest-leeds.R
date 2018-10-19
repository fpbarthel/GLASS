#######################################################
# Create manifest for GLASS LEEDS UK - samples
# Date: 2018.10.19
# Author: Kevin J.
#######################################################
# Local directory for github repo.
mybasedir = "/Users/johnsk/Documents/"
setwd(mybasedir)

# Files with information about fastq information and barcodes.
life_history_barcodes = "/Users/johnsk/Documents/Life-History/GLASS-WG/data/ref/glass_wg_aliquots_mapping_table.txt"
wxs_clinical = "/Users/johnsk/Documents/GLASS/data/manifest/glass_analysis-paired_bam_map_wes.clinic.txt"
wxs_bam_map = "/Users/johnsk/Documents/GLASS/data/manifest/glass_analysis-paired_bam_map_wes.variant.txt"
# wxs_bam_map = "/Users/johnsk/Documents/tmp/glass_analysis-paired_bam_map_wes.variant_20181017.txt"
glass_wxs_barcodes = "/Users/johnsk/Documents/tmp/glass_wxs_barcodes.txt"
wxs_revised_readgroups = "/Users/johnsk/Documents/tmp/glass_wxs_readgroups_revised_rgid.txt"

#######################################################