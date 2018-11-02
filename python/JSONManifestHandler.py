"""
JSON file
Manifest handler
"""

from python.ManifestHandler import ManifestHandler
import psycopg2
import psycopg2.extras
        
class JSONManifestHandler(ManifestHandler):
    
    def __init__(self, cases_json, samples_json, aliquots_json, files_json, readgroups_json, pairs_json):
        self.cases      = json.load(open(cases_json))
        self.samples    = json.load(open(samples_json))
        self.aliquots   = json.load(open(aliquots_json))
        self.files      = json.load(open(files_json))
        self.readgroups = json.load(open(readgroups_json))
        self.pairs      = json.load(open(pairs_json))

        ## CASES -> DICT
        self.cases_dict = build_dict(self.cases, "case_id")

        ## SAMPLES -> DICT
        self.samples_dict = build_dict(self.samples, "sample_id")

        ## ALIQUOTS -> DICT
        self.aliquots_dict = build_dict(self.aliquots, "aliquot_id")

        ## FILES -> DICT
        self.files_dict = build_dict(self.files, "file_uuid")

        ## Pair IDs are unique, PAIRS -> DICT
        self.pairs_dict = build_dict(self.pairs, "pair_id")

        ## Aliquot IDs and BAM files map 1:1
        self.ALIQUOT_TO_BAM_PATH = {}
        
        for file in self.files:
            if file["file_format"] == "BAM":
                self.ALIQUOT_TO_BAM_PATH[ file["aliquot_id"] ] = file["file_path"]
                
        ## Case to aliquots
        ## Dict of aliquots per case
        self.CASE_TO_ALIQUOT = {}
        for aliquot in ALIQUOTS:
            aliquot["case_id"] = self.samples_dict[ aliquot["sample_id"] ]["case_id"]  
            if aliquot["case_id"] not in self.CASE_TO_ALIQUOT:
                self.CASE_TO_ALIQUOT[ aliquot["case_id"] ] = [ aliquot["aliquot_id"] ]
            elif aliquot["aliquot_id"] not in self.CASE_TO_ALIQUOT[ aliquot["case_id"] ]:
                self.CASE_TO_ALIQUOT[ aliquot["case_id"] ].append(aliquot["aliquot_id"])

        ## Aliquots and RGIDs map 1:many
        self.ALIQUOT_TO_RGID = {}     
        self.ALIQUOT_TO_LEGACY_RGID = {}
        for readgroup in self.readgroups:
            if readgroup["aliquot_id"] not in self.ALIQUOT_TO_RGID:
                self.ALIQUOT_TO_RGID[ readgroup["aliquot_id"] ] = [ readgroup["readgroup_id"] ]
            else:
                self.ALIQUOT_TO_RGID[ readgroup["aliquot_id"] ].append(readgroup["readgroup_id"])
            if "legacy_readgroup_id" not in readgroup or len(readgroup["legacy_readgroup_id"]) == 0:
                continue
            if readgroup["aliquot_id"] not in self.ALIQUOT_TO_LEGACY_RGID:
                self.ALIQUOT_TO_LEGACY_RGID[ readgroup["aliquot_id"] ] = [ readgroup["legacy_readgroup_id"] ]
            else:
                self.ALIQUOT_TO_LEGACY_RGID[ readgroup["aliquot_id"] ].append(readgroup["legacy_readgroup_id"])

        ## Readgroup information and 
        ## Aliquots and RGIDs map 1:many
        ## RGIDs are unique within an aliquot
        ## Aliquot IDs and fastQ files map 1:many
        ## Because FQ files are also seperated by readgroup, create dictionary of FQ files here as well
        self.ALIQUOT_TO_READGROUP = {} 
        self.ALIQUOT_TO_FQ_PATH = {}
        for readgroup in self.readgroups:
            if readgroup["aliquot_id"] not in self.ALIQUOT_TO_READGROUP:
                self.ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ] = { readgroup["readgroup_id"] : readgroup }
            else:
                self.ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["readgroup_id"] ] = readgroup
            self.ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["readgroup_id"] ]["file_path"] = self.files_dict[ self.ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["readgroup_id"] ]["file_uuid"] ]["file_path"]
            self.ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["readgroup_id"] ]["file_format"] = self.files_dict[ self.ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["readgroup_id"] ]["file_uuid"] ]["file_format"]
            if self.ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["readgroup_id"] ]["file_format"] == "FQ":
                if readgroup["aliquot_id"] not in self.ALIQUOT_TO_FQ_PATH:
                    self.ALIQUOT_TO_FQ_PATH[ readgroup["aliquot_id"] ] = {}
                self.ALIQUOT_TO_FQ_PATH[ readgroup["aliquot_id"] ][ readgroup["readgroup_id"] ] = self.ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["readgroup_id"] ]["file_path"].split(",")

    ## IMPLEMENTATION PENDING

## END ##
