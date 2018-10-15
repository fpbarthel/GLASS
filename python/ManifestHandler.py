"""
ManifestHandler
"""

from python.glassfunc import locate

SOURCE_FASTQ_TYPE = "FASTQ"
SOURCE_BAM_TYPE = "uBAM"
ALIGNED_BAM_TYPE = "aligned BAM"

class ManifestHandler:
    """
    Class that defines the GLASS manifest
    It should not be implemented directly, but either as a JSONManifest using 
    JSON text files or a PostgreSQLManifest using a database connection
    """
    
    aliquots = {} ## Dictionary of aliquots, keys = aliquot_barcode
    readgroups = [] ## List of readgroups
 
    files = {} ## Dictionary of files, keys = file_name
    pairs = {} ## Dictionary of pairs, keys = pair_barcode
    
    files_readgroups = [] ## List of file to readgroup mappings

    selected_aliquots = set()
    selected_pairs = set()
    
    def __init__(self, source_file_basepath, aligned_file_basepath):
        
        self.locateFiles(source_file_basepath, aligned_file_basepath)
        
        self.selected_aliquots.update([f["aliquot_barcode"] for (file_name, f) in self.files.items() if len(f["file_path"]) > 0])

        for (pair_barcode, pair) in self.pairs.items():
            if(pair["tumor_barcode"] in self.selected_aliquots and pair["normal_barcode"] in self.selected_aliquots):
                self.selected_pairs.add(pair_barcode)
        
    def __str__(self):
        n_source_fastq = len([j["file_path"] for j in self.files.values() if j["file_format"] == SOURCE_FASTQ_TYPE]) 
        n_source_fastq_found = n_source_fastq - [j["file_path"] for j in self.files.values() if j["file_format"] == SOURCE_FASTQ_TYPE].count([])
        
        n_source_bam = len([j["file_path"] for j in self.files.values() if j["file_format"] == SOURCE_BAM_TYPE])
        n_source_bam_found = n_source_bam - [j["file_path"] for j in self.files.values() if j["file_format"] == SOURCE_BAM_TYPE].count([])
        
        n_aligned_bam = len([j["file_path"] for j in self.files.values() if j["file_format"] == ALIGNED_BAM_TYPE])
        n_aligned_bam_found = n_aligned_bam - [j["file_path"] for j in self.files.values() if j["file_format"] == ALIGNED_BAM_TYPE].count([])
        
        n_aliquots = len(self.aliquots)
        n_aliquots_selected = len(self.selected_aliquots)
        
        n_pairs = len(self.pairs)
        n_pairs_selected = len(self.selected_pairs)
        
        s = "Found {} of {} possible source FASTQ files.\n".format(n_source_fastq_found, n_source_fastq)
        s += "Found {} of {} possible source BAM files.\n".format(n_source_bam_found, n_source_bam)
        s += "Found {} of {} possible realigned BAM files.\n".format(n_aligned_bam_found, n_aligned_bam)
        s += "Selected {} of {} total aliquots.\n".format(n_aliquots_selected, n_aliquots)
        s += "Selected {} of {} possible pairs.\n".format(n_pairs_selected, n_pairs)
        return(s)
        
    def getSelectedAliquots(self):
        """
        Return a list of selected aliquots
        """
        return list(self.selected_aliquots)

    def getSelectedPairs(self):
        """
        Return a list of selected pairs
        """
        return list(self.selected_pairs)

    def getSelectedReadgroupsByAliquot(self):
        """
        Return a dictionary of readgroup_idtag (keys = aliquot_barcode) limited to selected aliquots
        """
        return {aliquot_barcode: self.getRGIDs(aliquot_barcode) for aliquot_barcode in self.getSelectedAliquots()}
    
    def locateFiles(self, source_file_basepath, aligned_file_basepath):
        """
        Locate FASTQ/BAM files
        """
        for (file_name, file) in self.files.items():
            if file["file_format"] == SOURCE_FASTQ_TYPE or file["file_format"] == SOURCE_BAM_TYPE:
                file["file_path"] = [f for f in locate(file_name, source_file_basepath)]
            elif file["file_format"] == ALIGNED_BAM_TYPE:
                file["file_path"] = [f for f in locate(file_name, aligned_file_basepath)]
    
    def initFiles(self):
        raise NotImplementedError("ManifestHandler should not be implemented directly")
        
    def initAliquots(self):
        raise NotImplementedError("ManifestHandler should not be implemented directly")
        
    def initReadgroups(self):
        raise NotImplementedError("ManifestHandler should not be implemented directly")
        
    def initPairs(self):
        raise NotImplementedError("ManifestHandler should not be implemented directly")
        
    def initFilesReadgroups(self):
        raise NotImplementedError("ManifestHandler should not be implemented directly")
        
    def getAllAliquots(self):
        """
        Returns a list of all aliquots
        """
        return list(self.aliquots.keys())
    
    def getSourceBAM(self, aliquot_barcode):
        """
        Returns a BAM filename (str) given an aliquot barcode (str)
        """     
        return [file_name for (file_name, f) in self.files.items() if f["aliquot_barcode"] == aliquot_barcode][0]

    def getFileFormatByAliquot(self, aliquot_barcode):
        """
        Returns a file format type (str) given an aliquot barcode (str)
        """     
        return [f["file_format"] for (file_name, f) in self.files.items() if f["aliquot_barcode"] == aliquot_barcode][0]
    
    def getAliquotsByCase(self, case_barcode):
        """
        Returns a list of aliquots given a case barcode
        """
        return [aliquot_barcode for (aliquot_barcode, al) in self.aliquots.items() if al["case_barcode"] == case_barcode]
    
    def getAliquotsByProject(self, case_project):
        """
        Returns a list of aliquots given a project name
        """
        return [aliquot_barcode for (aliquot_barcode, al) in self.aliquots.items() if al["case_project"] == case_project]

    def getPONAliquots(self):
        """
        Returns a list of all aliquots with sample_type = 'NB'
        """
        return [aliquot_barcode for (aliquot_barcode, al) in self.aliquots.items() if al["sample_type"] == "NB"]
    
    def getRGIDs(self, aliquot_barcode):
        """
        Returns a list of RGIDs given an aliquot barcode
        """
        return [rg["readgroup_idtag"] for rg in self.readgroups if rg["aliquot_barcode"] == aliquot_barcode]

    def getRGSampleTagByAliquot(self, aliquot_barcode):
        """
        Returns a list of RGIDs given an aliquot barcode
        """ 
        return [rg["readgroup_sample_id"] for rg in self.readgroups if rg["aliquot_barcode"] == aliquot_barcode][0]
    
    def getLegacyRGIDs(self, aliquot_barcode):  
        """
        Returns a list of legacy RGIDs given an aliquot barcode
        """       
        return [rg["readgroup_idtag_legacy"] for rg in self.readgroups if rg["aliquot_barcode"] == aliquot_barcode]

    def getRGTag(self, aliquot_barcode, readgroup_idtag, tag):
        """
        Returns the value (str) of a given RG tag (str) for a given aliquot barcode (str) and readgroup (str)
        """
        return [rg[tag] for rg in self.readgroups if (rg["aliquot_barcode"] == aliquot_barcode and rg["readgroup_idtag"] == readgroup_idtag)][0]
    
    def getAllReadgroups(self, limit_bam = False):
        """
        Returns all readgroups in BAM files
        """
        s = set()
        for rg in self.readgroups:
            if not limit_bam or (limit_bam and self.getFileFormatByAliquot(rg["aliquot_barcode"]) == "uBAM"):
                s.add(rg["readgroup_idtag"])
        return sorted(s)
    
    def getTumor(self, pair_barcode):
        """
        Returns a tumor aliquot ID given a pair ID
        """
        return [p["tumor_barcode"] for (barcode, p) in self.pairs.items() if barcode == pair_barcode][0]

    def getNormal(self, pair_id):
        """
        Returns a normal aliquot ID given a pair ID
        """
        return [p["normal_barcode"] for (barcode, p) in self.pairs.items() if barcode == pair_barcode][0]

    def getFiles(self):
        """
        Returns a dictionary (keys = file_name) of dictionaries containing all files and aliquots
        """     
        return self.files
    
    def getPairs(self):
        """
        Returns a dictionary (keys = pair_barcode) of dictionaries containing all pairs
        """
        return self.pairs
    
    def getRGIDsNotInAliquot(self, aliquot_barcode):
        """
        Returns a list of RGIDs NOT in a given aliquot barcode
        """
        s = set([fr["readgroup_idtag"] for fr in self.files_readgroups if fr["aliquot_barcode"] != aliquot_barcode])
        return sorted(s)

    def getFASTQ(self, aliquot_barcode, readgroup_idtag):
        """
        Returns a list of FASTQ filenames given an aliquot barcode and RGID tag
        """
        return [fr["file_name"] for fr in self.files_readgroups if fr["aliquot_barcode"] == aliquot_barcode and fr["readgroup_idtag"] == readgroup_idtag]

## END ##
