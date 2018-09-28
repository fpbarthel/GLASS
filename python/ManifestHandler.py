"""
ManifestHandler
"""

from python.glassfunc import locate

class ManifestHandler:
    """
    Class that defines the GLASS manifest
    It should not be implemented directly, but either as a JSONManifest using 
    JSON text files or a PostgreSQLManifest using a database connection
    """
    
    source_files = []
    aligned_files = []

    selected_aliquots = set()
    selected_pairs = set()
    selected_readgroups_by_aliquot = {}
    
    def __init__(self, source_file_basepath, aligned_file_basepath):
        self.source_files = self.locateSourceFiles(source_file_basepath)
        self.aligned_files = self.locateAlignedBAMFiles(aligned_file_basepath)
        
        self.selected_aliquots.update([i["aliquot_barcode"] for i in self.source_files if len(i["file_path"]) > 0])
        self.selected_aliquots.update([i["aliquot_barcode"] for i in self.aligned_files if len(i["file_path"]) > 0])
        
        for pair in self.getAllPairs():
            if(pair["tumor_barcode"] in self.selected_aliquots and pair["normal_barcode"] in self.selected_aliquots):
                self.selected_pairs.add(pair["pair_barcode"])

        for aliquot, readgroups in self.getAllReadgroupsByAliquot().items():
            if aliquot in self.selected_aliquots:
                self.selected_readgroups_by_aliquot[aliquot] = readgroups
        
    def __str__(self):
        n_source_fastq = len([j["file_path"] for j in self.source_files if j["file_format"] == "FASTQ"]) 
        n_source_fastq_found = n_source_fastq - [j["file_path"] for j in self.source_files if j["file_format"] == "FASTQ"].count([])
        
        n_source_bam = len([j["file_path"] for j in self.source_files if j["file_format"] == "uBAM"])
        n_source_bam_found = n_source_bam - [j["file_path"] for j in self.source_files if j["file_format"] == "uBAM"].count([])
        
        n_aligned_bam = len(self.aligned_files)
        n_aligned_bam_found = n_aligned_bam - [j["file_path"] for j in self.aligned_files].count([])
        
        n_aliquots = len(self.getAllAliquots())
        n_aliquots_selected = len(self.selected_aliquots)
        
        n_pairs = len(self.getAllPairs())
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
        Return a list of selected pairs
        """
        return self.selected_readgroups_by_aliquot

    def locateSourceFiles(self, source_file_basepath):
        """
        Locate raw/unaligned FASTQ/BAM files
        """
        
        sources_files = self.getAllFiles()
        for file in sources_files:
            file["file_path"] = [f for f in locate(file["file_name"], source_file_basepath)]
            
        return sources_files
        
    def locateAlignedBAMFiles(self, aligned_file_basepath):
        """
        Locate re-aligned BAM files
        """
        
        aligned_files = []
        for barcode in self.getAllAliquots():
            aligned_files.append({"aliquot_barcode" : barcode,
                                  "file_name" : "{}.realn.mdup.bqsr.bam".format(barcode),
                                  "file_format" : "aligned BAM",
                                  "file_path" : [f for f in locate("{}.realn.mdup.bqsr.bam".format(barcode), aligned_file_basepath)]})
        
        return aligned_files
    
    def getAllFiles(self):
        raise NotImplementedError("ManifestHandler should not be implemented directly")
        
    def getAllAliquots(self):
        raise NotImplementedError("ManifestHandler should not be implemented directly")
        
    def getAllPairs(self):
        raise NotImplementedError("ManifestHandler should not be implemented directly")

    def getAllReadgroupsByAliquot(self):
        raise NotImplementedError("ManifestHandler should not be implemented directly")

## END ##