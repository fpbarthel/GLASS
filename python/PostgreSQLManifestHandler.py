"""
PostgreSQL database
Manifest handler
"""

from python.ManifestHandler import ManifestHandler
import psycopg2
import psycopg2.extras
        
class PostgreSQLManifestHandler(ManifestHandler):
    __conn = None

    def __init__(self, user, password, database, host = "localhost", port = 5432, source_file_basepath = "data/", aligned_file_basepath = "results/align/bqsr"):
        try:
            self.__conn = psycopg2.connect(host = host, port = port, user = user, password = password, database = database)
        except(Exception, psycopg2.DatabaseError) as error:
            print(error)
            
        ManifestHandler.__init__(self, source_file_basepath, aligned_file_basepath)

    def __del__(self):
        try:
            self.__conn.close()
        except(Exception, psycopg2.DatabaseError) as error:
            print(error)

    def getCursor(self, cursor_factory = None):
        if cursor_factory is not None:
            return self.__conn.cursor(cursor_factory = cursor_factory)
        else:
            return self.__conn.cursor()
    
    def query(self, q, d = None, cursor_factory = None):
        """
        Runs a query
        Parameters:
            q Query to run
            d Query parameters to be supplied with cur.execute()
        """
        res = None
        try:
            cur = self.getCursor(cursor_factory)
            cur.execute(q, d)
            if cursor_factory is psycopg2.extras.RealDictCursor:
                res = cur.fetchall()
            else:
                res = [s[0] for s in cur.fetchall()] ## Returns first item of each tuple in a list of tuples
        except (Exception, psycopg2.DatabaseError) as error:
            res = error
            print(error)
        finally:
            cur.close()
            return res

    def getSourceBAM(self, aliquot_barcode):
        """
        Returns a BAM filename given an aliquot barcode
        """
        
        q = "SELECT f.file_name \
             FROM analysis.files AS f \
             WHERE f.aliquot_barcode = %s AND f.file_format = 'uBAM';"
        
        return self.query(q, (aliquot_barcode,))

    def getAliquotsByCase(self, case_barcode):
        """
        Returns a list of aliquots given a case barcode
        """
        
        q = "SELECT al.aliquot_barcode \
             FROM biospecimen.aliquots AS al \
                 INNER JOIN biospecimen.samples AS s ON al.sample_barcode = s.sample_barcode \
             WHERE s.case_barcode = %s;"
    
        return self.query(q, (case_barcode,))  
    
    def getAliquotsByProject(self, case_project):
        """
        Returns a list of aliquots given a project name
        """
        
        q = "SELECT al.aliquot_barcode \
             FROM biospecimen.aliquots AS al \
                INNER JOIN biospecimen.samples AS s ON al.sample_barcode = s.sample_barcode \
                INNER JOIN clinical.cases AS cl ON cl.case_barcode = s.case_barcode \
             WHERE cl.case_project = %s;"
    
        return self.query(q, (case_project,))  

    def getAllAliquots(self):
        """
        Returns a list of all aliquots
        """
        
        q = "SELECT aliquot_barcode \
             FROM biospecimen.aliquots;"
    
        return self.query(q)

    def getPONAliquots(self):
        """
        Returns a list of all aliquots
        """
        
        q = "SELECT al.aliquot_barcode \
             FROM biospecimen.aliquots AS al \
                INNER JOIN biospecimen.samples AS s ON al.sample_barcode = s.sample_barcode \
             WHERE s.sample_type = 'NB';"
    
        return self.query(q) 
    
    def getRGIDs(self, aliquot_barcode):
        """
        Returns a list of RGIDs given an aliquot barcode
        """
        q = "SELECT rg.readgroup_idtag \
             FROM biospecimen.readgroups AS rg \
             WHERE rg.aliquot_barcode = %s \
             ORDER BY rg.readgroup_idtag;"
        
        return self.query(q, (aliquot_barcode,))

    def getRGIDsNotInAliquot(self, aliquot_barcode):
        """
        Returns a list of RGIDs NOT in a given aliquot barcode
        """
        q = "SELECT rg.readgroup_idtag \
             FROM biospecimen.readgroups AS rg \
                INNER JOIN analysis.files AS f ON f.aliquot_barcode = rg.aliquot_barcode \
             WHERE rg.aliquot_barcode != %s AND f.file_format = 'uBAM' \
             ORDER BY rg.readgroup_idtag;"
        
        return self.query(q, (aliquot_barcode,))

    def getLegacyRGIDs(self, aliquot_barcode):  
        """
        Returns a list of legacy RGIDs given an aliquot barcode
        """
        q = "SELECT rg.readgroup_idtag_legacy \
             FROM biospecimen.readgroups AS rg \
             WHERE rg.aliquot_barcode = %s \
             ORDER BY rg.readgroup_idtag;"
        
        return self.query(q, (aliquot_barcode,))

    def getFASTQ(self, aliquot_barcode, readgroup_idtag):
        """
        Returns a list of FASTQ filenames given an aliquot barcode and RGID tag
        """
        q = "SELECT fr.file_name \
             FROM analysis.files_readgroups AS fr \
                 INNER JOIN biospecimen.readgroups AS rg ON fr.readgroup_idtag = rg.readgroup_idtag AND fr.readgroup_sample_id = rg.readgroup_sample_id \
             WHERE rg.aliquot_barcode = %s AND rg.readgroup_idtag = %s \
             ORDER BY fr.file_name;"
        
        return self.query(q, (aliquot_barcode, readgroup_idtag))

    def getRGTag(self, aliquot_barcode, readgroup_idtag, tag):
        """
        Returns a given RG tag for a given aliquot barcode and readgroup
        """
        q = "SELECT * \
             FROM biospecimen.readgroups AS rg \
             WHERE rg.aliquot_barcode = %s AND rg.readgroup_idtag = %s"

        try:
            return self.query(q, (aliquot_barcode, readgroup_idtag), cursor_factory = psycopg2.extras.RealDictCursor)[0][tag]
        except Exception as error:
            return error

    def getAllReadgroups(self, limit_bam = False):
        """
        Returns all readgroups in BAM files
        """
        q = "SELECT DISTINCT rg.readgroup_idtag \
             FROM biospecimen.readgroups AS rg \
                INNER JOIN analysis.files AS f ON f.aliquot_barcode = rg.aliquot_barcode"

        if limit_bam:
            q += " WHERE f.file_format = 'uBAM';"
        else:
            q += ";"

        return self.query(q)

    def getTumor(self, pair_id):
        """
        Returns a tumor aliquot ID given a pair ID
        """
        q = "SELECT pa.tumor_barcode \
             FROM analysis.pairs AS pa \
             WHERE pa.pair_barcode = %s;"

        return self.query(q, (pair_id, ))

    def getNormal(self, pair_id):
        """
        Returns a normal aliquot ID given a pair ID
        """
        q = "SELECT pa.normal_barcode \
             FROM analysis.pairs AS pa \
             WHERE pa.pair_barcode = %s;"

        return self.query(q, (pair_id, ))

    def getAllFiles(self):
        """
        Returns a list of dictionaries containing all files and aliquots
        """
        q = "SELECT aliquot_barcode, file_name, file_format \
             FROM analysis.files;"
        
        return self.query(q, cursor_factory = psycopg2.extras.RealDictCursor)
    
    def getAllPairs(self):
        """
        Returns a list of dictionaries containing all pairs
        """
        q = "SELECT pair_barcode, tumor_barcode, normal_barcode \
             FROM analysis.pairs;"
        
        return self.query(q, cursor_factory = psycopg2.extras.RealDictCursor)

    def getAllReadgroupsByAliquot(self):
        """
        Returns a dictionary with aliquots as keys and a list of rgids as value
        """
        q = "SELECT aliquot_barcode, readgroup_idtag \
             FROM biospecimen.readgroups;"
        
        res = self.query(q, cursor_factory = psycopg2.extras.RealDictCursor)

        d = {}
        for i in res:
            if i["aliquot_barcode"] in d:
                d[i["aliquot_barcode"]].append(i["readgroup_idtag"])
            else:
                d[i["aliquot_barcode"]] = [i["readgroup_idtag"]]

        return d

## END ##
