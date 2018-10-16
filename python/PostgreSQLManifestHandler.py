"""
PostgreSQL database
Manifest handler
"""

from python.ManifestHandler import ManifestHandler
from python.glassfunc import build_dict
import psycopg2
import psycopg2.extras
        
class PostgreSQLManifestHandler(ManifestHandler):
    __conn = None

    def __init__(self, user, password, database, host = "localhost", port = 5432, source_file_basepath = "data/", aligned_file_basepath = "results/align/bqsr", from_source = False):
        try:
            self.__conn = psycopg2.connect(host = host, port = port, user = user, password = password, database = database)
        except(Exception, psycopg2.DatabaseError) as error:
            print(error)
            
        self.initFiles()
        self.initAliquots()
        self.initReadgroups()
        self.initPairs()
        self.initFilesReadgroups()
            
        ManifestHandler.__init__(self, source_file_basepath, aligned_file_basepath, from_source)

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

    def initFiles(self):
        q = "SELECT f.aliquot_barcode, f.file_name, f.file_format, f.file_path \
             FROM analysis.files AS f;"
        
        res = self.query(q, cursor_factory = psycopg2.extras.RealDictCursor)

        self.files = build_dict(res, "file_name")
    
    def initAliquots(self):    
        q = "SELECT al.aliquot_barcode, s.sample_barcode, s.case_barcode, s.sample_type, cl.case_project \
             FROM biospecimen.aliquots AS al \
                 INNER JOIN biospecimen.samples AS s ON al.sample_barcode = s.sample_barcode \
                 INNER JOIN clinical.cases AS cl ON cl.case_barcode = s.case_barcode;"
    
        res = self.query(q, cursor_factory = psycopg2.extras.RealDictCursor)

        self.aliquots = build_dict(res, "aliquot_barcode")
        
    def initReadgroups(self):
        q = "SELECT rg.aliquot_barcode, rg.readgroup_idtag, rg.readgroup_sample_id, rg.readgroup_idtag_legacy, \
                rg.readgroup_platform, rg.readgroup_platform_unit, rg.readgroup_library, rg.readgroup_center, \
                rg.readgroup_timestamp \
             FROM biospecimen.readgroups AS rg \
             ORDER BY rg.readgroup_idtag;"
        
        res = self.query(q, cursor_factory = psycopg2.extras.RealDictCursor)

        self.readgroups = res
    
    def initPairs(self):
        q = "SELECT pair_barcode, tumor_barcode, normal_barcode \
             FROM analysis.pairs;"
        
        res = self.query(q, cursor_factory = psycopg2.extras.RealDictCursor)
        
        self.pairs = build_dict(res, "pair_barcode")
        
    def initFilesReadgroups(self):
        q = "SELECT fr.file_name, fr.readgroup_idtag, fr.readgroup_sample_id, rg.aliquot_barcode, f.file_format \
             FROM analysis.files_readgroups AS fr \
                 LEFT JOIN biospecimen.readgroups AS rg ON fr.readgroup_idtag = rg.readgroup_idtag AND fr.readgroup_sample_id = rg.readgroup_sample_id \
                 INNER JOIN analysis.files AS f ON fr.file_name = f.file_name \
             ORDER BY fr.file_name;"
        
        res = self.query(q, cursor_factory = psycopg2.extras.RealDictCursor)
        
        self.files_readgroups = res

## END ##
