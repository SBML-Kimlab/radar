#import subprocess
from radar import *
#from sero import dir_wgs
import radar
class db_statistics : 
    def __init__( self, strain ) :
        self.strain = strain
        pass
    @staticmethod
    def db_stat( strain ) :
        dir_db = radar.db_dir
        file_ntdb = dir_db + "radar.fasta"
        file_aadb = dir_db + "radar.faa"
        file_udb = dir_db + "radar.udb"
        #print dir_db
        if os.path.isfile( file_ntdb ) == False :
            print ( "Fatal error : No nucleotide database on ... %s" % dir_db )
            print ( "Please Download RADAR DB(nt) on website" )
        else :
            ntdb_open = open( file_ntdb )
            dic_ntdb = SeqIO.to_dict( SeqIO.parse( ntdb_open, "fasta" ) )
            print ( "Total Antibiotics resistance genes : %s" % len( dic_ntdb ) )
        if os.path.isfile( file_aadb ) == False :
            print ( "Fatal error : No peptide database on ... %s" % dir_db )
            print ( "Please Download RADAR DB(aa) on website" )
        else : print ( "The peptide database existed... %s" % file_aadb )