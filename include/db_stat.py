#import subprocess
from radar import *
#from sero import dir_wgs
import radar
class db_statistics : 
    def __init__( self, strain, db ) :
        self.strain = strain
        self.db = db
        pass
    @staticmethod
    def db_stat( strain, db ) :
        dir_db = radar.db_dir
        db0 = db.replace( ".faa", "" )
        if db0 == "BOARDS" : 
            file_ntdb = dir_db + "BOARDS.fasta"; file_aadb = dir_db + "BOARDS.faa"; file_udb = dir_db + "BOARDS.udb"
            if os.path.isfile( file_ntdb ) == False :
                print ( "Fatal error : No nucleotide database on ... %s" % dir_db )
                print ( "Please Download BOARDS DB(nt) on website" )
            else :
                ntdb_open = open( file_ntdb )
                dic_ntdb = SeqIO.to_dict( SeqIO.parse( ntdb_open, "fasta" ) )
                print ( "Total Antibiotics resistance genes : %s" % len( dic_ntdb ) )
            if os.path.isfile( file_aadb ) == False :
                print ( "Fatal error : No peptide database on ... %s" % dir_db )
                print ( "Please Download RADAR DB(aa) on website" )
            else : print ( "The peptide database existed... %s" % file_aadb )
        if db0 == "USER_DB" : 
            file_aadb = dir_db + "USER_DB.faa"; file_udb = dir_db + "USER_DB.udb"
            if os.path.isfile( file_aadb ) == False :
                print ( "Fatal error : No user-defined database on ... %s" % dir_db )
                print ( "Please Upload USER-DEFINED DATABASE..." )
            else :
                aadb_open = open( file_aadb )
                dic_aadb = SeqIO.to_dict( SeqIO.parse( aadb_open, "fasta" ) )
                print ( "Total included genes: %s" % len( dic_aadb ) )