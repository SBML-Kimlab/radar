#import subprocess
from radar import *
#from sero import dir_wgs
import radar
class prodigal : 
    def __init__( self, strain ) :
        self.strain = strain
        pass
    @staticmethod
    def prodigal_sc1( strain ) :
        dir_wgs = radar.dir_wgs
        dir_genome = radar.dir_genome
        file_wgs = dir_wgs + "*.fna" 
        lst_fna = glob.glob( file_wgs ) 
        lst_fna.sort()
        if lst_fna == [] :
            print ( "\n**********\nFatal Error : No Genome sequence File ...\n on %s\n**********\n" % dir_wgs )     count_fna = 0 
        for file_wgs in lst_fna :
            count_fna += 1
            fna_id = file_wgs.split( "/" )[ -1 ].replace( ".fna", "" )
            file_gbk = dir_genome + "%s.gbk" % fna_id 
            file_faa = file_gbk.replace( ".gbk", ".faa" )
            file_fna = file_gbk.replace( ".gbk", ".fna" )
            if os.path.exists( file_fna ) == True and os.path.getsize( file_fna ) >= 100 * 1024 :
                continue
            if count_fna % 10 == 0 :
                clear_output()
            print ( "Processing %s" %file_wgs )
            print ( "Processing %d wgs files" % count_fna )
            command = "prodigal -i %s -o %s -a %s -d %s" %( file_wgs, file_gbk, file_faa, file_fna ) 
            (exitstatus, outtext) = subprocess.getstatusoutput( "%s" % command )
            print ( outtext )