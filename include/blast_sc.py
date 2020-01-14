#import subprocess
from radar import *
#from sero import dir_wgs
import radar
class blastp : 
    def __init__( self, strain ) :
        global dir_db, program_dir, dir_genome, dir_blast, dir_blastp, dir_user, dir_result, file_ntdb, file_aadb, file_udb, file_usearch
        self.strain = strain
        dir_db = radar.db_dir
        program_dir = radar.program_dir
        dir_genome = radar.dir_genome
        dir_blast = radar.dir_blast
        dir_blastp = radar.dir_blastp
        dir_user = radar.dir_user
        dir_result = radar.dir_result
        file_ntdb = dir_db + "radar.fasta"
        file_aadb = dir_db + "radar.faa"
        file_udb = dir_db + "radar.udb"     
        file_usearch = program_dir + "usearch" 
        pass
    @staticmethod
    def make_udb( strain ) :
        #print dir_db
        #print program_dir
        #print dir_genome
        #print file_ntdb
        #print file_aadb
        #print file_udb
        #print file_usearch
        if os.path.isfile( file_usearch ) == False :
            print "Fatal error : Cannot find usearch on ... %s" % file_usearch
        if os.path.isfile( file_aadb ) == False :
            print "Fatal error : Cannot find the peptide database on ... %s" % dir_db
        else :
            command = "%s -makeudb_ublast %s -output %s" %( file_usearch, file_aadb, file_udb ) 
            print command
            (exitstatus, outtext) = commands.getstatusoutput( "%s" % command )
            print outtext
            print "Build Database completed...%s" % file_udb

    @staticmethod
    def blastp_run( strain ) :
        file_faa = dir_genome + "*.faa" 
        lst_faa = glob.glob( file_faa )
        lst_faa.sort()

        count_faa = 0 
        for file_faa in lst_faa :
            count_faa += 1
            faa_id = file_faa.split( "/" )[ -1 ].replace( ".faa", "" )
            file_out = dir_blastp + "%s.aln" % faa_id
            if os.path.exists( file_out ) == True and os.path.getsize( file_out ) > 10 * 1024 :
                continue
            if count_faa % 10 == 0 :
                clear_output()
            print "Processing %s" % file_faa
            print "Processing %d annotation files" % count_faa
            command = "%s -ublast %s -db %s -evalue 1e-9 -alnout %s" %( file_usearch, file_faa, file_udb, file_out )
            print command
            (exitstatus, outtext) = commands.getstatusoutput( "%s" % command )
            print outtext

    @staticmethod
    def blastp_parse( strain ) :
        def read_as_string( file0 ) :
            with open( file0, "r" ) as f :
                read = f.read()
            return read
        #print dir_user
        lst_aln = glob.glob( dir_blastp + "*.aln" )
        lst_aln.sort()
        count_file = 0 
        for file_aln in lst_aln :
            count_file += 1
            if count_file % 100 == 0 :
                clear_output()
                print "Processed %d files." % count_file
            file_tsv = file_aln.replace( ".aln", ".tsv" )
            if os.path.exists( file_tsv ) == True :
                continue
            output = read_as_string( file_aln )
            output = output.split( "\nQuery >" )
            with open( file_tsv, "w" ) as f :
                for output0 in output :
                    #if output0.startswith( "%s" %dir_user ) :
                    if output0.startswith( "%s" %file_usearch ) :
                        continue
                    output1 = output0.split( "\n" )
                    output2 = output1[ 2 ].split( )
                    gene_from = output1[ 0 ].split( " " )[ 0 ]
                    gene_to = output2[ 5 ]
                    score = output2[ 0 ]
                    evalue = output2[ 1 ]
                    identity = output2[ 2 ]

                    to_print = [ gene_from, gene_to, score, evalue, identity ]
                    to_print = [ str( m ) for m in to_print ]
                    to_print = ",".join( to_print )
                    f.write( "%s\r\n" % to_print )
            print "Wrote... %s" %file_tsv

    @staticmethod
    def merge_after_parse( strain ) :
        def file_readlines( file0 ) :
            f = open( file0, "r" )
            lines = f.read().splitlines()
            f.close()
            return lines
        print dir_blast
        count_file = 0 
        lst_identity = []
        for dir_ in os.listdir( dir_blast ) :
            file_i = dir_blast + dir_
            file_w = glob.glob( file_i + "/*.tsv" )
            if count_file == 0 :
                lst_files = list( file_w )
                count_file += 1 
            elif count_file > 0 :
                lst_files.extend( file_w )
                count_file += 1
        print ".tsv files: %d" % len( lst_files )
        for file0 in lst_files :
            lines = file_readlines( file0 )
            lst_identity0 = [ float( ( line.split( "," )[ 4 ] ).strip( "%" ) ) / 100 for line in lines ]
            lst_identity.extend( lst_identity0 )
        print len( lst_identity )

        dir_figure_wgs = dir_result + "wgs/figure/"
        dir_table_wgs = dir_result + "wgs/table/"
        os.system( "mkdir -p %s" % dir_figure_wgs )
        os.system( "mkdir -p %s" % dir_table_wgs )
        fig, ax = plt.subplots( 1, 1, figsize = ( 6, 6 ) )
        plt.hist( lst_identity, 25, normed = 1, facecolor = "green", alpha = 0.75, rwidth=0.9)
        ax.set_xlabel( "identintity" )
        ax.set_ylabel( "100K" )
        plt.draw()
        file_figure = dir_figure_wgs + "wgs_search_group.png" 
        #plt.axis( [ 0.1, 1.0, 0, 9 ] )
        #plt.grid(True)
        plt.savefig( file_figure )
        plt.savefig( file_figure.replace( ".png", ".pdf" ) )
        plt.show()
        
        file_blastp_merged = dir_table_wgs + "1+2_blastp_merged.tsv"
        file_blastp_merged0 = file_blastp_merged.replace( ".tsv", "0.tsv" )
        file_blastp_merged_esbl = file_blastp_merged.replace( ".tsv", "_esbl.tsv" )
        lst_esbl = [ "TEM-", "SHV-", "CTX-", "OXA-", "GES-", "PER-", "MCR-", "NDM-", "KPC-" ]
        dic_esbl = {}
        with open( file_blastp_merged, "w" ) as f :
            for file0 in lst_files :
                genome_id = file0.split( "/" )[ -1 ].replace( ".tsv", "" )
                if "_" in genome_id :
                    genome_id = genome_id.split( "_" )
                    genome_id = "%s_%s" % ( genome_id[ 0 ], genome_id[ 1 ] )
                lines = file_readlines( file0 )    
                for line in lines :
                    line0 = line.split( "," )
                    identity = float( line0[ 4 ].strip( "%" ) ) / 100
                    if identity >= 0.95 :
                        to_print = [ genome_id, line0[ 0 ], line0[ 1 ], line0[ 2 ], line0[ 3 ], identity ]
                        to_print = [ str( m ) for m in to_print ]
                        to_print = "\t".join( to_print )
                        f.write( "%s\r\n" % to_print )
        print "Wrote... %s" % file_blastp_merged
        
        with open( file_blastp_merged0, "w" ) as f, open( file_blastp_merged_esbl, "w" ) as f0 :
            for line in file_readlines( file_blastp_merged ) :
                line0 = line.split( "\t" )
                to_print = list( line0 )
                gene_to = line0[ 2 ]
                to_print = [ str( m ) for m in to_print ]
                to_print = "\t".join( to_print )
                is_esbl = False
                for esbl in lst_esbl :
                    if esbl in gene_to :
                        if esbl not in dic_esbl :
                            dic_esbl[ esbl ] = 1
                        else :
                            dic_esbl[ esbl ] += 1
                        is_esbl = True
                        break
                if is_esbl == True :
                    f0.write( "%s\r\n" % to_print )
                f.write( "%s\r\n" % to_print )
        print "Wrote... %s" % file_blastp_merged0
        print "Wrote... %s" % file_blastp_merged_esbl
        print dic_esbl