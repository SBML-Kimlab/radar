#import subprocess
from radar import *
#from sero import dir_wgs
import circ_contents
import radar
class circos : 
    def __init__( self, strain ) :
        global dir_vis, dir_wgs, dir_genome, dir_blastp, dir_circos, dir_user
        global file_circos
    	#global dir_db, program_dir, dir_genome, dir_blast, dir_blastp, dir_user, dir_result, file_ntdb, file_aadb, file_udb, file_usearch
    	self.strain = strain
        dir_vis = radar.dir_vis
        dir_wgs = radar.dir_wgs
        dir_genome = radar.dir_genome
        dir_blastp = radar.dir_blastp
        dir_circos = radar.dir_circos
        dir_user = radar.dir_user
        file_circos = radar.file_circos
    	#dir_db = radar.db_dir
    	#program_dir = radar.program_dir
    	#dir_genome = radar.dir_genome
    	#dir_blast = radar.dir_blast
    	#dir_blastp = radar.dir_blastp
    	#dir_user = radar.dir_user
    	#dir_result = radar.dir_result
    	#file_ntdb = dir_db + "radar.fasta"
    	#file_aadb = dir_db + "radar.faa"
    	#file_udb = dir_db + "radar.udb"    	
    	#file_usearch = program_dir + "usearch" 
        pass
    @staticmethod
    def make_config( strain ) :
        def read_as_string( file0 ) :
            with open( file0, "r" ) as f :
                read = f.read()
            return read
        print dir_vis
        circos_contents = circ_contents.contents_circos()
        dir_conf = dir_vis + "conf/"
        dir_contig = dir_vis + "circos/"
        dir_txt = dir_vis + "txt/"
        dir_output = dir_vis + "output/"
        os.system( "mkdir -p %s" % dir_conf )
        os.system( "mkdir -p %s" % dir_contig )
        os.system( "mkdir -p %s" % dir_txt )
        os.system( "mkdir -p %s" % dir_output )
        print "Make directory... %s" % dir_conf
        print "Make directory... %s" % dir_contig
        print "Make directory... %s" % dir_txt
        print "Make directory... %s" % dir_output
        fix_dir = os.listdir( dir_wgs )
        for fix_x in fix_dir :
            file_name = fix_x.replace( ".fna", "" )
            contig_file = dir_contig + file_name
            if not os.path.exists( contig_file ) :
                with open( contig_file, "w" ) as f :
                    for contig_parse in SeqIO.parse( dir_wgs + fix_x, "fasta" ) :
                        des = contig_parse.description.split( " " )[ 0 ]
                        des_len = len( contig_parse.seq )
                        f.write( "chr - "+ des + " " + des + " 1 " + str( des_len ) + " black" + "\r\n" )
                print "Wrote... %s" % contig_file
            else :
                print "Contig file already existed "
                pass
            txt_file0 = dir_txt + "under_95_hit_%s" % file_name
            txt_file1 = dir_txt + "to_95_hit_%s" % file_name
            gc_txt_file = dir_txt + "gc_skew_%s" % file_name
            forward_file = dir_txt + "foward_anno_%s" % file_name
            reverse_file = dir_txt + "reverse_anno_%s" % file_name
            file_conf = dir_conf + "circos_chro_%s.conf" % file_name
            
            output = read_as_string( dir_blastp + fix_x.replace( ".fna", ".aln" ) )
            #print output
            output = output.split( "\nQuery >" )
            with open( txt_file0, "w" ) as f0, open( txt_file1, "w" ) as f1 :
                for output0 in output :
                    #if output0.startswith( "%s" % dir_user ) :
                    if output0.startswith( "/home/sierray/sbl_nas/" ) :
                        continue
                    output1 = output0.split( "\n" )
                    #print dir_user
                    #print output1
                    output2 = output1[ 2 ].split( )
                    query = output1[ 0 ]
                    target = output2[ 5 ].split( "|" )[ -1 ]
                    identity = float( output2[ 2 ].strip( "%" ) ) / 100
                    query_id = query.split( "# " )[ 0 ].split( "_" )[ :-1 ]
                    query_id = "_".join( query_id )
                    gene_start = query.split( "# " )[ 1 ]
                    gene_end = query.split( "# " )[ 2 ]
                    if identity < 0.95 :
                        to_print = [ query_id, gene_start, gene_end, target ]
                        to_print = [ str( m ) for m in to_print ]
                        to_print = " ".join( to_print )
                        f0.write( to_print )
                        f0.write( "\r\n" )
                    elif identity >= 0.95 :
                        to_print = [ query_id, gene_start, gene_end, target ]
                        to_print = [ str( m ) for m in to_print ]
                        to_print = " ".join( to_print )
                        f1.write( to_print )
                        f1.write( "\r\n" )
            count = 0 
            gc_num = 0
            for parse0 in SeqIO.parse( dir_wgs + fix_x, "fasta" ) :
                count += len( parse0.seq )
                gc_num += parse0.seq.count("G") + parse0.seq.count("C")
            scale_lst = []
            for scale_x in SeqIO.parse( dir_wgs + fix_x, "fasta" ):
                rdian = len( scale_x.seq ) / float( count )
                temp_ = ( scale_x.id, rdian )
                scale_lst.append( temp_ )
            wgs_gc = round( float( gc_num ) / count, 3 )
            with open( gc_txt_file, "w" ) as f0, open( forward_file, "w" ) as f1, open( reverse_file, "w" ) as f2 :
                for anno_parse in SeqIO.parse( dir_genome + fix_x, "fasta" ) :
                    gene_id = anno_parse.id
                    gene_id = gene_id.split( "_" )
                    gene_id = gene_id[ :-1 ]
                    gene_id = "_".join(gene_id)
                    anno_gc = round( float( anno_parse.description.split( " #" )[ 4 ].split( ";" )[ -1 ].split( "=" )[ -1 ] ) , 3 )
                    gc = ( round( ( anno_gc - wgs_gc ) * 100 ,3 ) )
                    start = anno_parse.description.split( " #" )[ 1 ]
                    end = anno_parse.description.split( " #" )[ 2 ]
                    anno_strand = anno_parse.description.split( " #" )[ 3 ]
                    gc_print = [ gene_id, start, end, gc ]
                    gc_print = [ str( print_n ) for print_n in gc_print ]
                    gc_print = " ".join( gc_print )
                    f0.write( gc_print )
                    f0.write( "\r\n" )
                    if anno_strand == " 1" :
                        to_print = [ gene_id, start, end ]
                        to_print = [ str( print_m ) for print_m in to_print ]
                        to_print = " ".join( to_print )
                        f1.write( to_print )
                        f1.write( "\r\n" )
                    else :
                        to_print = [ gene_id, start, end ]
                        to_print = [ str( print_o ) for print_o in to_print ]
                        to_print = " ".join( to_print )
                        f2.write( to_print )
                        f2.write( "\r\n" )
            with open( file_conf, "w") as g :
                g.write( circos_contents.return_contents_1() % contig_file )
                num_0 = 0
                for num in scale_lst :
                    g.write( "%s=%s" %( num[ 0 ],str( num[ 1 ] ) +"rn" ) )
                    num_0 += 1
                    if num_0 == len( scale_lst ) :
                        break
                    g.write( ", ")
                g.write( "\r\n")
                g.write( circos_contents.return_contents_2() )
                for num1 in scale_lst :
                    g.write( circos_contents.return_contents_3() %( num1[ 0 ], num1[ 0 ], "%01.1f" ) )
                g.write( circos_contents.return_contents_4() %( txt_file0, txt_file1, forward_file,  reverse_file, gc_txt_file ) )
                g.write( circos_contents.return_contents_5() %( dir_circos + "circos-0.69-6", dir_circos + "circos-0.69-6", dir_circos + "circos-0.69-6" ))
            print "Wrote... %s " % txt_file0
            print "Wrote... %s " % txt_file1
            print "Wrote... %s " % gc_txt_file
            print "Wrote... %s " % forward_file
            print "Wrote... %s " % reverse_file
            print "Wrote... %s " % file_conf
            print "Wrote Process End " + fix_x
    @staticmethod
    def run_circos( strain ) :
        print dir_user
        print dir_vis
        #print dir_conf
        print dir_wgs
        count_circos = 0
        dir_conf = dir_vis + "conf/"
        fix_dir = os.listdir( dir_wgs )
        for file_conf in fix_dir :
            file_conf_0 = file_conf.replace( ".fna", "" )
            file_conf_ = dir_conf + "circos_chro_" + file_conf.replace( "fna", "conf" )
            (exitstatus, outtext) = commands.getstatusoutput( "%s " % file_circos + "-conf %s" %file_conf_ )
            print outtext
            results = dir_user + "/script/circos"
            os.rename( results + ".png", results + file_conf_0 + ".png" )
            os.rename( results + ".svg", results + file_conf_0 + ".svg" )
            shutil.move( results + file_conf_0 + ".png", dir_vis + "output/" + file_conf_0 + ".png" )
            shutil.move( results + file_conf_0 + ".svg", dir_vis + "output/" + file_conf_0 + ".svg" )

    	#print dir_db
    	#print program_dir