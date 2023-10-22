from radar import *
import radar
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from IPython.core.display import display, HTML

class blastp : 
    def __init__( self, strain, cutoff ) :
        global dir_db, include_dir, program_dir, dir_genome, dir_blastp, dir_user, dir_result, file_ntdb, file_aadb, file_udb, file_usearch
        self.strain = strain
        self.cutoff = cutoff
        dir_db = radar.db_dir
        include_dir = radar.include_dir
        program_dir = radar.program_dir
        dir_genome = radar.dir_genome        
        dir_blastp = radar.dir_blastp
        dir_user = radar.dir_user
        dir_result = radar.dir_result
        #file_ntdb = dir_db + "radar.fasta"
        #file_aadb = dir_db + "radar.faa"
        #file_udb = dir_db + "radar.udb"     
        file_usearch = program_dir + "usearch" 
        pass

    @staticmethod
    def blastp_parse( strain, cutoff ) :
        def read_as_string( file0 ) :
            with open( file0, "r" ) as f :
                read = f.read()
            return read

        def cal_coverage( query_align ) :
            spl_char = query_align.split( "-" )
            align_start = int( spl_char[ 0 ] )
            align_end = int( spl_char[ -1 ].split( "(" )[ 0 ] )
            un_align = int( spl_char[ -1 ].split( "(" )[ -1 ].split( ")" )[ 0 ] )
            align_query_len = float( align_end - align_start + 1 )
            total_query_len = float( align_end + un_align )
            coverage = align_query_len / total_query_len
            return coverage 

        def parse_mutations( msa_output ) :
            events = list()
            for i in range( 0, len( msa_output ), 3 ) :
                qry_line = msa_output[ i ]; tgt_line = msa_output[ i + 2 ]
                qry_seq = qry_line.split()[ 2 ]; tgt_seq = tgt_line.split()[ 2 ]
                qry_start = int( qry_line.split()[ 1 ] ); tgt_start = int( tgt_line.split()[ 1 ] )
                for j in range( len( qry_seq ) ) :
                    if qry_seq[ j ] == "-" and tgt_seq[ j ] != "-" :
                        events.append( f"{ tgt_seq[ j ] }{ qry_start + j }-" )
                    elif qry_seq[ j ] != "-" and tgt_seq[ j ] == "-" :
                        events.append( f"-{ qry_start + j }{ qry_seq[ j ] }" )
                    elif qry_seq[ j ] != "-" and tgt_seq[ j ] != "-" and qry_seq[ j ] != tgt_seq[ j ]:
                        events.append( f"{ tgt_seq[ j ] }{ qry_start + j }{ qry_seq[ j ] }" )
            return events

        def check_fod( hit_idx, lst_mutations ) : 
            #file_snp_tsv = "/mnt/d/pipeline/radar/dev/content/radar/include/snp_table_real.tsv"
            file_snp_tsv = include_dir + "snp_table_real.tsv"
            df_snp_table = pd.read_csv( file_snp_tsv, sep = "\t", header = None )
            df_snp_table.columns = [ "BOARDS_idx", "gene", "host", "genome", "blast_hit", "snp_details", "total_hit", "c_freq" ]
            filtered_row = df_snp_table[ df_snp_table[ "BOARDS_idx" ] == hit_idx ]
            
            dic_mut_cal = dict()
            for mutation in lst_mutations :            
                for idx, row in filtered_row.iterrows() :
                    snp_details_list = row[ "snp_details" ].split( ", " )
                    snp_details_list = [ x.replace( "]", "" ).replace( "[", "" ) for x in snp_details_list ]
                    updated_snp_details = list()
                    for snp_details in snp_details_list : 
                        if "," in snp_details : 
                            spl_mut = snp_details.split( "," )
                            ori_mut = spl_mut[ 0 ]; sub_mut0 = spl_mut[ 1: ]
                            updated_snp_details.append( ori_mut )
                            for x in sub_mut0 : 
                                sub_mut = "%s%s" %( ori_mut[ :-1 ], x )
                                updated_snp_details.append( sub_mut )
                        else : updated_snp_details.append( snp_details )
                    if mutation in updated_snp_details : 
                        if mutation not in dic_mut_cal.keys() :
                            dic_mut_cal[ mutation ] = { row[ "host" ] : { row[ "total_hit" ] : [ row[ "blast_hit" ] ] } }
                        else : 
                            if row[ "host" ] not in dic_mut_cal[ mutation ].keys() :
                                dic_mut_cal[ mutation ][ row[ "host" ] ] = { row[ "total_hit" ] : [ row[ "blast_hit" ] ] }
                            else : dic_mut_cal[ mutation ][ row[ "host" ] ][ row[ "total_hit" ] ].append( row[ "blast_hit" ] )
            return dic_mut_cal
        
        lst_aln = glob.glob( dir_blastp + "*.aln" )
        lst_aln.sort()
        count_file = 0 
        for file_aln in lst_aln :
            count_file += 1
            if count_file % 100 == 0 :
                clear_output()
                print ( "Processed %d files." % count_file )
            file_tsv0 = file_aln.replace( ".aln", "_all.tsv" )
            file_tsv1 = file_aln.replace( ".aln", "_cutoff.tsv" )
            
            if os.path.exists( file_tsv0 ) and os.path.exists( file_tsv1 ) == True :
                continue
            output = read_as_string( file_aln )
            output = output.split( "\nQuery >" )
            with open( file_tsv0, "w" ) as f, open( file_tsv1, "w" ) as g :
                f.write( "gene_id\ttarget\tscore\tevalue\tidentity\tcoverage\tdetection_summary\tsnp_event\tsnp_models\tsnp_model_idx\tfrequency of detection\n")
                g.write( "gene_id\ttarget\tscore\tevalue\tidentity\tcoverage\tdetection_summary\tsnp_event\tsnp_models\tsnp_model_idx\tfrequency of detection\n")
                for output0 in output : 
                    if output0.startswith( "program" ) : continue
                    output1 = output0.split( "\n" )
                    output2 = output1[ 2 ].split()
                    spl_msa0 = output0.split( "), Evalue " )[ 0 ] 
                    gene_from = output1[ 0 ].split( " " )[ 0 ]
                    gene_to = output2[ 5 ]
                    score = output2[ 0 ]
                    evalue = output2[ 1 ]
                    identity = float( output2[ 2 ].replace( "%", "" ) ) / 100
                    query_align = output2[ 3 ]
                    coverage = cal_coverage( query_align )        
                    if identity == 1.0 : 
                        detection_event = "Perfect"; snp_event = "no SNP"; snp_models = "N/A"; snp_model_idx = "N/A"; fod = "N/A"
                        to_print = [ gene_from, gene_to, score, evalue, identity, coverage, detection_event, snp_event, snp_models, snp_model_idx, fod ]
                    elif identity >= 0.9 and identity < 1.0 :     
                        spl_msa1 = ( "\nQry " + "\nQry ".join( spl_msa0.split( "\nQry ")[ 1: ] ) ).split( "\n" )[ :-1 ]
                        msa = list( filter( bool, spl_msa1 ) )
                        mutations = parse_mutations( msa )
                        out_mutations = "," .join( mutations ) 
                        detection_event = "Highly"; snp_event = out_mutations   
                        hit_boards_idx = int( gene_to.split( "|" )[ 0 ].replace( "_BOARDS", "" ) )
                        dic_mut = check_fod( hit_boards_idx, mutations )
                        if not dic_mut : 
                            snp_models = "Not detected"; snp_model_idx = "N/A"; fod = "N/A" 
                            to_print = [ gene_from, gene_to, score, evalue, identity, coverage, detection_event, snp_event, snp_models, snp_model_idx, fod ]
                        else :
                            snp_models = "Detected"; fod_list = list()        
                            snp_model_idx = hit_boards_idx
                            for mutation in mutations :
                                if mutation not in dic_mut :
                                    fod_list.append( "%s-NewSNP" % mutation )
                                else : 
                                    for host, total_hits in dic_mut[ mutation ].items() :
                                        for total_hit, blast_hit in total_hits.items() :
                                            total_blast_hits = sum( blast_hit ) 
                                            fod_val = round( float( total_blast_hits ) / float( total_hit ), 3 )
                                            fod_list.append( "%s-%s-%s" %( mutation, host, fod_val ) )
                                fod_list.append( "|" )
                            fod = ",".join( fod_list )                
                            to_print = [ gene_from, gene_to, score, evalue, identity, coverage, detection_event, snp_event, snp_models, snp_model_idx, fod ]            
                    elif identity < 0.9 : 
                        detection_event = "Poorly"; snp_event = "SNP not applicable"; snp_models = "N/A"; snp_model_idx = "N/A"; fod = "N/A"
                        to_print = [ gene_from, gene_to, score, evalue, identity, coverage, detection_event, snp_event, snp_models, snp_model_idx, fod ]
                    if float( to_print[ 4 ] ) >= cutoff : 
                        to_print1 = [ str( m ) for m in to_print ]
                        to_print1 = "\t".join( to_print1 ) 
                        g.write( "%s\n" %to_print1 )
                    to_print0 = [ str( m ) for m in to_print ]
                    to_print0 = "\t".join( to_print0 ) 
                    f.write( "%s\n" %( to_print0 ) )

    @staticmethod
    def merge_after_parse( strain, cutoff ) :
        def file_readlines( file0 ) :
            f = open( file0, "r" )
            lines = f.read().splitlines()[ 1: ]
            f.close()
            return lines
        
        count_file = 0 ; lst_identity = list()        
        dir_pr_blastp = os.path.dirname( os.path.dirname( dir_blastp ) )        
        for dir_ in os.listdir( dir_pr_blastp ) :
            file_i = "%s/%s" %( dir_pr_blastp, dir_ ); file_w = glob.glob( file_i + "/*_all.tsv" )
            if count_file == 0 :
                lst_files = list( file_w )
                count_file += 1 
            elif count_file > 0 :
                lst_files.extend( file_w )
                count_file += 1
        print ( ".tsv files: %d" % len( lst_files ) )
        for file0 in lst_files :
            lines = file_readlines( file0 )
            lst_identity0 = [ float( ( line.split( "\t" )[ 4 ] ) ) for line in lines ]
            lst_identity.extend( lst_identity0 )        

        dir_figure_wgs = dir_result + "wgs/figure/"
        dir_table_wgs = dir_result + "wgs/table/"
        os.system( "mkdir -p %s" % dir_figure_wgs )
        os.system( "mkdir -p %s" % dir_table_wgs )
        fig, ax = plt.subplots( 1, 1, figsize = ( 6, 6 ) )
        plt.hist( lst_identity, density = False, bins = 16,  facecolor = "gray", alpha = 0.75, rwidth = 0.9 )
        ax.set_xlabel( "P. identity" )
        #ax.set_ylabel( "100K" )
        plt.grid( True ) 
        plt.draw()
        file_figure = dir_figure_wgs + "wgs_search_group.png" 
        #plt.axis( [ 0.1, 1.0, 0, 9 ] )
        #plt.grid(True)
        plt.savefig( file_figure, dpi = 600 )
        plt.savefig( file_figure.replace( ".png", ".pdf" ), format = "pdf" )
        plt.savefig( file_figure.replace( ".png", ".svg" ), format = "svg" )
        plt.show()
        
        file_blastp_merged = dir_table_wgs + "blastp_merged.tsv"
        file_blastp_merged0 = file_blastp_merged.replace( ".tsv", "0.tsv" )
        file_blastp_merged_esbl = file_blastp_merged.replace( ".tsv", "_esbl.tsv" )
        lst_esbl = [ "TEM", "SHV", "CTX", "KPC", "VIM", "IMP", "NDM", "CMY", "Amp", "OXA", "GES", "PER" ]
        dic_esbl = {}
        data_list = list()
        with open( file_blastp_merged, "w" ) as f :
            f.write( "genome_id\tgene_id\ttarget\tscore\tevalue\tp.identity\tcoverage\tdetection_summary\tsnp_detect\tsnp_model\tsnp_model_idx\n" )
            for file0 in lst_files :
                genome_id = file0.split( "/" )[ -1 ].replace( ".tsv", "" )
                if "_" in genome_id :
                    genome_id = genome_id.split( "_" )
                    genome_id = "%s_%s" % ( genome_id[ 0 ], genome_id[ 1 ] )
                lines = file_readlines( file0 )    
                for line in lines :
                    line0 = line.split( "\t" )
                    if float( line0[ 4 ] ) >= cutoff : 
                        if line0[ 7 ] != "no SNP" and line0[ 7 ] != "SNP not applicable" :
                            to_print = [ genome_id, line0[ 0 ], line0[ 1 ], line0[ 2 ], line0[ 3 ], line0[ 4 ], line0[ 5 ], line0[ 6 ], "SNP detected", line0[ 8 ], line0[ 9 ] ]
                            if line0[ 8 ] == "Detected" : 
                                model_index = int( line0[ 9 ] )
                                base_url = "https://sbml.unist.ac.kr/psp/"
                                data_list.append( {
                                    "genome_id" : genome_id,
                                    "gene_id" : line0[ 0 ],
                                    "snp_detect" : line0[ 8 ],
                                    "SNPmodel_index" : line0[ 9 ],
                                    "ref_model": f'<a href="{base_url + line0[9] + "_BOARDS/ref_model/ranked_0.pdb"}">Download</a>' if line0[9] else None,
                                    "snp_model(AF2)": f'<a href="{base_url + line0[9] + "_BOARDS/mut_model/af2/0/ranked_0.pdb"}">Download</a>' if line0[9] else None,
                                    "snp_model(RTF)": f'<a href="{base_url + line0[9] + "_BOARDS/mut_model/rtf/0/model_1.lddt.pdb"}">Download</a>'if line0[9] else None
                                } )
                        else : 
                            to_print = [ genome_id, line0[ 0 ], line0[ 1 ], line0[ 2 ], line0[ 3 ], line0[ 4 ], line0[ 5 ], line0[ 6 ], line0[ 7 ], line0[ 8 ], line0[ 9 ] ]
                        to_print = [ str( m ) for m in to_print ]
                        to_print = "\t".join( to_print )
                        f.write( "%s\n" %to_print )
        print ( "Wrote... %s" % file_blastp_merged )
        df = pd.DataFrame(data_list)
        display(HTML(df.to_html(escape=False, index=False)))
        #display(df[['genome_id', 'gene_id', 'snp_detect', 'ref_model', 'snp_model(AF2)', 'snp_model(RTF)']])
        
        with open( file_blastp_merged0, "w" ) as f, open( file_blastp_merged_esbl, "w" ) as f0 :
            f.write( "genome_id\tgene_id\ttarget\tscore\tevalue\tp.identity\tcoverage\tdetection_summary\tsnp_detect\n" )
            f0.write( "genome_id\tgene_id\ttarget\tscore\tevalue\tp.identity\tcoverage\tdetection_summary\tsnp_detect\n" )
            for line in file_readlines( file_blastp_merged ) :
                line0 = line.split( "\t" )
                to_print = list( line0 )
                gene_to = line0[ 2 ]
                to_print =  [ str( m ) for m in to_print ]
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
        print ( "Wrote... %s" % file_blastp_merged0 )
        print ( "Wrote... %s" % file_blastp_merged_esbl )
        print ( dic_esbl )

    @staticmethod
    def wgs_snp_landscape( strain, cutoff ) :
        def snp_landscape( file_blastp, out_html ) :
            esbl_subclasses = [ "TEM", "SHV", "CTX", "KPC", "VIM", "IMP", "NDM", "CMY", "amp", "OXA", "GES", "PER" ]
            non_esbl_subclasses = [ "van", "mec", "MCR", "others" ]
            all_subclasses = esbl_subclasses + non_esbl_subclasses
            df_blastp = pd.read_csv( file_blastp, sep = "\t", header = 0 )
            filtered_df = df_blastp[ ~df_blastp[ "snp_event" ].isin( [ "no SNP", "SNP not applicable" ] ) ]
            dic_hits_with_snps = { index: row.to_dict() for index, row in filtered_df.iterrows() }
            data = { "Hits" : [ value[ "gene_id" ] for value in dic_hits_with_snps.values() ] }

            for subclass in all_subclasses :
                data[ subclass ] = [ 0 ] * len( data[ "Hits" ] )
            char_replace = [ "Escherichia_coli_", "bla" ]
            for index, ( k, v ) in enumerate( dic_hits_with_snps.items() ) :
                temp_target = v[ "target" ].split( "|" )[ -1 ].split( "_[" )[ 0 ]
                for char in char_replace :
                    temp_target = temp_target.replace( char, "" )
                target = temp_target.upper()
                matched_subclass = next( ( subclass for subclass in all_subclasses if target.startswith( subclass.upper() ) ), None )
                if matched_subclass:
                    subclass = matched_subclass
                else: subclass = "others"
                data[ subclass][ index ] = float( v[ "identity" ] * v[ "coverage" ] )
            df = pd.DataFrame( data ) 
            heatmap_data = df.set_index( "Hits" ).T
            hover_texts = list()
            for index in heatmap_data.index :
                column_texts = list()
                for column in heatmap_data.columns :
                    current_value = heatmap_data[ column ][ index ]
                    gene_info = None
                    for k, v in dic_hits_with_snps.items() :
                        if v[ "gene_id" ] == column : 
                            gene_info = v  
                    if index.upper() in gene_info[ "target" ].split( "|" )[ -1 ].split( "[" )[ 0 ].upper() : 
                        if gene_info : 
                            snp_models = gene_info.get( "snp_models" )
                            fod = str( gene_info.get( "frequency of detection" ) ).replace( ",|,", "<br>" ).replace( ",|", "" )
                    elif index == "others" and current_value != 0 : 
                        if gene_info :
                            snp_models = gene_info.get( "snp_models" )
                            fod = str( gene_info.get( "frequency of detection" ) ).replace( ",|,", "<br>" ).replace( ",|", "" )
                    else : snp_models = "Not detected"; fod = "N/A"
                    formatted_string = (
                        f"gene_id: {column}<br>"
                        f"subclass: {index}<br>"
                        f"score(identity*coverage): {current_value}<br>"
                        f"snp model: {snp_models}<br>"
                        f"frequency of detection: {fod}"
                    )
                    column_texts.append( formatted_string )
                hover_texts.append( column_texts )
            fig = make_subplots( rows = 1, cols = 1 )
            heatmap = go.Heatmap(
                z = heatmap_data.values,
                x = heatmap_data.columns,
                y = heatmap_data.index,
                colorscale = [ ( 0, "white" ), ( 1, "darkgray" ) ],
                text = hover_texts,
                hoverinfo = "text" )
            fig.add_trace( heatmap, row = 1, col = 1 )
            fig.update_layout(
                height = 800,
                width = 800,
                xaxis_nticks = len( heatmap_data.columns ), 
                xaxis = dict( side = "bottom" ), 
                yaxis_title = "Subclasses",
                xaxis_title = "Hits",
                coloraxis_colorbar = dict( len = 0.5, yanchor = "top", y = 0.5 ),
                title_text = "SNP landscape in WGS",
                title_x = 0.5,
                title_font = dict( color = "black" ) )                
            #fig.show()
            fig.write_html( out_html, auto_open = False )

        path_snp_landscape = dir_result + "wgs/figure/snp_out"
        os.system( "mkdir -p %s" %path_snp_landscape )
        lst_file_blastp = glob.glob( dir_blastp + "*_cutoff.tsv" )
        for file_blastp in lst_file_blastp :
            out_html_name = file_blastp.split( "/" )[ -1 ].replace( "_cutoff.tsv", "_snp.html" )
            out_html = "%s/%s"%( path_snp_landscape, out_html_name )
            snp_landscape( file_blastp, out_html )
