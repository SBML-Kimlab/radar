from radar import *
import radar
from Bio import SeqIO
import pandas as pd, numpy as np
import pickle as pkl
import torch, esm
import umap
import toytree, toyplot.png, toyplot.svg, toyplot.html
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import plotly.graph_objects as go
import plotly.express as px
import plotly.figure_factory as ff
from plotly.offline import plot

class cluster : 
    def __init__( self, strain, cutoff ) :
        global dir_db, include_dir, program_dir, dir_genome, dir_blastp, dir_cluster, dir_user, dir_result, file_ntdb, file_aadb, file_udb, file_usearch
        self.strain = strain
        self.cutoff = cutoff
        dir_db = radar.db_dir
        program_dir = radar.program_dir
        include_dir = radar.include_dir
        dir_genome = radar.dir_genome        
        dir_blastp = radar.dir_blastp
        dir_cluster = radar.dir_cluster
        dir_user = radar.dir_user
        dir_result = radar.dir_result
        file_usearch = program_dir + "usearch" 
        pass

    @staticmethod
    def cluster_hit( strain, cutoff ) :
        def read_fasta( file_fasta ) :
            dic_fasta = dict()
            for record in SeqIO.parse( open( file_fasta, "r" ), "fasta" ) :
                seq_id = record.description
                seq = record.seq
                dic_fasta[ seq_id ] = seq
            return dic_fasta

        def read_fasta_hit( file_fasta ) :
            dic_fasta = dict()
            for record in SeqIO.parse( open( file_fasta, "r" ), "fasta" ) :
                seq_id = record.description.split()[ 0 ]
                seq = record.seq
                genome_id = file_fasta.split( "/" )[ -1 ]
                if genome_id not in dic_fasta.keys() :
                    dic_fasta[ genome_id ] = { seq_id : seq }
                else : dic_fasta[ genome_id ].update( { seq_id : seq } )
            return dic_fasta

        def classification_amr( target ) :
            lst_esbl = [ "TEM", "SHV", "CTX", "KPC", "IMP", "VIM", "NDM", "Amp", "CMY", "OXA", "GES", "PER" ]
            lst_non_esbl = [ "van", "mec", "MCR", "others" ]
            found = False
            for esbl in lst_esbl : 
                if esbl in target :
                    out_esbl = "ESBL:%s" %( esbl )
                    return out_esbl
                    found = True
                    break
            if not found :
                for non_esbl in lst_non_esbl :
                    if non_esbl in target :
                        out_non_esbl = "non-ESBL:%s" %( non_esbl )
                        return out_non_esbl
                        found = True
                        break
                if not found :
                    out_non_esbl1 = "non-ESBL:others" 
                    return out_non_esbl1

        def cluster_division( file_cluster ) :
            out_esbl_cluster = file_cluster.replace( ".faa" ,"_esbl.faa" )
            out_non_esbl_cluster = file_cluster.replace( ".faa", "_non_esbl.faa" )
            def read_cluster( file_cluster ) :
                dic_cluster = dict()
                for record in SeqIO.parse( open( file_cluster, "r" ), "fasta" ) :
                    seq_id = record.description
                    seq = record.seq
                    dic_cluster[ seq_id ] = seq
                return dic_cluster
            read_cluster0 = read_cluster( file_cluster )
            with open( out_esbl_cluster, "w" ) as f, open( out_non_esbl_cluster, "w" ) as g :
                for k, v in read_cluster0.items() :
                    if "|ESBL:" in k : 
                        f.write( ">%s\n%s\n" %( k, v ) )
                    else : g.write( ">%s\n%s\n" %( k, v ) )

        def embedding_esm_model( file_hit, in_pkl, out_pkl ) :
            with open( in_pkl, "rb" ) as f :
                ref_data = pkl.load( f )
                f.close()
            model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
            layers_len = len( model.layers)
            batch_converter = alphabet.get_batch_converter()
            model.eval()
            whole_tokens = list(); whole_labels = list()
            data = dict()
            model.cuda()
            for record in SeqIO.parse( file_hit, "fasta" ) :
                temp = [ ( record.description, str(record.seq.split( "*" )[ 0 ] ) ) ]
                batch_labels, batch_strs, batch_tokens = batch_converter( temp )
                batch_tokens = batch_tokens.cuda()
                with torch.no_grad() :
                    results = model( batch_tokens, repr_layers = [ layers_len ], return_contacts = False )
                tokens = results[ "representations" ][ layers_len ]
                tokens = tokens.squeeze( 0 ).cpu().numpy()
                whole_tokens.append( np.mean( tokens, axis = 0 ) )
                whole_labels.append( batch_labels[ 0 ] )
            for i in range( len( whole_labels ) ) :
                data[ whole_labels[ i ] ] = whole_tokens[ i ]
            updated_data = { **ref_data, **data }
            with open( out_pkl, "wb" ) as g :
                pkl.dump( updated_data, g )
                g.close()

        def processing_umap( in_pkl, num ) :
            with open( in_pkl, "rb" ) as f :
                data = pkl.load( f )
                f.close()
            whole_labels = list(); whole_tokens = list()
            for key in data.keys() :
                whole_labels.append( key )
                whole_tokens.append( data[ key ] )
            whole_tokens = np.array( whole_tokens )
            reducer = umap.UMAP( n_neighbors = num, min_dist = 0.7 )
            embedding = reducer.fit_transform( whole_tokens )
            print( embedding.shape )
            return embedding

        def clustering_visualization( in_faa, embedding_array, title_text, out_html ) :
            subclass = list(); category = list(); description = list()
            anno = list(); seq = list(); 
            snp_detect = list(); detection_mode = list()
            color_map = { "TEM" : "#F81C2B", "CTX" : "#AD4772", "SHV" : "#D3314E", "KPC" : "#885C95",
            "IMP" : "#9D9634", "VIM" : "#E4C228", "NDM" : "#EDFC14",
            "AmpC" : "#B0663E", "CMY" : "#D3833C", "OXA" : "#D037F6",
            "GES" : "#00F7DA", "PER" : "#1BC0F7", "detection_hit" : "#000000" }

            for record in SeqIO.parse( in_faa, "fasta" ) :
                spl_des = record.description.split( "|" )
                description.append( str( record.description ) )
                seq.append( str( record.seq ) )
                if len( spl_des ) > 6 : #ref
                    category.append( spl_des[ 0 ] )
                    subclass.append( spl_des[ 1 ].split( ":" )[ -1 ] )
                    if "_[" in spl_des[ -1 ] :
                        anno.append( spl_des[ -1 ].split( "_[" )[ 0 ] )
                    else : anno.append( spl_des[ -1 ] )
                    snp_detect.append( "BOARDS_ref" )
                    detection_mode.append( "None" )
                else : #hit
                    #temp_subclass = "%s_HIT" %( spl_des[ 1 ].split( ":" )[ -1 ] )
                    subclass.append( "detection_hit" )
                    category.append( "HIT" )
                    anno.append( spl_des[ 2 ] )
                    snp_detect.append( spl_des[ 3 ] )
                    detection_mode.append( "|".join( spl_des[ 4: ] ) )

            df = pd.DataFrame( embedding_array, columns = [ "x", "y" ] )
            df[ "subclass" ] = subclass; df[ "allele" ] = anno; df[ "category" ] = category;
            df[ "snp_detect" ] = snp_detect; df[ "detection_mode(pidentity|coverage)" ] = detection_mode; df[ "seq" ] = seq
            df[ "description" ] = description
            #display( df )
            fig_scatter = px.scatter( df, x = "x", y = "y", color = "subclass", color_discrete_map = color_map, hover_data = [ "allele", "category", "snp_detect", "detection_mode(pidentity|coverage)" ] )
            fig_scatter.update_layout( title_text = title_text, title_x = 0.5, title_font = dict( color = "black" ) )
            fig_scatter.update_layout( legend = dict( x = 1, y = 0.5, xanchor = "left",) )
            #fig_scatter.show()
            plot( fig_scatter, filename = out_html, auto_open = False )
            return df

        def compute_cenroids( df ) :
            df_filtered = df[ df[ "subclass" ] != "detection_hit" ]
            centroid_df = df_filtered.groupby( "subclass" )[ [ "x", "y" ] ].mean().reset_index()
            centroid_indices = list()
            for idx, row in centroid_df.iterrows() :
                subclass_data = df[ df[ "subclass" ] == row[ "subclass" ] ]
                closest_index = ( ( subclass_data[ "x" ] - row[ "x" ] ) ** 2 + ( subclass_data[ "y" ] - row[ "y" ] ) ** 2 ).idxmin()
                centroid_indices.append( closest_index )
            centroid_df[ "closest_point_index" ] = centroid_indices
            return centroid_df

        def extract_anno_seq_from_centroids( df, centroids ) :
            anno_seq_data = list()
            for idx, row in centroids.iterrows() : 
                closest_index = row[ "closest_point_index" ] 
                anno = df.loc[ closest_index, "description" ]
                sequence = df.loc[ closest_index, "seq" ] 
                anno_seq_data.append( ( anno, sequence ) )
            return anno_seq_data

        def write_centroid( data, dic_hit, out_faa ) :
            with open( out_faa, "w" ) as f :
                for anno, seq in data : 
                    spl_anno = anno.split( "|" ) 
                    char_replace = [ "'", "(", ")", ":", " " ]
                    temp_anno = spl_anno[ -1 ]
                    for char in char_replace :
                        temp_anno = spl_anno[ -1 ].replace( char, "_" )
                    centroid_subclass = "REF|%s|%s|%s" % (spl_anno[0], spl_anno[1], temp_anno)
                    f.write( ">%s\n%s\n" %( centroid_subclass, seq ) )
                for k, v in dic_hit.items() :
                    f.write( ">%s\n%s\n" %( k.split( "|" )[ 2 ], v ) )

        def vis_phylogeny( file_dnd, title_text, out_phylo ) :
            treestring = open( file_dnd ).read().strip().replace( "\n", "" )
            tre = toytree.tree( treestring )
            canvas = toyplot.Canvas( width = 1000, height = 700 )
            axes = canvas.cartesian()
            axes.label.text = title_text
            axes.label.style = { "font-size": "15px", "fill": "black", "-toyplot-anchor-shift": "15px" }
            tre.draw( axes = axes, tip_labels_align = True, scalebar = True )
            #canvas, axes, mark = tre.draw( tip_labels_align = True, scalebar = True, width = 1000, height = 500 )
            canvas.style = { "background-color": "white" }
            toyplot.png.render( canvas, out_phylo )
            toyplot.svg.render( canvas, out_phylo.replace( ".png", ".svg" ) )
            with open( out_phylo.replace( ".png", ".html" ), "wb" ) as f :
                toyplot.html.render( canvas, f )

        file_hit = dir_result + "wgs/table/blastp_merged.tsv"
        df_hit = pd.read_csv( file_hit, sep = "\t", header = 0 ) 
        lst_blastp_results = glob.glob( dir_blastp + "/*_cutoff.tsv" )
        in_esbl_pkl = dir_db + "BOARDS_cluster_esbl.pkl" 
        in_non_esbl_pkl = dir_db + "BOARDS_cluster_non_esbl.pkl"

        for blastp_tsv in lst_blastp_results :
            genome = blastp_tsv.split( "/" )[ -1 ].replace( "_cutoff.tsv", "" )
            path_cluster = "%s%s" %( dir_cluster , genome )
            path_cluster_result = dir_result + "wgs/figure/cluster/%s" %genome 
            os.system( "mkdir -p %s" % path_cluster )
            os.system( "mkdir -p %s" % path_cluster_result )
            print ( genome )
            boards_esbl_cluster = dir_db + "BOARDS_cluster_esbl.faa" 
            boards_non_esbl_cluster = dir_db + "BOARDS_cluster_non_esbl.faa" 
            out_hit0 = path_cluster + "/hit_cluster.faa" 
            esbl_cluster = path_cluster + "/hit_cluster_esbl.faa" 
            non_esbl_cluster = path_cluster + "/hit_cluster_non_esbl.faa" 
            esbl_cluster_merged = path_cluster + "/hit_cluster_esbl_merged.faa" 
            non_esbl_cluster_merged = path_cluster + "/hit_cluster_non_esbl_merged.faa"
            
            out_esbl_pkl = path_cluster + "/hit_cluster_esbl.pkl" 
            out_non_esbl_pkl = path_cluster + "/hit_cluster_non_esbl.pkl"
            file_centroid_esbl = path_cluster + "/centroid_esbl.faa" 
            file_centroid_non_esbl = path_cluster + "/centroid_non_esbl.faa"

            out_esbl_cluster = path_cluster_result + "/esbl_hit_cluster.html"
            out_non_esbl_cluster = path_cluster_result + "/non_esbl_hit_cluster.html"

            faa = "%s/%s.faa" %( dir_genome, genome )
            dic_faa = read_fasta_hit( faa )
            with open( out_hit0, "w" ) as f :
                for idx0, row0 in df_hit.iterrows() :
                    matching_keys = [ key for key in dic_faa if row0[ "genome_id" ] in key ]
                    if matching_keys :
                        relevant_key = matching_keys[ 0 ]
                        if row0[ "gene_id" ] in dic_faa[ relevant_key ] :
                            value = dic_faa[ relevant_key ][ row0[ "gene_id" ] ]
                            BOARDS_idx = row0[ "target" ].split( "|" )[ 0 ]; pidentity = row0[ "p.identity" ]; coverage = row0[ "coverage" ]
                            genome_gene_id = "%s_%s" %( row0[ "genome_id" ], row0[ "gene_id" ] )
                            detection_mode = row0[ "detection_summary" ]
                            snp_detect = row0[ "snp_detect" ]
                            out0 = classification_amr( row0[ "target" ] )
                            annotation = "%s|%s|%s|%s|%s(%s|%s)" %( BOARDS_idx, out0, genome_gene_id, snp_detect, detection_mode, pidentity, coverage )
                            f.write( ">%s\n%s\n" %( annotation, value ) )
            cluster_division( out_hit0 )
            os.system( "cat %s %s > %s" %( boards_esbl_cluster, esbl_cluster, esbl_cluster_merged ) )
            os.system( "cat %s %s > %s" %( boards_non_esbl_cluster, non_esbl_cluster, non_esbl_cluster_merged ) )

            dic_esbl_hit = read_fasta( esbl_cluster )
            dic_non_esbl_hit = read_fasta( non_esbl_cluster )
            embedding_esm_model( esbl_cluster, in_esbl_pkl, out_esbl_pkl )
            embedding_esm_model( non_esbl_cluster, in_non_esbl_pkl, out_non_esbl_pkl )
            embedding_array_esbl = processing_umap( out_esbl_pkl, 5 )
            embedding_array_non_esbl = processing_umap( out_non_esbl_pkl, 5 )
            out_text0 = "HIT in ESBL CLUSTER"; out_text1 = "HIT in NON-ESBL CLUSTER"
            df_esbl = clustering_visualization( esbl_cluster_merged, embedding_array_esbl, out_text0, out_esbl_cluster )
            df_non_esbl = clustering_visualization( non_esbl_cluster_merged, embedding_array_non_esbl, out_text1, out_non_esbl_cluster )
            esbl_centroids = compute_cenroids( df_esbl )
            non_esbl_centroids = compute_cenroids( df_non_esbl )
            esbl_anno_seq_data = extract_anno_seq_from_centroids( df_esbl, esbl_centroids )
            non_esbl_anno_seq_data = extract_anno_seq_from_centroids( df_non_esbl, non_esbl_centroids )
            write_centroid( esbl_anno_seq_data, dic_esbl_hit, file_centroid_esbl )
            write_centroid( non_esbl_anno_seq_data, dic_non_esbl_hit, file_centroid_non_esbl )

            command0 = "clustalw -INFILE=%s" %( file_centroid_esbl )
            command1 = "clustalw -INFILE=%s" %( file_centroid_non_esbl )
            out_esbl_phylo = path_cluster_result + "/esbl_hit_phylo.png" 
            out_non_esbl_phylo = path_cluster_result + "/non_esbl_hit_phylo.png" 
            ( exitstatus, outtext ) = subprocess.getstatusoutput( "%s" %command0 )
            ( exitstatus, outtext ) = subprocess.getstatusoutput( "%s" %command1 )
            os.system( "rm %s" %file_centroid_esbl.replace( ".faa", ".aln" ) )
            os.system( "rm %s" %file_centroid_non_esbl.replace( ".faa", ".aln" ) )
            out_text2 = "Phylogenetic Tree (ESBL subclasses VS HIT)"
            out_text3 = "Phylogenetic Tree (Non-ESBL subclasses VS HIT)"
            vis_phylogeny( file_centroid_esbl.replace( ".faa", ".dnd" ), out_text2, out_esbl_phylo )
            vis_phylogeny( file_centroid_non_esbl.replace( ".faa", ".dnd" ), out_text3, out_non_esbl_phylo )