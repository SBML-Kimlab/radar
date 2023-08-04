from radar import *
import radar
import html_layouts
import glob
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
import base64

class report : 
    def __init__( self, strain, cutoff ) :
        global dir_db, include_dir, program_dir, dir_genome, dir_vis, dir_blastp, dir_cluster, dir_user, dir_result, file_ntdb, file_aadb, file_udb, file_usearch
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
        dir_vis = radar.dir_vis
        file_usearch = program_dir + "usearch" 
        pass

    @staticmethod
    def results_report( strain, cutoff ) :
        def read_fasta( file_fasta ) :
            dic_fasta = dict()
            for record in SeqIO.parse( open( file_fasta, "r" ), "fasta" ) :
                seq_id = record.description
                seq = record.seq
                dic_fasta[ seq_id ] = seq
            return dic_fasta

        def genome_html( image_path, html_file, title_text = "", width = 500, height = 500, table_data = None ) :
            with open( image_path, "rb" ) as image : 
                encoded_image = base64.b64encode(image.read()).decode('utf-8')
            img_tag = f'<img src="data:image/png;base64,{encoded_image}" alt="Image" title="{title_text}" width="{width}" height="{height}">'
            title_tag = f'<p class="image-title">{ title_text }</p>'
            table_html = ''
            if table_data:
                table_rows = []
                for row in table_data:
                    cells = ''.join([f'<td>{cell}</td>' for cell in row])
                    table_rows.append(f'<tr>{cells}</tr>')
                table_html = f'<div class="table-container"><table>{"".join(table_rows)}</table></div>'

            # 이미지와 표를 합친 div
            div_tag = f'<div class="image-container">{title_tag}{img_tag}{table_html}</div>'

            style_tag = layout_instance.return_layout_1()
            with open( html_file, "w" ) as f : 
                if f.tell() == 0 :
                    f.write( style_tag )
                f.write( div_tag ) 
        
        def html_merged( html_content, in_html ) :
            with open( in_html ) as in_file :
                in_html_content = in_file.read()
                in_html_content = in_html_content.replace("</body>", "").replace("<body>", "")
                combined_content = html_content + in_html_content
                return combined_content


        path_output = dir_vis + "output"
        path_report = dir_result + "wgs/report"
        os.system( "mkdir -p %s" %path_report )
        layout_instance = html_layouts.html_layout()
        lst_genome = glob.glob( path_output + "/*.png" )
        for temp_genome in lst_genome :
            radar_genome = temp_genome.split( "/" )[ -1 ].replace( ".png", "" )
            genome_report = path_report + "/%s.html" %radar_genome
            ab_img = "%s/%s" %( dir_user, temp_genome )
            #genome_report = "%s/%s/%s.html" %( dir_user, path_report, radar_genome )
            title_text = "GENOME local alignment ideogram"
            genome_data = "%s%s.faa" %( dir_genome, radar_genome )
            blast_all_data = "%s%s_all.tsv" %( dir_blastp, radar_genome )
            blast_cutoff_data = "%s%s_cutoff.tsv" %( dir_blastp, radar_genome ) 
            len_orf = len( read_fasta( genome_data ) )
            df_blast_all = pd.read_csv( blast_all_data, sep = "\t", header = 0 )
            df_blast_cutoff = pd.read_csv( blast_cutoff_data, sep ="\t", header = 0 )
            df_blast_95cut = df_blast_all[ df_blast_all.apply( lambda x: x.identity > 0.95, axis = 1 ) ]
            data = [ [ "GENOME ID", radar_genome ],
            [ "ORF", len_orf ], [ "CUTOFF", "{:.0f}%".format( cutoff * 100 ) ], [ "HIT COUNT (ALL/CUTOFF)", "%s / %s" %( len( df_blast_all ), len( df_blast_cutoff ) ) ],
            [ "HIGHLY AMR HIT", len( df_blast_95cut ) ] ]
            genome_html( ab_img, genome_report, title_text, table_data = data )
            with open( genome_report, "r" ) as f :
                html_content_in_memory = f.read()
            
            snp_html = "%swgs/figure/snp_out/%s_snp.html" %( dir_result, radar_genome )
            cluster_esbl_html = "%swgs/figure/cluster/%s/esbl_hit_cluster.html" %( dir_result, radar_genome )
            cluster_non_esbl_html = "%swgs/figure/cluster/%s/non_esbl_hit_cluster.html" %( dir_result, radar_genome )
            phylo_esbl_html = "%swgs/figure/cluster/%s/esbl_hit_phylo.html" %( dir_result, radar_genome )
            phylo_non_esbl_html = "%swgs/figure/cluster/%s/non_esbl_hit_phylo.html" %( dir_result, radar_genome )

            html_files_to_merge = [ snp_html, cluster_esbl_html, phylo_esbl_html, cluster_non_esbl_html, phylo_non_esbl_html ]
            for in_html in html_files_to_merge:
                html_content_in_memory = html_merged( html_content_in_memory, in_html )

            with open( genome_report, "w" ) as combined_file:
                combined_file.write( html_content_in_memory )
