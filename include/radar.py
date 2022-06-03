#for local
import os
import sys
import glob
import math
import datetime
import shutil
#import commands
import subprocess
import os.path as path
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib import rc
from IPython.display import clear_output
from Bio import SeqIO
import blast_sc0, blast_sc1, visual_ideo, circ_contents
from prodigal_sc import prodigal
from db_stat import db_statistics

class amr :
	def __init__( self, strain, db ) :
		#Built-in: include, database, program, script.
		print ( "Pipeline initiating..." ) 
		global dir_user, main_dir, program_dir, dir_wgs, dir_genome, db_dir, dir_blast, dir_blastp, dir_result, dir_vis, dir_circos
		global file_usearch, file_diamond, file_circos
		self.strain = strain
		self.db = db 
		
		dir_user = path.abspath( path.join( os.getcwd() ) )
		os.chdir( dir_user ) 
		#abs_piv = dir_user0.split( "/radar" )[ 0 ]
		main_dir = "pipeline/"
		program_dir = "program/"
		db_dir = "database/"
		dir_wgs = main_dir + "antibiotics/genome/sequence/1.fna/" + strain + "/"
		dir_genome = main_dir + "antibiotics/genome/annotation/2.anno/" + strain + "/"
		dir_blast = main_dir + "antibiotics/genome/blastp/3.align/" + strain + "/"
		dir_blastp = dir_blast + strain + "/"
		dir_result = main_dir + "antibiotics/result/"
		dir_vis = main_dir + "antibiotics/visual/" + strain + "/"
		dir_circos = program_dir + "circos/"
		file_circos = dir_circos + "circos-0.69-6/bin/circos"
		os.system( "mkdir -p %s" % main_dir )
		os.system( "mkdir -p %s" % dir_wgs ) 
		os.system( "mkdir -p %s" % dir_genome )
		os.system( "mkdir -p %s" % dir_blastp )
		os.system( "mkdir -p %s" % dir_result )
		os.system( "mkdir -p %s" % dir_vis )

	class method( ) :
		def __init__( self ) :
			pass
		@staticmethod
		def prodigal( strain ) :
			prodigal.prodigal_sc1( strain )
		@staticmethod
		def db_statistics( strain, db ) :
			db_statistics.db_stat( strain, db )
		@staticmethod
		class blast_method( ) :
			def __init__( self ) :
				pass
			@staticmethod
			def udb_making( strain, db ) :
				call_class1_1 = blast_sc0.blastp( strain, db ) 
				call_class1_1.make_udb( strain, db )
			@staticmethod
			def blastp_run( strain, db ) :
				call_class1_2 = blast_sc0.blastp( strain, db )
				call_class1_2.blastp_run( strain, db )
		class blast_parse_method( ) :
			def __init__( self ) :
				pass
			@staticmethod
			def blastp_parse( strain, cutoff ) :
				call_class1_3 = blast_sc1.blastp( strain, cutoff )
				call_class1_3.blastp_parse( strain, cutoff )
			@staticmethod
			def blastp_merge( strain, cutoff ) :
				call_class1_4 = blast_sc1.blastp( strain, cutoff )
				call_class1_4.merge_after_parse( strain, cutoff )
		@staticmethod
		class visual_method( ) :
			def __init__( self ) :
				pass
			@staticmethod
			def make_config_files( strain ) :
				call_class2_1 = visual_ideo.circos( strain )
				call_class2_1.make_config( strain )
			@staticmethod
			def run_circos( strain ) :
				call_class2_2 = visual_ideo.circos( strain )
				call_class2_2.run_circos( strain )

