class contents_circos :
    def __init__( self ) :
        pass
    @staticmethod
    def return_contents_1( ) :
        contents_1 = """karyotype = %s
<ideogram>

<spacing>
default = 0.005r
</spacing>

radius = 0.75r
thickness = 5p
fill = yes

</ideogram>

show_ticks = yes
show_tick_labels = yes

<ticks>

chromosomes_display_default     = no
chromosomes_units               = 1000000
chromosomes_scale               = """
        return contents_1
    @staticmethod
    def return_contents_2( ) :
        contents_2 = """chromosomes_radius              = 0.75r
        
radius = 1r
color = black
thickness = 2p
multiplier = 1e-6
format = %01.1f
suffix = " Mb"\
"""
        return contents_2
    @staticmethod
    def return_contents_3( ) :
        contents_3 = """\r\n<tick>
chromosomes = %s
spacing = 10000u
size = 5p
thickness = 2p
color = black
show_label = no
</tick>

<tick>
chromosomes = %s
spacing = 100000u
size = 90p
thickness = 2p
color = black
show_label = yes
label_size = 30p
label_offset = 10p
format = %s
</tick>\
"""
        return contents_3
    @staticmethod
    def return_contents_4( ) :
        contents_4 = """\r\n</ticks>


<plots>

<plot>
type = tile
file = %s
r1 = 0.74r
r0 = 0.65r
layers    = 1
margin    = 0.2u
thickness = 100
padding   = 0
orientation      = out
stroke_thickness = 0.01
stroke_color     = boards_col0
color            = boards_col0

</plot>
<plot>
type = tile
file = %s
r1 = 0.65r
r0 = 0.56r
layers    = 1
margin    = 0.2u
thickness = 100
padding   = 0
orientation      = out
stroke_thickness = 0.01
stroke_color     = boards_col1
color            = boards_col1

</plot>

<plot>
type = tile
file = %s
r1 = 1.04r
r0 = 1.01r
layers    = 1
margin    = 0.2u
thickness = 50
padding   = 0
orientation      = out
stroke_thickness = 0.01
stroke_color     = forward_col
color            = forward_col
</plot>

<plot>
type = tile
file = %s
r1 = 0.98r
r0 = 0.95r
layers    = 1
margin    = 0.2u
thickness = 50
padding   = 0
orientation      = out
stroke_thickness = 0.01
stroke_color     = reverse_col
color            = reverse_col
</plot>

<plot>
type             = histogram
file             = %s
r1               = 0.93r
r0               = 0.75r
stroke_type      = outline
thickness        = 0.1
color            = vdgrey
extend_bin       = no
<rules>
<rule>
condition  = var(value) <= 0
fill_color = vdgrey
</rule>
<rule>
condition  = var(value) > 0
fill_color = spectral-9-div-1
</rule>
</rules>
</plot>
</plots>
"""
        return contents_4
    @staticmethod
    def return_contents_5( ) :
        contents_5 ="""\r\n<image>
<<include %s/etc/image.conf>>
</image>

<<include %s/etc/colors_fonts_patterns.conf>> 

<<include %s/etc/housekeeping.conf>>
"""
        return contents_5