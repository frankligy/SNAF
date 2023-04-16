from matplotlib import cm
import pandas as pd
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_hex, to_rgb, to_rgba
from matplotlib import colors
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl
import os,sys

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# ASCI escape codes
def prRed(skk): return "\033[91m {}\033[00m" .format(skk)
def prGreen(skk): return "\033[92m {}\033[00m" .format(skk)
def prYellow(skk): return "\033[93m {}\033[00m" .format(skk)
def prLightPurple(skk): return "\033[94m {}\033[00m" .format(skk)
def prPurple(skk): return "\033[95m {}\033[00m" .format(skk)
def prCyan(skk): return "\033[96m {}\033[00m" .format(skk)
def prLightGray(skk): return "\033[97m {}\033[00m" .format(skk)
def prBlack(skk): return "\033[98m {}\033[00m" .format(skk)

# test_discrete_look
def generate_block(color_list,name):
    '''
    Given a list of color (each item is a hex code), visualize them side by side. See example.
    '''
    n = len(color_list)
    strip = np.empty(shape=(1,256),dtype='<U7')
    splitted = np.array_split(np.arange(strip.shape[1]),n)
    for i,c in enumerate(color_list):
        strip[:,splitted[i]] = c
    block = np.repeat(strip,10,axis=0)
    block_rgb = hex2_to_rgb3(block)
    fig,ax = plt.subplots()
    ax.imshow(block_rgb)
    ax.axis('off')
    ax.set_title('{}'.format(name))
    plt.savefig('{}_block.pdf'.format(name),bbox_inches='tight')
    plt.close()

# test_cmap_look
def generate_gradient(cmap,name):
    '''
    Given a continuous cmap, visualize them. See example.
    '''
    import numpy as np
    import matplotlib.pyplot as plt

    gradient = np.linspace(0, 1, 256).reshape(1,-1)
    gradient = np.repeat(gradient,10,axis=0)

    fig,ax = plt.subplots()
    ax.imshow(gradient,cmap=cmap)
    ax.axis('off')
    ax.set_title('{}'.format(name))
    plt.savefig('{}_gradient.pdf'.format(name),bbox_inches='tight')
    plt.close()

# background greyed colormap 
def bg_greyed_cmap(cmap_str):
    '''
    set 0 value as lightgrey, which will render better effect on umap

    :param cmap_str: string, any valid matplotlib colormap string

    :return: colormap object

    Examples::

        # normal cmap
        sc.pl.umap(sctri.adata,color='CD4',cmap='viridis')
        plt.savefig('normal.pdf',bbox_inches='tight')
        plt.close()

        # bg_greyed cmap
        sc.pl.umap(sctri.adata,color='CD4',cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
        plt.savefig('bg_greyed.pdf',bbox_inches='tight')
        plt.close()

    .. image:: ./_static/normal.png
        :height: 300px
        :width: 300px
        :align: left
        :target: target     

    .. image:: ./_static/bg_greyed.png
        :height: 300px
        :width: 300px
        :align: right
        :target: target    
    '''
    # give a matplotlib cmap str, for instance, 'viridis' or 'YlOrRd'
    cmap = copy.copy(cm.get_cmap(cmap_str))
    cmap.set_under('lightgrey')
    return cmap


# hex color 2d array, to (M,N,3) RGB array, used in imshow (plot_long_heatmap)
def hex2_to_rgb3(hex2):
    '''
    convert a hex color 2d array to (M,N,3) RGB array, very useful in ``ax.imshow``
    '''
    rgb3 = np.empty([hex2.shape[0],hex2.shape[1],3])
    for i in range(hex2.shape[0]):
        for j in range(hex2.shape[1]):
            hex_ = hex2[i][j]
            rgb_ = to_rgb(hex_)
            rgb3[i,j,:] = rgb_ 
    return rgb3



# 256 to [0,1]
def inter_from_256(x):
    return np.interp(x=x,xp=[0,255],fp=[0,1])

# [0,1] to 256
def infer_to_256(x):
    return int(np.interp(x=x,xp=[0,1],fp=[0,255]))


# choose colors
def retrieve_pretty_colors(name):
    '''
    retrieve pretty customized colors (discrete)

    :param name: string, valid value 'icgs2', 'shap'

    :return: list, each item is hex code

    Examples::

        generate_block(color_list = retrieve_pretty_colors('icgs2'),name='icgs2')
        generate_block(color_list = retrieve_pretty_colors('shap'),name='shap')

    .. image:: ./_static/colors.png
        :height: 100px
        :width: 550px
        :align: center
        :target: target      

    '''
    if name == 'icgs2':
        return _pub_icgs2
    elif name == 'shap':
        return _pub_shap


def retrieve_pretty_cmap(name):
    '''
    retrieve pretty customized colormap

    :param name: string, valid value 'altanalyze', 'shap', 'scphere'

    :return: cmap object

    Examples::

        generate_gradient(cmap=retrieve_pretty_cmap('shap'),name='shap')
        generate_gradient(cmap=retrieve_pretty_cmap('altanalyze'),name='altanalyze')
        generate_gradient(cmap=retrieve_pretty_cmap('scphere'),name='scphere')

    .. image:: ./_static/cmap.png
        :height: 250px
        :width: 550px
        :align: center
        :target: target   

    '''
    if name == 'altanalyze':
        return _ywb_cmap
    elif name == 'shap':
        return _pwb_cmap
    elif name == 'scphere':
        return _scphere_cmap

def pick_n_colors(n):
    '''
    a very handy and abstract function, pick n colors in hex code that guarantee decent contrast.

    1. n <=10, use tab10
    2. 10 < n <= 20, use tab20 
    3. 20 < n <= 28, use zeileis (take from scanpy)
    4. 28 < n <= 102, use godsnot (take from scanpy)
    5. n > 102, use jet cmap (no guarantee for obvious contrast)

    :param n: int, how many colors are needed

    :return: list, each item is a hex code.

    Examples::

        generate_block(color_list = pick_n_colors(10),name='tab10')
        generate_block(color_list = pick_n_colors(20),name='tab20')
        generate_block(color_list = pick_n_colors(28),name='zeileis')
        generate_block(color_list = pick_n_colors(102),name='godsnot')
        generate_block(color_list = pick_n_colors(120),name='jet')

    .. image:: ./_static/pick_n_colors.png
        :height: 300px
        :width: 550px
        :align: center
        :target: target      


    '''
    if n <= 10:
        _colors = [to_hex(color) for color in cm.get_cmap('tab10').colors[:n]]
    elif n > 10 and n <= 20:
        _colors = [to_hex(color) for color in cm.get_cmap('tab20').colors[:n]]
    elif n > 20 and n <= 28:
        _colors = _zeileis_28[:n]
    elif n > 28 and n <= 102:
        _colors = _godsnot_102[:n]
    elif n > 102:
        # _colors = [to_hex(cm.jet(round(i))) for i in np.linspace(0,255,n)]   # old way
        _colors = np.random.choice(r433,size=n,replace=False)
    return _colors

def colors_for_set(setlist):  # a list without redundancy
    '''
    given a set of items, based on how many unique item it has, pick the n color

    :param setlist: list without redundant items.

    :return: dictionary, {each item: hex code}

    Exmaples::

        cmap_dict = colors_for_set(['batch1','batch2])
        # {'batch1': '#1f77b4', 'batch2': '#ff7f0e'}

    '''
    length = len(setlist)
    _colors = pick_n_colors(n=length)
    cmap = pd.Series(index=setlist,data=_colors).to_dict()
    return cmap


def build_custom_continuous_cmap(*rgb_list):
    '''
    Generating any custom continuous colormap, user should supply a list of (R,G,B) color taking the value from [0,255], because this is
    the format the adobe color will output for you. 

    Examples::

        test_cmap = build_custom_continuous_cmap([64,57,144],[112,198,162],[230,241,146],[253,219,127],[244,109,69],[169,23,69])
        fig,ax = plt.subplots()
        fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(),cmap=diverge_cmap),ax=ax)

    .. image:: ./_static/custom_continuous_cmap.png
        :height: 400px
        :width: 550px
        :align: center
        :target: target     

    '''
    all_red = []
    all_green = []
    all_blue = []
    for rgb in rgb_list:
        all_red.append(rgb[0])
        all_green.append(rgb[1])
        all_blue.append(rgb[2])
    # build each section
    n_section = len(all_red) - 1
    red = tuple([(1/n_section*i,inter_from_256(v),inter_from_256(v)) for i,v in enumerate(all_red)])
    green = tuple([(1/n_section*i,inter_from_256(v),inter_from_256(v)) for i,v in enumerate(all_green)])
    blue = tuple([(1/n_section*i,inter_from_256(v),inter_from_256(v)) for i,v in enumerate(all_blue)])
    cdict = {'red':red,'green':green,'blue':blue}
    new_cmap = colors.LinearSegmentedColormap('new_cmap',segmentdata=cdict)
    return new_cmap

def build_custom_divergent_cmap(hex_left,hex_right):
    '''
    User supplies two arbitrary hex code for the vmin and vmax color values, then it will build a divergent cmap centers at pure white.

    Examples::

        diverge_cmap = build_custom_divergent_cmap('#21EBDB','#F0AA5F')
        fig,ax = plt.subplots()
        fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(),cmap=diverge_cmap),ax=ax)

    .. image:: ./_static/custom_divergent_cmap.png
        :height: 400px
        :width: 550px
        :align: center
        :target: target        

    '''
    left_rgb = colors.to_rgb(hex_left)
    right_rgb = colors.to_rgb(hex_right)
    # build each section
    n_section = 2
    red = ((0,left_rgb[0],left_rgb[0]),(0.5,1,1),(1,right_rgb[0],right_rgb[0]))
    green = ((0,left_rgb[1],left_rgb[1]), (0.5, 1, 1), (1, right_rgb[1], right_rgb[1]))
    blue = ((0,left_rgb[2],left_rgb[2]), (0.5, 1, 1), (1, right_rgb[2], right_rgb[2]))
    cdict = {'red':red,'green':green,'blue':blue}
    new_cmap = colors.LinearSegmentedColormap('new_cmap',segmentdata=cdict)
    return new_cmap


# zeileis_28 was taken from scanpy: https://github.com/theislab/scanpy/blob/master/scanpy/plotting/palettes.py
# and they noted the original source as below:
# https://graphicdesign.stackexchange.com/questions/3682/where-can-i-find-a-large-palette-set-of-contrasting-colors-for-coloring-many-d
# update 1
# orig reference http://epub.wu.ac.at/1692/1/document.pdf


_zeileis_28 = [
    "#023fa5",
    "#7d87b9",
    "#bec1d4",
    "#d6bcc0",
    "#bb7784",
    "#8e063b",
    "#4a6fe3",
    "#8595e1",
    "#b5bbe3",
    "#e6afb9",
    "#e07b91",
    "#d33f6a",
    "#11c638",
    "#8dd593",
    "#c6dec7",
    "#ead3c6",
    "#f0b98d",
    "#ef9708",
    "#0fcfc0",
    "#9cded6",
    "#d5eae7",
    "#f3e1eb",
    "#f6c4e1",
    "#f79cd4",
    # these last ones were added:
    '#7f7f7f',
    "#c7c7c7",
    "#1CE6FF",
    "#336600",
]

# godsnot_102 was taken from scanpy: https://github.com/theislab/scanpy/blob/master/scanpy/plotting/palettes.py
# the author noted the original source as below:
# from http://godsnotwheregodsnot.blogspot.de/2012/09/color-distribution-methodology.html
_godsnot_102 = [
    # "#000000",  # remove the black, as often, we have black colored annotation
    "#FFFF00",
    "#1CE6FF",
    "#FF34FF",
    "#FF4A46",
    "#008941",
    "#006FA6",
    "#A30059",
    "#FFDBE5",
    "#7A4900",
    "#0000A6",
    "#63FFAC",
    "#B79762",
    "#004D43",
    "#8FB0FF",
    "#997D87",
    "#5A0007",
    "#809693",
    "#6A3A4C",
    "#1B4400",
    "#4FC601",
    "#3B5DFF",
    "#4A3B53",
    "#FF2F80",
    "#61615A",
    "#BA0900",
    "#6B7900",
    "#00C2A0",
    "#FFAA92",
    "#FF90C9",
    "#B903AA",
    "#D16100",
    "#DDEFFF",
    "#000035",
    "#7B4F4B",
    "#A1C299",
    "#300018",
    "#0AA6D8",
    "#013349",
    "#00846F",
    "#372101",
    "#FFB500",
    "#C2FFED",
    "#A079BF",
    "#CC0744",
    "#C0B9B2",
    "#C2FF99",
    "#001E09",
    "#00489C",
    "#6F0062",
    "#0CBD66",
    "#EEC3FF",
    "#456D75",
    "#B77B68",
    "#7A87A1",
    "#788D66",
    "#885578",
    "#FAD09F",
    "#FF8A9A",
    "#D157A0",
    "#BEC459",
    "#456648",
    "#0086ED",
    "#886F4C",
    "#34362D",
    "#B4A8BD",
    "#00A6AA",
    "#452C2C",
    "#636375",
    "#A3C8C9",
    "#FF913F",
    "#938A81",
    "#575329",
    "#00FECF",
    "#B05B6F",
    "#8CD0FF",
    "#3B9700",
    "#04F757",
    "#C8A1A1",
    "#1E6E00",
    "#7900D7",
    "#A77500",
    "#6367A9",
    "#A05837",
    "#6B002C",
    "#772600",
    "#D790FF",
    "#9B9700",
    "#549E79",
    "#FFF69F",
    "#201625",
    "#72418F",
    "#BC23FF",
    "#99ADC0",
    "#3A2465",
    "#922329",
    "#5B4534",
    "#FDE8DC",
    "#404E55",
    "#0089A3",
    "#CB7E98",
    "#A4E804",
    "#324E72",
]

# r433 is generated in R using the code below
'''
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
library('gplots')
hex_vector = c()
for (c in color) {
  h = col2hex(c)
  hex_vector = append(hex_vector,h)
}

hex_matrix = t(as.matrix(hex_vector))
write.table(hex_matrix,'433colorhex.txt',sep='\t',row.names=F,col.names=F)
'''
r433 = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),'433colorhex.txt'),sep='\t',header=None).iloc[0,:].tolist()



_pub_icgs2 = [
    '#F26D6D',  # red
    '#BF9004',  # brown
    '#62BF04',  # blue
    '#2BB7EC',  # cyan
    '#A38BFD',  # purple
    '#F263DA',  # pink
]

_pub_shap = [
    '#F2075D',   # red
    '#158BFB',    # blue
]



'''
below stores the nice cmap I encoutered in my research
'''

# Nathan's Yellow-blue schema
# yellow-blue colormap
cdict = {
    'red':((0.0,0.0,0.0),
           (0.5,0.0,0.0),
           (1.0,1.0,1.0)),
    'green':((0.0,0.8,0.8),
             (0.5,0.0,0.0),
             (1.0,1.0,1.0)),
    'blue':((0.0,1.0,1.0),
            (0.5,0.0,0.0),
            (1.0,0.0,0.0))
}

_ywb_cmap = LinearSegmentedColormap('yellow_blue',segmentdata=cdict)



# SHAP pink-blue schema
cdict = {'red':((0.0,0.0,0.0),
                (1.0,1.0,1.0)),
         'green':((0.0,0.5,0.5),
                  (0.73,0.0,0.0),
                  (1.0,0.0,0.0)),
         'blue':((0.0,1.0,1.0),
                 (1.0,0.0,0.0))}
_pwb_cmap = LinearSegmentedColormap('shap', segmentdata=cdict)


# scPhere confusion matrix schema
cdict = {'red':((0.0,0.43,0.43),   # red chrome is 0.43 around (both left and right) 0.0 point, then increase to
                (0.45,1.0,1.0),    # 1.0 around 0.45 point, finally arrive
                (1.0,0.95,0.95)),  # 0.95 around 1.0 point
         'green':((0.0,0.61,0.61),
                  (0.45,1.0,1.0),
                  (1.0,0.27,0.27)),
          'blue':((0.0,0.85,0.85),
                  (0.4,0.96,0.96),
                  (1.0,0.18,0.18))}
_scphere_cmap = LinearSegmentedColormap('scphere', segmentdata=cdict)







    
