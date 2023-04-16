#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import os,sys
import subprocess
from tqdm import tqdm
from sklearn.preprocessing import MinMaxScaler

class Bed_File_Builder():
    '''
    Example1:

    track name=pairedReads description="Clone Paired Reads" useScore=1
    chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
    chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399, 0,3601

    Example2:

    browser position chr7:127471196-127495720
    browser hide all
    track name="ItemRGBDemo" description="Item RGB demonstration" visibility=2 itemRgb="On"
    chr7    127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
    chr7    127472363  127473530  Pos2  0  +  127472363  127473530  255,0,0
    chr7    127473530  127474697  Pos3  0  +  127473530  127474697  255,0,0
    chr7    127474697  127475864  Pos4  0  +  127474697  127475864  255,0,0
    chr7    127475864  127477031  Neg1  0  -  127475864  127477031  0,0,255
    chr7    127477031  127478198  Neg2  0  -  127477031  127478198  0,0,255
    chr7    127478198  127479365  Neg3  0  -  127478198  127479365  0,0,255
    chr7    127479365  127480532  Pos5  0  +  127479365  127480532  255,0,0
    chr7    127480532  127481699  Neg4  0  -  127480532  127481699  0,0,255

    Example3:

    browser position chr7:127471196-127495720
    browser hide all
    track name="ColorByStrandDemo" description="Color by strand demonstration" visibility=2 colorByStrand="255,0,0 0,0,255"
    chr7    127471196  127472363  Pos1  0  +
    chr7    127472363  127473530  Pos2  0  +
    chr7    127473530  127474697  Pos3  0  +
    chr7    127474697  127475864  Pos4  0  +
    chr7    127475864  127477031  Neg1  0  -
    chr7    127477031  127478198  Neg2  0  -
    chr7    127478198  127479365  Neg3  0  -
    chr7    127479365  127480532  Pos5  0  +
    chr7    127480532  127481699  Neg4  0  -
    '''
    def __init__(self,chrom,chromStart,chromEnd,name=[None],score=[None],strand=[None],thickStart=[None],thickEnd=[None],itemRgb=[None],blockCount=[None],blockSizes=[None],blockStarts=[None],sep='\t',track_name='track_name',track_description='track_description',track_useScore=None,track_visibility=None,track_itemRgb=None):
        '''
        the arguments are either a list, or list of tuples
        '''

        self.chrom = chrom
        self.chromStart = chromStart
        self.chromEnd = chromEnd
        self.name = name
        self.score = score
        self.strand = strand
        self.thickStart = thickStart
        self.thickEnd = thickEnd
        self.itemRgb = itemRgb
        self.blockCount = blockCount
        self.blockSizes = blockSizes
        self.blockStarts = blockStarts

        self.stream = ''
        self.sep = sep

        self.track_name = track_name
        self.track_description = track_description
        self.track_useScore = track_useScore
        self.track_visibility = track_visibility
        self.track_itemRgb = track_itemRgb

    def build(self):
        self.populate_none()
        self.preprocess()
        self.build_header()
        self.build_tracks()

    def populate_none(self):
        length = len(self.chrom)
        optional = [self.name,self.score,self.strand,self.thickStart,self.thickEnd,self.itemRgb,self.blockCount,self.blockSizes,self.blockStarts]
        for field in optional:
            if field == [None]:
                field *= length
            

    def preprocess(self):
        self.score = np.squeeze(MinMaxScaler(feature_range=(0,1000)).fit_transform(np.array(self.score).reshape(-1,1))).tolist()

    def build_header(self):
        header = 'track name="{}" description="{}" '.format(self.track_name,self.track_description)
        if self.track_useScore is not None:
            header += 'useScore={} '.format(self.track_useScore)
        if self.track_visibility is not None:
            header += 'visibility={} '.format(self.track_visibility)
        if self.track_itemRgb is not None:
            header += 'itemRgb="{}" '.format(self.track_itemRgb)
        header = header.rstrip(' ') + '\n'
        self.stream += header

    def build_tracks(self):
        tracks = ''
        for chrom,chromStart,chromEnd,name,score,strand,thickStart,thickEnd,itemRgb,blockCount,blockSizes,blockStarts in zip(self.chrom,self.chromStart,self.chromEnd,self.name,self.score,self.strand,self.thickStart,self.thickEnd,self.itemRgb,self.blockCount,self.blockSizes,self.blockStarts):
            required = [chrom,chromStart,chromEnd]
            optional = [name,score,strand,thickStart,thickEnd,itemRgb,blockCount,blockSizes,blockStarts]
            for field in required:
                tracks += '{}{}'.format(field,self.sep)
            for field in optional:
                if field is not None:
                    tracks += '{}{}'.format(field,self.sep)
            tracks += '\n'
        self.stream += tracks.rstrip('\n')
    

# main program
bed = pd.read_csv('/data/salomonis2/NCI-R01/Signatures/encode_eclip_bed/raw/HepG2/idr_bed/PCBP2.bed',sep='\t',header=None)
chrom = bed[0].tolist()
chromStart = bed[1].tolist()
chromEnd = bed[2].tolist()
name = bed[3].tolist()
score = bed[11].tolist()
strand = bed[5].tolist()
aars_bfb = Bed_File_Builder(chrom=chrom,chromStart=chromStart,chromEnd=chromEnd,name=name,score=score,strand=strand,track_useScore=1)
aars_bfb.build()
with open('PCBP2_HepG2.bed','w') as f:
    f.write(aars_bfb.stream)




            


