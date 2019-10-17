#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 10:34:36 2019

@author: brianne
"""

#!/usr/bin/env python
# coding: utf-8


from nilearn.plotting import plot_anat
get_ipython().run_line_magic('matplotlib', 'inline')
import os.path as op
import pathlib
import json
from nipype.interfaces.fsl import (ExtractROI, ImageMaths, Threshold, fslinfo)
from nipype.interfaces.spm import (Smooth, Coregister, NewSegment, Normalize, Realign)
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.algorithms.rapidart import ArtifactDetect
from nipype.algorithms.confounds import (TSNR)
from nipype import Workflow, SelectFiles, Node

from pathlib import Path
from shutil import copy2 # copy2 includes metadata
from collections import defaultdict
import matplotlib.pyplot as plt

# Homegrown packages
import sys
sys.path
sys.path.append('/home/brianne/tools/python_dev')
from genius import file_gui

# Locate the scans
files = file_gui()
#print ([method for method in dir(files) if callable(getattr(files, method))]) #dbg
files = files.find_files()


# String-based templates
templates = {'anat': 'vol001/t1/0*nii',
            'func': 'vol001/tap*/00*nii'}

#templates = {'anat': 'openneuro/ds000254-me/sub-1*/anat/sub*nii',
#            'func': 'openneuro/ds000254-me/sub-1*/func/sub*nii'}

# Create SelectFiles node
sf = Node(SelectFiles(templates),
         name='selectfiles')

#Specify the unchanging pieces of the path
sf.inputs.base_directory = proj_dir

filelist=sf.run().outputs
filelist.func






# Magic! Build the default dict with lambda to achieve a callable nested dict easily.
img_dict=lambda:defaultdict(img_dict)
results = img_dict()

for f in filelist.func:
    print(f)
    tsnr=TSNR()
    tsnr.inputs.in_file=f
    scanDir=op.dirname(f)
    scan=op.basename(f).split('.')[0]
#    src = str(Path(scanDir,scan + '.nii'))
#    dest = str(Path(scanDir, "cp_" + scan + '.nii'))
#    copy2(src, dest) # To make sure you don't have to run the normalization again.

    #print(scanDir,name)
    tsnr.inputs.tsnr_file=str(Path(scanDir,"tsnr_"+ scan + '_tsnr.nii'))
    tsnr.inputs.stddev_file=str(Path(scanDir,"tsnr_" + scan + '_std.nii'))
    tsnr.inputs.mean_file=str(Path(scanDir,"tsnr_"+scan+'_mean.nii'))
    tsnr.inputs.detrended_file=str(Path(scanDir,"tsnr_" + scan + '_detrended.nii'))
    res = tsnr.run().outputs

    for key in ['stddev_file','mean_file','tsnr_file']:
        results[f][key] = eval('res.'+key) # eval is BAD practice in python, but to build this dictionary, the variable must be called dynamically



for k,v in results.items():
    tsnr_map=v['tsnr_file']
    comp_file=v['mean_file']
    plot_anat(comp_file, title=comp_file, display_mode='z', dim=-1, draw_cross=False)
    plot_anat(tsnr_map, title=tsnr_map, display_mode='z', dim=-1, draw_cross=False)





