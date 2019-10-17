#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 13:21:36 2019

@author: brianne
"""
## From the fmri_spm tutorial
from __future__ import print_function
from builtins import str
from builtins import range

import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe
import nipype.algoritms.rapidart as ra
import nipype.algorithms.modelgen as model
import nipype.interfaces.matlab as mlab
from nipype.interfaces import spm

import os

# Specify where the data live
data_dir=
subj_list = []# Use my generic finder
timept_list = [1,2]
output_dir = 



subj_info = dict(
        func=[['subj_id',[run1, run2],'timept']],
        struct=[needed?]])
infosource = pe.Node(interface=util.IdentityInterface(fields=['subj_id']),name='infosource')
infosource.iterables('subj_id',subj_list,'timept',timept_list)


l2_wf = pe.Workflow(name,'level2_wf')
l2_wf.base_dir = output_dir

# Select the files
### Get the pre-subtracted con files.
templates={'func':'{subj_id}/smoothed_scan/*{subj_id}*{task}*{timept}/_fwhm_{kernel}/sw*.nii')
sf = Node(SelectFiles(templates,
                     base_directory=op.abspath(data_dir),
                     #raise_on_empty=False,
                     sort_filelist = True),
         #run_without_submitting=True, # Use this option to get around accidentally having too much processes accessing the same working directory b/c bigkahuna is too fast/good
         name='selectfiles')


# Define the stats
# OneSampleTTestDesign - creates one sample T-Test Design
onesamplettestdes = Node(OneSampleTTestDesign(),
                         name="onesampttestdes")
twosamplettest = Node(TwoSampleTTestDesign(),
                      name='twosample')


# EstimateModel - estimates the model
level2estimate = Node(EstimateModel(estimation_method={'Classical': 1}),
                      name="level2estimate")

# Get the pre-defined comparisons

cont1 = ('Inc. Activation', 'T',['Dieters','Exercisers'],[.5 .5])
cont2 = ('Exercise > diet','T',['Dieters','Exercisers'],[-1 1])
cont3 = ('Exercise < diet','T',['Dieters','Exercisers'],[1 -1])


# EstimateModel - estimates the model
level2estimate = Node(EstimateModel(estimation_method={'Classical': 1}),
                      name="level2estimate")

#Node: EstimateContrast - to estimate the first level contrasts we define later
l2conestimate = pe.Node(interface = spm.EstimateContrast(), name="contrastestimate")
l2conestimate.inputs.contrasts = [cont1, cont2, cont3]

# Threshold - thresholds contrasts
level2thresh = Node(Threshold(contrast_index=1,
                              use_topo_fdr=False,
                              use_fwe_correction=False,
                              extent_threshold=5,
                              height_threshold=0.005,
                              height_threshold_type='p-value',
                              extent_fdr_p_threshold=0.05),
                    name="level2thresh")

# Connect the nodes
l2_wf.connect([(infosource, selectfiles, [('contrast_id', 'contrast_id'),
                                               ('fwhm_id', 'fwhm_id')]),
                    (selectfiles, onesamplettestdes, [('cons', 'in_files')]),
                    (onesamplettestdes, level2estimate, [('spm_mat_file',
                                                          'spm_mat_file')]),
                    (level2estimate, level2conestimate, [('spm_mat_file',
                                                          'spm_mat_file'),
                                                         ('beta_images',
                                                          'beta_images'),
                                                         ('residual_image',
                                                          'residual_image')]),
                    (level2conestimate, level2thresh, [('spm_mat_file',
                                                        'spm_mat_file'),
                                                       ('spmT_images',
                                                        'stat_image'),
                                                       ])])