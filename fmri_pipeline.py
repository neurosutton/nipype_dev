# -*- coding: utf-8 -*-
"""
Initial attempt at using nipype modules to construct a parallel pipeline to the MATLAB-based standardized pipeline.

Version: 0.1
@author: neurosutton
Date: June 2019
"""
from __future__ import unicode_literals

from glob import glob
import os
import pandas as pd
import numpy as np
import pylab as plt
import json

#  Basic IO for nipype
from nipype import Node, Workflow, SelectFiles, IdentityInterface

# Load the SPM-related modules
import nipype.interfaces.spm as spm 

# Import artifact detection
from nipype.algorithms.rapidart import ArtifactDetect

# Import plotting tools
from nilearn.plotting import plot_anat

import nipype.interfaces.matlab as mlab  

class setup_spm_analysis:
    def __init__(self, proj_dir=None, task_names=[], analysis_param_file=None, scan_param_dir='scan_params'):
        self.proj_dir = proj_dir # The top-level, working directory.
        self.task_names = task_names
        self.scan_param_dir = scan_param_dir
        self.analysis_param_file = analysis_param_file
        self.opts = type('', (), {})()
        if self.analysis_param_file:
            self.select_options()
        
    def __repr__(self):
        return
    
    def __str__(self):
        return
    
    def select_options(self, fwhm=2., art=None, wm_csf_reg=None, smooth_fwhm= [8,8,8]):
        """If the analysis has already been performed or the analyst sets up the parameters beforehand, select_options will load the settings."""
        
        if self.analysis_param_file:
            with open(os.path.join(self.proj_dir,os.path.basename(self.analysis_param_file))) as json_file:
                var_dict = json.load(json_file)
                for k,v in var_dict.items():
                    k=v # Provide the value for fwhm, if fwhm was defined in the file.

        self.opts.fwhm = fwhm # input is either the argument to the method OR the value reset by the var_dict[key]
        self.opts.wm_csf_reg = wm_csf_reg
        self.opts.smooth_fwhm = smooth_fwhm
        self.opts.art = art

        self.analysis_param_file = os.path.join(self.proj_dir,self.task_names[0]+"_analysis_params.json")   
        with open(self.analysis_param_file,'w') as json_file:
            json.dump(vars(self.opts),json_file)
    
        if self.opts.art:
            self.set_art_params()
        
    def set_scan_params(self):
        import get_scan_info
        [self.conditions, self.onsets, self.durations, self.tr, self.t1_dims] = get_scan_info(self.base_dir, task_names, self.scan_param_dir)
        
    def get_vox_dims(volume):
        import nibabel as nb
        from nipype.utils import NUMPY_MMAP
        if isinstance(volume, list):
            volume = volume[0]
        nii = nb.load(volume, mmap=NUMPY_MMAP)
        hdr = nii.header
        voxdims = hdr.get_zooms()
        return [float(voxdims[0]), float(voxdims[1]), float(voxdims[2])]
            
    def set_art_params(self, zintensity=3, norm_threshold=2, abs_mvmt=3):
        if self.analysis:
            self.analysis = type('', (), {})()
        self.analysis.art = Node(ArtifactDetect(norm_threshold=2,
                                  zintensity_threshold=2,
                                  mask_type='spm_global',
                                  parameter_source='FSL',
                                  use_differences=[True, False],
                                  save_plot=True,
                                  plot_type='svg'),
                   name="art")
 
    def analysis_steps(self):
        self.analysis = type('', (), {})()
        # Get files
        subj_list = [subj.split('_')[:-1] for subj in next(os.walk(self.proj_dir))[1]]
        # TODO limit the subj_list to those without sw processed files.
        
        # for parallelization by subject, use idnetityInterface
        self.analysis.infosource = Node(IdentityInterface(fields=['subj_id','task']),name="infosource")
        self.analysis.infosource.iterables = [('subject_id', subj_list),('task',self.task_names)]

        templates={'anat': '{subj_id}/t1/{subj_id}_t1*.nii',
                  'func': '{subj_id}/{task}*/{subj_id}_{task}*.nii'}
        self.analysis.sf = Node(SelectFiles(templates),
                  name='selectfiles')
        self.analysis.sf.inputs.base_directory = self.proj_dir
        
        
        # Realign
        self.analysis.realign = Node(spm.Realign(register_to_mean=True,
                               fwhm=self.opts.fwhm), name='realign')
        
        # Coregister
        self.analysis.coreg = Node(spm.Coregister(), name='coregistration')
        # Normalize
        self.analysis.norm12 = Node(spm.Normalize12(bias_regularization=1e-05,
                                                    affine_regularization_type='mni'), name='normalize')

        #Smooth
        self.analysis.smooth = Node(spm.Smooth(), name='smooth')
        #smooth.inputs.in_files = 'functional.nii'
        self.analysis.smooth.inputs.fwhm = self.opts.smooth_fwhm
    
    def analysis_connections(self, test='no'):
        wf = Workflow(name='spm_preproc',base_dir=self.proj_dir)
        # The way wf.connect args seem to be constructed (based on the flat graph outputs) is input_node, output_node, ([in_arg,matching_out_arg), (in_arg, matching_out_arg)].
        # The in_arg, out_arg are named fields, not actual file names.
        wf.connect([(self.analysis.infosource, self.analysis.sf, [('subj_id','subj_id'),('task','task')]), # Matched on multiple criteria
                    (self.analysis.sf,self.analysis.realign, [('func','in_files')]), # '"func" b/c the template for the fMRI scans for selectfiles was "func"
                    (self.analysis.realign, self.analysis.coreg, [('mean_image', 'source'), ('realigned_files',
                                                      'apply_to_files')]),
                    (self.analysis.sf, self.analysis.coreg, [('anat','target')]), # Subject T1 and functional to be coregistered
                    (self.analysis.sf, self.analysis.norm12, [('anat','image_to_align')]), #Subject t1 to be segmented and normalized per SPM12
                    (self.analysis.coreg, self.analysis.norm12, [('coregistered_files','apply_to_files')] ),
                    (self.analysis.norm12, self.analysis.smooth, [('normalized_files','in_files')])])
    
        if hasattr(self.analysis, 'art'):
            wf.connect([(self.analysis.realign, self.analysis.art, [('realignment_parameters', 'realignment_parameters')])])

        if test=='yes':
            # Create preproc output graph
            wf.write_graph(graph2use='colored', format='png', simple_form=True)
        else:
            wf.run('MultiProc', plugin_args={'n_procs': 4})
        
    def plot_images(self):
        f = plt.figure(figsize=(12, 4))
        for i, img in enumerate([]):
            f.add_subplots(1,4,i+1)
            plot_slice(img)
            
            # Plot the motion paramters
        par = np.loadtxt('/output/work_preproc/_subject_id_07/mcflirt/'
                         'asub-07_ses-test_task-fingerfootlips_bold_roi_mcf.nii.gz.par')
        fig, axes = plt.subplots(2, 1, figsize=(15, 5))
        axes[0].set_ylabel('rotation (radians)')
        axes[0].plot(par[0:, :3])
        axes[1].plot(par[0:, 3:])
        axes[1].set_xlabel('time (TR)')
        axes[1].set_ylabel('translation (mm)');
        
        # Showing the artifact detection output
        from IPython.display import SVG
        SVG(filename='/output/work_preproc/_subject_id_07/art/'
            'plot.asub-07_ses-test_task-fingerfootlips_bold_roi_mcf.svg')

test = setup_spm_analysis('/data/images/exobk',['fp'])
test.select_options()
test.analysis_steps()
test.analysis_connections(test='yes')
