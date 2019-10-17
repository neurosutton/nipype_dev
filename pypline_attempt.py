#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 09:23:21 2019

@author: brianne
"""

#!/usr/bin/env python
# coding: utf-8




# Get the Node and Workflow object
from nipype import Node, Workflow, pipeline
import nipype.interfaces.utility as util  # utility (Needed?)

# Specify which SPM to use (useful for the SPM8 comparison testing)
from nipype.interfaces.matlab import MatlabCommand as mlabcmd
mlabcmd.set_default_paths('/usr/local/MATLAB/tools/spm12')

# Use nipype's version of collecting and inputing files.
from nipype import SelectFiles, DataSink, config
#config.enable_debug_mode()
config.set('execution', 'stop_on_first_crash', 'true') # Doesn't mean the whole pipeline will run properly if set to false, but can run through a couple times and hopefully hit the stragglers. Sometimes, it's because the scan doesn' exist, but the template is not flexible enough to catch it.
#config.set('execution', 'keep_inputs', 'true')
#config.set('execution', 'keep_unnecessary_files','true')
config.set('execution', 'hash_method', 'timestamp')
#config.set('execution', 'poll_sleep_duration','3')

import os
import glob
import os.path as op

# For now, hard code some of the paths that I can use to test the pipeline.
home_dir = op.abspath('/data/analysis/brianne/exobk') # Will get from tkinter eventually
paradigm = 'fp' # food pics
output_dir = op.abspath(op.join(home_dir, paradigm))  # Nipype may have trouble getting the correct path, if "abspath" is not used.
con_dir = op.join(home_dir,'ktl_fp_analysis/fasted_pre')
data_dir = op.abspath('/data/images/exobk') # Same potential issue of needing "abspath"


# #### Tip 1: Nipype needs "abspath" in some of its definitions or it will default to finding scans in expanded "~" folders

# Tip 15: If a pkl file won't update, try restarting the kernel.

# My way of getting the iterable list of subjects to pass to the workflow
def get_subject_list(search_directory, prefix):
    subjs=[]
    possible_subjs = glob.glob(os.path.join(search_directory,prefix+'*')) # change pattern to only * for most studies. Getting extraneoushits on EXOBK, b/c some participants withdrew.
    for subj in possible_subjs:
        subj = os.path.basename(subj)
        name = subj.split("_")[0]
    #   name = '_'.join(subj.split("_")[0:3]) # For more specificity
        subjs.append(name)
    subjs = sorted(list(set(subjs)))  # Preserve only the unique matches for the workflow to use
    return subjs

def paradigm_info(design_text_file):
    """Create the model design parameters based on a csv design file.
    The expected column titles for the csv are 'name' for the specific task (e.g., 'fp_run1'),'conditions' (e.g.,'highCal'), 'onsets',and 'durations'"""
    from nipype.interfaces.base import Bunch
    import pandas as pd
    try:
        df = pd.DataFrame(pd.read_csv(design_text_file))
        cols = ['onsets','durations']
        for col in cols:
            for ix,row in enumerate(df[col]):
                if ',' in row:
                    df.loc[ix,col] = [int(x) for x in str(row).strip('[]').strip().split(',')]
                else:
                    df.loc[ix,col] = [int(x) for x in str(row).strip('[]').strip().split()] 
        par_info={}
        
        for task in sorted(set(df['name'])):
            tmp = df.loc[df['name']==task,:] # limiting to the applicable data
            conditions = (tmp['conditions'].to_list())
            onsets = (tmp['onsets'].to_list())
            durations = (tmp['durations'].to_list())
            par_info[task] = Bunch(conditions = conditions, 
                    onsets = onsets,
                    durations = durations)
            
        return par_info
    except:
        # TODO test to see what types of errors may be common (lists in the durations, no column name, etc.)
        print ('Error in paradigm info.')


def contrast_info(contrast_text_file, con_level=1):
    """Create inputs for first- and second-level analyses using a csv definition file.
    Expected column titles are 'con_level' (1 or 2), title' ('low_greaterThan_baseline'),'con_type' ('T'),'conditions' (['high', 'low', 'crazy', 'baseline']),'matrix' ([0 1 0 -1].
    The returned list is appropriate for nipype contrast estimate inputs."""
    
    
# Tip 13   """The error states that the F-test must be “a tuple of the form: (a unicode string, ‘F’, a list of items which are a tuple of the form: (a unicode string, ‘T’, a list of items which are a unicode string, a list of items which are a float) or a tuple of the form: (a unicode string, ‘T’, a list of items which are a unicode string, a list of items which are a float, a list of items which are a float))”. I didn’t realize, until closer inspection, that the F-test tuple was building a matrix of contrasts, just like you would see if setting up the test in the SPM GUI (e.g., [1 0 0; 0 1 0; 0 0 1])."""
    
# Tip 14   The length of the conditions and weights must be the same. Either zero-pad or only list the included condition columns in the conditions list.
    
    import pandas as pd
    import re
    from ast import literal_eval
    try:
        df = pd.DataFrame(pd.read_csv(contrast_text_file))
        df = df.loc[df['con_level']==con_level,['title','con_type','conditions','matrix']]
        df = df.replace({'tcon':'T','fcon':'F','\'':''},regex=True)
        
        dfF = df.loc[df['con_type']=='F',:].reset_index(drop=True)
        dfF = dfF.replace({'T':"'T'"},regex=True)
        dfF.dropna(axis=1,inplace=True) # Since matrix is not defined for F-tests, eliminate the column from possible tuple conversion
        #for ix, row in enumerate(dfF['conditions']):
            #if ',' in row:
                # TODO F-testing does not work. Figure out how to import pieces systematically to build a list of tuples.
                #contrasts.append([(dfF['title'],dfF['con_type'], dfF['conditions'].values())])
                #dfF.loc[ix,'conditions'] = [tuple(x) for x in re.findall(r'\(([^\)]+)\)',dfF.loc[ix,'conditions'])]
        
        dfT = df.loc[df['con_type']=='T',:].reset_index(drop=True)
        dfT = dfT.replace({'\[':'',']':''},regex=True)
        for ix, row in enumerate(dfT['conditions']):
            if ',' in row:
                dfT.loc[ix,'conditions'] = [x for x in str(row).strip('[]').strip().split(',')]
            else:
                dfT.loc[ix,'conditions'] = [x for x in str(row).strip('[]').strip().split()]
        if not dfT['matrix'].str.contains(',').any():
            for ix,row in enumerate(dfT['matrix']):
                dfT.loc[ix,'matrix'] = [float(x) for x in str(row).strip('[]').strip().split()]
        contrasts = [x for x in dfT.itertuples(index=False, name=None)]
        #contrasts.append(tuple([x for x in dfF.itertuples(index=False, name=None)]))
        return contrasts
    except:
        # TODO test for the common exceptions and add useful error messages.
        print('Trouble in contrast info gathering')
        


timepts_dict = {'1':{'fmri':[1,2],'t1':[1]},'2':{'fmri':[3,4],'t1':[3]}}

for timept,scan_dict in timepts_dict.items():
    #Start creating the preprocessing workflow
    workflow_name ='spm_preproc' # Can be whatever description you desire, but the name is used for the output path (so make it useful)
    preproc_wf = Workflow(name=workflow_name,base_dir=op.abspath(op.join(output_dir)))
    
    # Make sure that there is a place to put the output, b/c nipype won't create the structure for you and the code will error out.
    if not os.path.isdir(os.path.join(output_dir,workflow_name)):
        os.mkdir(os.path.join(output_dir,workflow_name))
    else:
        print(os.path.join(output_dir,workflow_name))   
    
    infosource = Node(util.IdentityInterface(fields=['subj_id',
                                                     'task', 'timept','t1']), name='infosource')
    subjs = get_subject_list(data_dir, 'exo')
    print(subjs[-1:])
    tasks = ["fp_run1","fp_run2"] # TODO turn into paradigm object field
    timepts = scan_dict['fmri'] # access the embedded dictionary, rather than iterate over every t1 and fmri scan combo
    t1_timepts = scan_dict['t1']
    infosource.iterables = [("subj_id",subjs[-20:]),("task", tasks),('timept',timepts),('t1',t1_timepts)]
    # Tip 15: iterables run full combinations. Be sneaky with using a dictionary of limited iterables, if you are only trying to select combinations. For example, having the T1 be iterable with timepoints [1,1,3,3] and the fmri scans be iterable with [1,2,3,4] will create 16 results, instead of the intended 2 combinations (wher 1 went with 1, t1=1 went with fmri=2, t1=3 went with fmri=3,etc.)
    
    # Create the templates for the scan locations
    # The values from the template dictionary are added onto the end of the base_directory, when the templates key is invoked by the Workflow coneection
    templates = {'struct':'{subj_id}*{t1}*/t1/{subj_id}*t1*.nii',
                'func':'{subj_id}*_{timept}_*/{task}/{subj_id}*_????????.nii'}
    
    # Setup the options for the SelectFiles command. SelectFiles seems to be a bit more straightforward than DataGrabber and still invokes glob...
    sf = Node(SelectFiles(templates,
                         base_directory=op.abspath(data_dir),
                         raise_on_empty=False, # Don't worry about the potential filepath matches built off the templates that don't actually exist
                         sort_filelist = True),
             run_without_submitting=True, # Use this option to get around accidentally having too much processes accessing the same working directory b/c bigkahuna is too fast/good
             name='selectfiles')
    
    # Input the iterables for the experiment(s) to be analyzed.
    # N.B. usually these fields are set up in the sf definition as "sf.inputs.task", "sf.inputs.subj_id", etc. depending what you defined as your fields in the templates.
    # However, if you are using iterables to maximize the pipeline potential, the values are set in the sf.iterables structure.
    # sf.run().outputs.get() will show the {} fields from the templates dictionary as undefined, as iterables are only incorporated through MapNode (and iterfield) or during a Workflow process
    
    
    
    # #### Tip 2: Create the output folder, so that nipype will not error out on that technicality.
    # 
    # #### Tip 3: The I/O (and pipeline/iterable) arguments are basically setting up structures in the object. Therefore, "name" is a required argument, b/c it creates the {variable}.name field. Other arguments must also be named, but not necessarily by specific keyword. e.g., "timepts" is not a kw, but the {variable}.iterable will create {variable}.timept with the specified values.
    
    # ## Set up the basic pipeline
    
    
    
    
    # TODO add AC-PC check  RuntimeError("".join(result['traceback']))
    
    
    # Tip 11: Add any arguments to the Node as a field in the object (e.g., coreg). Nipype doesn't seem to register inputs in the Node definition... or at least the _report.rst makes it seem like the arguments may not be communicated.
    # Tip 12: Look at the _report/_report.rst and pyscript_{node_variable} to figure out what ran.
        
    import nipype.interfaces.spm as spm
    import nipype.algorithms.rapidart as ra  # artifact detection
    # SPM Realign (Estimate & Reslice - Reslice only the mean image)
    realign = Node(interface=spm.Realign(),
                   name="realign")
    realign.inputs.jobtype = 'estwrite'
    realign.inputs.register_to_mean=True
    realign.inputs.write_which= [0,1] # first field is for "reslice all"; 2nd is "reslice mean"
    
    # ART evaluation of realignment parameters
    art = Node(interface=ra.ArtifactDetect(),name='art')
    art.inputs.norm_threshold = 1
    art.inputs.zintensity_threshold = 5
    art.inputs.mask_type = 'spm_global'
    art.inputs.parameter_source = 'SPM'
    art.inputs.use_differences = [True, False]  # Reflects use difference for motion, not intensity
    art.inputs.use_norm = True
    
    # SPM Coregister (Estimate Only)
    coreg = Node(interface=spm.Coregister(), name="coregister")
    coreg.inputs.jobtype='estimate'
    
    # SPM Normalise (Estimate & Reslice)
    norm12 = Node(interface=spm.Normalize12(),name="normalize")
    #spm.Normalize12().input_spec.gm_output_type = [True, False, False] # Should be modulated normalised
    
    # Smooth
    smooth = Node(interface=spm.Smooth(),name="smooth")
    fwhmlist = [8]
    smooth.iterables = ('fwhm',fwhmlist) # Iterables seem to follow the convention of variable.input[key] = value of the iterable.
    
    
    # ## Set up the file management through ds to be able to wipe out unnecessary files?
    #  - There are all sorts of files labelled under the type of processing that produced them in each subject's folder, so I am unclear why this is a useful management step. It simply adds a "ds" folder next to the processing folders.
    
    
    
    
    # Collect key output files in a useful place
    ds = Node(interface=DataSink(),name="ds")
    ds.inputs.base_directory = op.abspath(output_dir)
    ds.substitutions = [('_subj_id',''),('_task',''),('_t1_','')]  
    
    
    # ## Connect the nodes.
    # #### Tip 4: Think of the connections like the "dependency" fields in SPM's batch editor.
    
    # TODO Add first-level estimations
    
    preproc_wf.connect([
        (infosource,sf,[('subj_id','subj_id'),('task','task'),('t1','t1'),('timept','timept')]),
        (infosource, ds, [('subj_id','container')]),
        (sf, realign, [('func', 'in_files')]),
        (realign, ds,[('realignment_parameters','motion.@param')]),
        (sf, coreg, [('struct', 'target')]),
        (realign, coreg, [('mean_image', 'source'),('modified_in_files', 'apply_to_files')]),  #Using modified_in_files b/c the estwrite jobtype is invoked with 'which_write'=[0,1], yielding no resliced files and foiling ensuing steps.
        (art, ds, [('outlier_files', 'art.@outliers'),('statistic_files', 'art.@stats')]),
        (sf, norm12, [('struct', 'image_to_align')]),
        (coreg, norm12, [('coregistered_files', 'apply_to_files')]),
        (norm12, smooth, [('normalized_files', 'in_files')]),
        (smooth, ds, [('smoothed_files','smoothed_scan')]),
        (realign, art, [('realignment_parameters', 'realignment_parameters')]), # Not sure if the realignment parameters that are sent to ds will forcibly include art-identified outliers or not.
        (norm12, art, [('normalized_files', 'realigned_files')]) 
    ])
    
    # Tip 6: Make sure that you don't call the datasink folder and the workflow the same thing. (e.g., "preproc") If you do, you overwrite the workflow as a simple object.
    # Tip 7: Commas must separate every tuple in the connect statement or you will get a generic "tuple" not callable error.
        
    # ## Check the flow of the connections.
    # - Are the pieces connected?
    # - Are they connected how you think you meant them to be?
    
    
    
    
    # Create preproc output graph
    preproc_wf.write_graph(graph2use='flat', format='png', simple_form=False)
    
    # Visualize the graph
    from PIL import Image
    img = Image.open(op.join(output_dir, 'spm_preproc/graph_detailed.png'))
    #img.show()
    
    
    # ## Run the pipeline (in parallel).
    # - Be prepared for A LOT of output (and probably crash reports).
    # 
    # #### Tip 5: Nipype tries to keep track of where the processing was sucessfully completed for the parameters set in the pipeline. If the process crashes, there should be a json with the completed steps stored somewhere in the participant's file.
    
    #preproc_wf.run() # Works
    #preproc_wf.run('PBS', plugin_args={'qsub_args': '-q many'})
#    try:
#        preproc_wf.run(plugin='MultiProc', plugin_args={'n_procs': 8})  #Doesn't work, b/c files are not copied into a realign folder (which is also not made [nodes.py lines 480ish]) and the processing stalls out.
#    except RuntimeError as rErr:
#        print(rErr)
#        print('Continuing with first levels as much as possible. Check the crashfile.')

# MultiProc didn't work until I stepped through _send_procs_to_workers with a breakpoint... needs a bunch of files or more sleep time???

# Tip 7: 'updatehash' is specifically to rerun with new working directory, not new analysis. Updatehash will just rename things, but doesn't run anything.

# MultiProc might not work in non-IPython environments. Suposition based on Notter's comment that, " For IPython the environment needs to be configured for distributed operation. Details are available at Using Nipype Plugins." And that MutilProc is from IPython

# Also seems that competing with Korey's processes for CPU is an issue. Is that related to timeout? Basically, if he is running stuff, MultiProc doesn't work correctly. If he is not running stuff, MultiProc is able to handle 20 cores and not error out unpredictably at some selectfiles or some realign or some other stage.


print("******* Part II: First-level estimates *****")
contrast_list = contrast_info(glob.glob(os.path.join(data_dir,'contrasts*.csv'))[0])  
info_lvl1 = Node(util.IdentityInterface(fields=['subj_id',
                                                 'task', 'timept',
                                                 'kernel',
                                                 'contrasts'], contrasts=contrast_list), name='info_lvl1')

subjs = get_subject_list(output_dir, 'exo')
print(len(subjs))
tasks = ["fp_run1","fp_run2"] # TODO turn into paradigm object field
timepts = [1,2,3,4]
kernels = fwhmlist
info_lvl1.iterables = [("subj_id",subjs[-20:]),("task", tasks),('timept',timepts[0:2]),('kernel',kernels)]

# Specify first level generically
from nipype.algorithms import modelgen
modelspec = Node(interface=modelgen.SpecifySPMModel(), name='modelspec')
modelspec.inputs.input_units = 'secs'
modelspec.inputs.time_repetition = 2 # TODO change to read from file header to make more flexible
modelspec.inputs.high_pass_filter_cutoff = 128
if len(tasks)>1:
    modelspec.inputs.concatenate_runs = False

# Load the parameters for the model design
templates_func = {}
par_info = paradigm_info(glob.glob(os.path.join(data_dir,'experimental_params*.csv'))[0]) # Returns dictionary with parameters

for t,task in enumerate(tasks):
    if not modelspec.inputs.subject_info:
        modelspec.inputs.subject_info = par_info[task]
    else:
        modelspec.inputs.subject_info.append(par_info[task])  # Rather than grabbing all the tasks associated with the study, just use the ones called early on in the script.
        
templates={'func':'{subj_id}/smoothed_scan/*{subj_id}*fp_run*{timept}/_fwhm_{kernel}/sw*.nii',
           'rp':'{subj_id}/motion/*{subj_id}*fp_run*{timept}/rp*txt'}
# Tip 16: templates create lists on their own, but do not accept lists for iterable input (i.e. don't specify run 1 and run 2). However, they do not append lists. SO if you iterate through tasks, only the final task will be passed on to the pipeline.
    
# Specify the functional runs
sf_func = Node(interface = SelectFiles(templates,
                     base_directory=output_dir,
                     sort_filelist = True),
         name='sf_level1files')
#sf_func.iterables = ('kernel',fwhmlist) 

# Specify the first level node
level1design = Node(interface=spm.Level1Design(), name="level1design")
level1design.inputs.timing_units = modelspec.inputs.output_units
level1design.inputs.interscan_interval = modelspec.inputs.time_repetition
level1design.inputs.bases = {'hrf': {'derivs': [0, 0]}}

level1estimate = Node(interface=spm.EstimateModel(), name="level1estimate")
level1estimate.inputs.estimation_method = {'Classical': 1}


# Setup the contrast estimation process
contrastestimate = Node(
    interface=spm.EstimateContrast(), name="contrastestimate")
contrastestimate.overwrite = True

l1pipeline = Workflow(name='level1')
l1pipeline.base_dir = output_dir


ds.inputs.base_directory = op.abspath(output_dir)
ds.substitutions = [('_subj_id',''),('_task',''),('_t1_',''), ('_kernel_8','')]  
   
l1pipeline.connect([
        (info_lvl1,sf_func,[('subj_id','subj_id'),('task','task'),('timept','timept'),('kernel','kernel')]),
        (sf_func,modelspec,[('func','functional_runs'),('rp','realignment_parameters')]),
        (modelspec,level1design,[('session_info','session_info')]),
        (info_lvl1, contrastestimate,[('contrasts','contrasts')]),
        (level1design,level1estimate,[('spm_mat_file','spm_mat_file')]),
        (level1estimate,contrastestimate,[('spm_mat_file','spm_mat_file'),
                                          ('beta_images','beta_images'),
                                          ('residual_image','residual_image')])
       ,(info_lvl1, ds, [('subj_id','container')]),
        (contrastestimate,ds, [('spm_mat_file','spm'),
                               ('spmT_images', '1stLevel.@T'),
                               ('con_images', '1stLevel.@con'),
                               ('spmF_images', '1stLevel.@F'),
                               ('ess_images', '1stLevel.@ess')])
])


# Create level 1 output graph
l1pipeline.write_graph(graph2use='flat', format='png', simple_form=False)

# Visualize the graph
img = Image.open(op.join(output_dir, 'level1/graph_detailed.png'))
#img.show()

l1pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 6})

# Figure out if steady state was reached for collection or if dummy volumes are still in the file
import nibabel as nb
import matplotlib.pyplot as plt
#print(func_file)
#plt.plot(nb.load(func_file[0]).get_fdata()[32,32,15,:]) # Pick a random voxel and plot raw signal across time.







