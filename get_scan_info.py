#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Helper function to identify scans, conditions, TR, and dims for nipype setups.

Version: 0.1
@author: neurosutton
Date: June 2019
"""

def get_scan_info(base_dir, task_name, param_dir='scan_params'):
    """Get info for a given study

    Parameters
    ----------
    base_dir : string
        Path to base directory of the dataset
    task_name : string
        Which task to process
    param_dir : string
        Directory name where generic, study json, xlsx, and text files are stored.


    Returns
    -------
    conditions : list of str
        Condition names
    onsets : list of ints
        Onsets corresponding to conditions
    durations : list of ints
        Durations associated with onsets/conditions
    TR : float
        Repetition time
    t1_dims: list of ints
        Anatomical dimensions to use during normalization.
        
    """

    from glob import glob
    import os.path as op
    import numpy as np
    import pandas as pd  # Will be used to modify and load conditions.
    
    condition_info = []
    cond_file = op.join(glob(base_dir, param_dir, task_name + '_params*xlsx')[0])
    
 
    json_info = os.path.join(base_dir, subject_id, 'BOLD', 'task%03d_run%03d' %
                             (task_name,
                              run_ids[task_name - 1][0]), 'bold_scaninfo.json')
    if os.path.exists(json_info):
        import json
        with open(json_info, 'rt') as fp:
            data = json.load(fp)
            TR = data['global']['const']['RepetitionTime'] / 1000.
            
    return  conditions, onsets, durations, tr, t1_dims