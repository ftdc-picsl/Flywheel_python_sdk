import flywheel
import pandas as pd
import numpy as np
import datetime
import pytz
import re
import os
import pathlib

fw = flywheel.Client()

def get_bids_nifti(acq):
    bids_niftis = [f for f in acq.files if ('info' in f.keys() and 'BIDS' in f.info.keys() and "nii" in f.name)]
    if len(bids_niftis) > 0:
        return(bids_niftis.pop())
    else:
        return(None)

def get_t1_file(sess):
    '''
        Function to pick the most desirable T1 file out of several in a session.
    '''
    #is_t1 = [any(['T1' in f.classification['Measurement'] for f in a.files \
    #    if 'Measurement' in f.classification.keys()]) for a in sess.acquisitions()]
    #t1_acq = [a for (a, v) in zip(sess.acquisitions(), is_t1) if v]
    
    t1_acq = []
    acqlist = sess.acquisitions()
    for acq in acqlist:
        if any(['T1' in f.classification['Measurement'] for f in acq.files \
            if 'Measurement' in f.classification.keys()]):
                t1_acq.append(acq)
    
    t1_file = None
    
    for acq in t1_acq:
        lab = acq.label.lower()
        if ("vnav" in lab) and ("moco" in lab) and ("rms" in lab) and not ("nd" in lab):
            t1_file = get_bids_nifti(acq)
            return(t1_file)
    
    for acq in t1_acq:
        lab = acq.label.lower()
        if ("vnav" in lab) and ("rms" in lab) and not ("nd" in lab):
            t1_file = get_bids_nifti(acq)
            return(t1_file)
    
    for acq in t1_acq:
        lab = acq.label.lower()
        if ("ax" in lab) and ("mprage" in lab):
            t1_file = get_bids_nifti(acq)
            return(t1_file)
    
    for acq in t1_acq:
        lab = acq.label.lower()
        if ("sag" in lab) and ("mprage" in lab):
            t1_file = get_bids_nifti(acq)
            return(t1_file)
    
    return(t1_file)

def run_fmriprep(subjectLabel, sessionLabel, group = 'pennftdcenter', projectLabel = 'HUP6', ignore = '', t1_file = None):
    projectPath = '{}/{}'.format(group, projectLabel)
    proj = fw.lookup(projectPath)
    fmriprep = fw.lookup('gears/fmriprep-fwheudiconv')
    
    # Destination session.
    sess = fw.lookup('/'.join([group,projectLabel,str(subjectLabel),str(sessionLabel)]))
    
    # Get the index of the FreeSurfer license file among the project files.
    findex = [(i,f.name) for i,f in enumerate(proj.files) if 'license.txt' in f.name][0][0]
    now = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M")
    analysis_label = '{}_{}_{}_{}_{}'.format(str(subjectLabel), 
        str(sessionLabel), fmriprep.gear.name, fmriprep.gear.version, now)
    
    # Find the T1 file.
    if t1_file is None:
        t1_file = get_t1_file(sess)
    
    # If there's more than one T1w image, do what?
    inputs = {
        'freesurfer_license': proj.files[findex],
        't1w_anatomy': t1_file
    }
    
    # If there are no fieldmaps, do what?
    config = {
        'anat_only': False,
        'aroma_melodic_dimensionality': -200,
        'bold2t1w_dof': 6,
        'cifti_output': 'None',
        'dummy_scans': 0,
        'fd_spike_threshold': 0.5,
        'fmap_bspline': False,
        'fmap_no_demean': False,
        'force_bbr': False,
        'force_no_bbr': False,
        'force_syn': False,
        'fs_no_reconall': False,
        'ignore': ignore,
        'intermediate_files': '',
        'intermediate_folders': '',
        'longitudinal': False,
        'low_mem': False,
        'medial_surface_nan': False,
        'no_submm_recon': False,
        'no_track': False,
        'output_spaces': 'anat fsnative fsaverage MNI152NLin2009cAsym:res-2',
        'return_all_components': False,
        'save_intermediate_work': False,
        'save_outputs': True,
        'sge-cpu': '8',
        'sge-ram': '64G',
        'sge-short': False,
        'singularity-debug': False,
        'singularity-writable': False,
        'skip_bids_validation': False,
        'skull_strip_fixed_seed': False,
        'skull_strip_template': 'OASIS30ANTs',
        'sloppy_mode': False,
        't2s_coreg': False,
        'task_id': '',
        'timeout': 2,
        'use_all_sessions': False,
        'use_aroma': False,
        'use_syn_sdc': False
        }
    
    result = None
    
    try:
        anz_id = fmriprep.run(config = config, analysis_label = analysis_label, 
            tags = None, destination = fw.get(sess.id), inputs = inputs)
        result = anz_id
    except Exception as e:
        print(e)
        result = e
    
    return result