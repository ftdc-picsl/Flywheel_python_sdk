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
        is_t1 = [any(['T1' in f.classification['Measurement'] for f in a.files \
            if 'Measurement' in f.classification.keys()]) for a in sess.acquisitions()]
        t1_acq = [a for (a, v) in zip(sess.acquisitions(), is_t1) if v]
        t1_acq = t1_acq.pop()
        t1_file = [f for f in t1_acq.files if ('info' in f.keys() and 'BIDS' in f.info.keys())].pop()
    
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