import flywheel
import pandas as pd
import numpy as np
import datetime
import pytz
import re
import os
import pathlib

fw = flywheel.Client()

def list_proj(projectPath):
    project = fw.lookup(projectPath)
    subjects = project.subjects()
    df = pd.DataFrame()
    for subject in subjects:
        sessions = subject.sessions()
        for sess in sessions:
            acq = sess.acquisitions()
            tmp = pd.DataFrame.from_dict({ 'subject': [subject['label']], 
                'firstname': [subject['firstname']], 'lastname': [subject['lastname']],
                'created': [sess['created']], 'label': [sess['label']],
                'acq': [[a['label'] for a in acq]], 'n_acq': [len(acq)] })
            df = pd.concat([df,tmp], axis = 0)
    df.index = range(df.shape[0])
    return(df)

def rename_sessions(metadata, project):
    md = metadata
    md['sesspath'] = 'pennftdcenter/' + project + '/' + md.subject + '/' + md.label
    for i in range(md.shape[0]):
        sess = fw.lookup(md.sesspath[i])
        sess.update(label = md.new_label[i])

def get_bids_nifti(acq):
    '''
        Returns a BIDS-format NIfTI T1 image from an acquisition.
    '''
    bids_niftis = [f for f in acq.files if ('info' in f.keys() and 'BIDS' in f.info.keys() and "nii" in f.name)]
    if len(bids_niftis) > 0:
        return(bids_niftis.pop())
    else:
        return(None)

def get_t1_file(sess):
    '''
        Function to pick the most desirable T1 file out of several in a session. Very, very FTDC-specific.
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

def get_latest_fmriprep(session, stateType = ['complete'], outputType = 'analysis'):
    '''
    Finds the latest available successful fmriprep job in a session.
    Written by Azeez Adebimpe.
    '''
        
    if session.analyses:
        
        timezone = pytz.timezone("UTC")
        init_date = datetime.datetime(2018, 1, 1)
        latest_date = timezone.localize(init_date)
        
        latest_run = None
        
        for i in session.analyses:
            gear_name = i.gear_info['name']
            state = i.job.state
            date = i.created
            if 'fmriprep' in gear_name and date > latest_date and state in stateType:
                latest_date = date
                latest_run = i
        
        if latest_run is not None:
            if outputType == 'analysis':
                return(latest_run)
            elif outputType == 'job':
                return(latest_run.job)
            elif outputType == 'file':
                fmriprep_out = [x for x in latest_run.files if ('zip' in x.name and not 'html' in x.name)].pop()
                return(fmriprep_out)
        else:
            return None
        
    else:
        return None

def get_latest_analysis(session, matchString = 'fmriprep', stateType = ['complete'], outputType = 'analysis'):
    '''
    Finds the latest available job from a particular gear (which user indicates by matchString) in a session.
    Based on Azeez Adebimpe's code.
    '''
        
    if session.analyses:
        
        timezone = pytz.timezone("UTC")
        # Just need an arbitrary start date for searching for analyses.
        # Since we adopted Flywheel in 2019-2020, 1/1/2018 will do.
        init_date = datetime.datetime(2018, 1, 1)
        latest_date = timezone.localize(init_date)
        
        latest_run = None
        
        for i in session.analyses:
            gear_name = i.gear_info['name']
            state = i.job.state
            date = i.created
            if matchString in gear_name and date > latest_date and state in stateType:
                latest_date = date
                latest_run = i
        
        if latest_run is not None:
            if outputType == 'analysis':
                return(latest_run)
            elif outputType == 'job':
                return(latest_run.job)
            elif outputType == 'file':
                zipFile = [x for x in latest_run.files if ('zip' in x.name and not 'html' in x.name)].pop()
                return(zipFile)
        else:
            return None
        
    else:
        return None

def find_modality_files(projectPath, modality):
    ''' Function for finding all of the files of a given modality (e.g., "bold")
    in a project. Returns a Pandas dataframe with key metadata.
    '''
    project = fw.lookup(projectPath)
#    acq_list = []
    df = pd.DataFrame()
    sessions = project.sessions()
    for sess in sessions:
        for acq in sess.acquisitions():
            for f in acq.files:
                if 'nii' in f.name:
                    tmp = pd.DataFrame.from_dict({ 
                        'subject': [sess.subject.label], 
                        'session': [sess.label],
                        'acquisition_label': [acq.label],
                        'acquisition_id': [acq.id],
                        'file': f.name, 
                        'classification': [str(f.classification)] })
                    df = pd.concat([df,tmp], axis = 0)                    
#                if modality in f.classification['
#                if ('BIDS' in f.info.keys() and 
#                    isinstance(f.info['BIDS'], dict) and
#                    modality in f.info['BIDS']['filename'] and
#                    'nii' in f.name) or
#                    ('classification' in f.keys() and
#                    'Modality' in f.classification.keys() and
#                    'nii' in f.name):
#                    acq_list.append(fw.get(acq.id))
    df.index = range(df.shape[0])
    return(df)

def run_dcm2niix(projectLabel, subjectLabel, sessionLabel, group = 'pennftdcenter'):
    fw = flywheel.Client()
    g = fw.lookup('gears/dcm2niix')
    s = fw.lookup('/'.join([group,projectLabel,str(subjectLabel),str(sessionLabel)]))
    config = {
        'bids_sidecar': 'y'
        }
    # Create a list of acquisitions to process.
    input_list = []    
    for acq in s.acquisitions():
        fnames = [f.name for f in acq.files]
        if not (any(['.nii.gz' in fn for fn in fnames]) or 'setter' in acq.label.lower()):
            input_list.append(acq)
    
    results = None
    
    if len(input_list) > 0:
        # Propose the batch
        proposal = g.propose_batch(input_list, config = config)
        try:
            results = proposal.run()
        except Exception as e:
            print(e)
            results = e    
    
    return results

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

def run_fw_fmriprep(projectLabel, subjectLabel, sessionLabel, group = 'pennftdcenter', gearName = 'bids-fmriprep', ignore = '', t1_file = None):
    projectPath = '{}/{}'.format(group, projectLabel)
    proj = fw.lookup(projectPath)
    fmriprep = fw.lookup('gears/' + gearName)
    
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
        'skip-bids-validation': False,
        'anat-only': False,
        'error-on-aroma-warnings': False,
        'verbose': 'vv',
        'ignore': '',
        'longitudinal': False,
        'output-spaces': 'anat fsnative fsaverage MNI152NLin2009cAsym:res-2',
        'bold2t1w-init': 'register',
        'bold2t1w-dof': 6,
        'force-bbr': False,
        'force-no-bbr': False,
        'medial-surface-nan': False,
        'dummy-scans': 0,
        'use-aroma': False,
        'aroma-melodic-dimensionality': -200,
        'return-all-components': False,
        'fd-spike-threshold': 0.5,
        'dvars-spike-threshold': 1.5,
        'skull-strip-template': 'OASIS30ANTs',
        'skull-strip-fixed-seed': False,
        'skull-strip-t1w': force,
        'fmap-bspline': False,
        'fmap-no-demean': False,
        'use-syn-sdc': False,
        'force-syn': False,
        'no-submm-recon': False,
        'cifti-output': False,
        'fs-no-reconall': False,
        'resource-monitor': False,
        'reports-only': False,
        'write-graph': True,
        'stop-on-first-crash': False,
        'notrack': False,
        'debug': '',
        'gear-log-level': INFO,
        'gear-run-bids-validation': False,
        'gear-save-intermediate-output': False,
        'gear-intermediate-files': '',
        'gear-intermediate-folders': '',
        'gear-dry-run': False,
        'gear-keep-output': True,
        'gear-keep-fsaverage': False
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

def run_qsiprep(projectLabel, subjectLabel, sessionLabel, group = 'pennftdcenter', gearName = 'qsiprep-fw', 
    recon_spec = None, t1w_anatomy = None, **kwargs):
    '''
        Function to run qsiprep. Must specify project, subject, and session labels. Variables recon_spec and
        t1w_anatomy are None by default but can be set with Flywheel file objects. Additional keyword arguments
        (**kwargs) can be specified to modify qsiprep config values.
    '''
    projectPath = '{}/{}'.format(group, projectLabel)
    proj = fw.lookup(projectPath)
    g = fw.lookup('gears/' + gearName)
    
    # Destination session.
    sess = fw.lookup('/'.join([group,projectLabel,str(subjectLabel),str(sessionLabel)]))
    
    # Get the index of the FreeSurfer license file among the project files.
    findex = [(i,f.name) for i,f in enumerate(proj.files) if 'license.txt' in f.name][0][0]
    now = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M")
    analysis_label = '{}_{}_{}_{}_{}'.format(str(subjectLabel), 
        str(sessionLabel), g.gear.name, g.gear.version, now)
 

    inputs = {
        'freesurfer_license': proj.files[findex],
    }
    if recon_spec is not None:
        inputs['recon_spec'] = recon_spec
    
    if t1w_anatomy is not None:
        inputs['t1w_anatomy'] = t1w_anatomy
   
    config = {
        'b0_motion_corr_to': 'iterative',
        'b0_threshold': 100,
        'b0_to_t1w_transform': 'Rigid',
        'combine_all_dwis': True,
        'denoise_before_combining': True,
        'distortion_group_merge': 'none',
        'do_reconall': False,
        'dwi_denoise_window': 5,
        'dwi_no_biascorr': False,
        'fmap_bspline': False,
        'fmap_no_demean': False,
        'force_spatial_normalization': True,
        'force_syn': False,
        'hmc_model': 'eddy',
        'hmc_transform': 'Affine',
        'ignore': '',
        'impute_slice_threshold': 0,
        'intermediate_files': '',
        'intermediate_folders': '',
        'intramodal_template_iters': 0,
        'intramodal_template_transform': 'BSplineSyN',
        'longitudinal': False,
        'no_b0_harmonization': False,
        'notrack': False,
        'output_resolution': 1.7,
        'output_space': 'T1w',
        'recon_builtin': '',
        'save_intermediate_work': False,
        'save_partial_outputs': False,
        'shoreline_iters': 2,
        'skip_bids_validation': False,
        'skull_strip_fixed_seed': False,
        'skull_strip_template': 'OASIS',
        'sloppy': False,
        'template': 'MNI152NLin2009cAsym',
        'timeout': 2,
        'unringing_method': 'mrdegibbs',
        'use_all_sessions': False,
        'use_syn_sdc': False,
        'write_local_bvecs': False
    }
    
    for key, value in kwargs.items():
        config[key] = value
    
    result = None
    
    try:
        anz_id = g.run(config = config, analysis_label = analysis_label, 
            tags = None, destination = fw.get(sess.id), inputs = inputs)
        result = anz_id
    except Exception as e:
        print(e)
        result = e
    
    return result

def run_xcp(projectLabel, subjectLabel, sessionLabel, fmriprep = None, group = 'pennftdcenter', designFile = 'fc-36p_despike.dsn'):
    projectPath = '{}/{}'.format(group, projectLabel)
    proj = fw.lookup(projectPath)
    
    # xcpEngine design file: configures the modules to be run.
    dsn = [x for x in proj.files if designFile in x.name][0]
    dsn = proj.get_file(dsn.name)
    
    # Destination session.
    sess = fw.lookup('/'.join([group,projectLabel,str(subjectLabel),str(sessionLabel)]))
    
    # Get the latest fmriprep path.
    if fmriprep is None:
        fmriprep = get_latest_fmriprep(sess, outputType = "file")
    
    if fmriprep is not None:
        # Get the xcpEngine gear.
        xcp = fw.lookup('gears/xcpengine-fw')
        
        # Find all input preproc_bold images in the fmriprep output
        boldfiles = get_zip_member(zipfile = fmriprep, regexp_member = ".*space-T1w_desc-preproc_bold.nii.gz$",)
        
        # Create a single-subject cohort file on the fly and attach to the session.
        boldfiles['id0'] = subjectLabel.replace("-", "").replace(".","").replace("_","")
        boldfiles['id1'] = sessionLabel.replace("-", "").replace(".","").replace("_","")
        boldfiles['img'] = ""
        for i in boldfiles.index:
            boldfiles.loc[i,"img"] = "/".join(boldfiles.loc[i,"ZipMember"].split("/")[2:])
        cohort = boldfiles[["id0","id1","img"]]
        cohort.to_csv("cohortfile.csv", index = False)
        sess_files = [f.name for f in sess.files]
        if "cohortfile.csv" in sess_files:
            sess.delete_file("cohortfile.csv")
        sess.upload_file("cohortfile.csv")

        # Create a label for the analysis.
        now = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M")
        analysis_label = '{}_{}_{}_{}_{}'.format(str(subjectLabel), 
            str(sessionLabel), xcp.gear.name, xcp.gear.version, now)
        
        myconfig = {
            'analysis_type': 'fc',
            'session': sessionLabel,
            'space': 'T1w'
        }
        
        myinput = {
            'designfile': dsn,
            'fmriprepdir': fmriprep,
        }
        
        # Actually run the thing.
        result = None
        
        try:
            anz_id = xcp.run(analysis_label = analysis_label, destination = sess, 
                inputs = myinput, config = myconfig)
            result = anz_id
        except Exception as e:
            print(e)
            result = e
    
    else:
        result = None
    
    return result

def get_zip_member(zipfile = None, group = 'pennftdcenter', project = None, subject = None, session = None, regexp_zip = "", regexp_anz = "", regexp_member = ".*space-T1w_desc-preproc_bold.nii.gz$", outPath = "./", download = False):
    '''
    Function to find and download specific files from an analysis zip file.
    Relies on matching of user-specified regular expressions for both the 
    analysis and the zip file member to be downloaded.

    project: name of Flywheel project (string)
    subject: subject ID/label (string or numeric)
    session: session ID/label (string or numeric)
    regexp_anz: pattern used to identify relevant analyses
    regexp_zip: if, for efficiency's sake, you want to add an additional filter
        on which zip files get searched, use a regular expression here. Defaults
        to a blank string, meaning all zip files in the candidate analyses will
        be searched.
    regexp_member: pattern used to identify zip file members to download
    outPath: path to which the download will be saved. If a directory is given
        without a base file name, the file name as it exists in the zip file
        will be used.
    '''
    
    results = pd.DataFrame()
    
    if zipfile is None:
        sesspath = '{}/{}/{}/{}'.format(group, project, str(subject), str(session))
        sess = fw.get(fw.lookup(sesspath).id)
    
        if sess.analyses is not None:
            # Search through session for analyses that match regexp_anz. If there aren't
            # any, return None and/or a warning message.
            # Is fw.resolve faster than fw.lookup or fw.get?
            anzlist = [a.id for a in sess.analyses if re.match(regexp_anz, a.label)]
        
            if len(anzlist) > 0:
            # For each matching analysis, search by regexp_member through attached files.
                for anz in anzlist:
                    anz = fw.get(anz)
                    flist = [f for f in anz.files if re.match(regexp_zip, 
                        f.name) and 'zip' in f.name]
                    if len(flist) > 0:
                        for zipfile in flist:
                            search_output = file_gopher(zipfile, regexp_member, download, outPath)
                            if search_output is not None:
                                pd.concat([results, search_output], axis = 0)
    else:
        search_output = file_gopher(zipfile, regexp_member, download, outPath)
        if search_output is not None:
            results = pd.concat([results, search_output], axis = 0)
    return(results)

def file_gopher(zipFile = None, regexp_zip_member = ".*nii.gz$", download = False, outPath = "."):
    zip_info = zipFile.get_zip_info()
    member_matches = [zm.path for zm in zip_info.members if \
        re.match(regexp_zip_member, zm.path)]
    outList = []
    output_df = None
    if len(member_matches) > 0:
        output_df = pd.DataFrame()
        for zm in member_matches:
            # If the output directory doesn't already
            # exist, create it.
            if download:
                if not os.path.isdir(outPath):
                    pathlib.Path(outPath).mkdir(parents=True,
                        exist_ok=True)
                fpath = '{}/{}'.format(outPath,os.path.basename(zm))
                zipFile.download_zip_member(zm, fpath)
            else:
                fpath = None
            tmp = (zipFile.name, zm, fpath)
            outList.append(tmp)
    output_df = pd.DataFrame(outList, columns = ['ZipFile','ZipMember','OutputFile'])
    return(output_df)

def fix_job_id(id):
    if type(id) is str:
        return(id[:(len(id)-1)] + chr(ord(str(id[-1]))+1))
    else:
        return(None)

def bids_info_from_file(fwFileEntry, keyName):
    '''
    Function called by get_bids_value for accessing BIDS dictionary fields in a
    NIfTI file object. Not necessary to call this function directly--call
    get_bids_value instead.
    
    fwFileEntry: a Flywheel file entry object.
    keyName: a string corresponding to a key in the BIDS dictionary.
    '''
    val = None
    if 'info' in fwFileEntry.keys():
        if 'BIDS' in fwFileEntry.info.keys():
            if keyName in fwFileEntry.info['BIDS'].keys():
                val = fwFileEntry.info['BIDS'][keyName]
                return(val)

def get_bids_value(input, keyName):
    '''
    Function for accessing a BIDS field in a Flywheel acquisition or file object.
    If you supply an acquisition object, the function will look for a NIfTI
    file object inside the acquisition.
    
    input: a Flywheel acquisition or file object.
    keyName: a string corresponding to a key in the BIDS dictionary.
    '''
    val = None
    if type(input) == flywheel.models.file_entry.FileEntry:
        val = bids_info_from_file(input, keyName)
    elif type(input) == flywheel.models.acquisition.Acquisition:
        f = [f for f in input.files if '.nii' in f.name]
        if len(f) > 0:
            f = f.pop()
            val = bids_info_from_file(f, keyName)
    return(val)
    
def info_from_file(fwFileEntry, keyName):
    '''
    Function called by get_bids_value for accessing info dictionary fields in a
    NIfTI file object. Not necessary to call this function directly--call
    get_info_value instead.
    
    fwFileEntry: a Flywheel file entry object.
    keyName: a string corresponding to a key in the BIDS dictionary.
    '''
    val = None
    if 'info' in fwFileEntry.keys():
        if keyName in fwFileEntry.info.keys():
            val = fwFileEntry.info[keyName]
            return(val)

def get_info_value(input, keyName):
    '''
    Function for accessing an info field in a Flywheel acquisition or file object.
    If you supply an acquisition object, the function will look for a NIfTI
    file object inside the acquisition.
    
    input: a Flywheel acquisition or file object.
    keyName: a string corresponding to a key in the info dictionary.
    '''
    val = None
    if type(input) == flywheel.models.file_entry.FileEntry:
        val = info_from_file(input, keyName)
    elif type(input) == flywheel.models.acquisition.Acquisition:
        f = [f for f in input.files if '.nii' in f.name]
        if len(f) > 0:
            f = f.pop()
            val = info_from_file(f, keyName)
    return(val)
