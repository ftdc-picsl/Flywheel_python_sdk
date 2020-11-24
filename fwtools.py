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

def get_latest_fmriprep(session):
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
            if 'fmriprep' in gear_name and date > latest_date and state =='complete':
                latest_date = date
                latest_run = i
        
        if latest_run is not None:
            fmriprep_out = [x for x in latest_run.files if 'fmriprep' in x.name][1]
            fmriprep_out
            return(fmriprep_out)
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

def run_fmriprep(projectLabel, subjectLabel, sessionLabel, group = 'pennftdcenter', ignore = '', t1_file = None):
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

def run_xcp(projectLabel, subjectLabel, sessionLabel, group = 'pennftdcenter', designFile = 'fc-36p_despike.dsn'):
    projectPath = '{}/{}'.format(group, projectLabel)
    proj = fw.lookup(projectPath)
    
    # xcpEngine design file: configures the moduules to be run.
    dsn = [x for x in proj.files if designFile in x.name][0]
    dsn = proj.get_file(dsn.name)
    
    # Destination session.
    sess = fw.lookup('/'.join([group,projectLabel,str(subjectLabel),str(sessionLabel)]))
    
    # Get the latest fmriprep path.
    fmriprep = get_latest_fmriprep(sess)
    
    if fmriprep is not None:
        # Get the xcpEngine gear.
        xcp = fw.lookup('gears/xcpengine-fw')
        
        # Create a label for the analysis.
        now = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M")
        analysis_label = '{}_{}_{}_{}_{}'.format(str(subjectLabel), 
            str(sessionLabel), xcp.gear.name, xcp.gear.version, now)
        
        myconfig = {
            'analysis_type': 'fc'
        }
        
        myinput = {
            'fmriprepdir': fmriprep,
            'designfile': dsn,
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


def get_zip_member(project, subject, session, regexp_anz, regexp_member, outPath, regexp_zip = "", group = 'pennftdcenter'):
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
    
    sesspath = '{}/{}/{}/{}'.format(group, project, str(subject), str(session))
    sess = fw.get(fw.lookup(sesspath).id)
    results = []
    
    if sess.analyses is not None:
        # Search through session for analyses that match regexp_anz. If there aren't
        # any, return None and/or a warning message.
        # Is fw.resolve faster than fw.lookup or fw.get?
        anzlist = [a.id for a in sess.analyses if re.match(regexp_anz, a.label)]
        
        if len(anzlist) > 0:
        # For each matching analysis, search by regexp_member through attached files.
            for anz in anzlist:
                anz = fw.get(anz)
                flist = [f.name for f in anz.files if re.match(regexp_zip, 
                    f.name) and 'zip' in f.name]
                if len(flist) > 0:
                    for f in flist:
                        zip_info = anz.get_file_zip_info(f)
                        member_matches = [m.path for m in zip_info.members if \
                            re.match(regexp_member, m.path)]
                        if len(member_matches) > 0:
                            for zm in member_matches:
#                                if os.path.isfile(outPath):
#                                    outd = os.path.dirname(outPath)
#                                    if not os.path.isdir(outd):
#                                        pathlib.Path(outd).mkdir(parents=True,
#                                            exist_ok=True)
#                                    print(str(subject),str(session),zm)
#                                    anz.download_file_zip_member(f, zm, outPath)
#                                    results.append((anz.label,f,zm,outPath))
#                                elif os.path.isdir(outPath):
                                # If the output directory doesn't already
                                # exist, create it.
                                if not os.path.isdir(outPath):
                                    pathlib.Path(outPath).mkdir(parents=True,
                                        exist_ok=True)
                                fpath = '{}/{}'.format(outPath,
                                    os.path.basename(zm))
                                print(str(subject),str(session),zm)
                                anz.download_file_zip_member(f, zm, fpath)
                                results.append((anz.label,f,zm,fpath))
    if results:
        results = pd.DataFrame(results, columns = ['analysis','file',
            'zip_member','destination'])
        return(results)
    else:
        return(None)

