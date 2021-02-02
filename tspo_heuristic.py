'''
    Heuristic file for BIDS classification of Columbia TSPO/MRI data. The main
    function, infotodict, is defined below and takes sequence information from
    dicom headers; you can see which information is extracted by running
    fw-heudiconv-tabulate on the session directory, which writes the sequence
    info to a tsv file that can subsequently be read in as a Pandas dataframe.
    Each row of seqinfo corresponds to one series/acquisition in an imaging 
    session.
'''

import datetime
import numpy as np
import pandas
from collections import defaultdict

def create_key(template, outtype=('nii.gz',), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return template, outtype, annotation_classes

# Anat
## MPRAGE
t1w = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_run-{item}_T1w')
## FLAIR
flair = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_run-{item}_FLAIR')

# DWI
dti = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-32dir_dwi')

def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where
        allowed template fields - follow python string module:
        item: index within category
        subject: participant id
        seqitem: run number during scanning
        subindex: sub index within group
    """

    info = {
        t1w: [], flair: [], dti: []
    }

    for s in seqinfo:
        protocol = s.protocol_name.lower()
        mydatetime = datetime.datetime.strptime(s.date, '%Y-%m-%dT%H:%M:%S.%f')
        if "mprage" in protocol:
            info[t1w].append(s.series_id)
        elif "flair" in protocol:
            info[flair].append(s.series_id)
        elif "dti" in protocol:
            info[dti].append(s.series_id)
        else:
            print("Series not recognized!: ", protocol, s.dcm_dir_name)
            
    return info
    
def ReplaceSession(sesname):
    return sesname.replace("-", "").replace(".","").replace("_","")

def ReplaceSubject(subjname):
    return subjname.replace("-", "").replace(".","").replace("_","")

MetadataExtras = {
    # No metadata so far.
}       

IntendedFor = {
    # No intentions so far.
}
