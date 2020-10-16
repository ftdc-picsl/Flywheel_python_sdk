import datetime
def create_key(template, outtype=('nii.gz',), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return template, outtype, annotation_classes

# anatomical images
t1w = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_T1w')
t1wn = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-norm_T1w')
t2w = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_T2w')
t2wn = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-norm_T2w')

# field maps
fm_ap1 = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-AP_epi')
fm_pa1 = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-PA_epi')

# functional scans
nback_ap = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-WM_acq-AP_bold')
nback_pa = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-WM_acq-PA_bold')
gamb_ap = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-gambling_acq-AP_bold')
gamb_pa = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-gambling_acq-PA_bold')
#rest_ap = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_acq-AP_bold')
#rest_pa = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_acq-PA_bold')

#ASL
asl = create_key('sub-{subject}/{session}/perf/sub-{subject}_{session}_task-rest_asl')
asl_mz = create_key('sub-{subject}/{session}/perf/sub-{subject}_{session}_task-rest_m0scan')
asl_mp = create_key('sub-{subject}/{session}/perf/sub-{subject}_{session}_task-rest_deltam')
           
#Diffusion scans
dti_98dir_ap = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-98dir_dir-AP_dwi')
dti_98dir_pa = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-98dir_dir-PA_dwi')
dti_99dir_ap = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-99dir_dir-AP_dwi')
dti_99dir_pa = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-99dir_dir-PA_dwi')

#sbref images
wm_ap_sbref = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-WM_acq-AP_sbref')
wm_pa_sbref = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-WM_acq-PA_sbref')
gamb_ap_sbref = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-gambling_acq-AP_sbref')
gamb_pa_sbref = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-gambling_acq-PA_sbref')

dti_98dir_ap_sbref = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-98dir_dir-AP_sbref')
dti_98dir_pa_sbref = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-98dir_dir-PA_sbref')
dti_99dir_ap_sbref = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-99dir_dir-AP_sbref')
dti_99dir_pa_sbref = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-99dir_dir-PA_sbref')


from collections import defaultdict
def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where
        allowed template fields - follow python string module:
        item: index within category
        subject: participant id
        seqitem: run number during scanning
        subindex: sub index within group"""

    info = {
        t1w: [], t1wn: [], t2w: [], t2wn: [], fm_ap1: [], fm_pa1: [], nback_ap: [],                 
        nback_pa: [], gamb_ap: [], gamb_pa: [], wm_ap_sbref: [], wm_pa_sbref: [],   
        gamb_ap_sbref: [], gamb_pa_sbref: [], #rest_ap: [], #rest_pa: [], 
        dti_98dir_ap: [], dti_98dir_pa: [], dti_99dir_ap:         
        [], dti_99dir_pa: [], asl: [], asl_mz: [], asl_mp: [], 
        dti_98dir_ap_sbref: [], dti_98dir_pa_sbref:      
        [], dti_99dir_ap_sbref: [], dti_99dir_pa_sbref: [],
    }

    fmap_times_pa = {}
    fmap_times_ap = {}
    
    for s in seqinfo:
        protocol = s.protocol_name.lower()
        if "t1w" in protocol and 'NORM' in s.image_type:
            info[t1wn].append(s.series_id)
        elif "t1w" in protocol:
            info[t1w].append(s.series_id)
        elif "t2w" in protocol and 'NORM' in s.image_type:
            info[t2wn].append(s.series_id)
        elif "t2w" in protocol:
            info[t2w].append(s.series_id)
        elif "tfMRI_WM_AP_SBRef" in s.series_description:
            info[wm_ap_sbref].append(s.series_id)
            tfmri_time = s.date
        elif "tfMRI_WM_PA_SBRef" in s.series_description:
            info[wm_pa_sbref].append(s.series_id)
        elif "tfMRI_GAMBLING_AP_SBRef" in s.series_description:
            info[gamb_ap_sbref].append(s.series_id)
        elif "tfMRI_GAMBLING_PA_SBRef" in s.series_description:
            info[gamb_pa_sbref].append(s.series_id)
        elif "gambling_ap" in protocol:
            info[gamb_ap].append(s.series_id)
        elif "gambling_pa" in protocol:
            info[gamb_pa].append(s.series_id)
            gambling_time = s.date
        elif "wm_ap" in protocol:
            info[nback_ap].append(s.series_id)
        elif "wm_pa" in protocol:
            info[nback_pa].append(s.series_id)       
        elif "spinechofieldmap_pa" in protocol:
            fmap_times_pa[s.date] = s
        elif "spinechofieldmap_ap" in protocol:
            fmap_times_ap[s.date] = s
        elif s.series_description.endswith("_M0"):
            info[asl_mz].append(s.series_id)
        elif s.series_description.endswith("_MeanPerf"):
            info[asl_mp].append(s.series_id)
        elif s.series_description.endswith("_ASL"):
            info[asl].append(s.series_id)
        elif "dMRI_dir98_AP_SBRef" in s.series_description:
            info[dti_98dir_ap_sbref].append(s.series_id)
        elif "dMRI_dir98_PA_SBRef" in s.series_description:
            info[dti_98dir_pa_sbref].append(s.series_id)
        elif "dMRI_dir99_AP_SBRef" in s.series_description:
            info[dti_99dir_ap_sbref].append(s.series_id)
        elif "dMRI_dir99_PA_SBRef" in s.series_description:
            info[dti_99dir_pa_sbref].append(s.series_id)
        elif "dmri_dir98_ap" in protocol:
            info[dti_98dir_ap].append(s.series_id)
        elif "dmri_dir98_pa" in protocol:
            info[dti_98dir_pa].append(s.series_id)
        elif "dmri_dir99_ap" in protocol:
            info[dti_99dir_ap].append(s.series_id)
        elif "dmri_dir99_pa" in protocol:
            info[dti_99dir_pa].append(s.series_id)
        #elif "rfmri_rest_ap" in protocol:
            #info[rest_ap].append(s.series_id)
        #elif "rfmri_rest_pa" in protocol:
            #info[rest_pa].append(s.series_id)
        else:
            print("Series not recognized!: ", protocol, s.dcm_dir_name)
            
    start_time = datetime.datetime.strptime(gambling_time, '%Y-%m-%dT%H:%M:%S.%f')
    end_time = datetime.datetime.strptime(tfmri_time, '%Y-%m-%dT%H:%M:%S.%f')
    
    ap_fmaps = []
    for k, v in fmap_times_ap.items():
        fmap_time = datetime.datetime.strptime(k, '%Y-%m-%dT%H:%M:%S.%f')
        if fmap_time < end_time and fmap_time > start_time:
            info[fm_ap1]. append(v.series_id)
    
    pa_fmaps = []
    for k, v in fmap_times_pa.items():
        fmap_time = datetime.datetime.strptime(k, '%Y-%m-%dT%H:%M:%S.%f')
        if fmap_time < end_time and fmap_time > start_time:
            info[fm_pa1]. append(v.series_id)
    
    return info
    
    def ReplaceSession(sesname):
        return sesname[:10].replace("-", '')


