'''
    Heuristic file for BIDS classification of legacy data collected by the Penn 
    Frontotemporal Degeneration Center on the HUP 3T MRI scanner. The main
    function, infotodict, is defined below and takes sequence information from
    dicom headers; you can see which information is extracted by running
    fw-heudiconv-tabulate on the session directory, which writes the sequence
    info to a tsv file that can subsequently be read in as a Pandas dataframe.
    Each row of seqinfo corresponds to one series/acquisition in an imaging 
    session.
    
    Edits by JSP, 4/21/2021:
        * Stop curating some sequences, leave as non-BIDS:
            * vNav ND images
            * Moco series
        * Sort the series by timestamp first (or some other chronological sort key).
        * After one pass through all the series, check for keys that have >1 value.
            * For these cases, create 1 key for each value, distinguished by run-{index}
            * Add 1 value to each key.
            * Delete the original key.
        * Alternatively, create add_series function that:
            * Checks whether a key already exists in the info dict
            * If so, increment the run number.
            * If not, add the key to the dict and add the series.
'''

# Changes, 2/17/2021:
# Put derivatives into derivatives folder, e.g., for ADC & FA images created from DTI series.
# ISSS BOLD fMRI: see 101162 20160525-1400

import datetime
import numpy as np

def create_key(template, outtype=('nii.gz',), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return template, outtype, annotation_classes

# Localizers and scouts.
locz = create_key('sub-{subject}/{session}/localizer/sub-{subject}_{session}_localizer{item}')
t2w_locz = create_key('sub-{subject}/{session}/localizer/sub-{subject}_{session}_acq-T2w_localizer{item}')
vess_scout = create_key('sub-{subject}/{session}/localizer/sub-{subject}_{session}_acq-vessel_localizer{item}')

# T1-weighted anatomical images
t1w = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_T1w')
t1w_ax = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_T1w')
t1w_sag = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-sag_T1w')
t1w_3d = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-3D_T1w')
t1w_norm = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_rec-norm_T1w')
t1w_gradwarp = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-sag_rec-gradwarp_T1w')
t1w_body = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_coil-body_T1w')
t1w_grappa = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-grappa_T1w')
t1w_vnav_moco = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-vnavmoco_T1w')
t1w_vnav_pass = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-vnavpass_T1w')

# T2-weighted anatomical images.
t2w = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_T2w')
t2w_norm = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_rec-norm_T2w')
t2w_tse_gradwarp = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-tse_rec-gradwarp_T2w')
t2w_gradwarp = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_rec-gradwarp_T2w')
t2w_gradwarp_phase = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_part-phase_rec-gradwarp_T2w')
t2w_hippo = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-hippo_T2w')
t2w_space = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-space_T2w')
t2w_sag = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-sag_T2w')
t2w_vnav_pass = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-vnavpass_T2w')

# T2 star:
t2star = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_T2star')
t2star_acpc = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-acpc_T2star')

# FLAIR
flair_3d = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-3D_FLAIR')
flair_3d_gradwarp = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-3D_rec-gradwarp_FLAIR')
flair_ax = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-axial_FLAIR')
flair_ax_gradwarp = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-axial_rec-gradwarp_FLAIR')
flair_sag = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-sag_FLAIR')

# SWI (susceptibility-weighted imaging)
swi_mag = create_key('sub-{subject}/{session}/swi/sub-{subject}_{session}_part-mag_GRE')
swi_ph = create_key('sub-{subject}/{session}/swi/sub-{subject}_{session}_part-phase_GRE')
swi_combo = create_key('sub-{subject}/{session}/swi/sub-{subject}_{session}_swi')
swi_minip = create_key('sub-{subject}/{session}/swi/sub-{subject}_{session}_minIP')

# Proton density
pd = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_PD')

# Combined proton density/T2 scan
pdt2 = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_PDT2')

# Angiography
tof = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_angio')

# Field maps
fm_phasediff1 = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_phasediff')
fm_mag1 = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_magnitude{item}')
fm_phasediff2 = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_phasediff')
fm_mag2 = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_magnitude{item}')
fm_phasediff3 = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_phasediff')
fm_mag3 = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_magnitude{item}')
fm_phasediff4 = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_phasediff')
fm_mag4 = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_magnitude{item}')
fm_phasediff5 = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_phasediff')
fm_mag5 = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_magnitude{item}')
fm_phasediff6 = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_phasediff')
fm_mag6 = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_magnitude{item}')
fm_ap1 = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_dir-AP_epi')
fm_pa1 = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_dir-PA_epi')
fm_gre = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_fieldmap')

# BOLD fMRI
bold_mb = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_acq-MB_bold')
rest_bold = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_bold')
rest_ap = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_dir-AP_bold')
rest_pa = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_dir-PA_bold')
rest_ap_sbref = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_dir-AP_sbref')
rest_pa_sbref = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_dir-PA_sbref')
pace = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_acq-pace_bold')

# Diffusion-weighted (white matter) imaging
dti_12dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-12dir_dwi')
dti_30dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_dwi')
dti_30dir_ap = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_dir-AP_dwi')
dti_30dir_pa = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_dir-PA_dwi')
dti_32dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-32dir_dwi')
dti_34dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_dwi')
dti_52dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-52dir_dwi')
dti_55dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-55dir_dwi')
dti_62dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-62dir_dwi')
dwi_abcd = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-96dir_dwi')
dwi_117dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-117dir_dwi')
noddi_b2000 =  create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-NODDIB2000_dwi')
noddi_b700 =  create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-NODDIB700_dwi')
noddi_b300 =  create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-NODDIB300_dwi')

#ASL
asl = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_asl')
asl_mz = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_m0scan')
asl_mp = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_cbf')
pcasl_3d = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-3D_asl')
pcasl_3d_mz = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-3D_m0scan')
pcasl_3d_mp = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-3D_cbf')
# Should we use the acq entity like this? How specific should filenames strive to be?
pasl = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-pasl_asl')
pasl_mp = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-pasl_cbf')
casl = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-casl_asl')
casl_mz = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-casl_m0scan')
fairest = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-fairest_asl')
fairest_mz = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-fairest_m0scan')

from collections import defaultdict
def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where
        allowed template fields - follow python string module:
        index: index within category
        subject: participant id
        seqindex: run number during scanning
        subindex: sub index within group
    """

    info = {
        locz: [], t2w_locz: [], vess_scout: [], 
        t1w: [], t1w_ax: [], t1w_sag: [],
        t1w_3d: [], t1w_norm: [], t1w_gradwarp: [], t1w_body: [], t1w_grappa: [],
        t1w_vnav_moco: [], t1w_vnav_pass: [],
        t2w: [], t2w_norm: [], t2w_gradwarp: [], t2w_gradwarp_phase: [],
        t2w_tse_gradwarp: [], t2w_space: [], t2w_sag: [], t2w_vnav_pass: [],
        t2w_hippo: [], 
        t2star: [], t2star_acpc: [],
        flair_3d: [], flair_3d_gradwarp: [],
        flair_ax: [], flair_ax_gradwarp: [],
        flair_sag: [],
        swi_mag: [], swi_ph: [], swi_combo: [], swi_minip: [],
        pd: [], pdt2: [], tof: [],
        fm_phasediff1: [], fm_mag1: [], 
        fm_phasediff2: [], fm_mag2: [],
        fm_phasediff3: [], fm_mag3: [], 
        fm_phasediff4: [], fm_mag4: [], 
        fm_phasediff5: [], fm_mag5: [],
        fm_phasediff6: [], fm_mag6: [],
        fm_ap1: [], fm_pa1: [],
        fm_gre: [],
        bold_mb: [], rest_bold: [],
        rest_ap: [], rest_pa: [],
        rest_ap_sbref: [], rest_pa_sbref: [],
        pace: [],
        dti_12dir: [],
        dti_30dir: [], 
        dti_34dir: [],
        dti_52dir: [],
        dti_55dir: [],
        dti_62dir: [],
        dwi_117dir: [], 
        dti_30dir_ap: [], dti_30dir_pa: [],
        dwi_abcd: [],
        noddi_b2000: [], noddi_b700: [], noddi_b300: [],
        asl: [], asl_mz: [], asl_mp: [],
        pcasl_3d: [], pcasl_3d_mz: [], pcasl_3d_mp: [],
        pasl: [], pasl_mp: [], 
        casl: [], casl_mz: [],
        fairest: [], fairest_mz: []
    }
    
    for s in seqinfo:
        protocol = s.protocol_name.lower()
        mydatetime = datetime.datetime.strptime(s.date, '%Y-%m-%dT%H:%M:%S.%f')
        if "localizer" in protocol:
            info[locz].append(s.series_id)
        elif 'aahead_scout' in protocol:
            info[locz].append(s.series_id)        
        elif 'vessel_scout' in protocol:
            info[vess_scout].append(s.series_id)      
        elif "t1_3d" in protocol:
            info[t1w_3d].append(s.series_id)
        elif protocol == '3dt1' and not "ND" in s.series_description:
            info[t1w_3d].append(s.series_id)
        elif protocol == 'mprage' and s.series_description == "MPRAGE":
            info[t1w_3d].append(s.series_id)
        elif "t1w" in protocol and 'NORM' in s.image_type and 'vnav' not in protocol:
            info[t1w_norm].append(s.series_id)
        elif "t1_mprage_iso" in protocol and "DIS3D" in s.image_type:
            info[t1w_gradwarp].append(s.series_id)
        elif "sag_t1_mprage" in protocol and "DIS3D" in s.image_type:
            info[t1w_gradwarp].append(s.series_id)
        elif "sag_t1_mprage" in protocol:
            info[t1w_sag].append(s.series_id)
        elif "t1w" in protocol and 'vnav' not in protocol:
            info[t1w].append(s.series_id)
        elif "t1_mpr_ax_mprage" in protocol and "PRIMARY" in s.image_type:
            info[t1w_ax].append(s.series_id)
        elif "mprage_bodycoil" in protocol and "PRIMARY" in s.image_type:
            info[t1w_body].append(s.series_id)
        elif "mprage_grappa" in protocol and "PRIMARY" in s.image_type:
            info[t1w_grappa].append(s.series_id)
        elif "sagittal_mp-rage" in protocol and "PRIMARY" in s.image_type:
            info[t1w_sag].append(s.series_id)
        elif "accelerated_sagittal_mprage" in protocol and "PRIMARY" in s.image_type:
            info[t1w_sag].append(s.series_id)
        elif "t1_se_sag" in protocol and "PRIMARY" in s.image_type:
            info[t1w_sag].append(s.series_id)
        elif "mprage_sag_iso" in protocol and "PRIMARY" in s.image_type:
            info[t1w_sag].append(s.series_id)
        elif "T1w_MPR_vNav_moco_RMS" in s.series_description:
            info[t1w_vnav_moco].append(s.series_id)
        elif "MPRAGE_sag_moco3" in s.series_description and "Moco3d1_ns" in s.sequence_name:
            info[t1w_vnav_moco].append(s.series_id)
        elif "T1w_MPR_vNav_passive_RMS" in s.series_description:
            info[t1w_vnav_pass].append(s.series_id)
        elif "T1w_MPR_vNav_passive_(sectra)_RMS" in s.series_description:
            info[t1w_vnav_pass].append(s.series_id)
        elif "t2_tse_cor8channel" in protocol:
            info[t2w_hippo].append(s.series_id)
        elif "t2_tse_obl_448_2mm" in protocol:
            info[t2w_hippo].append(s.series_id)
        elif "highreshippocampus" in protocol:
            info[t2w_hippo].append(s.series_id)
        elif "t2_tse_coronal_hippo" in protocol:
            info[t2w_hippo].append(s.series_id)
        elif "t2_cor_high_res" in protocol:
            info[t2w_hippo].append(s.series_id)
        elif "t2_tse_tra8channel" in protocol:
            info[t2w_locz].append(s.series_id)
        elif "t2_2d_0.4x0.4x1.2mm_180flip" in protocol and not "ND" in s.series_description:
            info[t2w_hippo].append(s.series_id)
        elif "t2w" in protocol and 'NORM' in s.image_type and 'vnav' not in protocol:
            info[t2w_norm].append(s.series_id)
        elif "t2w" in protocol and 'vnav' not in protocol:
            info[t2w].append(s.series_id)
        elif "t2_axial_space" in protocol:
            info[t2w_space].append(s.series_id)
        elif "t2_sag" in protocol:
            info[t2w_sag].append(s.series_id)
        elif "T2w_SPC_vNav_passive" in s.series_description and not "_ND" in s.series_description:
            info[t2w_vnav_pass].append(s.series_id)
        elif "t2_tse_axial" in protocol and "DIS2D" in s.image_type:
            info[t2w_tse_gradwarp].append(s.series_id)
        elif "axial_t2_acpc_angle" in protocol and "DIS2D" in s.image_type and "M" in s.image_type:
            info[t2w_gradwarp].append(s.series_id)
        elif "t2_fl2d_tra_hemo" in protocol and "DIS2D" in s.image_type and "M" in s.image_type:
            info[t2w_gradwarp].append(s.series_id)
        elif "t2_fl2d_tra_hemo" in protocol and "DIS2D" in s.image_type and "P" in s.image_type:
            info[t2w_gradwarp_phase].append(s.series_id)
        elif "axt2_star_acpc_angle" in protocol:
            info[t2star_acpc].append(s.series_id)
        elif "t2_star" in protocol:
            info[t2star].append(s.series_id)
        elif protocol == "t2":
            info[t2star].append(s.series_id)
        elif "flair_3d" in protocol and 'DIS3D' in s.image_type:
            info[flair_3d_gradwarp].append(s.series_id)
        elif "flair_3d" in protocol:
            info[flair_3d].append(s.series_id)
        elif "flair_axial" in protocol and ("DIS2D" in s.image_type or "DIS3D" in s.image_type):
            info[flair_ax_gradwarp].append(s.series_id)
        elif "flair_axial" in protocol:
            info[flair_ax].append(s.series_id)
        elif protocol == 'flair' and ("DIS2D" in s.image_type or "DIS3D" in s.image_type):
            info[flair_ax_gradwarp].append(s.series_id)
        elif protocol == 'flair':
            info[flair_ax].append(s.series_id)
        elif "flair_acpc_angle" in protocol and ("DIS2D" in s.image_type or "DIS3D" in s.image_type):
            info[flair_ax_gradwarp].append(s.series_id)
        elif "flair_acpc_angle" in protocol:
            info[flair_ax].append(s.series_id)
        elif "ax_flair" in protocol:
            info[flair_ax].append(s.series_id)
        elif "axial_flair" in protocol:
            info[flair_ax].append(s.series_id)
        elif "flair_sag" in protocol:
            info[flair_sag].append(s.series_id)
        elif "sagittal_3d_flair" in protocol:
            info[flair_sag].append(s.series_id)
        elif "swi" in protocol and s.series_description == 'Mag_Images':
            info[swi_mag].append(s.series_id)
        elif "swi" in protocol and s.series_description == 'Pha_Images':
            info[swi_ph].append(s.series_id)
        elif "swi" in protocol and s.series_description == 'SWI_Images':
            info[swi_combo].append(s.series_id)
        elif "swi" in protocol and s.series_description == 'mIP_Images(SW)':
            info[swi_minip].append(s.series_id)
        elif "pd_tse_tra" in protocol:
            info[pd].append(s.series_id)
        elif "axial_pd-t2_tse" in protocol:
            info[pdt2].append(s.series_id)
        elif protocol == "pdt2":
            info[pdt2].append(s.series_id)
        elif "tof" in protocol:
            info[tof].append(s.series_id)
        elif "ep2d_se_pcasl_m0" in protocol:
            info[asl_mz].append(s.series_id)
        elif "axial_3d_pasl_(eyes_open)" in protocol and "Perfusion_Weighted" in s.series_description:
            info[pasl_mp].append(s.series_id)
        elif "3d_pasl" in protocol and "Perfusion_Weighted" in s.series_description:
            info[pasl_mp].append(s.series_id)
        elif "3d_pasl" in protocol and "3D_PASL" in s.series_description:
            info[pasl].append(s.series_id)
        elif "3d_pcasl" in protocol and "Perfusion_Weighted" in s.series_description:
            info[pcasl_3d_mp].append(s.series_id)
        elif "3d_pcasl" in protocol and "3D_PCASL" in s.series_description:
            info[pcasl_3d].append(s.series_id)
        elif "ep2d_se_pcasl_phc_1500ms" in protocol and not "moco" in protocol:
            info[asl].append(s.series_id)
        elif "3dsp_pcasl_singleshot_2arms" in protocol:
            info[asl].append(s.series_id)
        elif "pcasl_resting_90mm" in protocol:
            info[asl].append(s.series_id)
        elif "ep2d_casl_ui_1500ms" in protocol and not "moco" in protocol:
            info[casl].append(s.series_id)
        elif "ep2d_casl_am_ui_1500ms" in protocol and not s.series_description == 'MoCoSeries':
            info[casl].append(s.series_id)
        elif "ep2d_casl_1500ms" in protocol and not s.series_description == 'MoCoSeries':
            info[casl].append(s.series_id)
        elif "ep2d_casl_1500ms" in protocol and not "MOSAIC" in s.image_type:
            info[casl_mz].append(s.series_id)
        elif "spiral_v20_hcp" in protocol and "M0" in s.series_description:
            info[asl_mz].append(s.series_id)
        elif "spiral_v20_hcp" in protocol and "MeanPerf" in s.series_description:
            info[asl_mp].append(s.series_id)
        elif "spiral_v20_hcp" in protocol:
            info[asl].append(s.series_id)
        elif "asl_3dspiral_4shot_2.5mm_1daccel_v20" in protocol and "M0" in s.series_description:
            info[asl_mz].append(s.series_id)
        elif "asl_3dspiral_4shot_2.5mm_1daccel_v20" in protocol and "MeanPerf" in s.series_description:
            info[asl_mp].append(s.series_id)
        elif "asl_3dspiral_4shot_2.5mm_1daccel_v20" in protocol:
            info[asl].append(s.series_id)
        elif "ep2d_fairest_ui_1500ms" in protocol and not "moco" in protocol:
            info[fairest].append(s.series_id)
        elif "ep2d_fairest_ui_m0" in protocol and not "moco" in protocol:
            info[fairest_mz].append(s.series_id)
        elif 'ep2d_diff_mddw_12' in protocol:
            info[dti_12dir].append(s.series_id)
        elif 'dti_12dir' in protocol:
            info[dti_12dir].append(s.series_id)
        elif 'ep2d_dti_30dir_t' in protocol:
            info[dti_62dir].append(s.series_id)
        elif protocol == 'dti_30dir_nodico_vox2_1000':
            info[dti_34dir].append(s.series_id)
        elif 'dti_30dir' in protocol and not "moco" in protocol:
            info[dti_30dir].append(s.series_id)
        elif s.series_description == 'DTI_A-P':
            info[dti_30dir_ap].append(s.series_id)
        elif s.series_description == 'DTI_P-A':
            info[dti_30dir_pa].append(s.series_id)
        elif s.series_description == 'DTI_A-P_BW2394':
            info[dti_30dir_ap].append(s.series_id)
        elif s.series_description == 'DTI_P-A_BW2394':
            info[dti_30dir_pa].append(s.series_id)
        elif 'dti_34dir' in protocol:
            info[dti_34dir].append(s.series_id)
        elif 'dti_34_dir' in protocol:
            info[dti_34dir].append(s.series_id)
        elif protocol == 'axial_dti' and 'ORIGINAL' in s.series_description:
            info[dti_52dir].append(s.series_id)
        elif 'dti_2x32' in protocol:
            info[dti_32dir].append(s.series_id)
        elif 'long_dti_p-a' in protocol:
            info[dti_55dir].append(s.series_id)
        elif 'axial_dti_a>>p' in protocol or 'axial_dti_p>>a' in protocol:
            info[dti_55dir].append(s.series_id)
        elif protocol == 'abcd_dmri':
            info[dwi_abcd].append(s.series_id)
        elif protocol == 'multishell_117dir':
            info[dwi_117dir].append(s.series_id)
        elif 'noddi_b_2000' in protocol:
            info[noddi_b2000].append(s.series_id)
        elif 'noddi_b_700' in protocol:
            info[noddi_b700].append(s.series_id)
        elif 'noddi_b_300' in protocol:
            info[noddi_b300].append(s.series_id)
        elif "spinechofieldmap_ap" in protocol:
            info[fm_ap1].append(s.series_id)
        elif "spinechofieldmap_pa" in protocol:
            info[fm_pa1].append(s.series_id)
        elif "b0map_onesizefitsall_v4" in protocol and "P" in s.image_type:
            info[fm_phasediff1].append(s.series_id)
        elif "b0map_onesizefitsall_v4" in protocol and "M" in s.image_type:
            info[fm_mag1].append(s.series_id)
        elif "b0map_siemens" in protocol and "P" in s.image_type:
            info[fm_phasediff2].append(s.series_id)
        elif "b0map_siemens" in protocol and "M" in s.image_type:
            info[fm_mag2].append(s.series_id)
        elif "b0map_v4" in protocol and "P" in s.image_type:
            info[fm_phasediff2].append(s.series_id)
        elif "b0map_v4" in protocol and "M" in s.image_type:
            info[fm_mag2].append(s.series_id)
        elif "b0_map" in protocol and "P" in s.image_type:
            info[fm_phasediff3].append(s.series_id)
        elif "b0_map" in protocol and "M" in s.image_type:
            info[fm_mag3].append(s.series_id)
        elif "fieldmap_gre" in protocol:
            info[fm_gre].append(s.series_id)
        elif "gre_field_mapping" in protocol:
            info[fm_gre].append(s.series_id)
        elif "field_mapping" in protocol and "P" in s.image_type:
            info[fm_phasediff4].append(s.series_id)
        elif "field_mapping" in protocol and "M" in s.image_type:
            info[fm_mag4].append(s.series_id)
        elif "b0map" in protocol and "P" in s.image_type:
            info[fm_phasediff5].append(s.series_id)
        elif "b0map" in protocol and "M" in s.image_type:
            info[fm_mag5].append(s.series_id)
        elif protocol == 'field_map' and "P" in s.image_type:
            info[fm_phasediff6].append(s.series_id)
        elif protocol == 'field_map' and "M" in s.image_type:
            info[fm_mag6].append(s.series_id)
        elif "axial mb rsfmri (eyes open)" in protocol:
            info[bold_mb].append(s.series_id)
        elif "axial_mb_rsfmri_(eyes_open)" in protocol:
            info[bold_mb].append(s.series_id)
        elif "resting_bold_124" in protocol:
            info[rest_bold].append(s.series_id)
        elif "long_rsfmri_p-a" in protocol and 'sbref' in s.series_description.lower():
            info[rest_pa_sbref].append(s.series_id)
        elif "long_rsfmri_p-a" in protocol:
            info[rest_pa].append(s.series_id)
        elif "long_rsfmri_a-p" in protocol and 'sbref' in s.series_description.lower():
            info[rest_ap_sbref].append(s.series_id)
        elif "long_rsfmri_a-p" in protocol:
            info[rest_ap].append(s.series_id)
        elif "rfmri_rest_ap" in protocol and 'sbref' in s.series_description.lower():
            info[rest_ap_sbref].append(s.series_id)
        elif "rfmri_rest_pa" in protocol and 'sbref' in s.series_description.lower():
            info[rest_pa_sbref].append(s.series_id)
        elif "rfmri_rest_ap" in protocol:
            info[rest_ap].append(s.series_id)
        elif "rfmri_rest_pa" in protocol:
            info[rest_pa].append(s.series_id)
        elif "ax_rsfmri_a>>p" in protocol:
            info[rest_ap].append(s.series_id)
        elif "ax_rsfmri_p>>a" in protocol:
            info[rest_pa].append(s.series_id)
        elif "ax rsfmri_a>>p" in protocol:
            info[rest_ap].append(s.series_id)
        elif "ax rsfmri_p>>a" in protocol:
            info[rest_pa].append(s.series_id)
        elif "ax_rsfmri" in protocol:
            info[rest_pa].append(s.series_id)
        elif "axial_rsfmri_(eyes_open)" in protocol:
            info[rest_pa].append(s.series_id)
        elif "bold_resting_2x2x2" in protocol:
            info[rest_bold].append(s.series_id)
        elif "ep2d_bold_restingze6min" in protocol:
            info[rest_bold].append(s.series_id)
        elif protocol == "fmri":
            info[rest_bold].append(s.series_id)
        elif "ep2d_max_pace" in protocol and not "MoCo" in s.series_description:
            info[pace].append(s.series_id)
        elif "ep2d_pace_max" in protocol and not "MoCo" in s.series_description:
            info[pace].append(s.series_id)
        else:
            print("Series not recognized!: ", protocol, s.dcm_dir_name)

    # Get timestamp info to use as a sort key.
    def get_date(series_info):
        return(datetime.datetime.strptime(series_info.date, '%Y-%m-%dT%H:%M:%S.%f'))


    # Before returning the info dictionary, 1) get rid of empty dict entries; and
    # 2) for entries that have more than one series, differentiate them by run-{index}.
    def update_key(series_key, runindex):
        series_name = series_key[0]
        s = series_name.split("_")
        nfields = len(s)
        s.insert(nfields-1, "run" + str(runindex))
        new_name = "_".join(s)
        return((new_name, series_key[1], series_key[2]))

    newdict = {}
    delkeys = []
    for k in info.keys():
        ids = info[k]
        unique_ids = list(set(ids))
        if len(unique_ids) > 1:
            series_list = [s for s in seqinfo if (s.series_id in unique_ids)]
            unique_tuples = list(set([(s.series_uid, get_date(s)) for s in series_list]))
            # Sort the series UIDs by time of acquisition.
            def sortfunc(mytuple):
                return(mytuple[1])
            unique_tuples.sort(key = sortfunc)
            uids = []
            for val in unique_tuples:
                if not val[0] in uids:
                    uids.append(val[0])
#            print("unique_tuples: ",unique_tuples)
            nseries = len(uids)
            newkeys = [update_key(k, i) for i in range(1, nseries + 1)]
#            print("newkeys: ", newkeys)
#            print("series_list: ", [(s.series_id, get_date(s)) for s in series_list])
            delkeys.append(k)
            for i in range(nseries):
                series_matches = [s for s in series_list if s.series_uid == uids[i]]
                newdict[newkeys[i]] = []
                for match in series_matches:
                    newdict[newkeys[i]].append(match.series_id)
    # Merge the two dictionaries.
    info.update(newdict)

    # Delete keys that were expanded on in the new dictionary.
    for k in delkeys:
        info.pop(k, None)

#    for k,v in info.items():
#        if len(info[k]) > 0:
#            for vals in v:
#                print(k,vals)

    return info

def ReplaceSession(sesname):
    return sesname[:13].replace("-", "x").replace(".","x").replace("_","x")

def ReplaceSubject(subjname):
    return subjname[:10].replace("-", "x").replace(".","x").replace("_","x")

MetadataExtras = {

   fm_phasediff1: {
       "EchoTime1": 0.00255,
       "EchoTime2": 0.00501
   },

   fm_phasediff2: {
       "EchoTime1": 0.00406,
       "EchoTime2": 0.00652
   },

   fm_phasediff3: {
       "EchoTime1": 0.00411,
       "EchoTime2": 0.00657
   },

   fm_phasediff4: {
       "EchoTime1": 0.00492,
       "EchoTime2": 0.00738
   },

   fm_phasediff5: {
       "EchoTime1": 0.00412,
       "EchoTime2": 0.00658
   },

   fm_phasediff6: {
       "EchoTime1": 0.00519,
       "EchoTime2": 0.00765
   }

}       

IntendedFor = {

    # B0 fieldmap
    fm_phasediff1: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],
    fm_mag1: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],

    fm_phasediff2: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],
    fm_mag2: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],

    fm_phasediff3: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],
    fm_mag3: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],

    fm_phasediff4: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz',
        '{session}/func/sub-{subject}_{session}_task-rest_acq-MB_bold.nii.gz'],
    fm_mag4: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz',
        '{session}/func/sub-{subject}_{session}_task-rest_acq-MB_bold.nii.gz'],

    fm_phasediff5: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],
    fm_mag5: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],

    fm_phasediff6: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],
    fm_mag6: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],

    fm_ap1: [ '{session}/func/{subject}_{session}_task-rest_dir-AP_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-AP_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-AP_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-AP_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-AP_run-1_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-AP_run-2_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_run-1_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_run-2_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz',
        '{session}/dwi/{subject}_{session}_acq-96dir_dwi.nii.gz' ],
        
    fm_pa1: [ '{session}/func/{subject}_{session}_task-rest_dir-AP_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-AP_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-AP_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-AP_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-AP_run-1_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-AP_run-2_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_run-1_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_run-2_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz',
        '{session}/dwi/{subject}_{session}_acq-96dir_dwi.nii.gz' ],

    fm_gre: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_run-4_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_run-5_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_run-6_bold.nii.gz']

}

# Non-BIDS graveyard: a bunch of keys and conditions that were used in prior versions of
# this heuristic but that we've walked back.
#
# EPI distortion map for correcting eddy currents along phase-encoding
# direction in ABCD diffusion protocol. Don't think a P->A map was ever
# collected.
# Edit, JSP, 12/07/2020: the BIDS validator doesn't like this distortion map, so make
# it a non-BIDS file.
# dwi_distmap_ap: [ '{session}/dwi/{subject}_{session}_acq-96dir_dwi.nii.gz' ],
# dwi_distmap_pa: [ '{session}/dwi/{subject}_{session}_acq-96dir_dwi.nii.gz' ],
#t1w_3d_nd = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-3D_rec-ND_T1w')
#t1w_vnav_moco_nd = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-vnavmoco_rec-ND_T1w')
#t1w_vnav_pass_nd = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-vnavpass_rec-ND_T1w')
# Kinda making this up, but t2w_gap will include the t2_nex_1 and t2_tse_tra
# sequences--both axial T2 scans with only partial slice coverage due to large 
# gaps (4+ mm) between slices.
#t2w_gap = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-gap_T2w')
# T2 SPC vNavs
#t2w_vnav_pass_nd = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-vnavpass_rec-ND_T2w')
#t2w_hippo_nd = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-hippo_rec-ND_T2w')
#dti_30dir_run1 = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_run-1_dwi')
#dti_30dir_run2 = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_run-2_dwi')
#dti_30dir_run3 = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_run-3_dwi')
#dti_32dir_run1 = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-32dir_run-1_dwi')
#dti_32dir_run2 = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-32dir_run-2_dwi')
#dti_32dir_run3 = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-32dir_run-3_dwi')
#dti_34dir_run1 = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_run-1_dwi')
#dti_34dir_run2 = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_run-2_dwi')
#dti_34dir_run3 = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_run-3_dwi')
#dti_34dir_run1_moco = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_rec-moco_run-1_dwi')
#dti_34dir_run2_moco = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_rec-moco_run-2_dwi')
#dti_34dir_run3_moco = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_rec-moco_run-3_dwi')
#dti_34dir_moco = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_rec-moco_dwi')
#pace_moco = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_acq-pace_rec-moco_bold')
#asl_moco = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_rec-moco_asl')
#casl_moco = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-casl_rec-moco_asl')
#fairest_moco = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-fairest_rec-moco_asl')
#fairest_moco_mz = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-fairest_rec-moco_m0scan')
# Non-BIDS
# 'SWI_axial'
# 'ViSTa_3D_segEPI_2.5x2.5x2.5'
# 'GRE_3D_segEPI_2.5x2.5x2.5'
# 'ABCD_T1w_MPR_vNav_passive_(SECTRA)'
# 'T2w_SPC_vNav_passive'
#'ABCD_T1w_MPR_vNav_setter'
# 'MDDW_DTI_AXIAL_30_DIRECTION'
# ACR_SAG_T1: 2d image, some kind of localizer?
# ACR_AX_T1: 2d image, some kind of localizer?
# MPRAGE_NAVprotocol: these seem to just be the navigator images, not the T1s
# DWI: only 3 directional images plus one B0, not usable for white matter imaging
#        elif "t2_nex_1" in protocol:
#            info[t2w_gap].append(s.series_id)
#        elif "t2_tse_tra" in protocol:
#            info[t2w_gap].append(s.series_id)
#        elif s.series_description == '3DT1_ND':
#            info[t1w_3d_nd].append(s.series_id)
