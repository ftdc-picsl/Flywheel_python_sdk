'''
    Heuristic file for BIDS classification of legacy data collected by the Penn 
    Frontotemporal Degeneration Center on the HUP 3T MRI scanner. The main
    function, infotodict, is defined below and takes sequence information from
    dicom headers; you can see which information is extracted by running
    fw-heudiconv-tabulate on the session directory, which writes the sequence
    info to a tsv file that can subsequently be read in as a Pandas dataframe.
    Each row of seqinfo corresponds to one series/acquisition in an imaging 
    session.
'''

import datetime
import numpy as np

def create_key(template, outtype=('nii.gz',), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return template, outtype, annotation_classes

# Non-BIDS
# 'SWI_axial'
# 'ViSTa_3D_segEPI_2.5x2.5x2.5'
# 'GRE_3D_segEPI_2.5x2.5x2.5'
# 'ABCD_T1w_MPR_vNav_passive_(SECTRA)'
# 'T2w_SPC_vNav_passive'
# 'T2w_SPC_vNav_passive_(SECTRA)'
#'ABCD_T1w_MPR_vNav_setter'
#'ABCD_T1w_MPR_vNav_moco'
#'ABCD_T1w_MPR_vNav_passive'
# 'MDDW_DTI_AXIAL_30_DIRECTION'

# Localizers
# 'localizer','Localizer'
# 'AAHead_Scout'
# vessel_scout Vessel_scout_12.nii.gz
locz = create_key('sub-{subject}/{session}/localizer/sub-{subject}_{session}_run{item}_localizer')
t2w_locz = create_key('sub-{subject}/{session}/localizer/sub-{subject}_{session}_acq-T2w_localizer')
vess_scout = create_key('sub-{subject}/{session}/localizer/sub-{subject}_{session}_acq-vessel_localizer')

# anatomical images
#'t1_mpr_AX_MPRAGE'
# t1_3d_0.8x0.8x0.8 T1_3D_0.8x0.8x0.8_6.nii.gz
# 'MPRAGE_BodyCoil_ref'
# 'MPRAGE_GRAPPA2'
# 'MPRAGE_SAG_ISO'
# 'Sagittal_MP-Rage'
# 't1_se_sag'
t1w = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_T1w')
t1w_ax = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_T1w')
t1w_sag = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-sag_T1w')
t1w_3d = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-3D_T1w')
t1w_norm = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_rec-norm_T1w')
t1w_body = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-bodycoil_T1w')
t1w_grappa = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-grappa_T1w')
# T1w vNav
t1w_vnav_moco_nd = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-vnavmoco_rec-ND_T1w')
t1w_vnav_pass_nd = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-vnavpass_rec-ND_T1w')
t1w_vnav_moco = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-vnavmoco_T1w')
t1w_vnav_pass = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-vnavpass_T1w')

# T2-weighted anatomical images.
#'t2_tse_obl_448_2mm'
#'HighResHippocampus'
#'t2_tse_CORONAL_HIPPO'
#'t2_tse_tra8CHANNEL'
#
# 'T2_2D_0.4x0.4x1.2mm_180flip_aa_temp': t2w_hippo
# 'T2_2D_0.4x0.4x1.2mm_180flip_aa_temp_(SECTRA)': t2w_hippo
# 'T2_2D_0.4x0.4x1.2mm_180flip_aa_temp_repeat': t2w_hippo
# 'T2_COR_HIGH_RES': t2w_hippo
# 'T2_Axial_space': t2w_space
# 'T2_sag': t2w_sag
# 't2_nex_1': t2w_gap
# 't2_tse_tra': t2w_gap
t2w = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_T2w')
t2w_norm = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_rec-norm_T2w')
t2w_hippo = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-hippo_T2w')
t2w_space = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-space_T2w')
t2w_sag = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-sag_T2w')
# Kinda making this up, but t2w_gap will include the t2_nex_1 and t2_tse_tra
# sequences--both axial T2 scans with only partial slice coverage due to large 
# gaps (4+ mm) between slices.
t2w_gap = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-gap_T2w')
# T2 SPC vNavs
t2w_vnav_pass_nd = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-vnavpass_rec-ND_T2w')
t2w_vnav_pass = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-vnavpass_T2w')

# FLAIR
# 'AX_FLAIR_3mm'
# 'Axial_FLAIR_3mm'
# 'FLAIR_sag'
# 'T2_AXIAL_FLAIR'
flair_ax = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-axial_FLAIR')
flair_sag = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-sag_FLAIR')

# T2-star
t2star_ax = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-axial_T2star')

# Proton density
'pd_tse_tra'
pd = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_PD')

# Combined proton density/T2 scan
'Axial_PD-T2_TSE'
pdt2 = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_PDT2')

# Angiography
# fl_tof_neck fl_tof_neck_16.nii.gz
# tof_3d_multi-slab TOF_3D_multi-slab_25.nii.gz
# tof_3d_multi-slab TOF_3D_multi-slab_24.nii.gz
tof = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_angio')

# Field maps
#'b0map_v4'
#'B0map_onesizefitsall_v4'
#'SpinEchoFieldMap_AP'
# 'B0_map'
# 'fieldmap_GRE'
# 'SpinEchoFieldMap_PA'
# 'ABCD_dMRI_DistortionMap_AP'
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
fm_ap1 = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_dir-AP_epi')
fm_pa1 = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_dir-PA_epi')
fm_gre = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_fieldmap')

# BOLD fMRI
# 'Resting_bold_124'
#
# 'Ax_rsfMRI'
# 'Ax_rsfMRI_A>>P'
# 'Ax_rsfMRI_P>>A'
# 'BOLD_resting_2x2x2'
# 'ep2d_bold_restingZE6min'
# 'long_rsfMRI_P-A'
# 'rfMRI_REST_PA'
# 'ep2d_max_pace'
# 'ep2d_max_pace_151'
# 'ep2d_max_pace_151_MoCo'
# 'ep2d_max_pace_165'
# 'ep2d_max_pace_165_MoCo'
# 'ep2d_max_pace_205'
# 'ep2d_max_pace_205_MoCo'
# 'ep2d_max_pace_35'
# 'ep2d_max_pace_35_MoCo'
# 'ep2d_max_pace_45'
# 'ep2d_max_pace_45_MoCo'
# 'ep2d_max_pace_MoCo'
rest_bold = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_bold')
rest_ap = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_dir-AP_bold')
rest_pa = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_dir-PA_bold')
rest_ap_sbref = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_dir-AP_sbref')
rest_pa_sbref = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_dir-PA_sbref')
pace = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_acq-pace_run-{item}_bold')
pace_moco = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_acq-pace_rec-moco_run-{item}_bold')

# Diffusion-weighted (white matter) imaging
# DTI_30dir
# DTI_43dir
# 'Axial_DTI': dti_52dir
    # series_description
    # 'Axial_DTI_ADC'
    # 'Axial_DTI_FA'
    # 'Axial_DTI_ColFA'
    # 'Axial_DTI_TRACEW'
# 'Axial_DTI_A>>P': dti_34dir_trace, dti_34dir_fa
# 'DTI_2x32_35': dti_32dir
# 'DTI_2x32_36'
# 'DTI_30dir_noDiCo_vox2_1000': dti_34dir or dti_34dir_multirun
# 'DTI_30dir_noDiCo_vox2_1000_MoCo'
# 'ep2d_diff_MDDW_12': dti_12dir
# 'ep2d_DTI_30dir_T': dti_62dir
# 'long_DTI_P-A': dti_55dir
# 'no_AA_long_DTI_P-A':  dti_55dir
# 'NODDI_B_2000'
# 'NODDI_B_300'
# 'NODDI_B_700'
# MultiShell_117dir
dti_12dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-12dir_dwi')
dti_30dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_dwi')
dti_30dir_multirun = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_run-{item}_dwi')
dti_30dir_multirun_bval = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_run-{item}_dwi')
dti_30dir_multirun_bvec = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_run-{item}_dwi')
dti_32dir_multirun = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-32dir_run-{item}_dwi')
dti_32dir_multirun_bval = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-32dir_run-{item}_dwi')
dti_32dir_multirun_bvec = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-32dir_run-{item}_dwi')
dti_34dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_dwi')
dti_34dir_multirun = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_run-{item}_dwi')
dti_34dir_multirun_bval = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_run-{item}_dwi')
dti_34dir_multirun_bvec = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_run-{item}_dwi')
dti_34dir_moco = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_rec-moco_dwi')
dti_34dir_multirun_moco = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_rec-moco_run-{item}_dwi')
dti_34dir_multirun_moco_bval = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_rec-moco_run-{item}_dwi')
dti_34dir_multirun_moco_bvec = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_rec-moco_run-{item}_dwi')
dti_52dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-52dir_dwi')
dti_55dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-55dir_dwi')
dti_62dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-62dir_dwi')
dti_trace = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_desc-trace_dwi')
dwi_abcd = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-96dir_dwi')
dwi_117dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-117dir_dwi')
dwi_distmap_ap = create_key(
    'sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-dwi_dir-AP_epi')
dwi_distmap_pa = create_key(
    'sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-dwi_dir-PA_epi')
dti_adc = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_desc-adc_dwi')
dti_fa = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_desc-fa_dwi')
dti_exp = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_desc-exp_dwi')
noddi_b2000 =  create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-NODDIB2000_dwi')
noddi_b700 =  create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-NODDIB700_dwi')
noddi_b300 =  create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-NODDIB300_dwi')

#ASL
# Input based on ASL draft certification?
#'3DSP_pcasl_singleshot_2arms'
#'Axial_3D_PASL_(Eyes_Open)'
#'ep2d_se_pcasl_PHC_1500ms_MoCo'
#'ep2d_se_pcasl_PHC_1500ms'
#'3D_PASL'
#'ep2d_se_pcasl_M0'
#
# 'ep2d_casl_UI_1500ms'
# 'ep2d_casl_UI_1500ms_MoCo'
# 'PCASL_resting_90mm'
# 'SPIRAL_V20_HCP'
# 'ep2d_fairest_UI_1500ms'
# 'ep2d_fairest_UI_1500ms_MoCo'
# 'ep2d_fairest_UI_M0'
# 'ep2d_fairest_UI_M0_MoCo'
# asl_3dspiral_4shot_2.5mm_1daccel_v20 ASL_3DSpiral_4shot_2.5mm_1Daccel_V20_30.nii.gz
# asl_3dspiral_4shot_2.5mm_1daccel_v20 ASL_3DSpiral_4shot_2.5mm_1Daccel_V20_9.nii.gz
asl = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_asl')
asl_mz = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_m0scan')
asl_mp = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_cbf')
asl_moco = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_rec-moco_asl')
# Should we use the acq entity like this? How specific should filenames strive to be?
pasl = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-pasl_asl')
pasl_mp = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-pasl_cbf')
casl = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-casl_asl')
casl_moco = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-casl_rec-moco_asl')
fairest = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-fairest_asl')
fairest_moco = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-fairest_rec-moco_asl')
fairest_mz = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-fairest_m0scan')
fairest_moco_mz = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-fairest_rec-moco_m0scan')

from collections import defaultdict
def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where
        allowed template fields - follow python string module:
        item: index within category
        subject: participant id
        seqitem: run number during scanning
        subindex: sub index within group
    """

    info = {
        locz: [], t2w_locz: [], vess_locz: [], t1w: [], t1w_ax: [], t1w_sag: [],
        t1w_3d: [], t1w_norm: [],
        t1w_body: [], t1w_grappa: [], t1w_vnav_moco_nd: [], t1w_vnav_pass_nd: [], 
        t1w_vnav_moco: [], t1w_vnav_pass: [], t2w: [], t2w_norm: [],
        t2w_hippo: [], t2w_space: [], t2w_sag: [], t2w_gap: [], 
        t2w_vnav_pass_nd: [], t2w_vnav_pass: [], flair_ax: [], flair_sag: [], 
        t2star_ax: [], pd: [], pdt2: [], tof: [],
        fm_phasediff1: [], fm_mag1: [], fm_phasediff2: [], fm_mag2: [],
        fm_phasediff3: [], fm_mag3: [], fm_phasediff4: [], fm_mag4: [], 
        fm_phasediff5: [], fm_mag5: [],
        fm_ap1: [], fm_pa1: [], fm_gre: [],
        dwi_distmap_ap: [], dwi_distmap_pa: [],
        rest_bold: [], rest_ap: [], rest_pa: [], rest_ap_sbref: [],
        rest_pa_sbref: [], pace: [], pace_moco: [],
        dti_12dir: [], dti_30dir: [], dti_30dir_multirun: [], 
        dti_32dir_multirun: [], dti_34dir: [], dti_34dir_multirun: [], 
        dti_34dir_moco: [], dti_52dir: [], dti_55dir: [], dti_62dir: [], 
        dwi_abcd: [], dwi_117dir: [], dti_trace: [], dti_adc: [], dti_fa: [],
        dti_exp: [], noddi_b2000: [], noddi_b700: [], noddi_b300: [],
        asl: [], asl_mz: [], asl_mp: [], asl_moco: [], pasl: [], pasl_mp: [],
        casl: [], casl_moco: [], fairest: [], fairest_moco: [], fairest_mz: [],
        fairest_moco_mz: []
    }

    # Find how many DTI runs were performed by identifying unique series times.
    dwi_seq = [ s for s in seqinfo if ('dti' in s.protocol_name.lower() or 
        'diff_mddw' in s.protocol_name.lower() or 'noddi' in s.protocol_name.lower() or
        'dmri' in s.protocol_name.lower()) ]
    dwi_times = [ datetime.datetime.strptime(s.date, '%Y-%m-%dT%H:%M:%S.%f') for s in dwi_seq ]
    dwi_times = np.unique(dwi_times).tolist()
    dwi_times.sort()
    
    for s in seqinfo:
        protocol = s.protocol_name.lower()
        if "localizer" in protocol:
            info[locz].append(s.series_id)
        elif 'aahead_scout' in protocol:
            info[locz].append(s.series_id)        
        elif 'vessel_scout' in protocol:
            info[vess_locz].append(s.series_id)      
        elif "t1_3d" in protocol:
            info[t1w_3d].append(s.series_id)
        elif "t1w" in protocol and 'NORM' in s.image_type and 'vnav' not in protocol:
            info[t1w_norm].append(s.series_id)
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
        elif "T1w_MPR_vNav_moco_ND_RMS" in s.series_description:
            info[t1w_vnav_moco_nd].append(s.series_id)
        elif "T1w_MPR_vNav_moco_RMS" in s.series_description:
            info[t1w_vnav_moco].append(s.series_id)
        elif "T1w_MPR_vNav_passive_ND_RMS" in s.series_description:
            info[t1w_vnav_pass_nd].append(s.series_id)
        elif "T1w_MPR_vNav_passive_RMS" in s.series_description:
            info[t1w_vnav_pass].append(s.series_id)
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
        elif "t2_2d_0.4x0.4x1.2mm_180flip" in protocol:
            info[t2w_hippo].append(s.series_id)
        elif "t2w" in protocol and 'NORM' in s.image_type and 'vnav' not in protocol:
            info[t2w_norm].append(s.series_id)
        elif "t2w" in protocol and 'vnav' not in protocol:
            info[t2w].append(s.series_id)
        elif "t2_axial_space" in protocol:
            info[t2w_space].append(s.series_id)
        elif "t2_sag" in protocol:
            info[t2w_sag].append(s.series_id)
        elif "t2_nex_1" in protocol:
            info[t2w_gap].append(s.series_id)
        elif "t2_tse_tra" in protocol:
            info[t2w_gap].append(s.series_id)
        elif "T2w_SPC_vNav_passive" in s.series_description and "_ND" in s.series_description:
            info[t2w_vnav_pass_nd].append(s.series_id)
        elif "T2w_SPC_vNav_passive" in s.series_description:
            info[t2w_vnav_pass].append(s.series_id)
        elif "ax_flair" in protocol:
            info[flair_ax].append(s.series_id)
        elif "axial_flair" in protocol:
            info[flair_ax].append(s.series_id)
        elif "flair_sag" in protocol:
            info[flair_sag].append(s.series_id)
        elif "sagittal_3d_flair" in protocol:
            info[flair_sag].append(s.series_id)
        elif "axial_t2_star" in protocol:
            info[t2star_ax].append(s.series_id)
        elif "pd_tse_tra" in protocol:
            info[pd].append(s.series_id)
        elif "axial_pd-t2_tse" in protocol:
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
        elif "ep2d_se_pcasl_phc_1500ms_moco" in protocol:
            info[asl_moco].append(s.series_id)
        elif "ep2d_se_pcasl_phc_1500ms" in protocol:
            info[asl].append(s.series_id)
        elif "3dsp_pcasl_singleshot_2arms" in protocol:
            info[asl].append(s.series_id)
        elif "pcasl_resting_90mm" in protocol:
            info[asl].append(s.series_id)
        elif "ep2d_casl_ui_1500ms_moco" in protocol:
            info[casl_moco].append(s.series_id)
        elif "ep2d_casl_ui_1500ms" in protocol:
            info[casl].append(s.series_id)
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
        elif "ep2d_fairest_ui_1500ms_moco" in protocol:
            info[fairest_moco].append(s.series_id)
        elif "ep2d_fairest_ui_m0_moco" in protocol:
            info[fairest_moco_mz].append(s.series_id)
        elif "ep2d_fairest_ui_1500ms" in protocol:
            info[fairest].append(s.series_id)
        elif "ep2d_fairest_ui_m0" in protocol:
            info[fairest_mz].append(s.series_id)
        elif 'ep2d_diff_mddw_12' in protocol:
            info[dti_12dir].append(s.series_id)
        elif len(dwi_times) == 1 and 'dti_30dir' in protocol:
            info[dti_30dir].append(s.series_id)
        elif len(dwi_times) == 1 and 'dti_34dir' in protocol:
            info[dti_34dir].append(s.series_id)
        elif len(dwi_times) == 1 and protocol == 'dti_30dir_nodico_vox2_1000':
            info[dti_34dir].append(s.series_id)
        elif len(dwi_times) == 1 and protocol == 'dti_30dir_nodico_vox2_1000_moco':
            info[dti_34dir_moco].append(s.series_id)
        elif len(dwi_times) > 1 and 'dti_30dir' in protocol and 'bval' in s.dcm_dir_name:
            info[dti_30dir_multirun_bval].append(s.series_id)
        elif len(dwi_times) > 1 and 'dti_30dir' in protocol and 'bvec' in s.dcm_dir_name:
            info[dti_30dir_multirun_bvec].append(s.series_id)
        elif len(dwi_times) > 1 and 'dti_30dir' in protocol:
            info[dti_30dir_multirun].append(s.series_id)
        elif len(dwi_times) > 1 and 'dti_34dir' in protocol and 'bval' in s.dcm_dir_name:
            info[dti_34dir_multirun_bval].append(s.series_id)
        elif len(dwi_times) > 1 and 'dti_34dir' in protocol and 'bvec' in s.dcm_dir_name:
            info[dti_34dir_multirun_bvec].append(s.series_id)
        elif len(dwi_times) > 1 and 'dti_34dir' in protocol:
            info[dti_34dir_multirun].append(s.series_id)
        elif len(dwi_times) > 1 and protocol == 'dti_30dir_nodico_vox2_1000' and 'bval' in s.dcm_dir_name:
            info[dti_34dir_multirun_bval].append(s.series_id)
        elif len(dwi_times) > 1 and protocol == 'dti_30dir_nodico_vox2_1000' and 'bvec' in s.dcm_dir_name:
            info[dti_34dir_multirun_bvec].append(s.series_id)
        elif len(dwi_times) > 1 and protocol == 'dti_30dir_nodico_vox2_1000':
            info[dti_34dir_multirun].append(s.series_id)
        elif len(dwi_times) > 1 and protocol == 'dti_30dir_nodico_vox2_1000_moco' and 'bval' in s.dcm_dir_name:
            info[dti_34dir_multirun_moco_bval].append(s.series_id)
        elif len(dwi_times) > 1 and protocol == 'dti_30dir_nodico_vox2_1000_moco' and 'bvec' in s.dcm_dir_name:
            info[dti_34dir_multirun_moco_bvec].append(s.series_id)
        elif len(dwi_times) > 1 and protocol == 'dti_30dir_nodico_vox2_1000_moco':
            info[dti_34dir_multirun_moco].append(s.series_id)
        elif protocol == 'axial_dti' and 'ORIGINAL' in s.series_description:
            info[dti_52dir].append(s.series_id)
        elif 'dti_2x32' in protocol:
            info[dti_32dir_multirun].append(s.series_id)
        elif 'long_dti_p-a' in protocol:
            info[dti_55dir].append(s.series_id)
        elif 'axial_dti' in protocol:
            info[dti_55dir].append(s.series_id)
        elif 'ep2d_dti_30dir_t' in protocol:
            info[dti_62dir].append(s.series_id)
        elif protocol == 'abcd_dmri':
            info[dwi_abcd].append(s.series_id)
        elif protocol == 'multishell_117dir':
            info[dwi_117dir].append(s.series_id)
        elif 'dti' in protocol and 'ADC' in s.image_type:
            info[dti_adc].append(s.series_id)
        elif 'dti' in protocol and 'FA' in s.image_type:
            info[dti_fa].append(s.series_id)
        elif 'dti' in protocol and 'EXP' in s.image_type:
            info[dti_exp].append(s.series_id)
        elif 'dti' in protocol and 'TRACEW' in s.image_type:
            info[dti_trace].append(s.series_id)
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
        elif "b0map_v4" in protocol and "P" in s.image_type:
            info[fm_phasediff2].append(s.series_id)
        elif "b0map_v4" in protocol and "M" in s.image_type:
            info[fm_mag2].append(s.series_id)
        elif "b0_map" in protocol and "P" in s.image_type:
            info[fm_phasediff3].append(s.series_id)
        elif "b0_map" in protocol and "M" in s.image_type:
            info[fm_mag3].append(s.series_id)
        elif "field_mapping" in protocol and "P" in s.image_type:
            info[fm_phasediff4].append(s.series_id)
        elif "field_mapping" in protocol and "M" in s.image_type:
            info[fm_mag4].append(s.series_id)
        elif "b0map" in protocol and "P" in s.image_type:
            info[fm_phasediff5].append(s.series_id)
        elif "b0map" in protocol and "M" in s.image_type:
            info[fm_mag5].append(s.series_id)
        elif "fieldmap_gre" in protocol:
            info[fm_gre].append(s.series_id)
        elif "abcd_dmri_distortionmap_ap" in protocol:
            info[dwi_distmap_ap].append(s.series_id)
        elif "abcd_dmri_distortionmap_pa" in protocol:
            info[dwi_distmap_pa].append(s.series_id)
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
        elif "rfmri_rest_pa" in protocol:
            info[rest_pa].append(s.series_id)
        elif "ax_rsfmri_a>>p" in protocol:
            info[rest_ap].append(s.series_id)
        elif "ax_rsfmri_p>>a" in protocol:
            info[rest_pa].append(s.series_id)
        elif "ax_rsfmri" in protocol:
            info[rest_pa].append(s.series_id)
        elif "axial_rsfmri_(eyes_open)" in protocol:
            info[rest_pa].append(s.series_id)
        elif "bold_resting_2x2x2" in protocol:
            info[rest_bold].append(s.series_id)
        elif "ep2d_bold_restingze6min" in protocol:
            info[rest_bold].append(s.series_id)
        elif "ep2d_max_pace" in protocol and "moco" in protocol:
            info[pace_moco].append(s.series_id)
        elif "ep2d_max_pace" in protocol:
            info[pace].append(s.series_id)
        else:
            print("Series not recognized!: ", protocol, s.dcm_dir_name)
            
    return info
    
def ReplaceSession(sesname):
    return sesname[:13].replace("-", '')

def ReplaceSubject(subjname):
    return subjname[:10].replace("-", '')

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
   }

   fm_phasediff5: {
       "EchoTime1": 0.00412,
       "EchoTime2": 0.00658
   }

}       

IntendedFor = {

    # B0 fieldmap
    fm_phasediff1: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz' ],
    fm_mag1: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz' ],

    fm_phasediff2: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz' ],
    fm_mag2: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz' ],

    fm_phasediff3: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz' ],
    fm_mag3: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz' ],

    fm_phasediff4: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz' ],
    fm_mag4: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz' ],

    fm_phasediff5: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz' ],
    fm_mag5: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz' ],

    fm_ap1: [ '{session}/func/{subject}_{session}_task-rest_dir-AP_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_bold.nii.gz',
        '{session}/dwi/{subject}_{session}_acq-96dir_dwi.nii.gz' ],

#    fm_ap1: [ '{session}/func/{subject}_{session}_task-rest_dir-PA_bold.nii.gz',
#        '{session}/dwi/{subject}_{session}_acq-96dir_dwi.nii.gz' ],
        
    fm_pa1: [ '{session}/func/{subject}_{session}_task-rest_dir-AP_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_bold.nii.gz',
        '{session}/dwi/{subject}_{session}_acq-96dir_dwi.nii.gz' ],

#    fm_pa1: [ '{session}/func/{subject}_{session}_task-rest_dir-PA_bold.nii.gz',
#        '{session}/dwi/{subject}_{session}_acq-96dir_dwi.nii.gz' ],

    fm_gre: [ '{session}/func/{subject}_{session}_task-rest_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_run-{item}_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_rec-moco_run-{item}_bold.nii.gz' ],

    # EPI distortion map for correcting eddy currents along phase-encoding
    # direction in ABCD diffusion protocol. Don't think a P->A map was ever
    # collected.
    dwi_distmap_ap: [ '{session}/dwi/{subject}_{session}_acq-96dir_dwi.nii.gz' ],
    dwi_distmap_pa: [ '{session}/dwi/{subject}_{session}_acq-96dir_dwi.nii.gz' ],

}


