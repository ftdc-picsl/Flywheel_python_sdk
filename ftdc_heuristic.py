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

# Changes, 2/17/2021:
# Put derivatives into derivatives folder, e.g., for ADC & FA images created from DTI series.
# ISSS BOLD fMRI: see 101162 20160525-1400

import datetime
import numpy as np
import pandas

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
#'ABCD_T1w_MPR_vNav_setter'
# 'MDDW_DTI_AXIAL_30_DIRECTION'
# ACR_SAG_T1: 2d image, some kind of localizer?
# ACR_AX_T1: 2d image, some kind of localizer?
# MPRAGE_NAVprotocol: these seem to just be the navigator images, not the T1s
# DWI: only 3 directional images plus one B0, not usable for white matter imaging

# Localizers
# 'localizer','Localizer'
# 'AAHead_Scout'
# vessel_scout Vessel_scout_12.nii.gz
locz = create_key('sub-{subject}/{session}/localizer/sub-{subject}_{session}_run-{item}_localizer')
t2w_locz = create_key('sub-{subject}/{session}/localizer/sub-{subject}_{session}_acq-T2w_run-{item}_localizer')
vess_scout = create_key('sub-{subject}/{session}/localizer/sub-{subject}_{session}_acq-vessel_run-{item}_localizer')

# anatomical images
#'t1_mpr_AX_MPRAGE'
# t1_3d_0.8x0.8x0.8 T1_3D_0.8x0.8x0.8_6.nii.gz
# 'MPRAGE_BodyCoil_ref'
# 'MPRAGE_GRAPPA2'
# 'MPRAGE_SAG_ISO'
# 'Sagittal_MP-Rage'
# 't1_se_sag'
# MPRAGE_sag_moco3: this is actually a vNav sequence
# MPRAGE: a 3D T1 acquisition with 1.2 mm sagittal slices and 1.0 x 1.0 y-z in-plane resolution
# 3DT1
# t1_mpr_AX_MPRAGE
# 3DT1
t1w = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_run-{item}_T1w')
t1w_ax = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_run-{item}_T1w')
t1w_sag = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-sag_run-{item}_T1w')
t1w_3d = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-3D_run-{item}_T1w')
t1w_3d_nd = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-3D_rec-ND_run-{item}_T1w')
t1w_norm = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_rec-norm_run-{item}_T1w')
t1w_gradwarp = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-sag_rec-gradwarp_run-{item}_T1w')
t1w_body = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_coil-body_run-{item}_T1w')
t1w_grappa = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-grappa_run-{item}_T1w')
# T1w vNav
t1w_vnav_moco_nd = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-vnavmoco_rec-ND_run-{item}_T1w')
t1w_vnav_pass_nd = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-vnavpass_rec-ND_run-{item}_T1w')
t1w_vnav_moco = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-vnavmoco_run-{item}_T1w')
t1w_vnav_pass = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-vnavpass_run-{item}_T1w')

# T2-weighted anatomical images.
#'t2_tse_obl_448_2mm'
#'HighResHippocampus'
#'t2_tse_CORONAL_HIPPO'
#'t2_tse_tra8CHANNEL'
# t2_tse_cor8channel t2_tse_COR8CHANNEL_5.nii.gz
# 'T2_2D_0.4x0.4x1.2mm_180flip_aa_temp': t2w_hippo
# 'T2_2D_0.4x0.4x1.2mm_180flip_aa_temp_(SECTRA)': t2w_hippo
# 'T2_2D_0.4x0.4x1.2mm_180flip_aa_temp_repeat': t2w_hippo
# 'T2_COR_HIGH_RES': t2w_hippo
# 'T2_Axial_space': t2w_space
# 'T2_sag': t2w_sag
# 't2_nex_1': t2w_gap
# 't2_tse_tra': t2w_gap
# t2_fl2d_tra_hemo_6, t2w_gradwarp if "M" in image_type, t2w_gradwarp_phase if "P" in image_type
# t2_fl2d_tra_hemo, t2w_gradwarp if "M" in image_type, t2w_gradwarp_phase if "P" in image_type
# T2_TSE_AXIAL, t2w_tse_gradwarp
# t2_fl2d_tra_hemo_5, t2w_gradwarp if "M" in image_type, t2w_gradwarp_phase if "P" in image_type
# Axial_T2_ACPC_ANGLE, t2w_gradwarp
t2w = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_run-{item}_T2w')
t2w_norm = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_rec-norm_run-{item}_T2w')
t2w_tse_gradwarp = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-tse_rec-gradwarp_run-{item}_T2w')
t2w_gradwarp = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_rec-gradwarp_run-{item}_T2w')
t2w_gradwarp_phase = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_part-phase_rec-gradwarp_run-{item}_T2w')
t2w_hippo = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-hippo_run-{item}_T2w')
t2w_hippo_nd = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-hippo_rec-ND_run-{item}_T2w')
t2w_space = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-space_run-{item}_T2w')
t2w_sag = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-sag_run-{item}_T2w')
# Kinda making this up, but t2w_gap will include the t2_nex_1 and t2_tse_tra
# sequences--both axial T2 scans with only partial slice coverage due to large 
# gaps (4+ mm) between slices.
t2w_gap = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-gap_run-{item}_T2w')
# T2 SPC vNavs
t2w_vnav_pass_nd = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-vnavpass_rec-ND_run-{item}_T2w')
t2w_vnav_pass = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-vnavpass_run-{item}_T2w')

# T2 star:
# AXT2_STAR_ACPC_ANGLE
# T2_STAR
# t2 T2_5.nii.gz: the "T2" protocol is actually "T2*", but the asterisk gets stripped
t2star = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_run-{item}_T2star')
t2star_acpc = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-acpc_run-{item}_T2star')

# FLAIR
# 'AX_FLAIR_3mm'
# 'Axial_FLAIR_3mm'
# 'FLAIR_sag'
# 'T2_AXIAL_FLAIR'
# FLAIR_3D_1x1x1
# T2_FLAIR_AXIAL
# T2-FLAIR_ACPC_ANGLE
# flair FLAIR_4.nii.gz
# flair FLAIR_5.nii.gz

flair_3d = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-3D_run-{item}_FLAIR')
flair_3d_gradwarp = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-3D_rec-gradwarp_run-{item}_FLAIR')
flair_ax = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-axial_run-{item}_FLAIR')
flair_ax_gradwarp = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-axial_rec-gradwarp_run-{item}_FLAIR')
flair_sag = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-sag_run-{item}_FLAIR')

# SWI (susceptibility-weighted imaging)
swi_mag = create_key('sub-{subject}/{session}/swi/sub-{subject}_{session}_part-mag_run-{item}_GRE')
swi_ph = create_key('sub-{subject}/{session}/swi/sub-{subject}_{session}_part-phase_run-{item}_GRE')
swi_combo = create_key('sub-{subject}/{session}/swi/sub-{subject}_{session}_run-{item}_swi')
swi_minip = create_key('sub-{subject}/{session}/swi/sub-{subject}_{session}_run-{item}_minIP')

# Proton density
#'pd_tse_tra'
pd = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_run-{item}_PD')

# Combined proton density/T2 scan
#'Axial_PD-T2_TSE'
# pdt2 PDT2_6_e1.nii.gz
pdt2 = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_run-{item}_PDT2')

# Angiography
# fl_tof_neck fl_tof_neck_16.nii.gz
# tof_3d_multi-slab TOF_3D_multi-slab_25.nii.gz
# tof_3d_multi-slab TOF_3D_multi-slab_24.nii.gz
tof = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_run-{item}_angio')

# Field maps
#'b0map_v4'
#'B0map_onesizefitsall_v4'
#'SpinEchoFieldMap_AP'
# 'B0_map'
# 'fieldmap_GRE'
# 'SpinEchoFieldMap_PA'
# 'ABCD_dMRI_DistortionMap_AP'
# field_map field_map_10_ph.nii.gz
# field_map field_map_9_e0.nii.gz
# field_map field_map_9_e1.nii.gz

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
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_run-{item}_dir-AP_epi')
fm_pa1 = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_run-{item}_dir-PA_epi')
fm_gre = create_key(
    'sub-{subject}/{session}/fmap/sub-{subject}_{session}_run-{item}_fieldmap')

# BOLD fMRI
# 'Resting_bold_124'
#
# 'Ax_rsfMRI'
# 'Ax_rsfMRI_A>>P'
# 'Ax_rsfMRI_P>>A'
# Ax rsfMRI A>>P, Ax rsfMRI P>>A
# 'BOLD_resting_2x2x2'
# 'ep2d_bold_restingZE6min'
# 'long_rsfMRI_P-A'
# rfMRI_REST_AP
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
# ep2d_pace_max_1 ep2d_pace_max_1_3.nii.gz
# ep2d_pace_max_1_moco ep2d_pace_max_1_MoCo_4.nii.gz
# ep2d_pace_max_165 ep2d_pace_max_165_5.nii.gz
# ep2d_pace_max_165_moco ep2d_pace_max_165_MoCo_6.nii.gz
# ep2d_pace_max_205 ep2d_pace_max_205_7.nii.gz
# ep2d_pace_max_205_moco ep2d_pace_max_205_MoCo_8.nii.gz
# fmri fMRI_12.nii.gz
# fmri fMRI_13.nii.gz
# Axial MB rsfMRI (Eyes Open)
bold_mb = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_acq-MB_run-{item}_bold')
rest_bold = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_run-{item}_bold')
rest_ap = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_dir-AP_run-{item}_bold')
rest_pa = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_dir-PA_run-{item}_bold')
rest_ap_sbref = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_dir-AP_run-{item}_sbref')
rest_pa_sbref = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_dir-PA_run-{item}_sbref')
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
# ep2d_diff_3scan_trace_p2, dti_trace if "TRACEW" or "ORIGINAL" in image_type, or dti_adc if "ADC" in image_type
# DTI_P-A, dti_30dir_pa
# DTI_A-P, dti_30dir_ap
# DTI_P-A_BW2394, dti_30dir_pa
# DTI_A-P_BW2394, dti_30dir_ap
# DTI_34_DIR
# DTI_12dir
#        elif "abcd_dmri_distortionmap_ap" in protocol:
#            info[dwi_distmap_ap].append(s.series_id)
#        elif "abcd_dmri_distortionmap_pa" in protocol:
#            info[dwi_distmap_pa].append(s.series_id)
dti_12dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-12dir_dwi')
dti_30dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_dwi')
dti_30dir_ap = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_dir-AP_dwi')
dti_30dir_pa = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_dir-PA_dwi')
dti_30dir_run1 = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_run-1_dwi')
dti_30dir_run2 = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_run-2_dwi')
dti_30dir_run3 = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_run-3_dwi')
dti_32dir_run1 = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-32dir_run-1_dwi')
dti_32dir_run2 = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-32dir_run-2_dwi')
dti_32dir_run3 = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-32dir_run-3_dwi')
dti_34dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_dwi')
dti_34dir_run1 = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_run-1_dwi')
dti_34dir_run2 = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_run-2_dwi')
dti_34dir_run3 = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_run-3_dwi')
dti_34dir_moco = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_rec-moco_dwi')
dti_34dir_run1_moco = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_rec-moco_run-1_dwi')
dti_34dir_run2_moco = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_rec-moco_run-2_dwi')
dti_34dir_run3_moco = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_rec-moco_run-3_dwi')
dti_52dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-52dir_dwi')
dti_55dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-55dir_dwi')
dti_62dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-62dir_dwi')
dwi_abcd = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-96dir_dwi')
dwi_117dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-117dir_dwi')
noddi_b2000 =  create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-NODDIB2000_dwi')
noddi_b700 =  create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-NODDIB700_dwi')
noddi_b300 =  create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-NODDIB300_dwi')

# Leave ABCD DWI distortion maps as non-BIDS.
#dwi_distmap_ap = create_key(
#    'sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-dwi_dir-AP_epi')
#dwi_distmap_pa = create_key(
#    'sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-dwi_dir-PA_epi')
# Leave derived images (trace, ADC, FA, colFA, and exp) as non-BIDS.
#dti_trace = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_desc-trace_dwi')
#dti_adc = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_desc-adc_dwi')
#dti_fa = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_desc-fa_dwi')
#dti_colfa = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_desc-colfa_dwi')
#dti_exp = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_desc-exp_dwi')
#dti_trace_55dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_desc-trace_acq-55dir_dwi')
#dti_adc_55dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_desc-adc_acq-55dir_dwi')
#dti_fa_55dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_desc-fa_acq-55dir_dwi')
#dti_colfa_55dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_desc-colfa_acq-55dir_dwi')
#dti_exp_55dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_desc-exp_acq-55dir_dwi')
# Leave derived files as non-BIDS.
#        dti_trace: [], dti_adc: [], dti_fa: [], dti_colfa: [], dti_exp: [],
#        dti_trace_55dir: [], dti_adc_55dir: [], dti_fa_55dir: [], dti_colfa_55dir: [], dti_exp_55dir: [],
#        elif 'long_dti_p-a' in protocol and 'ADC' in s.image_type:
#            info[dti_adc_55dir].append(s.series_id)
#        elif 'long_dti_p-a' in protocol and 'ColFA' in s.series_description:
#            info[dti_colfa_55dir].append(s.series_id)
#        elif 'long_dti_p-a' in protocol and 'FA' in s.image_type:
#            info[dti_fa_55dir].append(s.series_id)
#        elif 'long_dti_p-a' in protocol and 'EXP' in s.image_type:
#            info[dti_exp_55dir].append(s.series_id)
#        elif 'long_dti_p-a' in protocol and 'TRACEW' in s.image_type:
#            info[dti_trace_55dir].append(s.series_id)
#        elif ('axial_dti_a>>p' in protocol or 'axial_dti_p>>a' in protocol) and 'ADC' in s.image_type:
#            info[dti_adc_55dir].append(s.series_id)
#        elif ('axial_dti_a>>p' in protocol or 'axial_dti_p>>a' in protocol) and 'ColFA' in s.series_description:
#            info[dti_colfa_55dir].append(s.series_id)
#        elif ('axial_dti_a>>p' in protocol or 'axial_dti_p>>a' in protocol) and 'FA' in s.image_type:
#            info[dti_fa_55dir].append(s.series_id)
#        elif ('axial_dti_a>>p' in protocol or 'axial_dti_p>>a' in protocol) and 'EXP' in s.image_type:
#            info[dti_exp_55dir].append(s.series_id)
#        elif ('axial_dti_a>>p' in protocol or 'axial_dti_p>>a' in protocol) and 'TRACEW' in s.image_type:
#            info[dti_trace_55dir].append(s.series_id)
#        elif 'dti' in protocol and 'ADC' in s.image_type:
#            info[dti_adc].append(s.series_id)
#        elif 'dti' in protocol and 'ColFA' in s.series_description:
#            info[dti_colfa].append(s.series_id)
#        elif 'dti' in protocol and 'FA' in s.image_type:
#            info[dti_fa].append(s.series_id)
#        elif 'dti' in protocol and 'EXP' in s.image_type:
#            info[dti_exp].append(s.series_id)
#        elif 'dti' in protocol and 'TRACEW' in s.image_type:
#            info[dti_trace].append(s.series_id)
#        elif 'ep2d_diff_3scan_trace' in protocol and 'TRACEW' in s.image_type:
#            info[dti_trace].append(s.series_id)
#        elif 'ep2d_diff_3scan_trace' in protocol and 'ORIGINAL' in s.image_type:
#            info[dti_trace].append(s.series_id)
#        elif 'ep2d_diff_3scan_trace' in protocol and 'ADC' in s.image_type:
#            info[dti_adc].append(s.series_id)
#dti_30dir_run1_bval = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_run-1_dwi')
#dti_30dir_run1_bvec = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_run-1_dwi')
#dti_30dir_run2_bval = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_run-2_dwi')
#dti_30dir_run2_bvec = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-30dir_run-2_dwi')
#dti_32dir_run1_bval = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-32dir_run-1_dwi')
#dti_32dir_run1_bvec = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-32dir_run-1_dwi')
#dti_34dir_run1_bval = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_run-1_dwi')
#dti_34dir_run1_bvec = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_run-1_dwi')
#dti_34dir_run1_moco_bval = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_rec-moco_run-1_dwi')
#dti_34dir_run1_moco_bvec = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-34dir_rec-moco_run-1_dwi')

#ASL
# Input based on ASL draft certification?
#'3DSP_pcasl_singleshot_2arms'
#'Axial_3D_PASL_(Eyes_Open)'
#'ep2d_se_pcasl_PHC_1500ms_MoCo'
#'ep2d_se_pcasl_PHC_1500ms'
#'3D_PASL'
#'ep2d_se_pcasl_M0'
# 3D_PCASL
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
# ep2d_casl_1500ms ep2d_casl_1500ms_6.nii.gz
# ep2d_casl_1500ms ep2d_casl_1500ms_8.nii.gz
# ep2d_casl_1500ms_moco ep2d_casl_1500ms_MoCo_7.nii.gz
# ep2d_casl_am_ui_1500ms ep2d_casl_AM_UI_1500ms_3.nii.gz
# ep2d_casl_am_ui_1500ms_moco ep2d_casl_AM_UI_1500ms_MoCo_4.nii.gz
asl = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_run-{item}_asl')
asl_mz = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_run-{item}_m0scan')
asl_mp = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_run-{item}_cbf')
asl_moco = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_rec-moco_run-{item}_asl')
pcasl_3d = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-3D_run-{item}_asl')
pcasl_3d_mz = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-3D_run-{item}_m0scan')
pcasl_3d_mp = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-3D_run-{item}_cbf')
# Should we use the acq entity like this? How specific should filenames strive to be?
pasl = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-pasl_run-{item}_asl')
pasl_mp = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-pasl_run-{item}_cbf')
casl = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-casl_run-{item}_asl')
casl_mz = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-casl_run-{item}_m0scan')
casl_moco = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-casl_rec-moco_run-{item}_asl')
fairest = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-fairest_run-{item}_asl')
fairest_moco = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-fairest_rec-moco_run-{item}_asl')
fairest_mz = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-fairest_run-{item}_m0scan')
fairest_moco_mz = create_key('sub-{subject}/{session}/asl/sub-{subject}_{session}_acq-fairest_rec-moco_run-{item}_m0scan')

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
        locz: [], t2w_locz: [], vess_scout: [], t1w: [], t1w_ax: [], t1w_sag: [],
        t1w_3d: [], t1w_3d_nd: [], t1w_norm: [], t1w_gradwarp: [],
        t1w_body: [], t1w_grappa: [], t1w_vnav_moco_nd: [], t1w_vnav_pass_nd: [], 
        t1w_vnav_moco: [], t1w_vnav_pass: [],
        t2w: [], t2w_norm: [], t2w_gradwarp: [], t2w_gradwarp_phase: [],
        t2w_tse_gradwarp: [], t2w_hippo: [], t2w_hippo_nd: [], t2w_space: [],
        t2w_sag: [], t2w_gap: [], t2w_vnav_pass_nd: [], t2w_vnav_pass: [],
        t2star: [], t2star_acpc: [],
        flair_3d: [], flair_3d_gradwarp: [], flair_ax: [], flair_ax_gradwarp: [],
        flair_sag: [], 
        swi_mag: [], swi_ph: [], swi_combo: [], swi_minip: [],
        pd: [], pdt2: [], tof: [],
        fm_phasediff1: [], fm_mag1: [], fm_phasediff2: [], fm_mag2: [],
        fm_phasediff3: [], fm_mag3: [], fm_phasediff4: [], fm_mag4: [], 
        fm_phasediff5: [], fm_mag5: [], fm_phasediff6: [], fm_mag6: [],
        fm_ap1: [], fm_pa1: [], fm_gre: [],
        bold_mb: [], rest_bold: [], rest_ap: [], rest_pa: [], rest_ap_sbref: [],
        rest_pa_sbref: [], pace: [], pace_moco: [],
        dti_12dir: [], dti_30dir: [], dti_30dir_ap: [], dti_30dir_pa: [],
        dti_30dir_run1: [], dti_30dir_run2: [], dti_30dir_run3: [],
        dti_32dir_run1: [], dti_32dir_run2: [], dti_32dir_run3: [],
        dti_34dir: [], dti_34dir_moco: [],
        dti_34dir_run1: [], dti_34dir_run2: [], dti_34dir_run3: [],
        dti_52dir: [], dti_55dir: [], dti_62dir: [], 
        dwi_abcd: [], dwi_117dir: [],
        noddi_b2000: [], noddi_b700: [], noddi_b300: [],
        asl: [], asl_mz: [], asl_mp: [], asl_moco: [], 
        pcasl_3d: [], pcasl_3d_mz: [], pcasl_3d_mp: [],
        pasl: [], pasl_mp: [], casl: [], casl_mz: [], casl_moco: [], 
        fairest: [], fairest_moco: [], fairest_mz: [],
        fairest_moco_mz: []
    }

# Find how many DTI runs were performed by identifying unique series times.
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
# DTI_P-A: 30-direction sequence + 2 B0 scans
    # Find all of the DWI sequences that were run more than once in the session.
    dwi_seq = [ s for s in seqinfo if ('dti' in s.protocol_name.lower() or 
        'diff_mddw' in s.protocol_name.lower() or 'noddi' in s.protocol_name.lower() or
        'dmri' in s.protocol_name.lower() or 'multishell_117dir' in s.protocol_name.lower()) ]
    dwi_names = [ s.protocol_name.lower() for s in dwi_seq ]
    n_scans = pandas.value_counts(dwi_names)
    dwi_multi_names = [n for n in n_scans.index if n_scans[n] > 1]
    if len(dwi_multi_names) > 0:
        dwi_multi = [ s for s in dwi_seq if s.protocol_name.lower() in dwi_multi_names ]
        dwi_times = [ datetime.datetime.strptime(s.date, '%Y-%m-%dT%H:%M:%S.%f') for s in dwi_multi ]
    else:
        dwi_times = [ datetime.datetime.strptime(s.date, '%Y-%m-%dT%H:%M:%S.%f') for s in dwi_seq ]    
    dwi_times = np.unique(dwi_times).tolist()
    dwi_times.sort()
    
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
        elif s.series_description == '3DT1_ND':
            info[t1w_3d_nd].append(s.series_id)
        elif protocol == '3dt1':
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
        elif "T1w_MPR_vNav_moco_ND_RMS" in s.series_description:
            info[t1w_vnav_moco_nd].append(s.series_id)
        elif "T1w_MPR_vNav_moco_RMS" in s.series_description:
            info[t1w_vnav_moco].append(s.series_id)
        elif "MPRAGE_sag_moco3" in s.series_description and "Moco3d1_ns" in s.sequence_name:
            info[t1w_vnav_moco].append(s.series_id)
        elif "T1w_MPR_vNav_passive_ND_RMS" in s.series_description:
            info[t1w_vnav_pass_nd].append(s.series_id)
        elif "T1w_MPR_vNav_passive_RMS" in s.series_description:
            info[t1w_vnav_pass].append(s.series_id)
        elif "T1w_MPR_vNav_passive_(sectra)_ND_RMS" in s.series_description:
            info[t1w_vnav_pass_nd].append(s.series_id)
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
        elif "t2_2d_0.4x0.4x1.2mm_180flip" in protocol and "ND" in s.series_description:
            info[t2w_hippo_nd].append(s.series_id)
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
        elif "T2w_SPC_vNav_passive" in s.series_description and "_ND" in s.series_description:
            info[t2w_vnav_pass_nd].append(s.series_id)
        elif "T2w_SPC_vNav_passive" in s.series_description:
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
        elif "ep2d_casl_am_ui_1500ms" in protocol and s.series_description == 'MoCoSeries':
            info[casl_moco].append(s.series_id)
        elif "ep2d_casl_am_ui_1500ms" in protocol:
            info[casl].append(s.series_id)
        elif "ep2d_casl_1500ms" in protocol and s.series_description == 'MoCoSeries':
            info[casl_moco].append(s.series_id)
        elif "ep2d_casl_1500ms" in protocol and not "MOSAIC" in s.image_type:
            info[casl_mz].append(s.series_id)
        elif "ep2d_casl_1500ms" in protocol:
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
        elif 'dti_12dir' in protocol:
            info[dti_12dir].append(s.series_id)
        elif 'ep2d_dti_30dir_t' in protocol:
            info[dti_62dir].append(s.series_id)
        elif protocol == 'dti_30dir_nodico_vox2_1000' and not protocol in dwi_multi_names:
            info[dti_34dir].append(s.series_id)
        elif protocol == 'dti_30dir_nodico_vox2_1000_moco' and not protocol in dwi_multi_names:
            info[dti_34dir_moco].append(s.series_id)
        elif 'dti_30dir' in protocol and not 'dti_30dir' in dwi_multi_names:
            info[dti_30dir].append(s.series_id)
        elif s.series_description == 'DTI_A-P':
            info[dti_30dir_ap].append(s.series_id)
        elif s.series_description == 'DTI_P-A':
            info[dti_30dir_pa].append(s.series_id)
        elif s.series_description == 'DTI_A-P_BW2394':
            info[dti_30dir_ap].append(s.series_id)
        elif s.series_description == 'DTI_P-A_BW2394':
            info[dti_30dir_pa].append(s.series_id)
        elif 'dti_34dir' in protocol and not 'dti_34dir' in dwi_multi_names:
            info[dti_34dir].append(s.series_id)
        elif 'dti_34_dir' in protocol and not 'dti_34_dir' in dwi_multi_names:
            info[dti_34dir].append(s.series_id)
        elif len(dwi_times) > 1 and mydatetime == dwi_times[0] and 'dti_30dir' in protocol and 'dti_30dir' in dwi_multi_names:
            info[dti_30dir_run1].append(s.series_id)
        elif len(dwi_times) > 1 and mydatetime == dwi_times[1] and 'dti_30dir' in protocol and 'dti_30dir' in dwi_multi_names:
            info[dti_30dir_run2].append(s.series_id)
        elif len(dwi_times) > 2 and mydatetime == dwi_times[2] and 'dti_30dir' in protocol and 'dti_30dir' in dwi_multi_names:
            info[dti_30dir_run3].append(s.series_id)
        elif len(dwi_times) > 1 and mydatetime == dwi_times[0] and 'dti_34dir' in protocol and 'dti_34dir' in dwi_multi_names:
            info[dti_34dir_run1].append(s.series_id)
        elif len(dwi_times) > 1 and mydatetime == dwi_times[1] and 'dti_34dir' in protocol and 'dti_34dir' in dwi_multi_names:
            info[dti_34dir_run2].append(s.series_id)
        elif len(dwi_times) > 2 and mydatetime == dwi_times[2] and 'dti_34dir' in protocol and 'dti_34dir' in dwi_multi_names:
            info[dti_34dir_run3].append(s.series_id)
        elif len(dwi_times) > 1 and mydatetime == dwi_times[0] and 'dti_34_dir' in protocol and 'dti_34_dir' in dwi_multi_names:
            info[dti_34dir_run1].append(s.series_id)
        elif len(dwi_times) > 1 and mydatetime == dwi_times[1] and 'dti_34_dir' in protocol and 'dti_34_dir' in dwi_multi_names:
            info[dti_34dir_run2].append(s.series_id)
        elif len(dwi_times) > 2 and mydatetime == dwi_times[2] and 'dti_34_dir' in protocol and 'dti_34_dir' in dwi_multi_names:
            info[dti_34dir_run3].append(s.series_id)
        elif len(dwi_times) > 1 and mydatetime == dwi_times[0] and protocol == 'dti_30dir_nodico_vox2_1000' and protocol in dwi_multi_names:
            info[dti_34dir_run1].append(s.series_id)
        elif len(dwi_times) > 1 and mydatetime == dwi_times[1] and protocol == 'dti_30dir_nodico_vox2_1000' and protocol in dwi_multi_names:
            info[dti_34dir_run2].append(s.series_id)
        elif len(dwi_times) > 2 and mydatetime == dwi_times[2] and protocol == 'dti_30dir_nodico_vox2_1000' and protocol in dwi_multi_names:
            info[dti_34dir_run3].append(s.series_id)
        elif len(dwi_times) > 1 and mydatetime == dwi_times[0] and protocol == 'dti_30dir_nodico_vox2_1000_moco' and protocol in dwi_multi_names:
            info[dti_34dir_run1_moco].append(s.series_id)
        elif len(dwi_times) > 1 and mydatetime == dwi_times[1] and protocol == 'dti_30dir_nodico_vox2_1000_moco' and protocol in dwi_multi_names:
            info[dti_34dir_run2_moco].append(s.series_id)
        elif len(dwi_times) > 2 and mydatetime == dwi_times[2] and protocol == 'dti_30dir_nodico_vox2_1000_moco' and protocol in dwi_multi_names:
            info[dti_34dir_run3_moco].append(s.series_id)
        elif protocol == 'axial_dti' and 'ORIGINAL' in s.series_description:
            info[dti_52dir].append(s.series_id)
        elif 'dti_2x32' in protocol:
            info[dti_32dir_run1].append(s.series_id)
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
        elif "ep2d_max_pace" in protocol and "MoCo" in s.series_description:
            info[pace_moco].append(s.series_id)
        elif "ep2d_max_pace" in protocol:
            info[pace].append(s.series_id)
        elif "ep2d_pace_max" in protocol and "MoCo" in s.series_description:
            info[pace_moco].append(s.series_id)
        elif "ep2d_pace_max" in protocol:
            info[pace].append(s.series_id)
        else:
            print("Series not recognized!: ", protocol, s.dcm_dir_name)
            
    return info
    
def ReplaceSession(sesname):
    return sesname[:13].replace("-", "").replace(".","").replace("_","")

def ReplaceSubject(subjname):
    return subjname[:10].replace("-", "").replace(".","").replace("_","")

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
    fm_phasediff1: [ '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],
    fm_mag1: [ '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],

    fm_phasediff2: [ '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],
    fm_mag2: [ '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],

    fm_phasediff3: [ '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],
    fm_mag3: [ '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],

    fm_phasediff4: [ '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz',
        '{session}/func/sub-{subject}_{session}_task-rest_acq-MB_run-{item}_bold.nii.gz'],
    fm_mag4: [ '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz',
        '{session}/func/sub-{subject}_{session}_task-rest_acq-MB_run-{item}_bold.nii.gz'],


    fm_phasediff5: [ '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],
    fm_mag5: [ '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],

    fm_phasediff6: [ '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],
    fm_mag6: [ '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz'],

    fm_ap1: [ '{session}/func/{subject}_{session}_task-rest_dir-AP_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-AP_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-AP_run-1_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-AP_run-2_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_run-1_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_run-2_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz',
        '{session}/dwi/{subject}_{session}_acq-96dir_dwi.nii.gz' ],
        
    fm_pa1: [ '{session}/func/{subject}_{session}_task-rest_dir-AP_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-AP_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-AP_run-1_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-AP_run-2_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_run-1_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_dir-PA_run-2_sbref.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz',
        '{session}/dwi/{subject}_{session}_acq-96dir_dwi.nii.gz' ],

    fm_gre: [ '{session}/func/{subject}_{session}_task-rest_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_run-4_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_rec-moco_run-1_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_rec-moco_run-2_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_rec-moco_run-3_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_run-4_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_rec-moco_run-4_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_run-5_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_rec-moco_run-5_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_run-6_bold.nii.gz',
        '{session}/func/{subject}_{session}_task-rest_acq-pace_rec-moco_run-6_bold.nii.gz'],

    # EPI distortion map for correcting eddy currents along phase-encoding
    # direction in ABCD diffusion protocol. Don't think a P->A map was ever
    # collected.
    # Edit, JSP, 12/07/2020: the BIDS validator doesn't like this distortion map, so make
    # it a non-BIDS file.
    # dwi_distmap_ap: [ '{session}/dwi/{subject}_{session}_acq-96dir_dwi.nii.gz' ],
    # dwi_distmap_pa: [ '{session}/dwi/{subject}_{session}_acq-96dir_dwi.nii.gz' ],

}
