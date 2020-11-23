#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 15:02:49 2020

@author: seburke
"""
import flywheel
import pandas as pd

fw=flywheel.Client()

group = 'pennftdcenter'
project = ''

#csv with column "subject" and column "session"
sessions=pd.read_csv("path/to/csv")

sesslist=[]
for i in range(0,len(sessions)-1):
    sess=fw.lookup('{}/{}/{}/{}'.format(group,project,sessions['subject'][i],sessions['session'][i]))
    sesslist.append(sess.id)


file_log=[]
run=[]
sub_log=[]
file_id_log=[]
acq_list=[]
sess_log=[]
run_sub=[]
run_sess=[]
for i in range(0,len(sesslist)-1):    
    session = fw.get(sesslist[i])
    acqs=fw.get_session_acquisitions(session.id)
    for j,a in enumerate(acqs):
      files=a.files
      types=[x.type for x in files]
      if('dicom' in types) and ('nifti' not in types):
          run.append(a.files[0].name)
          run_sub.append(fw.get(acqs[j]["parents"]["subject"]).label)
          run_sess.append(fw.get(acqs[j]["parents"]["session"]).label)
      else:
          acq_list.append(acqs[j].id)
          fw.get(acqs[0]["parents"]["subject"]).label
          sub_log.append(fw.get(acqs[0]["parents"]["subject"]).label)
          sess_log.append(fw.get(acqs[0]["parents"]["session"]).label)
          file_log.append(a.files[1].name)
          file_id_log.append(files[1]['_id'])

d = {'subject':sub_log,'session':sess_log,'file':file_log}
run_log = {'subject':run_sub, 'session':run_sess,'file':run}
success_log=pd.DataFrame(d)
fail_log=pd.DataFrame(run_log)

#########extraction of MR parameters, needs troubleshooting for empty dictionary fields.
###this command will list all keys in acquisition.file
#ll=pd.DataFrame(scan_of_interest['files'][0]['info'].keys())

PixelSpacing=[]
SliceThickness=[]
EchoTime=[]
Modality=[]
RepetitionTime=[]
Acquisition_matrix=[]
ImageType=[]
sub_list=[]
sess_log=[]
flip_angle=[]
scan_type=[]
zip_count=[]
for a in range(0,len(acq_list)-1):
    scan_of_interest=fw.get(acq_list[a])
    sub_list.append(fw.get(scan_of_interest["parents"]["subject"]).label)
    sess_log.append(fw.get(scan_of_interest["parents"]["session"]).label)
    scan_type.append(scan_of_interest['files'][0]['name'])
    EchoTime.append(scan_of_interest['files'][0]['info']['EchoTime'])
    Modality.append(scan_of_interest['files'][0]['info']['Modality'])
    RepetitionTime.append(scan_of_interest['files'][0]['info']['RepetitionTime'])
    Acquisition_matrix.append(scan_of_interest['files'][0]['info']['Rows'])
    ImageType.append(scan_of_interest['files'][0]['info']['ImageType'][0])
    flip_angle.append(scan_of_interest['files'][0]['info']['FlipAngle'])
    zip_count.append(scan_of_interest['files'][0]['zip_member_count'])
    if not bool(scan_of_interest['files'][0]['info']):
        PixelSpacing.append(["NA"])
    else:
        PixelSpacing.append(scan_of_interest['files'][0]['info']['PixelSpacing'])
    if not bool(scan_of_interest['files'][0]['info']):
        SliceThickness.append("NA")
    else:
        SliceThickness.append(scan_of_interest['files'][0]['info']['SliceThickness'])

e={'subject':sub_list,'session':sess_log,'scantype':scan_type,'EchoTime':EchoTime,'Modality':Modality,'RepetitionTime':RepetitionTime,'Acquisition_matrix':Acquisition_matrix,'ImageType':ImageType,'FlipAngle':flip_angle,'SliceThickness':SliceThickness,'PixelSpacing':PixelSpacing,'ZipCount':zip_count}
parameter_log=pd.DataFrame(e)


parameter_log.to_csv('/Users/seburke/projects/flywheel_migration/parameter_log.csv',index=False)
fail_log.to_csv('/Users/seburke/projects/flywheel_migration/failed_nii.csv',index=False)
success_log.to_csv('/Users/seburke/projects/flywheel_migration/success_nii.csv',index=False)
