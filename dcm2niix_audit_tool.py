#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 13:55:23 2020

@author: seburke
"""

import flywheel
import pandas as pd

fw=flywheel.Client()

group = 'pennftdcenter'
project = ''

sesslist=pd.read_csv("path/to/csv")


for i in range(0,len(sesslist)-1):
    sess=fw.lookup('{}/{}/{}/{}'.format(group,project,sesslist['subject'][i],sesslist['session'][i]))
    sesslist.append(sess.id)

file_log=[]
run=[]
sub_log=[]
zip_log=[]
file_id_log=[]
acq_list=[]
for i in range(0,2):
    session = fw.get(sesslist[i])
    acqs=fw.get_session_acquisitions(session.id)
    for j,a in enumerate(acqs):
      files=a.files
      types=[x.type for x in files]
      if('dicom' in types) and ('nifti' in types):
          acq_list.append(acqs[j].id)
          fw.get(acqs[0]["parents"]["subject"]).label
          sub_log.append(fw.get(acqs[0]["parents"]["subject"]).label)
          file_log.append(a.files[0].name)
          zip_log.append(a.files[0]["zip_member_count"])
          file_id_log.append(files[1]['_id'])
      else:    
          run.append(a.files[0].name)
                      
d = {'subject':sub_log,'file':file_log,'zip_size':zip_log}
success_log=pd.DataFrame(d)
fail_log=pd.DataFrame(run)


PixelSpacing=[]
SliceThickness=[]
EchoTime=[]
Modality=[]
RepetitionTime=[]
LargestImagePixelValue=[]    
ImageType=[]
sub_list=[]
sess_log=[]
for a in range(0,len(acq_list)-1):
    scan_of_interest=fw.get(acq_list[a])
    sub_list.append(fw.get(scan_of_interest["parents"]["subject"]).label)
    sess_log.append(fw.get(scan_of_interest["parents"]["session"]).label)
    PixelSpacing.append(scan_of_interest['files'][0]['info']['PixelSpacing'])
    SliceThickness.append(scan_of_interest['files'][0]['info']['SliceThickness'])
    EchoTime.append(scan_of_interest['files'][0]['info']['EchoTime'])
    Modality.append(scan_of_interest['files'][0]['info']['Modality'])
    RepetitionTime.append(scan_of_interest['files'][0]['info']['RepetitionTime'])
    LargestImagePixelValue.append(scan_of_interest['files'][0]['info']['LargestImagePixelValue'])    
    ImageType.append(scan_of_interest['files'][0]['info']['ImageType'][0])
    
#xx=pd.DataFrame(scan_of_interest['files'][0]['info'].keys())
e={'subject':sub_list,'session':sess_log,'PixelSpacing':PixelSpacing,'SliceThickness':SliceThickness,'EchoTime':EchoTime,'Modality':Modality,'RepetitionTime':RepetitionTime,'LargestImagePixelValue':LargestImagePixelValue,'ImageType':ImageType}

parameter_log=pd.DataFrame(e)
