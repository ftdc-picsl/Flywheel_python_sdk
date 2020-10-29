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

file_log=[]
run=[]
sub_log=[]
zip_log=[]
file_id_log=[]
for i in range(0,4):
    session = fw.get(sesslist[i])
    acqs=fw.get_session_acquisitions(session.id)
    for j,a in enumerate(acqs):
      files=a.files
      types=[x.type for x in files]
      if('dicom' in types) and ('nifti' in types):
          print(acqs[j].label)
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

