#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##this script will take a .csv file containing a column of "subjects_sessions" and return a boolean column to verify if subject_session exists in the defined project.  

import flywheel
import pandas as pd
import os
fw=flywheel.Client()

#create a project object and read in sublist of interest
group = 'pennftdcenter'
project_label = ''
project= fw.lookup('{}/{}'.format(group,project_label))
audit_subs=pd.read_csv('path/to/sublist')

#create a session level view
view = fw.View(columns='session')
df = fw.read_view_dataframe(view, project.id)

#subet view output to include only columns of interest
sublist = df[['subject.label','session.label','session.timestamp','subject.id', 'session.id','project.label', 'project.id']]

#format dataframe: convert ints to str for boolean operation and merge sub_session for easier indexing
proj_log=df[['subject.label','session.label','session.timestamp']].astype(str)
proj_log['sess_merge'] = proj_log[['subject.label','session.label']].apply(lambda x: '_'.join(x),axis = 1)

#add boolean column to original sublist dataframe and send to .csv
fly_log=proj_log['sess_merge'] 
audit_log=audit_subs['Session']
audit_subs["audit_log"]=audit_log.isin(fly_log)
audit_subs.to_csv('audit_log.csv', index=False)
