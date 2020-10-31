#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Version as of Thu Oct  8 17:42:07 2020

@author: isabelannwingert
"""

# logging into Flywheel & importing packages
import flywheel
import os # to access input list / excel sheet
import sys
import pandas as pd # to read in input list / excel sheet
fw = flywheel.Client() # logging into fw

# group id = group name
# fetching info container into variable
# fw.get() uses hash ids to print out info container
group_id = 'pennftdcenter'
group = fw.get(group_id)

# accessing particular project with prompt for user
# & printing only that hash id and label
# fetching project info container into project
input_project = input("Please type in the exact name of the project you are using: ")
print('searching selected project... \n ... \n')
projects_labels = []
projects_ids = []
for project in fw.projects():
    projects_labels.append(project.label)
    projects_ids.append(project.id)

projects_dict = dict(zip(projects_labels, projects_ids)) # for use if user wants to see label and ids

  
  
if input_project in projects_labels:
    print('%s: %s \n \n *****' % (input_project, projects_dict[input_project]))
    project = fw.get(projects_dict[input_project])
    project_label = project.label
    project_id = project.id
else:
    raise SystemExit('The project you typed in was not found in FlyWheel.' 
                     ' Please make sure you typed in the project name correctly,' 
                     ' or check your permissions in FlyWheel. Exiting...')


# Prompting user to input .csv file of session and subject labels
input_prompt = input("Do you have a list of subjects that you want to run the dcm2niix gear on? Please type 'yes' or press enter if no: ") # Cancel if no
if input_prompt == "yes":
    input_directory = input("Please type in which directory the list of your subjects exists (.csv): ")
    os.chdir(input_directory)
    input_list = input("Please type in the name of your list file: ")
    test = pd.read_csv(input_list, dtype = str)
    #test['Session'] = test['Session'].str[1:-1] # Session Labels may need quotation marks due to date-time interpretation (@Phil Cook ask IW about this if you see this)
    print('\n \n printing your inputted sessions under %s \n Column1 = Session hash-id Column2 = Session label Column 3 = Subject label...' % (input_project))
    # print(test)
    #test_dict = dict(zip(test.SubjectLabel, test.SessionLabel))
    subject_labels = test['SubjectLabel'].to_list()
    session_labels = test['SessionLabel'].to_list()
    subjects_ids = []
    sessions_ids = []
    for i in range(len(subject_labels)):
        #subject = fw.lookup('{}/{}/{}'.format(group_id, project_label, k))
        #print('%s : %s' % (subject.id, subject.label))
        #subjects_ids.append(subject.id)
        session = fw.lookup('{}/{}/{}/{}'.format(group_id, project.label, subject_labels[i], session_labels[i]))
        subject = fw.get(session.parents.subject)
        print('%s : %s : %s' % (session.id, session.label, subject.label))
        sessions_ids.append(session.id)          


# This was just done for writing purposes - 
# to find the gear label for fw-heudiconv that is FlyWheel specific
# this loop is useful for looking at the entire list of all gears on FlyWheel
# for gear in fw.get_all_gears():
    # print(gear.gear.name, ':', gear.gear.label)
    
# From the above loop, fw-heudiconv is the gear of interest
# FlyWheel has a lookup function that loads the container into a variable
fwheudiconv = fw.lookup('gears/fw-heudiconv')
print(fwheudiconv.gear.label)

# Importing Datetime for analysis label differentiation with timestamps
import datetime
now = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M")
for i in range(len(subject_labels)):
    analysis_label = '{}_{}_{}_{}_{}'.format(fwheudiconv.gear.name, 
                                   subject_labels[i], 
                                   session_labels[i],
                                   now,
                                   fwheudiconv.gear.version)
    print(analysis_label)

# The gear will look like this:
# gearname_subjectlabel_sessionlabel_acq_label_date_time
# Example: fw-heudiconv_119156_20170609-1051_2020-10-10_11:46_0.2.15_0.3.3
# print(analysis_label is for book-keeping)


config = {
    "action": "Curate",
    "default_heuristic": "",
    'do_whole_project': False,
    'dry_run': False,
    'extended_bids': True
    }

for index_num, file_obj in enumerate(project.files):
    print(index_num, ":", file_obj.name)
    
inputs = {'heuristic': project.files[3785]
    }

sessions_to_run = []
for i in sessions_ids:
    session = fw.get(i)
    sessions_to_run.append(session)
len(sessions_to_run)

# Potential prompt if running on CL: Do these analysis labels look correct?
analysis_ids = []
fails = []
for ses in sessions_to_run:
    try:
        _id = fwheudiconv.run(analysis_label=analysis_label,
                          config=config, inputs=inputs, destination=ses)
        analysis_ids.append(_id)
    except Exception as e:
        print(e)
        fails.append(ses)
