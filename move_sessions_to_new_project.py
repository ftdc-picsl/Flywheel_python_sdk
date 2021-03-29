import flywheel
import csv
import sys


if (len(sys.argv) == 1):
    usage = '''
Usage: {} [srcProj] [destProj] sessions.csv

  srcProj and destProj full path to project, eg pennftdcenter/srcProj

  sessions.csv should be CSV file with header Subject.Label,Session.Label

'''

    print(usage.format(sys.argv[0]))

    sys.exit(1)    
    

fw = flywheel.Client()

src_proj_path = sys.argv[1]
dest_proj_path = sys.argv[2]

src_proj = fw.lookup(src_proj_path)
dest_proj = fw.lookup(dest_proj_path)

with open(sys.argv[3], newline='') as csvfile:

    sess_reader = csv.reader(csvfile, delimiter=',', quotechar='"')

    line_count = 0

    for row in sess_reader:
        line_count += 1
        if line_count == 1:
            print(f'Column names are {", ".join(row)}')
            continue
        src_subj_label = row[0]
        src_sess_label = row[1]

        # First check that the session exists
        src_subj = None
        src_sess = None

        try:
            src_subj = fw.lookup(f'{src_proj_path}/{src_subj_label}')
            src_sess = fw.lookup(f'{src_proj_path}/{src_subj_label}/{src_sess_label}')
        except flywheel.rest.ApiException:
            print(f'Cannot find session {src_proj_path}/{src_subj_label}/{src_sess_label}')
            continue

        # Make sure session doesn't already exist in dest_proj
        try:
            sess = fw.lookup(f'{dest_proj_path}/{src_subj_label}/{src_sess_label}')
            print(f'Session {src_sess_label} already exists in {dest_proj_path}')
            continue
        except flywheel.rest.ApiException:
            pass

        # subject to move the session to in the dest project
        dest_subj = None
 
        try:
            dest_subj = fw.lookup(f'{dest_proj_path}/{src_subj_label}')
        except flywheel.rest.ApiException:
            # Get subject info from src
            dest_subj = dest_proj.add_subject(label=src_subj_label, sex=src_subj.sex)

        # Finally, move session by changing its subject to the dest_subj
        src_sess.update(subject={'_id':dest_subj.id})
        
            

