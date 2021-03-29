import flywheel

fw = flywheel.Client()

project = fw.lookup("pennftdcenter/HUP6")

sessions = project.sessions()

# Checks the file size of the first file of the first acquisition
# Doesn't specifically check for dicom files but any sensible zip, nii, or json should be
# bigger than 100 bytes
#

print('SubjectLabel,SessionLabel,Acquisitions,FirstAcqFileSize,FirstAcqOrSessionUploadDate')

for sess in sessions:

    line = {'subject': '', 'session': '', 'numAcq': 0, 'firstAcqFileSize': 0, 'uploadDate': 'NA'}

    acquisitions = sess.acquisitions()
    line['subject'] = sess.subject.label
    line['session'] = sess.label
    line['numAcq'] = len(acquisitions)

    d = sess.created

    line['uploadDate'] = '{}{:02}{:02}'.format(d.year, d.month, d.day)
    
    if (len(acquisitions) > 0):
      if (len(acquisitions[0].files) > 0):
          line['firstAcqFileSize'] = acquisitions[0].files[0].size
      else:
          line['firstAcqFileSize'] = 0
      d = acquisitions[0].created
      line['uploadDate'] = '{}{:02}{:02}'.format(d.year, d.month, d.day) 
    
    if (line['numAcq'] == 0 or int(line['firstAcqFileSize']) < 100):
        print('{subject},{session},{numAcq},{firstAcqFileSize},{uploadDate}'.format(**line))
