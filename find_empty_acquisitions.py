import flywheel

fw = flywheel.Client()

project = fw.lookup("pennftdcenter/HUP3TLegacyUIDFix")

sessions = project.sessions()

# print subject,session,numAcquisitions,firstAcqFileSize,uploadDate

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
    
    if (line['numAcq'] == 0 or line['firstAcqFileSize'] < 100):
        print('{subject},{session},{numAcq},{firstAcqFileSize},{uploadDate}'.format(**line))
