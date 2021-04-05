import flywheel

fw = flywheel.Client()

group = fw.lookup("pennftdcenter")

print('Project,NumSessions,DateCreated')

for proj in group.projects():

    line = {'project': '', 'numSessions': '', 'dateCreated': 'NA'}

    line['project'] = proj.label

    line['numSessions'] = len(proj.sessions())
 
    d = proj.created

    line['dateCreated'] = '{}-{:02}-{:02}'.format(d.year, d.month, d.day)
    
    print('{project},{numSessions},{dateCreated}'.format(**line))    
