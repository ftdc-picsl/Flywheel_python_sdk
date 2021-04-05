import flywheel

fw = flywheel.Client()

group = fw.lookup("pennftdcenter")

print('Project,NumSessions,DateCreated,AccountNumber')

for proj in group.projects():

    proj = proj.reload()
    
    line = {'project': '', 'numSessions': '', 'dateCreated': 'NA', 'accountNumber': 'None'}

    line['project'] = proj.label

    line['numSessions'] = len(proj.sessions())
 
    d = proj.created

    line['accountNumber'] = proj.info['ProjectFunding']['accountNumber']

    line['dateCreated'] = '{}-{:02}-{:02}'.format(d.year, d.month, d.day)
    
    print('{project},{numSessions},{dateCreated},{accountNumber}'.format(**line))    
