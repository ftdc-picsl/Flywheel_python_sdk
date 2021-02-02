import flywheel
import pytz

fw = flywheel.Client()

proj = fw.lookup('pennftdcenter/HUP6') # Use fw.lookup to get the Flywheel project object based on its path.
must_fix = [] # Create an empty list

# For-loop to find all of the sessions that have this meaningless session label.
for mysess in proj.sessions():
  if mysess.label == "BRAIN RESEARCH^GROSSMAN":
    must_fix.append(mysess)     

# Another for-loop to replace the session label with the informative one.
#for mysess in must_fix:
  # Create a new session label by piecing together year/month/day/hour/minute info from the session timestamp.
#  tstamp = mysess.timestamp
#  mylabel = '{}{}{}-{}{}'.format(tstamp.year,f'{tstamp.month:02}',f'{tstamp.day:02}',f'{tstamp.hour:02}',
#    f'{tstamp.minute:02}')
  # Update the session label using the update() method, whose input is a dictionary of the fields to be changed and
  # their new values.
#  mysess.update({'label': mylabel}) 
  # Not sure if this last line is necessary, but it reloads the session object with updates.
#  mysess.reload()

# Update, 12/10/2020: must adjust from UTC to US/Eastern timezone to get correct session label.
# Another for-loop to replace the session label with the informative one.
fix_list = []
for mysess in proj.sessions():
  # Create a new session label by piecing together year/month/day/hour/minute info from the session timestamp.
  tstamp = mysess.timestamp.astimezone(pytz.timezone("US/Eastern"))
  mylabel = '{}{}{}-{}{}'.format(tstamp.year,f'{tstamp.month:02}',f'{tstamp.day:02}',f'{tstamp.hour:02}',
    f'{tstamp.minute:02}')
  if not mylabel == mysess.label:
    anzlist = mysess.analyses
    if anzlist is not None:
      anz = [a.label for a in anzlist]
    else:
      anz = None
    fix_list.append((mysess.id,anz))

fixdf = pd.DataFrame(fix_list, columns = ['id','analyses'])

for i in range(fixdf.shape[0]):
  mysess = fw.get(fixdf.id.iloc[i])
  tstamp = mysess.timestamp.astimezone(pytz.timezone("US/Eastern"))
  mylabel = '{}{}{}-{}{}'.format(tstamp.year,f'{tstamp.month:02}',f'{tstamp.day:02}',f'{tstamp.hour:02}',
    f'{tstamp.minute:02}')
  mysess.update({'label': mylabel}) 
  # Not sure if this last line is necessary, but it reloads the session object with updates.
  mysess.reload()
  print(mysess.label, mysess.timestamp)