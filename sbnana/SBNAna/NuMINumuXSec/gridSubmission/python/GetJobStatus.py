import os

def GetJobStatusByUser(user):
  os.system('jobsub_q -G icarus --user %s &> /tmp/jobstatus_%s.log'%(user,user))
  return open('/tmp/jobstatus_%s.log'%(user)).readlines()

def GetJobStatusByJobID(lines_JobStatus, this_jobID):
  for line in lines_JobStatus:

    # '33667526.0@jobsub03.fnal.gov          jskim           06/13 19:53   0+00:30:20 H   0   0.0 test.sh_20220613_195324_2194810_0_1_wrap.sh \n'
    words = line.split()
    '''
    0: "33667526.0@jobsub03.fnal.gov"
    1: "jskim"
    2: "06/13"
    3: "19:53"
    4: "0+00:30:20"
    5: "H"
    '''
         
    if len(words)==0:
      break

    if words[0]==this_jobID:
      jobFlag = words[5]
      if jobFlag=="H":
        return "ERROR"
      elif jobFlag=="R":
        return "RUNNING::%s"%(words[4])
      elif jobFlag=="I":
        return "IDLE"
      else:
        return jobFlag
     
  return "FINISHED"
