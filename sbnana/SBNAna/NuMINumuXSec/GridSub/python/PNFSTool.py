import os
import Logger

class PNFSTool:

  def __init__(self, CreationDir):

    self.CreationDir = CreationDir

  def Open(self, fname):

    return open('%s/%s'%(self.CreationDir, fname), 'w')

  def Close(self, f, To):

    f.close()
    os.system('ifdh cp %s %s >> /dev/null 2>&1'%(f.name, To))
    os.system('rm %s'%(f.name))
