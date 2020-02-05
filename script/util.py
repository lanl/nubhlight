################################################################################
#                                                                              #
#  UTILITY FUNCTIONS                                                           #
#                                                                              #
################################################################################

import glob
import os

# COLORIZED OUTPUT
class color:
  BOLD    = '\033[1m'
  WARNING = '\033[1;31m'
  BLUE    = '\033[94m'
  NORMAL  = '\033[0m'
  YELLOW  = '\033[93m'

def get_files(PATH, NAME):                                                       
  return sorted(glob.glob(os.path.join(PATH,'') + NAME))

# PRINT ERROR MESSAGE
def warn(mesg):
  print(color.WARNING + "\n  ERROR: " + color.NORMAL + mesg + "\n")
def gentle_warn(mesg):
  print(color.YELLOW + "\n  WARNING: " + color.NORMAL + mesg + "\n")

# APPEND '/' TO PATH IF MISSING
def sanitize_path(path):
  return os.path.join(path, '')

# SEND OUTPUT TO LOG FILE AS WELL AS TERMINAL
def log_output(sys, logfile_name):
  import re
  f = open(logfile_name, 'w')
  class split(object):
    def __init__(self, *files):
      self.files = files
    def write(self, obj):
      n = 0
      ansi_escape = re.compile(r'\x1b[^m]*m')
      for f in self.files:
        if n > 0:
          f.write(ansi_escape.sub('', obj))
        else:
          f.write(obj)
        f.flush()
        n += 1
    def flush(self):
      for f in self.files:
        f.flush()
  sys.stdout = split(sys.stdout, f)
  sys.stderr = split(sys.stderr, f)

# CREATE DIRECTORY
def make_dir(path):
  if not os.path.exists(path):
    os.makedirs(path)

# CALL MAKE CLEAN IN A DIRECTORY
def make_clean(path):
  import subprocess
  if os.path.isfile(os.path.join(path,'makefile')):
    currdir = os.getcwd()
    os.chdir(path)
    print("\n  CALLING MAKE CLEAN")
    # suppress errors but don't actually care.
    with open(os.devnull,'w') as devnull:
      subprocess.call(['make','clean'],
                      stdout=devnull,
                      stderr=devnull)
    os.chdir(currdir)
    return

# CALL rm -rf ON RELATIVE PATHS ONLY
def safe_remove(path):
  import sys
  from subprocess import call

  # if path doesn't exist, do nothing
  if not os.path.exists(path):
    return
  
  # ONLY ALLOW RELATIVE PATHS
  if path[0] in ['/', '~']:
    warn("DIRECTORY\n"
         + "\t " + path + "\n"
         + "\t IS NOT A RELATIVE PATH!"
         + " DANGER OF DATA LOSS!\n"
         + "\t If you really meant to do this, remove the path manally.")
    sys.exit()
  
  call(['rm', '-rf', path])
  return

def parm_is_active(PARMS, key):
  if key in PARMS:
    if PARMS[key]:
      return True
  return False

def change_cparm(key, val, cfname):
  buffer = ''
  for line in open(cfname, 'r'):
    if key in line:
      buffer += 'bhl.config.set_cparm(\'' + key + '\', ' + str(val) + ')\n'
    else:
      buffer += line

  cfile = open(cfname, 'w')
  for line in buffer:
    cfile.write(line)
  cfile.close()

def change_rparm(key, val, rfname):
  buffer = ''
  for line in open(rfname, 'r'):
    if key in line:
      buffer += line.split(key)[0] + key + ' = ' + str(val) + '\n'
    else:
      buffer += line

  rfile = open(rfname, 'w')
  for line in buffer:
    rfile.write(line)
  rfile.close()

def write_rparm(pf, name):                                                       
  if name not in RPARMS:                                                         
    warn(name + ' NOT SET')                                                 
    sys.exit()                                                                   
                                                                                 
  datatype = RPARMS[name].datatype                                               
  default = RPARMS[name].value                                                   
  if datatype == 'integer':                                                      
    if default == None:                                                          
      pf.write('[int] ' + name + ' = \n')                                        
    else:                                                                        
      pf.write('[int] ' + name + ' = %d\n' % default)                            
  elif datatype == 'double':                                                     
    if default == None:                                                          
      pf.write('[dbl] ' + name + ' = \n')                                        
    else:                                                                        
      pf.write('[dbl] ' + name + ' = %e\n' % default)                            
  elif datatype == 'string':                                                     
    if default == None:                                                          
      pf.write('[str] ' + name + ' = \n')                                        
    else:                                                                        
      pf.write('[str] ' + name + ' = %s\n' % default)                            
  else:                                                                          
    print(name)                                                                  
    util.warn("DATATYPE " + str(datatype) + " NOT RECOGNIZED")                   
    sys.exit()                                                                   
                                                                                 
  del RPARMS[name]
