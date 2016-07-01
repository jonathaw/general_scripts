#!/usr/bin/env python3.5
import os
import sys
import subprocess
pwd = os.getcwd().split('/')
jonathaw = pwd.index('jonathaw')
correct_path = '/Volumes/labs/fleishman/jonathaw/' + '/'.join(pwd[jonathaw+1:]) + '/' + sys.argv[1]
print(correct_path)
cmd = 'ssh -X jonathaw@132.77.79.159 "/usr/local/bin/atom %s"' % correct_path
print(cmd)
subprocess.call(cmd)