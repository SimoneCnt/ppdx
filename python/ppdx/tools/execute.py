#!/usr/bin/env python3

import os

#import subprocess
#import shlex
#import logging
#log = logging.getLogger(__name__)

#def execute(cmd):
#    """
#        Execute a command on the terminal.
#    """
#    p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=-1)
#    stdout, stderr = p.communicate()
#    stdout = stdout.decode()
#    stderr = stderr.decode()
#    ret = p.returncode
#    if ret!=0:
#        print("WARN! <%s> returned %d" % (cmd, ret))
#        print("STDOUT:")
#        print(stdout)
#        print("STDERR:")
#        print(stderr)
#    return stdout, stderr, ret

def execute(cmd):
    ret = os.system(cmd)
    return ret

