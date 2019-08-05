#!/usr/bin/env python3

import subprocess

def execute(cmd):
    """
        Execute a command on the terminal.
    """
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    stdout = stdout.decode()
    stderr = stderr.decode()
    ret = p.wait()
    if ret!=0:
        print("WARN! <%s> returned %d" % (cmd, ret))
        print("STDOUT:")
        print(stdout)
        print("STDERR:")
        print(stderr)
    return stdout, stderr, ret

