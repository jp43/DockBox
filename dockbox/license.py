import sys
import os
import time
import subprocess

def wrap_command(cmd, prgm):

    exe = os.path.abspath(__file__)
    if exe[-1] == 'c':
        exe = exe[:-1] # get .py from .pyc

    cmd_s = cmd.split()
    exe_cmd = cmd_s[0]

    if prgm == 'moe':
        logfile = 'moebatch.log'
        # write eval command until license is found
        newcmd = """while true; do
  %(cmd)s &> %(logfile)s
  status=`python %(exe)s %(prgm)s %(logfile)s`
  if [ "$status" == "0" ]; then break; fi
  sleep 10s
done"""% locals()

    elif prgm == 'gold':
        logfile = 'gold.err'
        newcmd = """while true; do
  %(cmd)s > /dev/null
  status=`python %(exe)s %(prgm)s %(logfile)s`
  if [ "$status" == "0" ]; then break; fi
  sleep 10s
done"""% locals()

    elif prgm == 'schrodinger':
        if exe_cmd == 'glide':
            filename1 = cmd_s[-1]
        elif exe_cmd == 'prepwizard':
            filename1 = cmd_s[-2]
        elif exe_cmd == 'ifd':
            filename1 = cmd_s[-1]
        else:
            raise ValueError("Schrodinger's command %s not recognized!"%exe_cmd)
        splitext_0 = os.path.splitext(filename1)[0]
        suffix = os.path.basename(splitext_0)
        logfile = suffix + '.log'
        newcmd = """while true; do
  output=`%(cmd)s`
  jobid=`echo "$output" | sed -n -e 's/^.*JobId: //p'`
  status=`python %(exe)s %(prgm)s %(logfile)s $jobid`
  if [ "$status" == "0" ]; then break; fi
  sleep 10s
done"""% locals()

    return newcmd

def check_schrodinger_license(logfile, jobid):
    """Check if schrodinger exe had license issues, design to avoid retry every 60 sec"""

    status = 0
    is_job_done = False
    is_job_killed = False

    while True:
        # (A) check if the job is still running
        output = subprocess.check_output('jobcontrol -list', shell=True, executable='/bin/bash')
        if jobid in output:
            time.sleep(2) # sleep for 2 sec
        else:
            is_job_done = True # the job is done
        # (B) check if the job has license issues
        if not is_job_killed:
            with open(logfile) as logf:
                for line in logf:
                    if 'Licensed number of users already reached' in line:
                        output = subprocess.check_output('jobcontrol -killnooutput %s'%jobid, shell=True, executable='/bin/bash')
                        status = 1
                        is_job_killed = True
        if is_job_done:
            break
    return status

def check_moe_license(logfile):

    status = 0
    with open(logfile) as logf:
        for line in logf:
            if 'Licensed number of users already reached' in line:
                status = 1
    return status

def check_gold_license(logfile):

    status = 0
    if os.path.exists(logfile):
        with open(logfile) as logf:
            for line in logf:
                if 'Licensed number of users already reached' in line:
                    status = 1
    return status

def run(args):

    if len(args) < 2:
        raise ValueError("check_licence.py should have at least two arguments")

    # first argument should be the program name
    prgm = args[1]
    status = 0

    # second argument should be the log file where to look for warning/error messages
    logfile = args[2]

    if prgm == 'moe':
        status = check_moe_license(logfile)
    elif prgm == 'gold':
        status = check_gold_license(logfile)
    elif prgm == 'schrodinger':
        jobid = args[3] # for Schrodinger, an extra argument is expected (job ID)
        status = check_schrodinger_license(logfile, jobid)

    return status

if __name__ == '__main__':
    status = run(sys.argv)
    print(status)
