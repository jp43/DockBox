import sys
import os

def eval(cmd, prgm, redirect=None):

    exe = os.path.abspath(__file__)
    if exe[-1] == 'c':
        exe = exe[:-1] # get .py from .pyc
    tmpfile = 'tmp.out'

    if redirect:
        redstr = '&> ' + redirect
    else:
        redstr = '' 

    # write eval command until license is found
    newcmd = """while true; do
  %(cmd)s &> %(tmpfile)s
  status=`python %(exe)s %(prgm)s`
  if [ "$status" == "0" ]; then cat %(tmpfile)s %(redstr)s; rm -rf %(tmpfile)s; break; fi
  sleep 120s
done"""% locals()

    return newcmd

def run(args):
    prgm = args[1]
    with open('tmp.out') as logf:
        status = 0
        for line in logf:
            if prgm == 'schrodinger':
                if 'FATAL' in line and 'license key' in line:
                    status = 1
                    break
            elif prgm == 'moe':
                if 'Licensed number of users already reached' in line:
                    status = 1
                    break
            else:
                ValueError('program %s has no license to check'%prgm)
    print status

if __name__ == '__main__':
    run(sys.argv)
