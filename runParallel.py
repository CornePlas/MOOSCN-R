import multiprocessing
from subprocess import call
import sys

SCRIPT=sys.argv[1]

def callScript(dim):
    print "Starting with dimensions ", dim
    call("R --vanilla --slave --args " + str(dim) + " < " + SCRIPT, shell=True)

## Use 8 processes not to hog up all resources at the server, and
## for not running out of memory
count = multiprocessing.cpu_count()
pool = multiprocessing.Pool(processes=count)
tasks = range(1, 21) # dimensions in [3,20]
print "Making computational tests for dimensions: ", tasks
r = pool.map_async(callScript, tasks)
r.wait() # Wait on the results


    
