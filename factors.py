import numpy as np 
import sys, getopt
##  input parameters 
### use the arguments given by the user ##########################################
args=sys.argv[1:]   # list of arguments
try:
   opts, args = getopt.getopt(args,"h",["number="])#, "m=","incl_angle="])

except getopt.GetoptError, exc:                 # display right format to parse arguments
   print exc.msg
   sys.exit(2)  # exit if the format of arguments is not right
print opts
num= 0

#print args 
for opt, arg in opts:
        if opt in ('--number'):
                num = float(arg)


def factors(n):    
    return set(reduce(list.__add__, 
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))

print factors(num)
