import numpy
import sys

if(len(sys.argv)<4):
  print("syntax:",sys.argv[0],"[FLASK info infile] [bias, one for each galaxy/matter field] [outfile]")
  sys.exit(1)

t=numpy.genfromtxt(sys.argv[1], dtype=(int,int,float,float,int,float,float))
print("t: {}".format(t))

#print(sys.argv[1])
#biases_string = sys.argv[2:-1][0]
#bias_strings = biases_string.split(",")
#
##print("")
##print("")
##print(sys.argv[2:-1])
##print("")
##bias=[float(s) for s in sys.argv[2:-1]]
#bias=[float(s) for s in bias_strings]
##bias = []
##for s in sys.argv[2:-1]:
##  print(s)
##  print(float(s))
##  bias.append(float(s))
#
#bi=0
#
#for f in t:
#  if f[4]==1: # galaxy field
#   print("multiplying by bias",bias[bi],bi)
#   print(f)
#   f[3]*=bias[bi]
#   bi+=1
#   print(f)

numpy.savetxt(sys.argv[-1],t,fmt='%i    %i   %f   %f   %i   %f    %f   ')
