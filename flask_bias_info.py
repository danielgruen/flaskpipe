import numpy
import sys

if(len(sys.argv)<4):
  print("syntax:",sys.argv[0],"[FLASK info infile] [bias, one for each galaxy/matter field] [outfile]")
  sys.exit(1)

t=numpy.genfromtxt(sys.argv[1], dtype=(int,int,float,float,int,float,float))

bias=[float(s) for s in sys.argv[2:-1]]

bi=0

for f in t:
  if f[4]==1: # galaxy field
   print("multiplying by bias",bias[bi],bi)
   print(f)
   f[3]*=bias[bi]
   bi+=1
   print(f)

numpy.savetxt(sys.argv[-1],t,fmt='%i    %i   %f   %f   %i   %f    %f   ')
