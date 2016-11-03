#!/u/ki/dgruen/anaconda3/bin/python
import numpy
import sys

if(len(sys.argv)!=4):
  print("syntax:",sys.argv[0],"[FLASK info infile] [bias] [outfile]")
  sys.exit(1)

t=numpy.genfromtxt(sys.argv[1], dtype=(int,int,float,float,int,float,float))
bias=float(sys.argv[2])

for f in t:
  if f[4]==1:
    f[3] *= bias

numpy.savetxt(sys.argv[3],t,fmt='%i    %i   %f   %f   %i   %f    %f   ')
