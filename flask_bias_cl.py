#!/u/ki/dgruen/anaconda3/bin/python
import numpy
import sys

if(len(sys.argv)!=5):
  print("syntax:",sys.argv[0],"[FLASK info file] [FLASK Cl input prefix] [bias] [FLASK Cl output prefix]")
  sys.exit(1)

t=numpy.genfromtxt(sys.argv[1], dtype=(int,int,float,float,int,float,float))
bias=float(sys.argv[3])

for f1 in t:
 for f2 in t:
  if(f1[0]>f2[0] or (f1[0]==f2[0] and f1[1]>f2[1])):
   continue
  cl=numpy.genfromtxt(sys.argv[2]+"f"+str(f1[0])+"z"+str(f1[1])+"f"+str(f2[0])+"z"+str(f2[1])+".dat") 
  if f1[4]==1:
    cl[:,1] *= bias
  if f2[4]==1:
    cl[:,1] *= bias
  numpy.savetxt(sys.argv[4]+"f"+str(f1[0])+"z"+str(f1[1])+"f"+str(f2[0])+"z"+str(f2[1])+".dat",cl,fmt='%i %.11e')
