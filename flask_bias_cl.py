import numpy as np
import sys

if(len(sys.argv)<5):
  print("syntax:",sys.argv[0],"[FLASK info file] [FLASK Cl input prefix] [biases, one for each galaxy field] [FLASK Cl output prefix]")
  sys.exit(1)

t=np.genfromtxt(sys.argv[1], dtype=(int,int,float,float,int,float,float))

print("t: {}".format(t))

#biases_string = sys.argv[3:-1][0]
#bias_strings = biases_string.split(",")
#
##print("")
##print("")
##print(sys.argv[2:-1])
##print("")
##bias=[float(s) for s in sys.argv[2:-1]]
#bias=[float(s) for s in bias_strings]
#
#
#
#
#
##bias=[float(s) for s in sys.argv[3:-1]]
#
#bi=0
#bai=0
#ba=np.ones(len(t))
#
#for f in t:
#  if f[4]==1: # galaxy field
#   print(f)
#   ba[bai]=bias[bi]
#   bi+=1
#  bai+=1
#
#if(bi!=len(bias)):
#  print("incorrect number of bias parameters given for this set of fields. exiting")
#  sys.exit(1)
#
#print(ba)
#
#for f1,b1 in zip(t,ba):
# for f2,b2 in zip(t,ba):
for f1 in t:
 for f2 in t:
#  if f1[0]==10:
#   f1[0]=0
#  if f2[0]==10:
#   f2[0]=0
  if(f1[0]>f2[0] or (f1[0]==f2[0] and f1[1]>f2[1])):
#   if f1[0]==0:
#    f1[0]=10
#   if f2[0]==0:
#    f2[0]=10
   continue
#  if f1[0]==0:
#   f1[0]=10
#  if f2[0]==0:
#   f2[0]=10
  cl=np.genfromtxt(sys.argv[2]+"f"+str(f1[0])+"z"+str(f1[1])+"f"+str(f2[0])+"z"+str(f2[1])+".dat") 
#  print("multiplying Cl",f1[0],f1[1]," ",f2[0],f2[1],"by bias",b1,"*",b2)
#  cl[:,1] *= b1*b2
  np.savetxt(sys.argv[-1]+"f"+str(f1[0])+"z"+str(f1[1])+"f"+str(f2[0])+"z"+str(f2[1])+".dat",cl,fmt='%i %.11e')
