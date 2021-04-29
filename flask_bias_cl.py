import numpy as np
import sys
import copy

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
try:
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
   print("")
   print("")
   print("start")
   print("f1:")
   print(type(f1))
   print(f1)
   print("f2:")
   print(type(f2))
   print(f2)
   print("t:")
   print(type(t))
   print(t)
   print(sys.argv[2]+"f"+str(f1[0])+"z"+str(f1[1])+"f"+str(f2[0])+"z"+str(f2[1])+".dat")
   print("end")
   print("")
   print("")
   cl=np.genfromtxt(sys.argv[2]+"f"+str(f1[0])+"z"+str(f1[1])+"f"+str(f2[0])+"z"+str(f2[1])+".dat") 
#  print("multiplying Cl",f1[0],f1[1]," ",f2[0],f2[1],"by bias",b1,"*",b2)
#  cl[:,1] *= b1*b2
   np.savetxt(sys.argv[-1]+"f"+str(f1[0])+"z"+str(f1[1])+"f"+str(f2[0])+"z"+str(f2[1])+".dat",cl,fmt='%i %.11e')
except:
 for f1 in [t]:
  for f2 in [t]:
   f1_copy = str(copy.deepcopy(f1)).split(", ")
   f2_copy = str(copy.deepcopy(f2)).split(", ")
   f1 = [f1_copy[0].split("(")[1],f1_copy[1]]
   f2 = [f2_copy[0].split("(")[1],f2_copy[1]]
#  if f1[0]==10:
#   f1[0]=0
#  if f2[0]==10:
#   f2[0]=0
#   if f1[0]==0:
#    f1[0]=10
#   if f2[0]==0:
#    f2[0]=10
#  if f1[0]==0:
#   f1[0]=10
#  if f2[0]==0:
#   f2[0]=10
   print("")
   print("")
   print("start")
   print("f1:")
   print(type(f1))
   print(f1)
   print("f1[0]:")
   print(f1[0])
   print("f2[0]:")
   print(f2[0])
   print("f2:")
   print(type(f2))
   print(f2)
   print("t:")
   print(type(t))
   print(t)
   print(sys.argv[2]+"f"+str(f1[0])+"z"+str(f1[1])+"f"+str(f2[0])+"z"+str(f2[1])+".dat")
   print("end")
   print("")
   print("")

   cl=np.genfromtxt(sys.argv[2]+"f"+str(f1[0])+"z"+str(f1[1])+"f"+str(f2[0])+"z"+str(f2[1])+".dat")
#  print("multiplying Cl",f1[0],f1[1]," ",f2[0],f2[1],"by bias",b1,"*",b2)
#  cl[:,1] *= b1*b2

   print(sys.argv[-1]+"f"+str(f1[0])+"z"+str(f1[1])+"f"+str(f2[0])+"z"+str(f2[1])+".dat")
   print(cl)
   np.savetxt(sys.argv[-1]+"f"+str(f1[0])+"z"+str(f1[1])+"f"+str(f2[0])+"z"+str(f2[1])+".dat",cl,fmt='%i %.11e')
   #f1_copy = str(copy.deepcopy(t)).split(", ")
   #f1 = [f1_copy[0].split("(")[1], f1_copy[1], f1_copy[2], f1_copy[3], f1_copy[4], f1_copy[5], f1_copy[6]]
