#!/usr/bin/env python
import healpy as hp
import numpy as np
from math import pi
import sys

if(len(sys.argv)<3):
  print("syntax:",sys.argv[0]," [input maps...] [output map]")
  sys.exit(0)

out_map = hp.read_map(sys.argv[1])

for f in sys.argv[2:-1]:
  add_map = hp.read_map(f)
  out_map += add_map
  add_map=0

hp.write_map(sys.argv[-1],out_map)

