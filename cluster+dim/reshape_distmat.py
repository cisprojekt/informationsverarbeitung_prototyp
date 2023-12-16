#!/usr/bin/env python3

import random
from matplotlib import pyplot as plt
import sys, re
import numpy as np
import math

x = []
y = []




inputfile = "testfile_d.csv"
try:
  stream = open(inputfile)
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)


re_se_m = '(-*\d+.\d*),(-*\d+.\d*)'
n = 0
  
for line in stream:
  re_se = re.search(r'{}'.format(re_se_m), line)
  if re_se:
    x.append(float(re_se.group(1)))
    y.append(float(re_se.group(2)))
    n += 1

stream.close

outputfile2 = "distmat.csv"
try:
  stream = open(outputfile2, "w")
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)
for i in range(1,len(x)):
  for j in range(0,i):
    stream.write("{:.4f}\n".format(math.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)))
    #d = sqrt(x[]*x[]+y[]*y[])
stream.close

outputfile3 = "distmat_inv.csv"
try:
  stream = open(outputfile3, "w")
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)
for i in range(len(x)):
  for j in range(i+1,len(x)):
    stream.write("{:.4f}\n".format(math.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)))
    #d = sqrt(x[]*x[]+y[]*y[])
stream.close

outputfile4 = "distmat_full.csv"
try:
  stream = open(outputfile4, "w")
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)
for i in range(len(x)):
  for j in range(len(x)):
    stream.write("{:.4f}\n".format(math.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)))
    #d = sqrt(x[]*x[]+y[]*y[])
stream.close

'''
for a in range(len(x)):
  print('{}  {}'.format(x[a],y[a]))
'''
fig, ax = plt.subplots()
ax.scatter(x,y)
fig.savefig('orig_data.pdf')
