#!/usr/bin/env python3

import random
from matplotlib import pyplot as plt
import sys, re
import numpy as np
import math

dist = np.array([])

inputfile = "distmat_it_1.csv"
try:
  stream = open(inputfile)
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)


re_se_m = '(-*\d+.\d*)'
n = 0
  
for line in stream:
  re_se = re.search(r'{}'.format(re_se_m), line)
  if re_se:
    a = np.array([(float)(re_se.group(1))])
    dist = np.append(dist, a)
    n += 1

stream.close
npoints = (int) (0.5 + np.sqrt(0.25 + 2*n))

dist_new = np.zeros((npoints, npoints))

z = 0
for i in range(npoints-1):
  for j in range(i+1, npoints):
    dist_new[i][j] = dist[z]
    z += 1
    dist_new[j][i] = dist_new[i][j]
for k in range(npoints):
  dist_new[npoints-1][k] = dist_new[k][npoints-1]

outputfile2 = "distmat_full_it_1.csv"
try:
  stream = open(outputfile2, "w")
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)
for i in range(npoints):
  for j in range(npoints):
    stream.write("{:.4f}\n".format(dist_new[i][j]))
    #d = sqrt(x[]*x[]+y[]*y[])
stream.close
'''
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
'''
'''
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
'''
for a in range(len(x)):
  print('{}  {}'.format(x[a],y[a]))
'''
'''
fig, ax = plt.subplots()
ax.scatter(x,y)
fig.savefig('orig_data.pdf')
'''

dist = np.array([])


inputfile = "distmat_it_2.csv"
try:
  stream = open(inputfile)
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)


re_se_m = '(-*\d+.\d*)'
n = 0
  
for line in stream:
  re_se = re.search(r'{}'.format(re_se_m), line)
  if re_se:
    a = np.array([(float)(re_se.group(1))])
    dist = np.append(dist, a)
    n += 1

stream.close
npoints = (int) (0.5 + np.sqrt(0.25 + 2*n))

dist_new = np.zeros((npoints, npoints))

z = 0
for i in range(npoints-1):
  for j in range(i+1, npoints):
    dist_new[i][j] = dist[z]
    z += 1
    dist_new[j][i] = dist_new[i][j]
for k in range(npoints):
  dist_new[npoints-1][k] = dist_new[k][npoints-1]

outputfile2 = "distmat_full_it_2.csv"
try:
  stream = open(outputfile2, "w")
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)
for i in range(npoints):
  for j in range(npoints):
    stream.write("{:.4f}\n".format(dist_new[i][j]))
    #d = sqrt(x[]*x[]+y[]*y[])
stream.close

dist = np.array([])


inputfile = "distmat_it_3.csv"
try:
  stream = open(inputfile)
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)


re_se_m = '(-*\d+.\d*)'
n = 0
  
for line in stream:
  re_se = re.search(r'{}'.format(re_se_m), line)
  if re_se:
    a = np.array([(float)(re_se.group(1))])
    dist = np.append(dist, a)
    n += 1

stream.close
npoints = (int) (0.5 + np.sqrt(0.25 + 2*n))

dist_new = np.zeros((npoints, npoints))

z = 0
for i in range(npoints-1):
  for j in range(i+1, npoints):
    dist_new[i][j] = dist[z]
    z += 1
    dist_new[j][i] = dist_new[i][j]
for k in range(npoints):
  dist_new[npoints-1][k] = dist_new[k][npoints-1]

outputfile2 = "distmat_full_it_3.csv"
try:
  stream = open(outputfile2, "w")
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)
for i in range(npoints):
  for j in range(npoints):
    stream.write("{:.4f}\n".format(dist_new[i][j]))
    #d = sqrt(x[]*x[]+y[]*y[])
stream.close

dist = np.array([])


inputfile = "distmat_it_4.csv"
try:
  stream = open(inputfile)
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)


re_se_m = '(-*\d+.\d*)'
n = 0
  
for line in stream:
  re_se = re.search(r'{}'.format(re_se_m), line)
  if re_se:
    a = np.array([(float)(re_se.group(1))])
    dist = np.append(dist, a)
    n += 1

stream.close
npoints = (int) (0.5 + np.sqrt(0.25 + 2*n))

dist_new = np.zeros((npoints, npoints))

z = 0
for i in range(npoints-1):
  for j in range(i+1, npoints):
    dist_new[i][j] = dist[z]
    z += 1
    dist_new[j][i] = dist_new[i][j]
for k in range(npoints):
  dist_new[npoints-1][k] = dist_new[k][npoints-1]

outputfile2 = "distmat_full_it_4.csv"
try:
  stream = open(outputfile2, "w")
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)
for i in range(npoints):
  for j in range(npoints):
    stream.write("{:.4f}\n".format(dist_new[i][j]))
    #d = sqrt(x[]*x[]+y[]*y[])
stream.close

dist = np.array([])


inputfile = "distmat_it_5.csv"
try:
  stream = open(inputfile)
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)


re_se_m = '(-*\d+.\d*)'
n = 0
  
for line in stream:
  re_se = re.search(r'{}'.format(re_se_m), line)
  if re_se:
    a = np.array([(float)(re_se.group(1))])
    dist = np.append(dist, a)
    n += 1

stream.close
npoints = (int) (0.5 + np.sqrt(0.25 + 2*n))

dist_new = np.zeros((npoints, npoints))

z = 0
for i in range(npoints-1):
  for j in range(i+1, npoints):
    dist_new[i][j] = dist[z]
    z += 1
    dist_new[j][i] = dist_new[i][j]
for k in range(npoints):
  dist_new[npoints-1][k] = dist_new[k][npoints-1]

outputfile2 = "distmat_full_it_5.csv"
try:
  stream = open(outputfile2, "w")
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)
for i in range(npoints):
  for j in range(npoints):
    stream.write("{:.4f}\n".format(dist_new[i][j]))
    #d = sqrt(x[]*x[]+y[]*y[])
stream.close

dist = np.array([])


inputfile = "distmat_it_6.csv"
try:
  stream = open(inputfile)
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)


re_se_m = '(-*\d+.\d*)'
n = 0
  
for line in stream:
  re_se = re.search(r'{}'.format(re_se_m), line)
  if re_se:
    a = np.array([(float)(re_se.group(1))])
    dist = np.append(dist, a)
    n += 1

stream.close
npoints = (int) (0.5 + np.sqrt(0.25 + 2*n))

dist_new = np.zeros((npoints, npoints))

z = 0
for i in range(npoints-1):
  for j in range(i+1, npoints):
    dist_new[i][j] = dist[z]
    z += 1
    dist_new[j][i] = dist_new[i][j]
for k in range(npoints):
  dist_new[npoints-1][k] = dist_new[k][npoints-1]

outputfile2 = "distmat_full_it_6.csv"
try:
  stream = open(outputfile2, "w")
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)
for i in range(npoints):
  for j in range(npoints):
    stream.write("{:.4f}\n".format(dist_new[i][j]))
    #d = sqrt(x[]*x[]+y[]*y[])
stream.close

dist = np.array([])


inputfile = "distmat_it_7.csv"
try:
  stream = open(inputfile)
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)


re_se_m = '(-*\d+.\d*)'
n = 0
  
for line in stream:
  re_se = re.search(r'{}'.format(re_se_m), line)
  if re_se:
    a = np.array([(float)(re_se.group(1))])
    dist = np.append(dist, a)
    n += 1

stream.close
npoints = (int) (0.5 + np.sqrt(0.25 + 2*n))

dist_new = np.zeros((npoints, npoints))

z = 0
for i in range(npoints-1):
  for j in range(i+1, npoints):
    dist_new[i][j] = dist[z]
    z += 1
    dist_new[j][i] = dist_new[i][j]
for k in range(npoints):
  dist_new[npoints-1][k] = dist_new[k][npoints-1]

outputfile2 = "distmat_full_it_7.csv"
try:
  stream = open(outputfile2, "w")
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)
for i in range(npoints):
  for j in range(npoints):
    stream.write("{:.4f}\n".format(dist_new[i][j]))
    #d = sqrt(x[]*x[]+y[]*y[])
stream.close

dist = np.array([])


inputfile = "distmat_it_8.csv"
try:
  stream = open(inputfile)
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)


re_se_m = '(-*\d+.\d*)'
n = 0
  
for line in stream:
  re_se = re.search(r'{}'.format(re_se_m), line)
  if re_se:
    a = np.array([(float)(re_se.group(1))])
    dist = np.append(dist, a)
    n += 1

stream.close
npoints = (int) (0.5 + np.sqrt(0.25 + 2*n))

dist_new = np.zeros((npoints, npoints))

z = 0
for i in range(npoints-1):
  for j in range(i+1, npoints):
    dist_new[i][j] = dist[z]
    z += 1
    dist_new[j][i] = dist_new[i][j]
for k in range(npoints):
  dist_new[npoints-1][k] = dist_new[k][npoints-1]

outputfile2 = "distmat_full_it_8.csv"
try:
  stream = open(outputfile2, "w")
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)
for i in range(npoints):
  for j in range(npoints):
    stream.write("{:.4f}\n".format(dist_new[i][j]))
    #d = sqrt(x[]*x[]+y[]*y[])
stream.close

dist = np.array([])


inputfile = "distmat_it_9.csv"
try:
  stream = open(inputfile)
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)


re_se_m = '(-*\d+.\d*)'
n = 0
  
for line in stream:
  re_se = re.search(r'{}'.format(re_se_m), line)
  if re_se:
    a = np.array([(float)(re_se.group(1))])
    dist = np.append(dist, a)
    n += 1

stream.close
npoints = (int) (0.5 + np.sqrt(0.25 + 2*n))

dist_new = np.zeros((npoints, npoints))

z = 0
for i in range(npoints-1):
  for j in range(i+1, npoints):
    dist_new[i][j] = dist[z]
    z += 1
    dist_new[j][i] = dist_new[i][j]
for k in range(npoints):
  dist_new[npoints-1][k] = dist_new[k][npoints-1]

outputfile2 = "distmat_full_it_9.csv"
try:
  stream = open(outputfile2, "w")
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)
for i in range(npoints):
  for j in range(npoints):
    stream.write("{:.4f}\n".format(dist_new[i][j]))
    #d = sqrt(x[]*x[]+y[]*y[])
stream.close

dist = np.array([])


inputfile = "distmat_it_10.csv"
try:
  stream = open(inputfile)
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)


re_se_m = '(-*\d+.\d*)'
n = 0
  
for line in stream:
  re_se = re.search(r'{}'.format(re_se_m), line)
  if re_se:
    a = np.array([(float)(re_se.group(1))])
    dist = np.append(dist, a)
    n += 1

stream.close
npoints = (int) (0.5 + np.sqrt(0.25 + 2*n))

dist_new = np.zeros((npoints, npoints))

z = 0
for i in range(npoints-1):
  for j in range(i+1, npoints):
    dist_new[i][j] = dist[z]
    z += 1
    dist_new[j][i] = dist_new[i][j]
for k in range(npoints):
  dist_new[npoints-1][k] = dist_new[k][npoints-1]

outputfile2 = "distmat_full_it_10.csv"
try:
  stream = open(outputfile2, "w")
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)
for i in range(npoints):
  for j in range(npoints):
    stream.write("{:.4f}\n".format(dist_new[i][j]))
    #d = sqrt(x[]*x[]+y[]*y[])
stream.close

