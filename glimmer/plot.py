#!/usr/bin/env python3

import random
from matplotlib import pyplot as plt
import sys, re
import numpy as np

inputfile = sys.argv[1]
try:
  stream = open(inputfile)
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)

x = []
y = []
re_se_m = '(-*\d+.\d*),(-*\d+.\d*)'
n = 0
  
for line in stream:
  re_se = re.search(r'{}'.format(re_se_m), line)
  if re_se:
    x.append(float(re_se.group(1)))
    y.append(float(re_se.group(2)))
    n += 1

stream.close

for a in range(n):
  print('{}  {}'.format(x[a],y[a]))
fig, ax = plt.subplots()
ax.scatter(x,y)
fig.savefig('plot.pdf')
