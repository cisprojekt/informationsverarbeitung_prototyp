#!/usr/bin/env python3

import random
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys, re
import numpy as np


def save_image(filename): 
    
    p = PdfPages(filename) 
      
    fig_nums = plt.get_fignums()   
    figs = [plt.figure(n) for n in fig_nums] 
      
    for fig in figs:  
        
        fig.savefig(p, format='pdf')  
          
    p.close()   


plt.rcParams["figure.figsize"] = [7.00, 3.50] 
plt.rcParams["figure.autolayout"] = True

inputfile = "result_zoom_1.csv"
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
fig1 = plt.figure()
plt.plot(x,y,'bo')
#fig, ax = plt.subplots()
#ax.scatter(x,y)
#fig.savefig('plot_zoom_1.pdf')

inputfile = "result_zoom_2.csv"
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
fig2 = plt.figure()
plt.plot(x,y,'bo')
#fig.savefig('plot_zoom_2.pdf')

inputfile = "result_zoom_3.csv"
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
fig3 = plt.figure()
plt.plot(x,y,'bo')
#fig.savefig('plot_zoom_3.pdf')

inputfile = "result_zoom_4.csv"
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
fig4 = plt.figure()
plt.plot(x,y,'bo')
#fig.savefig('plot_zoom_4.pdf')

inputfile = "result_zoom_5.csv"
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
fig5 = plt.figure()
plt.plot(x,y,'bo')
#fig.savefig('plot_zoom_5.pdf')

inputfile = "result_zoom_6.csv"
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
fig6 = plt.figure()
plt.plot(x,y,'bo')
#fig.savefig('plot_zoom_6.pdf')

inputfile = "result_zoom_7.csv"
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
fig7 = plt.figure()
plt.plot(x,y,'bo')
#fig.savefig('plot_zoom_7.pdf')

inputfile = "result_zoom_8.csv"
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
fig8 = plt.figure()
plt.plot(x,y,'bo')
#fig.savefig('plot_zoom_8.pdf')

inputfile = "result_zoom_9.csv"
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
fig9 = plt.figure()
plt.plot(x,y,'bo')
#fig.savefig('plot_zoom_9.pdf')

inputfile = "result_zoom_10.csv"
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
fig10 = plt.figure()
plt.plot(x,y,'bo')
#fig.savefig('plot_zoom_10.pdf')

filename = "plot.pdf"
save_image(filename)

