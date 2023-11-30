#!/usr/bin/env python3

import random

random.seed()
for i in range(1000):
  x = random.random()*100
  y = random.random()*100
  print("{:.4f},{:.4f}".format(x,y))
