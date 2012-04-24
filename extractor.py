#!/usr/bin/python

import sys

i = 0
while True:
  line = sys.stdin.readline()
  if not line:
	break
  i = i + 1
  if i % 10 == 0:
	line = line.rstrip()
	cols = line.split()
	print ' '.join(cols[1:3])

