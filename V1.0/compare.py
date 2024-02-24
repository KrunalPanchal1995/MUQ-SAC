import numpy as np

g = open("locations","r").readlines()
f = open("progress","r").readlines()

for i in g:
	if i not in f:
		print(i)

