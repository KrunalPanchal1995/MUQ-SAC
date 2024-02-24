import numpy as np
import time
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib import style

style.use("fivethirtyeight")

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)

def animate(i):
	graph_data = open("samplefile.txt","r").readlines()
	xs = []
	ys = []
	for line in graph_data[-20:]:
		if len(line)>1:
			x, y = line.split(",")
			xs.append(x)
			ys.append(y)
			
	ax1.clear()
	ax1.plot(xs,ys)
	ax1.set_xlabel("Iteration number")
	ax1.set_ylabel("Objective function")
	
ani = animation.FuncAnimation(fig,animate, interval = 1000) #1000 ms the graph will update

plt.show()
