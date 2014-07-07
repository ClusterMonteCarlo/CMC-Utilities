def parseStoryDynamics(filename):
	data = []
	i = -1
	count = 0
	
	file = open(filename)
	line = file.readline() #can't use 'with open as', the iterator breaks for some reason
	
	while line:
		if line.split()[0] == 'system_time': #grab shit from the root node
			t = float(line.split()[2])
			data.append([t])
			i += 1
			while line.split()[0] != ')Dynamics':
				line = file.readline() #advance forward to end up dynamics
				if line.split()[0] == 'com_pos':
					xCOM,yCOM,zCOM = stringsToFloats(line.split()[2:])
				elif line.split()[0] == 'com_vel': 
					vxCOM,vyCOM,vzCOM = stringsToFloats(line.split()[2:])
		elif line.split()[0] == 'i':
			idNum = int(line.split()[2])
		elif line.split()[0] == 'N':
			if line.split()[2] == '2':
				binFlag = 10
				count = 2
			elif line.split()[2] == '3':
				binFlag = 50
				count = 3
			elif line.split()[2] == '4':
				binFlag = 90
				count = 4
			else:
				binFlag = 1
		elif line.split()[0] == 'm':
			m = float(line.split()[2])
		elif line.split()[0] == 'r':
			x,y,z = line.split()[2:]
			x = float(x) + xCOM
			y = float(y) + yCOM
			z = float(z) + zCOM
		elif line.split()[0] == 'v':
			vx,vy,vz = line.split()[2:]
			vx = float(vx) + vxCOM
			vy = float(vy) + vyCOM
			vz = float(vz) + vzCOM
		elif line.split()[0] == ')Dynamics':
			data[i].append([x,y,z,vx,vy,vz,binFlag,m,idNum])
			line = file.readline()
			while count != 0:
				if line.split()[0] == ')Dynamics':
					count -= 1
				line = file.readline()
		line = file.readline()
		
	file.close()
	return data	

def stringsToFloats(stringArray):
	for i in range(len(stringArray)): #probably a more pythonic way, but whatever
		stringArray[i] = float(stringArray[i])
	return stringArray

def parseStoryTotals(filename):
	data = []
	i = 0
	count = 0
	eKin = 0
	ePot = 0
	eDyn = 0
	from numpy import zeros
	eTot,ePot,eKin,virial,numSingle,numBinary,numTriple,numQuad = zeros(8) 
	lagRad = 0

	file = open(filename)
	line = file.readline() #can't use 'with open as', the iterator breaks for some reason
	
	while line.split()[0] != 'system_time':
		line = file.readline()
	t = float(line.split()[2])
	if t == 0.: #need to skip the initial print, if present
		line = file.readline()
		while line.split()[0] != 'system_time':
			line = file.readline()
		t = float(line.split()[2])
	data.append([t])
	line = file.readline()
	scale = 0
	numSingle = 0
	numBinary = 0
	numTriple = 0
	numQuad= 0
	
	while line: #C SYNTAX 4EVAR!!!
		if line.strip() == 'N = 1':
			numSingle += 1
		elif line.strip() == 'N = 2':
			numBinary += 1
		elif line.strip() == 'N = 3':
			numTriple += 1
		elif line.strip() == 'N = 4':
			numQuad += 1
		elif line.split()[0] == 'total_energy':
			eTot = float(line.split()[2])
		elif line.split()[0] == 'potential_energy':
			ePot = float(line.split()[2])
		elif line.split()[0] == 'kinetic_energy':
			eKin = float(line.split()[2])
		elif line.split()[0]  =='new_plummer_scale':
			scale = float(line.split()[2])
		elif line.split()[0]  =='friction_energy_lost':
			eDyn = float(line.split()[2])
		elif line.split()[0] == 'r_lagr':
			lagRad = line.split()[2:]
			lagRad = stringsToFloats(lagRad)
		elif line.split()[0] == 'system_time':
			virial = abs(2*eKin/(ePot))
			data[i].append([eTot,ePot,eKin,virial,eDyn,scale,numSingle,numBinary,numTriple,numQuad])
			data[i].append(lagRad)
			i += 1
			t = float(line.split()[2])
			data.append([t])
			numSingle = 0
			numBinary = 0
			numTriple = 0
			numQuad = 0
		line = file.readline()
		
	file.close()
	del data[-1]

	return data 

def plotNTotals(data):
	import matplotlib.pyplot as plt
	from numpy import array
	plt.ion()
	plt.figure()
	t = array([step[0] for step in data])
	totals = array([step[1][5:] for step in data])
	totalsTranspose = totals.T
	for count in totalsTranspose:
		plt.plot(t,count)
	plt.grid(True)
	plt.legend((r'$N_{total}$',r'$N_{bin}$',
	            r'$N_{trip}$',r'$N_{quad}$'),loc=6)
	plt.xlabel('Time (crossing times)')
	plt.ylabel('Number')

def plotLagRad(data):
	import matplotlib.pyplot as plt
	from numpy import array
	plt.ion()
	plt.figure()
	t = array([step[0] for step in data])
	lagRad = array([step[2] for step in data])
	lagRadTranspose = lagRad.T
	for radius in lagRadTranspose:
		plt.plot(t,radius)
	plt.grid(True)
	plt.xlabel('Time (crossing times)')
	plt.ylabel('Lagrange Radii')
	
def plotEnergies(data):
	import matplotlib.pyplot as plt
	from numpy import array
	plt.ion()
	plt.figure()
	t = array([step[0] for step in data])
	energies = array([step[1][:4] for step in data])
	eTot, ePot, eKin, virial = energies.T
	plt.plot(t,eTot)
	plt.plot(t,ePot)
	plt.plot(t,eKin)
	plt.plot(t,virial)
	plt.grid(True)
	plt.legend((r'$E_{total}$',r'$E_{potential}$',
	            r'$E_{kinetic}$',r'$Virial =2E_k/E_p$'),loc=0)
	plt.xlabel('Time (crossing times)')
	plt.ylabel('Energy (N-body Units)')

def rToXY(r):
	from numpy import arccos,cos,sin,pi
	from numpy.random import uniform
	sinTheta = sin(arccos(uniform(-1,1)))
	phi = uniform(0,2*pi)
	x = r*sinTheta*cos(phi)
	y = r*sinTheta*sin(phi)
	return x,y

def rToXYZ(r):
	from numpy import arccos,cos,sin,pi
	from numpy.random import uniform
	angle = arccos(uniform(-1,1))
	sinTheta = sin(angle)
	phi = uniform(0,2*pi)
	x = r*sinTheta*cos(phi)
	y = r*sinTheta*sin(phi)
	z = r*cos(angle)
	return x,y,z


def make2dDynamicsMovie(data, filename="output.mp4", Min=-20, Max=20, size=8,
fps=30, noShow=True, streaks=True):
	"""
	Make a 2D movie of the dynamics from kira, using the data from parseStoryDynamics
	filename - duh
	Min/Max  - size of the plot range
	size	 - size of the saved figures
	fps	  - duh
	noShow   - don't show the movie plots as they're being made
	"""
	from matplotlib import animation, pyplot as plt
	from numpy import loadtxt

	#mANDr = loadtxt('64k-20-0.02-1.snap0015.dat',usecols=(1,2))
	#points = [rToXY(p[1]) for p in mANDr]
	#x = [p[0] for p in points]
	#y = [p[1] for p in points]
	mThresh = 0.0002

	def animate(i):
		plt.cla()
		#scat = plt.scatter(x,y,s=0.05,color='black')
		t = data[i][0]
		scat = ax.scatter([p[0] for p in data[i][1:] if p[7] > mThresh], [p[1] for p in data[i][1:] if p[7] > mThresh],
			  s=[8*p[6] for p in data[i][1:] if p[7] > mThresh],
				alpha=1,color='r')
		if streaks:
			for j in range(1,min(i,10)):
				scat = ax.scatter([p[0] for p in data[i-j][1:] if p[7] > mThresh], [p[1] for p in data[i-j][1:] if p[7] > mThresh],
					alpha=(1.-j/10.)/4, s = 3*(1-j/10.),color='r')
		scat = ax.scatter([p[0] for p in data[i][1:] if p[7] < mThresh], [p[1] for p in data[i][1:] if p[7] < mThresh],
			  s=[1*p[6] for p in data[i][1:] if p[7] < mThresh],
				alpha=0.5,color='black')
		plt.title("Time = " + str(t))
		plt.xlim(Min,Max)
		plt.ylim(Min,Max)
		plt.grid(True)
		return scat 
	
	if noShow:
		plt.ioff()

	fig, ax = plt.subplots(figsize=(size,size))
	animate(0)

	anim = animation.FuncAnimation(fig,animate,frames=len(data),interval=10,blit=False)
	anim.save(filename,fps=fps,extra_args=['-vb','5M','-vcodec', 'mpeg4'])


def make3dDynamicsMovie(data, filename="output.mp4", Min=-20, Max=20, size=8, fps=30, noShow=True):
	"""
	Make a 3D movie of the dynamics from kira, using the data from parseStoryDynamics
	filename - duh
	Min/Max  - size of the plot range
	size	 - size of the saved figures
	fps	  - duh
	noShow   - don't show the movie plots as they're being made
	"""
	from mpl_toolkits.mplot3d import Axes3D
	from matplotlib import animation, pyplot as plt

	def animate(i):
		plt.cla()
		x = [p[0] for p in data[i][1:]]
		y = [p[1] for p in data[i][1:]]
		z = [p[2] for p in data[i][1:]]
		scat = ax.scatter3D(x,y,z,
				 s= [20*p[6] for p in data[i][1:]],
						  alpha=0.8)
		for j in range(1,min(i,50)):
			scat = ax.scatter3D([p[0] for p in data[i-j][1:]],
								[p[1] for p in data[i-j][1:]],
								[p[2] for p in data[i-j][1:]],
					alpha=(1-j/50.), s = 5*(1-j/50.))
		ax.set_xlim(Min,Max) 
		ax.set_ylim(Min,Max)
		ax.set_zlim(Min,Max)
		ax.grid("off")
		ax.set_axis_off()
		return scat

	if noShow:
		plt.ioff()

	fig = plt.figure(figsize=(size,size))
	#ax = fig.add_subplot(111, projection='3d')
	ax = Axes3D(fig)
	anim = animation.FuncAnimation(fig,animate,frames=len(data),interval=10,blit=False)
	anim.save(filename,fps=fps,extra_args=['-vb','5M','-vcodec', 'mpeg4'])

def make3dRotatingMovie(data, filename="output.mp4", size=8, fps=30, noShow=True):
	"""
	Make a 3D movie of the dynamics from kira, using the data from parseStoryDynamics
	data	 - a *single* frame of data from parseStoryDynamics
	filename - duh
	Min/Max  - size of the plot range
	size	 - size of the saved figures
	fps	  - duh
	noShow   - don't show the movie plots as they're being made
	"""
	from mpl_toolkits.mplot3d import Axes3D
	from matplotlib import animation, pyplot as plt
	from numpy import pi, sum, abs, array, median
	from itertools import combinations, product

	def animate(i):
		plt.cla()
		scat = ax.scatter(x,y,z, s=sizes,c=color,alpha=1.)
		ax.set_xlim(Min,Max) 
		ax.set_ylim(Min,Max)
		ax.set_zlim(Min,Max)
		a = i*360./250.
		ax.view_init(azim=a)
		ax.grid(False)
		r = [Min, Max]
		for s, e in combinations(array(list(product(r,r,r))), 2):
			if sum(abs(s-e)) == r[1]-r[0]:
					scat =  ax.plot3D(*zip(s,e), color="black",lw=0.5)
		ax.set_axis_off()
		return scat

	if noShow:
		plt.ioff()

	fig = plt.figure(figsize=(size,size))
	ax = fig.add_subplot(111, projection='3d')

	x = [p[0] for p in data[1:]]
	y = [p[1] for p in data[1:]]
	z = [p[2] for p in data[1:]]
	color = [p[6] for p in data[1:]]
	sizes = [5*p[6] for p in data[1:]]
	Min = -median(x+y+z)*100
	Max = -Min


	anim = animation.FuncAnimation(fig,animate,frames=25,interval=10,blit=False)
	anim.save(filename,fps=fps,extra_args=['-vb','10M','-vcodec', 'mpeg4'])
