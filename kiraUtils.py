def parseStoryDynamics(filename,scattering=False):
	"""
	Imports a Starlab Story for the dynamical information

	For a file with S snapshots and N particles, it yields an object of size

	story[S][N][9]

	Where the 8 elements of the array for each particle are:
		[x,y,z,vx,vy,vz,binFlag,m,id]
	
	If Scattering is true, it ignores the COM correction, since Fewbody doesn't
	care
	"""
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
			if scattering:
				count = 0
			while line.split()[0] != ')Dynamics':
				line = file.readline() #advance forward to end up dynamics
				if line.split()[0] == 'com_pos':
					xCOM,yCOM,zCOM = stringsToFloats(line.split()[2:])
				elif scattering and line.split()[0] == 'r':
					xCOM,yCOM,zCOM = stringsToFloats(line.split()[2:])
					xCOM,yCOM,zCOM = 0,0,0 
				elif line.split()[0] == 'com_vel': 
					vxCOM,vyCOM,vzCOM = stringsToFloats(line.split()[2:])
				elif scattering and line.split()[0] == 'v': 
					vxCOM,vyCOM,vzCOM = stringsToFloats(line.split()[2:])
		elif line.split()[0] == 'i':
			idNum = int(line.split()[2])
		elif line.split()[0] == 'N':
			if line.split()[2] == '2':
				binFlag = 2
				count = 2
			elif line.split()[2] == '3':
				binFlag = 3
				count = 3
			elif line.split()[2] == '4':
				binFlag = 4
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
			while count != 0 and len(line.split()) > 0:
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
	from numpy import zeros,array,hstack,isnan
	eTot,ePot,eKin,virial,numSingle,numBinary,numTriple,numQuad = zeros(8) 
	lagRad = 0

	file = open(filename)
	line = file.readline() #can't use 'with open as', the iterator breaks for some reason

	while line.split()[0] != ')Dynamics':
		line = file.readline() #advance forward to end up dynamics
		if line.split()[0] == 'com_pos':
			xCOM,yCOM,zCOM = stringsToFloats(line.split()[2:])
		elif line.split()[0] == 'com_vel': 
			vxCOM,vyCOM,vzCOM = stringsToFloats(line.split()[2:])
	
	while line.split()[0] != 'system_time':
		line = file.readline()
	t = float(line.split()[2])

	if t == 0.: #need to skip the initial print, if present
		line = file.readline()
		while line.split()[0] != 'system_time':
			line = file.readline()
		t = float(line.split()[2])
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
		elif line.split()[0] == 'com_pos':
			xCOM,yCOM,zCOM = stringsToFloats(line.split()[2:])
		elif line.split()[0] == 'com_vel': 
			vxCOM,vyCOM,vzCOM = stringsToFloats(line.split()[2:])
		elif line.split()[0] == 'system_time':
			virial = abs(2*eKin/(ePot))
			data.append(tuple(hstack([t,array([eTot,ePot,eKin,virial,eDyn,scale,numSingle,numBinary,numTriple,numQuad]),lagRad,xCOM,yCOM,zCOM,vxCOM,vyCOM,vzCOM])))
			t = float(line.split()[2])
			numSingle = 0
			numBinary = 0
			numTriple = 0
			numQuad = 0
		line = file.readline()
		
	file.close()
	del data[-1]

	i = 0
	while isnan(data[i][4]):
		i+=1

	out_array = array(data[i:], dtype=[('t', '<f8'), ('e', '<f8'),
	('pe', '<f8'), ('ke', '<f8'), ('vir', '<f8'), ('efric', '<f8'), ('scale',
	'<f8'), ('n', '<i8'), ('nb', '<i8'), ('nt', '<i8'), ('nq', '<i8'), ('lag90',
	'<f8'), ('lag80', '<f8'), ('lag70', '<f8'), ('lag60', '<f8'), ('lag50',
	'<f8'), ('lag40', '<f8'), ('lag30', '<f8'), ('lag20', '<f8'), ('lag10',
	'<f8'), ('xCOM', '<f8'), ('yCOM', '<f8'), ('zCOM', '<f8'), ('vxCOM', '<f8'),
	('vyCOM', '<f8'), ('vzCOM', '<f8')]) 
	return out_array


def plotNTotals(data):
	"""
	Plots number of singles, binaries, and higher over time
	Takes the data from parseStoryTotals
	"""
	import matplotlib.pyplot as plt
	from numpy import array
	plt.ion()
	plt.figure()
	t = data['t']
	plt.plot(t,data['n'])
	plt.plot(t,data['nb'])
	plt.plot(t,data['nt'])
	plt.plot(t,data['nq'])
	plt.grid(True)
	plt.legend((r'$N_{total}$',r'$N_{bin}$',
	            r'$N_{trip}$',r'$N_{quad}$'),loc=6)
	plt.xlabel('Time (crossing times)')
	plt.ylabel('Number')

def plotLagRad(data):
	"""
	Plots number Lagrange Radii over time 
	Takes the data from parseStoryTotals
	"""
	import matplotlib.pyplot as plt
	from numpy import array
	plt.ion()
	plt.figure()
	lagrads = ['lag90','lag80','lag70','lag60','lag50','lag40','lag30','lag20','lag10',]
	for lagrad in lagrads:
		plt.plot(data['t'],data[lagrad])
	plt.grid(True)
	plt.xlabel('Time (Crossing Times)')
	plt.ylabel('Lagrange Radii')
	
def plotEnergies(data):
	"""
	Plots various energy quantities over time 
	Takes the data from parseStoryTotals
	"""
	import matplotlib.pyplot as plt
	from numpy import array
	plt.ion()
	plt.figure()
	plt.plot(data['t'],data['e'])
	plt.plot(data['t'],data['pe'])
	plt.plot(data['t'],data['ke'])
	plt.plot(data['t'],data['vir'])
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

def pos_and_vel_to_xyz(r,vr,vt):
	from numpy import arccos,cos,sin,pi,sqrt
	from numpy.random import uniform
	theta = arccos(uniform(-1,1))
	phi = uniform(0,2*pi)
	angle = uniform(0,2*pi)

	x = r*sin(theta)*cos(phi)
	y = r*sin(theta)*sin(phi)
	z = r*cos(theta)

	v = sqrt(vr**2 + vt**2)
	thetaDot = cos(angle)*vt / r
	phiDot = sin(angle)*vt / (r * sin(theta))
	
	vx = vr*sin(theta)*cos(phi) + \
	     r*thetaDot*cos(theta)*cos(phi) - \
		 r*phiDot*sin(theta)*sin(phi)
	vy = vr*sin(theta)*sin(phi) + \
	     r*thetaDot*cos(theta)*sin(phi) + \
		 r*phiDot*sin(theta)*cos(phi)
	vz = vr*cos(theta) - r*thetaDot*sin(theta)
	return x,y,z,vx,vy,vz

def make2dDynamicsMovie(data, filename="output.mp4", Min=-20, Max=20, size=8,
fps=30, point_size = 100000, noShow=True, streaks=False, num_frames = 0,COM=False):
	"""
	Make a 2D movie of the dynamics from kira, using the data from parseStoryDynamics
	filename - duh
	Min/Max  - size of the plot range
	size	 - size of the saved figures
	fps	  - duh
	point_size - how big the points on the figure are (default 100000)
	noShow   - don't show the movie plots as they're being made
	streaks - put those crappy streaks behind the stars as 
	num_frames - how many frames to draw (default is all)
	COM - put a + at the center of mass
	"""
	from matplotlib import animation, pyplot as plt
	from numpy import loadtxt,sqrt
	from numpy.linalg import norm

	def animate(i,data,ax):
		plt.cla()
		t = data[i][0]
		scat = ax.scatter([p[0] for p in data[i][1:] if norm(p[:3]) < Max and p[6] == 1],
						  [p[1] for p in data[i][1:] if norm(p[:3]) < Max and p[6] == 1],
				  s=[point_size*p[7] for p in data[i][1:] if norm(p[:3]) < Max and p[6] == 1],
				    alpha=1,color='r')
		scat = ax.scatter([p[0] for p in data[i][1:] if norm(p[:3]) < Max and
		p[7] >= 0.00051],
						  [p[1] for p in data[i][1:] if norm(p[:3]) < Max and
						  p[7] >= 0.00051],
				  s=[2*point_size*p[7] for p in data[i][1:] if norm(p[:3]) < Max
				  and p[7] >= 0.00051],
					alpha=1,color='b')
		scat = ax.scatter([p[0] for p in data[i][1:] if norm(p[:3]) < Max and p[6] > 2],
						  [p[1] for p in data[i][1:] if norm(p[:3]) < Max and p[6] > 2],
				  s=[4*point_size*p[7] for p in data[i][1:] if norm(p[:3]) < Max and p[6] > 2],
					alpha=1,color='g')
		if streaks:
			for j in range(1,min(i,10)):
				scat = ax.scatter([p[0] for p in data[i-j][1:]], [p[1] for p in data[i-j][1:]],
					alpha=(1.-j/10.)/4, s = 3*(1-j/10.),color='r')
		if COM:
			X_COM = sum([p[0]*p[7] for p in data[i][1:]]) / sum([p[7] for p in data[i][1:]])
			Y_COM = sum([p[1]*p[7] for p in data[i][1:]]) / sum([p[7] for p in data[i][1:]])
			scat = ax.scatter(X_COM,Y_COM,marker='x',color='black',s=100)

			X_COV = sum([p[3]*p[7] for p in data[i][1:]]) / sum([p[7] for p in data[i][1:]])
			Y_COV = sum([p[4]*p[7] for p in data[i][1:]]) / sum([p[7] for p in data[i][1:]])
			scat = ax.scatter(X_COV,Y_COV,marker='x',color='r',s=100)
		plt.title("Time = " + str(t)[:5])
		plt.xlim(Min,Max)
		plt.ylim(Min,Max)
		plt.grid(True)
		return scat 
	
	if noShow:
		plt.ioff()

	fig, ax = plt.subplots(figsize=(size,size))


	if num_frames == 0:
		num_frames = len(data) - 1

	if num_frames > len(data):
		print "More frames than snapshots!"
		return

	anim = animation.FuncAnimation(fig,animate,frames=num_frames,fargs=[data,ax],blit=False)
	anim.save(filename,writer='ffmpeg',fps=20,dpi=300)


#anim.save(filename, writer=None, fps=None, dpi=None, codec=None, bitrate=None,extra_args=None, metadata=None, extra_anim=None, savefig_kwargs=None)


def gauplot(centers, radiuses, xr=None, yr=None):
	import matplotlib.pyplot as plt,numpy as np
	nx, ny = 1000.,1000.
	xgrid, ygrid = np.mgrid[xr[0]:xr[1]:(xr[1]-xr[0])/nx,yr[0]:yr[1]:(yr[1]-yr[0])/ny]
	im = xgrid*0 + np.nan
	xs = np.array([np.nan])
	ys = np.array([np.nan])
	fis = np.concatenate((np.linspace(-np.pi,np.pi,100), [np.nan]) )
	cmap = plt.cm.gray
	cmap.set_bad('white')
	thresh = 3
	for curcen,currad in zip(centers,radiuses):
			curim=(((xgrid-curcen[0])**2+(ygrid-curcen[1])**2)**.5)/currad*thresh
			im[curim<thresh]=np.exp(-.5*curim**2)[curim<thresh]
			xs = np.append(xs, curcen[0] + currad * np.cos(fis))
			ys = np.append(ys, curcen[1] + currad * np.sin(fis))
	plt.imshow(im.T, cmap=cmap, extent=xr+yr)
	plt.plot(xs, ys, 'r-')

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
	 

def computeLagrad(data,lagrad = [1,2,5,10,25,50,75,90]):
	"""
	Compute the lagrange radii for the percentages given by lagrad
	from the given data (story dynamics) 
	"""
	from numpy import array, sum, argsort
	from numpy.linalg import norm

	output = []

	for snap in data:
		time = snap[0]
		particles = snap[1:]
		mANDr = array([[p[7],norm(p[:3])] for p in particles])
		mass = sum(mANDr.T[0])
		args = argsort(mANDr.T[1])
		mANDr = mANDr[args]
		radii = array(lagrad)*mass/100. 

		l = 0
		cumMass = 0
		LagRad = []
		for particle in mANDr:
			cumMass += particle[0]
			if cumMass >= radii[l]:
				LagRad.append(particle[1])
				l += 1
				if l == len(radii):
					break
		
		output.append([time,LagRad])

	return output

def plotComputedLagrad(LagRad,lagrad = [1,2,5,10,25,50,75,90]):
	"""
	Plots the lagrange radii, given by lagrad
	Takes output from computeLagrad
	"""
	import matplotlib.pyplot as plt
	from numpy import array
	times = []
	radii = []
	for things in LagRad:
		times.append(things[0])
		radii.append(things[1])

	times = array(times) # * 0.11440499869028585 #  0.057407700000000006
	radii = array(radii)

	for rad in radii.T:
		plt.plot(times,rad)

	plt.legend(lagrad)

	
def computeAngularMomentum(data,threshold):
	"""
	Computes the angular momentum support of the innermost 'threshold' number of
	particles

	Returns [support, L], where support is an array of the form [time, radii at
	threshold, T/W, beta] and L is the angular momentum vectors
	"""
	things = []
	angMom = []
	from numpy import array, dot
	from numpy.linalg import norm

	for snapshot in data:
		posBHs = array([p[:3] for p in snapshot[1:]])
		velBHs = array([p[3:6] for p in snapshot[1:]])
		masses = array([p[7] for p in snapshot[1:]])
		radii = array([norm(p) for p in posBHs])
		posBHs = posBHs[radii.argsort()]
		velBHs = velBHs[radii.argsort()]
		masses = masses[radii.argsort()]

		vr2Ave = 0
		vt2Ave = 0

		posBH = posBHs[:threshold]
		velBH = velBHs[:threshold]
		m = masses[:threshold]

		for i in range(len(posBH)):
			rHat = posBH[i] / norm(posBH[i])
			vRad = dot(rHat,velBH[i])
			vTan = norm(velBH[i] - vRad*rHat)
			vr2Ave += vRad**2
			vt2Ave += vTan**2
		vr2Ave /= len(posBH)
		vt2Ave /= len(posBH)

		W = 0
		L = zeros(3)
		I = 0

		for i in range(len(posBH)):
			for j in range(len(posBH)):
				if i == j:
					continue
				W += m[i]*m[j]/norm( posBH[i] - posBH[j])

		L = zeros(3)
		for i in range(len(posBH)):
			L += m[i]*cross(posBH[i],velBH[i])

		I = 0
		for i in range(len(posBH)):
			I += m[i]*dot(posBH[i],posBH[i])

		things.append(array([snapshot[0],radii[threshold-1],norm(L)**2/(I*W),1-vt2Ave/vr2Ave]))
		angMom.append(L)

	return array(things), array(angMom)
	

def make2dScatteringMovie(data, filename="output.mp4", Min=-20, Max=20, size=8,
fps=30, point_size = 100000, noShow=True, streaks=False, num_frames = 0):
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
		x = [p[0] for p in data[i][1:] if p[7] < .25]
		y = [p[1] for p in data[i][1:] if p[7] < .25]
		z = [p[2] for p in data[i][1:] if p[7] < .25]
		particles = len(data[i][1:])
		scat = ax.plot(x,y,z,lw=0,zdir='z',
				 markersize=20, markeredgecolor='black',
						  alpha=1.0,color='cyan',marker='*',zorder=5)
		x = [p[0] for p in data[i][1:] if p[7] > .25 and p[7] < 0.85]
		y = [p[1] for p in data[i][1:] if p[7] > .25 and p[7] < 0.85]
		z = [p[2] for p in data[i][1:] if p[7] > .25 and p[7] < 0.85]
		particles = len(data[i][1:])
		scat = ax.plot(x,y,z,lw=0,zdir='z',
				 markersize=20, markeredgecolor='C1',
						  alpha=1.0,marker="*",c='C1',zorder=5)
		x = [p[0] for p in data[i][1:] if p[7] > .85]
		y = [p[1] for p in data[i][1:] if p[7] > .85]
		z = [p[2] for p in data[i][1:] if p[7] > 0.85]
		particles = len(data[i][1:])
		scat = ax.plot(x,y,z,lw=0,zdir='z',
				 markersize=0.030,markeredgecolor='black',
						  alpha=1.0,marker="*",c='C3',zorder=5)

		if streaks:
			for par in range(len(data[0][1:])):
				X = []
				Y = []
				Z = []
				last_idx = 0
				for j in range(i+1):
					if len(data[j][1:]) < len(data[0][1:]) and last_idx == 0:
						last_idx = j-1
						if data[last_idx][par+1][7] < .25:
							C = 'cyan'
						elif data[last_idx][par+1][7] > .25: 
							C = 'C1'
						else:
							C = 'C3'
						AL = 0.751
						ax.plot(X,Y,Z,color=C,alpha=AL,lw=1.5,zorder=2)
						X = []
						Y = []
						Z = []
					if par >= len(data[j][1:]):
						break
					X.append(data[j][par+1][0])
					Y.append(data[j][par+1][1])
					Z.append(data[j][par+1][2])
						
				if par >= len(data[i][1:]):
					continue
				if data[i][par+1][7] < .25:
					C = 'cyan'
					wid = 1.5
				elif data[i][par+1][7] > .25: 
					C = 'C1'
					wid = 1.5
				else:
					C = 'C3'
					wid = 2.5
#				if data[i][par+1][8] == 2:
#					AL = 0.0
#				else:
				AL = 0.751
				ax.plot(X,Y,Z,color=C,alpha=AL,lw=wid,zorder=2)

		#ax.set_xlim(0.37,0.3705) 
		#ax.set_ylim(-0.09295,-0.09306)
		#ax.set_zlim(-0.01828,-0.01831)
		#ax.set_xlim(0.378,0.380) 
		#ax.set_ylim(-0.08985,-0.09916)
		#ax.set_zlim(-0.01818,-0.02241)
		ax.set_xlim(-30,40) 
		ax.set_ylim(Min/1.,Max/1.)
		ax.set_zlim(Min/1.,Max/1.)
		ax.view_init(azim=0,elev=0)
		ax.grid(False)
		ax.set_axis_off()
		ax.set_frame_on(False)
		ax.set_facecolor('black')
		for axis in ax.w_xaxis, ax.w_yaxis, ax.w_zaxis: 
		  for elt in axis.get_ticklines() + axis.get_ticklabels(): 
			  elt.set_visible(False) 

		  axis.pane.set_visible(False) 
		  axis.gridlines.set_visible(False) 
		  axis.line.set_visible(False) 
		return scat

	if noShow:
		plt.ioff()

	if num_frames == 0:
		num_frames = len(data) - 1

	if num_frames > len(data):
		print "More frames than snapshots!"
		return

	fig = plt.figure(figsize=(size,10))
	fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
	ax = fig.add_subplot(111, projection='3d')

	anim = animation.FuncAnimation(fig,animate,frames=num_frames,interval=50,blit=False)
	anim.save(filename,fps=fps,extra_args=['-codec:v','libx264','-profile:v','high','-preset','slow','-b:v', '25M','-pix_fmt','yuv420p'])

	thing = animate(num_frames)
	plt.savefig("/Users/carl/Desktop/test.png")

	#anim.save(filename,fps=fps,extra_args=['-vb','5M','-vcodec','mpeg'])
