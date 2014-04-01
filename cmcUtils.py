def plotLagRad(filename = '', gamma = 0, N = 0,Color='r'):
    """ plots the Lagrangian Radii for the given filename
    Filename -- obvious (defaults to *.lagrad.rad if in right folder)
    gamma and N -- for conversion to N-body time"""
    from os import listdir
    from numpy import loadtxt, log
    from matplotlib.pyplot import plot,yscale,ylim,grid,xlabel,ylabel,text,title
    if filename == '':
        filename = [name for name in listdir('.') if name.endswith('.lagrad.dat')][0]
    colsToImport = (0,9,14,22,24,28,29,30,32,33,34,35,36,37,38,39)
    data = loadtxt(filename,unpack=True,usecols=colsToImport)
    xName = 'Time (FP Units)'
    yName = 'Lagrange Radii'
    radii = 'Lagrange Radii (%): 0.1,0.2,0.5,1,2,5,10,20,30,40,50,60,70,80,90'
    if not gamma == 0 and not N == 0:
        data[0] = data[0] * N / log(gamma*N)
        xName = 'Time (N-Body Units)'
    for index in range(1,len(colsToImport)):
        plot(data[0],data[index],color=Color)
    yscale('log')
    grid(True) 
    xlabel(xName)
    ylabel(yName)
    #text(0.1,0.5*ylim()[1],radii)
    title(radii)

def find_nearest(array,value):
    from numpy import abs,argmin
    return (abs(array-value)).argmin()

def plotEnergy(filename='',gamma=0,N=0,quant=['etot']):
    """ plots the energies to the below quantities 
    Filename -- obvious (defaults to *.dyn.dat if in right folder)
    gamma and N -- for conversion to N-body time
    quant -- various energies to plot
        etot       - Total energy(default)
        KE         - Kinetic Energy
        PE         - Potential Energy
        Etot_int   - Total internal energy
        Etot_bin   - Total binary binding energy
        Etot_cenma - Total Center Mass Energy
        Eesc       - Total escaped Energy
        Ebesc      - Total escaped Binary Energies
        Eintesc    - Total escaped internal energy
        Eoops      - Total energy lost during regularization
        EoT        - Total + Eoops"""
    from os import listdir
    from numpy import loadtxt, log
    from matplotlib.pyplot import figure,plot,legend,ylim,grid,xlabel,ylabel,text,title
    if filename == '':
        filename = [name for name in listdir('.') if name.endswith('.dyn.dat')][0]
    colsToImport = tuple([0])
    legendText = []
    xName = 'Time (FP Units)'
    yName = 'Energy (E_TOT(0) = 0.25)'
    if not gamma == 0 and not N == 0:
        data[0] = data[0] * N / log(gamma*N)
        xName = 'Time (N-Body Units)'
    if 'etot' in quant:
        colsToImport += tuple([9])
        legendText.append('Total')
    if 'KE' in quant:
        colsToImport += tuple([10])
        legendText.append('KE')
    if 'PE' in quant:
        colsToImport += tuple([11])
        legendText.append('PE')
    if 'Etot_int' in quant:
        colsToImport += tuple([12])
        legendText.append('Etot_internal')
    if 'Etot_bin' in quant:
        colsToImport += tuple([13])
        legendText.append('Etot_binary')
    if 'Etot_centma' in quant:
        colsToImport += tuple([14])
        legendText.append('Etot_centerMass')
    if 'Eesc' in quant:
        colsToImport += tuple([15])
        legendText.append('E_escape')
    if 'Ebesc' in quant:
        colsToImport += tuple([16])
        legendText.append('E_binEscape')
    if 'Eintesc' in quant:
        colsToImport += tuple([17])
        legendText.append('E_internalEscape')
    if 'Eoops' in quant:
        colsToImport += tuple([18])
        legendText.append('E_oops')
    if 'EoT' in quant:
        colsToImport += tuple([19])
        legendText.append('E_oops + E_total')
    data = loadtxt(filename,unpack=True,usecols=colsToImport)
    figure()
    for index in range(1,len(colsToImport)):
        plot(data[0],data[index])
    grid(True) 
    xlabel(xName)
    ylabel(yName)
    legend(legendText,loc='best')
    title('Cluster Energy')


def plotEnergyPerR(filename='',movie=False,index=0,gamma=0,N=0):
    from os import listdir,system
    from subprocess import check_output
    from numpy import loadtxt,nonzero,min,max
    from matplotlib.pyplot import subplot,subplots_adjust,figure,xlabel,ylabel,title,plot,xscale,yscale,xlim,ylim,grid,scatter,savefig,suptitle,ioff,text 

    #Get the filenames of the snapshots
    names = listdir('.')
    snaps = [snap for snap in names if 'snap' in snap] #lol
    snaps.sort()
    filebase = snaps[0].split('.snap')[0]
    if movie:
        system('rm *png') #clears out any already existing files
        if '.gz' in snaps[0] or '.gz' in snaps[-1]:
            print 'Snapshots haven\'t been unzipped; it might take a while'
            system('gunzip *.dat.gz') #there are like a million better ways to do this
            print 'Done with that'
        names = listdir('.') #redo this
        snaps = [snap for snap in names if 'snap' in snap] #lol
        snaps.sort()
    else:
        if '.gz' in snaps[index]:
            system('gunzip ' + snaps[index]) 
            names = listdir('.')
            snaps = [snap for snap in names if 'snap' in snap] #lol
            snaps.sort()

    #Find times from top of snapshots
    times = []    
    if movie:
        gen = [snap for snap in snaps]
    else:
        gen = [snaps[index]]

    for snap in gen:
        headCommand = 'head -n 1 ' + snap
        time = check_output(headCommand.split()) #obtuse, but much faster than python
        times.append(time.split()[1])

    #Import the Lagrange Radii
    colsToImport = (0,9,14,22,24,28,29,30,32,33,34,35,36,37,38,39)
    lagRad = loadtxt(filebase + '.lagrad.dat',unpack=True,usecols=colsToImport)
    radii = 'Lagrange Radii (%): 0.1,0.2,0.5,1,2,5,10,20,30,40,50,60,70,80,90'
    minRadius = min(lagRad[1])/10.
    maxRadius = max(lagRad[-1])*10
    endTime = max(lagRad[0])

    if movie:
        ioff()
        print "Generating actual plots now"
    i = 0
    for snap in gen:
        data = loadtxt(snap,usecols=(1,2,5,7,23,24),unpack=True) #m,r,energy,binary,m1,m2,a

        #calculate the binding energies of the binaries
        binInd = nonzero(data[3])
        binR = data[1][binInd]
        eta = data[5][binInd] 

        #and then plot things
        figure(figsize=(8,10))
        subplot(311)
        subplots_adjust(hspace=0.35,bottom=0.05,top=0.98)
        plot(data[1],data[2],markevery=100)
        xscale('log')
        yscale('symlog')
        xlim(minRadius,maxRadius)
            #xlim(1e-4,1e4)
        ylim(-100,1.1)
        Limits = xlim()
        grid(True)
        xlabel('Radius')
        ylabel('Energy')
        subplot(312)
        scatter(binR,eta)
        xlabel('Radius')
        ylabel('Binary Hardness ($\eta$)')
        grid(True)
        xlim(Limits)
        if movie:
            ylim(-100,1.1)
        xscale('log')
        yscale('symlog')
        #suptitle('Energy at ' + times[i] + ' (FP Time)',fontsize=18)
        subplot(313)
        yscale('log')
        ylim(Limits)
        xlim(0,endTime)
        timeStop = find_nearest(lagRad[0],float(times[i].split('=')[1]))
        for index in range(1,len(colsToImport)):
            plot(lagRad[0][:timeStop],lagRad[index][:timeStop],color='r')
        title(radii) 
        xlabel('Time')
        ylabel('Radius')
        grid(True)
        if movie:
            savefig(snap[:-3] + 'png')
        i += 1

    if movie:
        print "And finally for the movie"
        system("rm *mp4")
        ffmpegCommand = "ffmpeg -r 30 -i " + filebase + ".snap%04d.png -r 30 -vb 5M -vcodec mpeg4 " + filebase + ".mp4"
        #system(ffmpegCommand)
        #system("rm *.png")


def escapedBH(filename = '', gamma = 0, N = 0, nBH = 65, thresh = 0.98094568, Color='r'):
    """ plots the escaped heavy particles with mass above threshold from the .esc.dat file 
    (it kinda assumes everything is a multiple of thresh)
    Filename -- obvious (defaults to *.esc.dat if in right folder)
    gamma and N -- for conversion to N-body time
    nBH -- number of initial black holes"""
    from os import listdir
    from numpy import loadtxt, log, transpose
    from matplotlib.pyplot import plot,yscale,ylim,grid,xlabel,ylabel,text,title
    if filename == '':
        filename = [name for name in listdir('.') if name.endswith('.esc.dat')][0]
    data = loadtxt(filename,usecols=(1,2))
    new = [[dat[0],round(dat[1]/thresh/2/10)] for dat in data if dat[1] > thresh]
    for dat in new:
        nBH += -dat[1]
        dat[1] = nBH
    new = transpose(new)
    if not gamma == 0 and not N == 0:
        new[0] = N * new[0] / log( gamma * N )
    plot(new[0],new[1],color=Color)


def rToXY(r):
    from numpy import arccos,cos,sin,pi
    from numpy.random import uniform
    sinTheta = sin(arccos(uniform(-1,1)))
    phi = uniform(0,2*pi)
    x = r*sinTheta*cos(phi)
    y = r*sinTheta*sin(phi)
    return x,y

def makeDynamicsMovie(filename='output.mp4',noShow=True,size=8,fps=30):
    from matplotlib import animation, pyplot as plt
    from os import listdir,system
    from numpy import loadtxt

    names = listdir('.')
    snaps = [snap for snap in names if 'snap' in snap] #lol
    snaps.sort()

    if '.gz' in snaps[0] or '.gz' in snaps[-1]:
        print 'Snapshots haven\'t been unzipped; it might take a while'
        system('gunzip *.dat.gz') #there are like a million better ways to do this
        print 'Done with that'

    def animate(i):
        r,m = loadtxt(snaps[i],usecols=(3,2),unpack=True)
        points = [rToXY(p) for p in r]
        x = [p[0] for p in points]
        y = [p[1] for p in points]
        scat = plt.scatter(x,y,s=1)
        return scat

    if noShow:
        plt.ioff()

    fig, ax = plt.subplots(figsize=(size,size))
    animate(0)

    anim = animation.FuncAnimation(fig,animate,frames=len(snaps)/100,interval=10,blit=False)
    anim.save(filename,fps=fps)#,extra_args=['-vb','5M','-vcodec', 'mpeg4'])




