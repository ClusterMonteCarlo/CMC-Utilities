def parseStoryDynamics(filename):
    data = []
    i = -1
    count = 0
    
    file = open(filename)
    line = file.readline() #can't use 'with open as', the iterator breaks for some reason
    
    while line:
        if line.split()[0] == 'system_time':
            t = float(line.split()[2])
            data.append([t])
            i += 1
            while line.split()[0] != ')Dynamics':
                line = file.readline() #advance forward to end up dynamics
        elif line.split()[0] == 'N':
            if line.split()[2] == '2':
                binFlag = 10
                count = 2
            elif line.split()[2] == '3':
                binFlag = 50
                count = 3
            else:
                binFlag = 1
        elif line.split()[0] == 'r':
            x,y,z = line.split()[2:]
            x = float(x)
            y = float(y)
            z = float(z)
        elif line.split()[0] == 'v':
            vx,vy,vz = line.split()[2:]
            vx = float(vx)
            vy = float(vy)
            vz = float(vz)
        elif line.split()[0] == ')Dynamics':
            data[i].append([x,y,z,vx,vy,vz,binFlag])
            while count != 0:
                if line.split()[0] != ')Dynamics':
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
    
    file = open(filename)
    line = file.readline() #can't use 'with open as', the iterator breaks for some reason
    
    while line.split()[0] != 'system_time':
        line = file.readline()
    t = float(line.split()[2])
    data.append(t)
    line = file.readline()
    numSingle = 0
    numBinary = 0
    numTriple = 0
    
    while line: #C SYNTAX 4EVAR!!!
        if line.strip() == 'N = 1':
            numSingle += 1
        elif line.strip() == 'N = 2':
            numBinary += 1
        elif line.strip() == 'N = 3':
            numTriple += 1
        elif line.split()[0] == 'total_energy':
            eTot = float(line.split()[2])
        elif line.split()[0] == 'potential_energy':
            ePot = float(line.split()[2])
        elif line.split()[0] == 'kinetic_energy':
            eKin = float(line.split()[2])
        elif line.split()[0] == 'r_lagr':
            lagRad = line.split()[2:]
            lagRad = stringsToFloats(lagRad)
        elif line.split()[0] == 'system_time':
            virial = abs(2*eKin/ePot)
            data[i].append(eTot,ePot,eKin,virial,numSingle,numBinary,numTriple)
            data[i].append(lagRad)
            i += 1
            t = float(line.split()[2])
            data.append([t])
            numSingle = 0
            numBinary = 0
            numTriple = 0
        line = file.readline()
        
    file.close()

    return data 
