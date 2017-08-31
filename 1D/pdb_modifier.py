#! /usr/bin/env python
import time

inputFile=open('pre_WAT1.pdb','r')
outputFile=open('WAT1.pdb','w')
resID = 1
atomID =0
for line in inputFile:
    line = line.strip()
    columns = line.split()
    if columns[0] == 'REMARK' or columns[0] =='CRYST1':
        outputFile.write(line)
        outputFile.write('\n')
    if columns[0] == 'TER':
        resID+=1
    if columns[0]=='ATOM':
        atom     = columns[0]
        atomID+=1
        atomIndex= '{:>7}'.format(atomID)
        atomName = '{:<3}'.format(columns[2])
        resName  = '{:>6}'.format(columns[3])
        if len(resName) > 6:
            residueID = resName[5:]
            x        = '{:>12}'.format(columns[4])
            y        = '{:>8}'.format(columns[5])
            z        = '{:>8}'.format(columns[6])
            if float(z) < 271:
                print z
                z = str('{:>8}'.format('%6.3f' % (float(z) + 69)))
            beta     = '{:>6}'.format(columns[7])
            occupancy= '{:>6}'.format(columns[8])
        else:
            residueID= columns[4]
            x        = '{:>12}'.format(columns[5])
            y        = '{:>8}'.format(columns[6])
            z        = '{:>8}'.format(columns[7])
            z = z
            if float(z) < 271:
                z = str('{:>8}'.format('%6.3f' % (float(z) + 69)))
            beta     = '{:>6}'.format(columns[8])
            occupancy= '{:>6}'.format(columns[9])
        #atomSym  = '{:>6}'.format(columns[10])
        outputFile.write(atom)
        outputFile.write(atomIndex)
        outputFile.write('  ')
        outputFile.write(atomName)
        outputFile.write(resName[:5])
        outputFile.write('{:>4}'.format(residueID))
        outputFile.write(x)
        outputFile.write(y)
        outputFile.write(z)
        outputFile.write(beta)
        outputFile.write(occupancy)
        outputFile.write('\n')

outputFile.write('TER')

inputFile.close()
outputFile.close()