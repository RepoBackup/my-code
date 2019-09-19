# -*- coding: utf-8 -*-

loc = 'C:/Users/edanb/Desktop'
filename = 'tempPOS.vasp'

axis = 1 # 0, 1, 2 = x, y, z
l = 43.2661589986
coordType = 'Cartesian'
  
with open(loc + '/' + filename, 'r') as f, open(loc + '/SD.vasp', 'w') as g:

    for line in range(5): 
        g.write(f.readline())
   
    atoms = f.readline(); g.write(atoms); atoms = atoms.split()
    nums = f.readline(); g.write(nums); nums = nums.split()
    g.write('SD\n')
    g.write(f.readline())

    for s in range(len(nums)):
        
        lines = [[float(value) for value in f.readline().split()] \
                  for line in range(int(nums[s]))]
        
        for line in sorted(lines, key=lambda line: line[axis]):
            SD = ' F F F' if .17 * l < line[axis] < .36 * l else ' T T T'
            g.write('%18.9f %18.9f %18.9f' % tuple(line) + SD + ' \n')
            
#    trigger = 0
#    for line in f.read().splitlines():
#        if trigger == 0:
#            if re.match(coordType, line):
#                trigger = 1
#                g.write('SD\n')
#            SD = ''
#        else:
#            if axis == 'x':
#                ax = float(line.split()[0])
#            elif axis == 'y':
#                ax = float(line.split()[1])
#            else:
#                ax = float(line.split()[2])
#            if .13 * l < ax < .35 * l: # or .61 * l < ax < .89 * l:
#                SD = ' F F F'
#            else:
#                SD = ' T T T'
#        g.write(line + SD + '\n')