from vpython import *
import numpy as np
from hw6_histogram import *

N = 200
m, size = 4E-3/6E23, 31E-12*10 # He atom are 10 times bigger for easiear collision but not too big for accuracy
L = ((24.4E-3/(6E23))*N)**(1/3.0)/2 + size # 2L is the cubic container's original length, width, and height
k, T = 1.38E-23, 298.0 # Boltzmann Constant and initial temperature
t, dt = 0, 3E-13
vrms = (2*k*1.5*T/m)**0.5 # the initial root mean square velocity
stage = 0 # stage number
atoms = [] # list to store atoms
delta_p = 0 #for recording momentum

# histogram setting
deltav = 50. # slotwidth for v histogram
vdist = graph(x=800, y=0, ymax = N*deltav/1000.,width=500, height=300, xtitle='v', ytitle='dN', align = 'left')
theory_low_T = gcurve(color=color.cyan) # for plot of the curve for the atom speed distribution
dv = 10.
for v in arange(0.,4201.+dv,dv): # theoretical speed distribution
    theory_low_T.plot(pos=(v,(deltav/dv)*N*4.*pi*((m/(2.*pi*k*T))**1.5)*exp((-0.5*m*v**2)/(k*T))*(v**2)*dv))
observation = ghistogram(graph = vdist, bins=arange(0.,4200.,deltav), color=color.red) # for the simulation speed distribution

#initialization
scene = canvas(width=500, height=500, background=vector(0.2,0.2,0), align = 'left')
container = box(length = 2*L, height = 2*L, width = 2*L, opacity=0.2, color = color.yellow )
p_a, v_a = np.zeros((N,3)), np.zeros((N,3)) # particle position array and particle velocity array, N particles and 3 for x, y, z
for i in range(N):
    p_a[i] = [2 * L*random() - L, 2 * L*random() - L, 2 * L*random() - L] # particle is initially random positioned in container
    if i== N-1: # the last atom is with yellow color and leaves a trail
        atom = sphere(pos=vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]), radius = size, color=color.yellow, make_trail = True, retain = 50)
    else: # other atoms are with random color and leaves no trail
        atom = sphere(pos=vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]), radius = size, color=vector(random(), random(), random()))
    ra = pi*random()
    rb = 2*pi*random()
    v_a[i] = [vrms*sin(ra)*cos(rb), vrms*sin(ra)*sin(rb), vrms*cos(ra)] # particle initially same speed but random direction
    atoms.append(atom)
    
def vcollision(a1p, a2p, a1v,a2v): # the function for handling velocity after collisions between two atoms
    v1prime = a1v - (a1p - a2p) * sum((a1v-a2v)*(a1p-a2p)) / sum((a1p-a2p)**2)
    v2prime = a2v - (a2p - a1p) * sum((a2v-a1v)*(a2p-a1p)) / sum((a2p-a1p)**2)
    return v1prime, v2prime

#event triggered
def keyinput(event):
    global stage
    if event.key == 'n':
        stage+=1
scene.bind('keydown', keyinput)

runs = 0
while stage == 0:
    t += dt
    runs+=1
    rate(10000)

    p_a += v_a*dt # calculate new positions for all atoms
    for i in range(N): atoms[i].pos = vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]) # to display atoms at new positions
    if stage != 1 : observation.plot(data = np.sqrt(np.sum(np.square(v_a),-1))) ## freeze histogram for stage != 1
    
    ### find collisions between pairs of atoms, and handle their collisions
    r_array = p_a-p_a[:,np.newaxis] # array for vector from one atom to another atom for all pairs of atoms
    rmag = np.sqrt(np.sum(np.square(r_array),-1)) # distance array between atoms for all pairs of atoms
    hit = np.less_equal(rmag,2*size)-np.identity(N) # if smaller than 2*size meaning these two atoms might hit each other
    hitlist = np.sort(np.nonzero(hit.flat)[0]).tolist() # change hit to a list
    for ij in hitlist: # i,j encoded as i*Natoms+j
        i, j = divmod(ij,N) # atom pair, i-th and j-th atoms, hit each other
        hitlist.remove(j*N+i) # remove j,i pair from list to avoid handling the collision twice
        if sum((p_a[i]-p_a[j])*(v_a[i]-v_a[j])) < 0 : # only handling collision if two atoms are approaching each other
            v_a[i], v_a[j] = vcollision(p_a[i], p_a[j], v_a[i], v_a[j]) # handle collision
            
    #find collisions between the atoms and the walls, and handle their elastic collisions
    for i in range(N):
        if abs(p_a[i][0]) >= L - size and p_a[i][0]*v_a[i][0] > 0 :
            v_a[i][0] = - v_a[i][0]
            delta_p += 2*m*abs(v_a[i][0])
        if abs(p_a[i][1]) >= L - size and p_a[i][1]*v_a[i][1] > 0 :
            v_a[i][1] = - v_a[i][1]
            delta_p += 2*m*abs(v_a[i][1])
        if abs(p_a[i][2]) >= L - size and p_a[i][2]*v_a[i][2] > 0 :
            v_a[i][2] = - v_a[i][2]
            delta_p += 2*m*abs(v_a[i][2])

    #hw request 1
    if runs%1000 == 0:
        ke = 0 #total translational ke
        for i in range(N):
            ke+=0.5*m*(v_a[i][0]**2+v_a[i][1]**2+v_a[i][0]**2)
        SurfaceArea = 2*(container.length*container.width + container.width*container.height + container.height*container.length)
        Volume  = container.length*container.width*container.height
        print("temperature:     "+str(ke / (3*N*k/2)))
        print("pressure:            "+str( delta_p/((1000*dt)*SurfaceArea) ) )
        print("Volume:             "+str(Volume))
        print("PV:                      "+str(delta_p/((1000*dt)*SurfaceArea) * Volume))
        print("NkT:                   "+str(2/3*ke))
        print("p*(V**gamma)"+str(   delta_p/((1000*dt)*(SurfaceArea))   * (Volume)**(5/3)  ))
        delta_p = 0

#=====================================      stage 1    ========================================
v_w =  L / (20000.0*dt)
delta_p = 0
runs = 0
container.length = 2*L
container.width = 2*L
container.height = 2*L
while stage == 1:
    t += dt
    runs += 1
    rate(10000)

    #walls are moving closer
    container.length -= 2*v_w*dt
    if(container.length < L):
        stage+=1
    
    p_a += v_a*dt # calculate new positions for all atoms
    for i in range(N): atoms[i].pos = vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]) # to display atoms at new positions
    
    ### find collisions between pairs of atoms, and handle their collisions
    r_array = p_a-p_a[:,np.newaxis] # array for vector from one atom to another atom for all pairs of atoms
    rmag = np.sqrt(np.sum(np.square(r_array),-1)) # distance array between atoms for all pairs of atoms
    hit = np.less_equal(rmag,2*size)-np.identity(N) # if smaller than 2*size meaning these two atoms might hit each other
    hitlist = np.sort(np.nonzero(hit.flat)[0]).tolist() # change hit to a list
    for ij in hitlist: # i,j encoded as i*Natoms+j
        i, j = divmod(ij,N) # atom pair, i-th and j-th atoms, hit each other
        hitlist.remove(j*N+i) # remove j,i pair from list to avoid handling the collision twice
        if sum((p_a[i]-p_a[j])*(v_a[i]-v_a[j])) < 0 : # only handling collision if two atoms are approaching each other
            v_a[i], v_a[j] = vcollision(p_a[i], p_a[j], v_a[i], v_a[j]) # handle collision
            
    #find collisions between the atoms and the walls, and handle their elastic collisions
    for i in range(N):
        if abs(p_a[i][0]) >= container.length/2 - size and p_a[i][0]*v_a[i][0] > 0 :
            if(v_a[i][0] > 0): v_a[i][0] = - v_a[i][0] - 2 * v_w
            else: v_a[i][0] = - v_a[i][0] + 2 * v_w
            delta_p += 2*m*(abs(v_a[i][0]) + 2*v_w)
        if abs(p_a[i][1]) >= L - size and p_a[i][1]*v_a[i][1] > 0 :
            v_a[i][1] = - v_a[i][1]
            delta_p += 2*m*abs(v_a[i][1])
        if abs(p_a[i][2]) >= L - size and p_a[i][2]*v_a[i][2] > 0 :
            v_a[i][2] = - v_a[i][2]
            delta_p += 2*m*abs(v_a[i][2])

    #hw request
    if runs%1000 == 0:
        ke = 0 #total translational ke
        for i in range(N):
            ke+=0.5*m*(v_a[i][0]**2+v_a[i][1]**2+v_a[i][0]**2)
        SurfaceArea = 2*(container.length*container.width + container.width*container.height + container.height*container.length)
        Volume  = container.length*container.width*container.height
        print("temperature:     "+str(ke / (3*N*k/2)))
        print("pressure:            "+str( delta_p/((1000*dt)*SurfaceArea) ) )
        print("Volume:             "+str(Volume))
        print("PV:                      "+str(delta_p/((1000*dt)*SurfaceArea) * Volume))
        print("NkT:                   "+str(2/3*ke))
        print("p*(V**gamma)"+str(   delta_p/((1000*dt)*(SurfaceArea))   * (Volume)**(5/3)  ))
        delta_p = 0

#=====================================      stage 2   ========================================
ke = 0 #total translational ke
for i in range(N):
    ke+=0.5*m*(v_a[i][0]**2+v_a[i][1]**2+v_a[i][0]**2)
T = ke / (3*N*k/2)

#new graph
theory_high_T = gcurve(color=color.green)
for v in arange(0.,4201.+dv,dv): # theoretical speed distribution
    theory_high_T.plot(pos=(v,(deltav/dv)*N*4.*pi*((m/(2.*pi*k*T))**1.5)*exp((-0.5*m*v**2)/(k*T))*(v**2)*dv))
observation2 = ghistogram(graph = vdist, bins=arange(0.,4200.,deltav), color=color.blue) # for the simulation speed distribution


delta_p = 0
runs = 0
while stage == 2:
    t += dt
    runs += 1
    rate(10000)

    #walls are moving closer
    
    p_a += v_a*dt # calculate new positions for all atoms
    for i in range(N): atoms[i].pos = vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]) # to display atoms at new positions
    observation2.plot(data = np.sqrt(np.sum(np.square(v_a),-1))) 
    
    ### find collisions between pairs of atoms, and handle their collisions
    r_array = p_a-p_a[:,np.newaxis] # array for vector from one atom to another atom for all pairs of atoms
    rmag = np.sqrt(np.sum(np.square(r_array),-1)) # distance array between atoms for all pairs of atoms
    hit = np.less_equal(rmag,2*size)-np.identity(N) # if smaller than 2*size meaning these two atoms might hit each other
    hitlist = np.sort(np.nonzero(hit.flat)[0]).tolist() # change hit to a list
    for ij in hitlist: # i,j encoded as i*Natoms+j
        i, j = divmod(ij,N) # atom pair, i-th and j-th atoms, hit each other
        hitlist.remove(j*N+i) # remove j,i pair from list to avoid handling the collision twice
        if sum((p_a[i]-p_a[j])*(v_a[i]-v_a[j])) < 0 : # only handling collision if two atoms are approaching each other
            v_a[i], v_a[j] = vcollision(p_a[i], p_a[j], v_a[i], v_a[j]) # handle collision
            
    #find collisions between the atoms and the walls, and handle their elastic collisions
    for i in range(N):
        if abs(p_a[i][0]) >= container.length/2 - size and p_a[i][0]*v_a[i][0] > 0 :
            v_a[i][0] = - v_a[i][0]
            delta_p += 2*m*abs(v_a[i][0])
        if abs(p_a[i][1]) >= L - size and p_a[i][1]*v_a[i][1] > 0 :
            v_a[i][1] = - v_a[i][1]
            delta_p += 2*m*abs(v_a[i][1])
        if abs(p_a[i][2]) >= L - size and p_a[i][2]*v_a[i][2] > 0 :
            v_a[i][2] = - v_a[i][2]
            delta_p += 2*m*abs(v_a[i][2])

    #hw request
    if runs%1000 == 0:
        ke = 0 #total translational ke
        for i in range(N):
            ke+=0.5*m*(v_a[i][0]**2+v_a[i][1]**2+v_a[i][0]**2)
        SurfaceArea = 2*(container.length*container.width + container.width*container.height + container.height*container.length)
        Volume  = container.length*container.width*container.height
        print("temperature:     "+str(ke / (3*N*k/2)))
        print("pressure:            "+str( delta_p/((1000*dt)*SurfaceArea) ) )
        print("Volume:             "+str(Volume))
        print("PV:                      "+str(delta_p/((1000*dt)*SurfaceArea) * Volume))
        print("NkT:                   "+str(2/3*ke))
        print("p*(V**gamma)"+str(   delta_p/((1000*dt)*(SurfaceArea))   * (Volume)**(5/3)  ))
        delta_p = 0

        
#=====================================      stage 3   ========================================
container.length = 2*L
while stage == 3:
    t += dt
    runs+=1
    rate(10000)

    p_a += v_a*dt # calculate new positions for all atoms
    for i in range(N): atoms[i].pos = vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]) # to display atoms at new positions
    observation2.plot(data = np.sqrt(np.sum(np.square(v_a),-1)))
    
    ### find collisions between pairs of atoms, and handle their collisions
    r_array = p_a-p_a[:,np.newaxis] # array for vector from one atom to another atom for all pairs of atoms
    rmag = np.sqrt(np.sum(np.square(r_array),-1)) # distance array between atoms for all pairs of atoms
    hit = np.less_equal(rmag,2*size)-np.identity(N) # if smaller than 2*size meaning these two atoms might hit each other
    hitlist = np.sort(np.nonzero(hit.flat)[0]).tolist() # change hit to a list
    for ij in hitlist: # i,j encoded as i*Natoms+j
        i, j = divmod(ij,N) # atom pair, i-th and j-th atoms, hit each other
        hitlist.remove(j*N+i) # remove j,i pair from list to avoid handling the collision twice
        if sum((p_a[i]-p_a[j])*(v_a[i]-v_a[j])) < 0 : # only handling collision if two atoms are approaching each other
            v_a[i], v_a[j] = vcollision(p_a[i], p_a[j], v_a[i], v_a[j]) # handle collision
            
    #find collisions between the atoms and the walls, and handle their elastic collisions
    for i in range(N):
        if abs(p_a[i][0]) >= L - size and p_a[i][0]*v_a[i][0] > 0 :
            v_a[i][0] = - v_a[i][0]
            delta_p += 2*m*abs(v_a[i][0])
        if abs(p_a[i][1]) >= L - size and p_a[i][1]*v_a[i][1] > 0 :
            v_a[i][1] = - v_a[i][1]
            delta_p += 2*m*abs(v_a[i][1])
        if abs(p_a[i][2]) >= L - size and p_a[i][2]*v_a[i][2] > 0 :
            v_a[i][2] = - v_a[i][2]
            delta_p += 2*m*abs(v_a[i][2])

    #hw request 1
    if runs%1000 == 0:
        ke = 0 #total translational ke
        for i in range(N):
            ke+=0.5*m*(v_a[i][0]**2+v_a[i][1]**2+v_a[i][0]**2)
        SurfaceArea = 2*(container.length*container.width + container.width*container.height + container.height*container.length)
        Volume  = container.length*container.width*container.height
        print("temperature:     "+str(ke / (3*N*k/2)))
        print("pressure:            "+str( delta_p/((1000*dt)*SurfaceArea) ) )
        print("Volume:             "+str(Volume))
        print("PV:                      "+str(delta_p/((1000*dt)*SurfaceArea) * Volume))
        print("NkT:                   "+str(2/3*ke))
        print("p*(V**gamma)"+str(   delta_p/((1000*dt)*(SurfaceArea))   * (Volume)**(5/3)  ))
        delta_p = 0
