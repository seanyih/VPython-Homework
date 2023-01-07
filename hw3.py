from vpython import *

G=6.6743E-11

mass = {'earth': 5.97E24, 'moon': 7.36E22, 'sun':1.99E30}
radius = {'earth': 6.371E6*10, 'moon': 1.317E6*10, 'sun':6.95E8*10} #10 times larger for better view 
earth_orbit = {'r': 1.495E11, 'v': 2.9783E4}
moon_orbit = {'r': 3.84E8, 'v': 1.022E3}
theta = 5.145*pi/180.0



scene = canvas(width=800, height=800, center = vec(0,400,0), background=vec(0.5,0.5,0), range = 4E8)
floor = box(length=30, height=0.01, width=10, color=color.blue)


earth_com = -moon_orbit['r']*mass['moon']/(mass['earth']+mass['moon'])
moon_com = moon_orbit['r']*mass['earth']/(mass['earth']+mass['moon'])
sun = sphere(radius=radius['sun'],color=color.red)
sun.pos= vector(0,0,0)
sun.m =mass['sun']

earth = sphere(radius=radius['earth'], texture={'file':textures.earth})
earth.pos = vector(earth_com*cos(theta),-earth_com*sin(theta),0)+vector(earth_orbit['r'], 0, 0)
earth.m = mass['earth']
earth.v = vector(0,0,moon_orbit['v']*mass['moon']/mass['earth'])+vector(0, 0, -earth_orbit['v'])


moon = sphere(radius=radius['moon'])
moon.pos = vector(moon_com*cos(theta),-moon_com*sin(theta),0)+vector(earth_orbit['r'], 0, 0)
moon.m = mass['moon']
moon.v = vector(0,0,-moon_orbit['v'])+vector(0, 0, -earth_orbit['v'])

scene.light = []
local_light(pos=vector(0,0,0))


c1 = arrow(color = color.white, shaftwidth=1000000) #normal vec of lunar orbit
c2 = arrow(color = color.yellow, shaftwidth=500000) #normal vec of earth's orbit
c2.axis = vector(0,moon_orbit['r']*0.5,0)
t=0
dt = 3600*6 #a half day
record = 0
recorded = False
while True:
    t+=dt
    rate(1000)
    
    moon.a = -norm(moon.pos)*G*mass['sun']/(mag2(moon.pos))+norm(earth.pos-moon.pos)*G*earth.m/(mag2(moon.pos-earth.pos))
    earth.a = -norm(earth.pos)*G*mass['sun']/(mag2(earth.pos))+norm(moon.pos-earth.pos)*G*moon.m/(mag2(moon.pos-earth.pos))

    moon.v = moon.v + moon.a * dt
    moon.pos = moon.pos + moon.v * dt
    
    earth.v = earth.v + earth.a * dt
    earth.pos = earth.pos + earth.v * dt

    scene.center = earth.pos
    c1.pos = earth.pos
    c2.pos = earth.pos
    c1.axis = 0.3*moon_orbit['r']*cross(norm(moon.pos-earth.pos),norm(moon.v-earth.v))

    
    if(c1.axis.x>0 and c1.axis.z < 0 and recorded == False): #quadrant IV --> record ONCE
        recorded = True
        if(record != 0 and (t-record)/60/60/24 > 365):
            print(str((t-record)/60/60/24)+" days")
        record = t
    elif(c1.axis.x<0): # quadrant II, III --> reset
        recorded = False

    
    
    
    
    

