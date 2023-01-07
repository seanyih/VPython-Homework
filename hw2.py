from vpython import *
g = vec(0,-9.8,0)
size, m = 0.2, 1
L, k = 2, 150000
N=2
scene = canvas(width=500, height=500, center=vec(0, 0.2, 0), background=vec(0.5,0.5,0))
ceiling = box(length=5, height=0.005, width=5, color=color.blue)
# setting some shit
balls =[]
for i in range(5):
    balls.append(sphere(pos = vec(0.4*(i-2), -2 , 0), radius = size, color=color.red))
    balls[i].v = vec(0,0,0)
    balls[i].m = m

for i in range(N):
    balls[i].pos = vec(-0.44441+0.4*(i-2), -1.95, 0)

springs = []
for i in range(5):
    springs.append(cylinder(radius=0.005))# default pos = vec(0, 0, 0)
    springs[i].pos = (vec(0.4*(i-2),0, 0))
    
#plot 1
ieplot = graph(width = 400, align = 'top')
kefunc = gcurve(graph = ieplot,color=color.blue, width =2)
pefunc = gcurve(graph = ieplot, color=color.red, width =2)

#plot 2
avgke = 0
avgpe = 0
avgeplot = graph(width = 400, align = 'top')
avgkefunc = gcurve(graph =avgeplot,color=color.blue, width =2)
avgpefunc = gcurve(graph = avgeplot, color=color.red, width =2)
    
t=0
dt = 0.0001
while True:
    rate(5000)
    t += dt
    ke=0
    pe=0
    for i in range(5):
        springs[i].axis = balls[i].pos - springs[i].pos #spring extended from endpoint to ball
        spring_force = - k * (mag(springs[i].axis) - L) * springs[i].axis.norm() # to get spring force vector
        balls[i].a = g + spring_force / m # ball acceleration = - g in y + spring force /m
        balls[i].v += balls[i].a *dt
        balls[i].pos += balls[i].v *dt
        if (i!=4 and mag(balls[i].pos - balls[i+1].pos) <= 0.4 and dot(balls[i].pos-balls[i+1].pos, balls[i].v-balls[i+1].v) <= 0) :
            temp = balls[i].v
            balls[i].v = balls[i+1].v
            balls[i+1].v = temp
        #graph
        ke+=0.5*m*mag(balls[i].v)
        pe+=m*(balls[i].pos.y+2)*9.8
    kefunc.plot(pos = (t,ke))
    pefunc.plot(pos = (t,pe))
    avgke = (avgke*(t-dt) + ke*dt)/t
    avgpe = (avgpe*(t-dt) + pe*dt)/t
    avgkefunc.plot(pos = (t,avgke))
    avgpefunc.plot(pos = (t,avgpe))
    
        
