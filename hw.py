from vpython import *
g=9.8 # g = 9.8 m/s^2
size = 0.25 # ball radius = 0.25 m
height = 15.0 # ball center initial height = 15 m
C_drag = 0.9
t=0 # for graphing
touchcount = 0 #count how many time did the ball touch the ground
prevx = -15 #pos record for counting total distance
prevy = size #pos record for counting total distance
arcl = 0 #for counting total distance
highest = 0

#set bg
scene = canvas(width=500, height=500, center =vec(-5,height/2,0), background=vec(0.5,0.5,0))
floor = box(length=30, height=0.01, width=10, color=color.blue)
ball = sphere(radius = size, color=color.red, make_trail = True)

#plot
vtplot = graph(width = 400, align = 'right')
func = gcurve(graph = vtplot,color=color.blue, width =4)

#set arrow
c1 = arrow(color = color.yellow, shaftwidth=0.2)
#initial conditions
ball.pos = vec(-15, size, 0)
ball.v = vec(20*cos(pi/4), 20*sin(pi/4), 0)

dt = 0.001 # time step
while touchcount < 3: # until the ball hit the ground 3 times
    rate(1000) # run 1000 times per real second
    ball.v += vec(0, -g, 0) * dt - C_drag*ball.v*dt
    ball.pos += ball.v*dt
    #arrow
    c1.pos = ball.pos
    c1.axis = 0.5*ball.v
    #graph
    func.plot(pos = (t,sqrt(ball.v.x**2+ball.v.y**2)))
    t+=1
    #arc length counting
    arcl += sqrt((ball.pos.x-prevx)**2+(ball.pos.y-prevy)**2)
    prevx = ball.pos.x
    prevy = ball.pos.y
    if ball.pos.y > highest:
        highest = ball.pos.y
    if ball.pos.y <= size and ball.v.y <= 0:
        touchcount += 1
        ball.v.y = -ball.v.y

msg1 = text(text = "total distance  " +str(arcl) + "m", pos = vec(5, 15, 0))
msg2 = text(text = "displacement  " +str(ball.pos.x+15) + "m", pos = vec(5, 10, 0))
msg3 = text(text = "largest height  " +str(highest) + "m", pos = vec(5, 5, 0))




