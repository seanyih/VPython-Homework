import numpy as np
from vpython import *
A, N= 0.10, 50
m, k, d = 0.1, 10.0, 0.4 

#no animation
class Obj:
    pass

balls = [Obj for i in range(N)]
springs = [Obj for i in range(N - 1)]


#c = curve([vector(i*d, 1.0, 0) for i in range(N)], color=color.black)

t, dt = 0, 0.0003

g = graph(width=1200, height=600, align='left', xtitle="k",
          ytitle="omega")
p = gcurve(graph=g, color=color.blue, width=2)

for n in arange(1,N/2-1,1.0):
    Unit_K= 2 * pi/(N*d)
    Wavevector = n * Unit_K
    phase = Wavevector * arange(N) * d
    ball_pos, ball_orig, ball_v, spring_len = np.arange(N)*d + A*np.sin(phase), np.arange(N)*d, np.zeros(N), np.ones(N)*d
    waves = 0
    t = 0
    while waves < 10:
        #rate(1000)
        t += dt
        #ball_pos[0] = A * sin(omega * t ) #4
        spring_len[:-1] = ball_pos[1:] - ball_pos[:-1]
        ball_v[1:] += (spring_len[1:] - spring_len[:-1]) * k / m * dt
        #periodically boundary condition
        spring_len[N - 1] = ball_pos[0] - ball_pos[N - 1] + N * d
        ball_v[0] += (spring_len[0] - spring_len[N - 1]) * k / m * dt
        #
        if ball_pos[0] * (ball_pos[0] + ball_v[0] * dt) < 0 and t > 5 * dt:
            waves += 0.5
        ball_pos += ball_v*dt
        #for i in range(N): balls[i].pos.x = ball_pos[i] #3
        #for i in range(N-1): #3
        #    springs[i].pos = balls[i].pos #3
        #    springs[i].axis = balls[i+1].pos - balls[i].pos #3
        
        #ball_disp = ball_pos - ball_orig
        
        
    T = t / waves #average period
    p.plot(Wavevector, 2.0 * pi / T)
