from vpython import *
import numpy as np

fd = 120  # 120Hz
T = 1/fd
w = 2*pi *fd
# (Your Parameters here)
#Given: 
C = 20E-6
L = 0.2
R = 30

#init:
i = 0
v = 0
E = 0

Q = 0
i_old = 0
v_old = 0
E0 = 0

task = 0

phi_i = 0
phi_v = 0

t = 0
dt = 1.0 / (fd * 1000)  # 1000 simulation points per cycle

scene1 = graph(align='left', xtitle='t', ytitle='i (A) blue, v (100V) red,',
               background=vector(0.2, 0.6, 0.2))
scene2 = graph(align='left', xtitle='t', ytitle='Energy (J)', background=vector(0.2, 0.6, 0.2))


i_t = gcurve(color=color.blue, graph=scene1)
v_t = gcurve(color=color.red, graph=scene1)
E_t = gcurve(color=color.red, graph=scene2)

# (Your program here)
while t <= 20 * T:
    if t < 12 * T: #Power ON
        v = 36 * sin(w * t) #power
        i_theoretical = 0.40156 * sin(w* t - 70 / 180 * pi)
        
        if v * v_old <= 0 and v_old - v > 0: 
            phi_v = w*t
        
        if task == 0 and t >= 9 * T:
                print("i                               =", abs(i))
                print("i_theoretical value =", abs(i_theoretical))
                print("===============================")
                phi = (phi_v - phi_i)
                phi_theoretical = np.arctan(-(w * L - 1 / w / C) / R) #=Arg(Z)
                print("phi                                =", phi)
                print("phi_theoretical value =", phi_theoretical)
                print("===============================")
                task = 1
                
    else: #Power OFF
        if task == 1:
            E0 = E #E0 = E(t=12T)
            task = 2
        elif task == 2:
            if E <= E0 / 10:
                print("Decay time:", t)
                task = -1


    i = (v + L * i_old / dt - Q / C) / (R + dt / C + L / dt)

    if i_old * i <= 0 and i_old - i > 0: 
        phi_i = w*t

    Q += i * dt
    i_old = i
    v_old = v
    t += dt
    E = Q * Q / C / 2 + L * i * i / 2

    i_t.plot(pos=(t, i))
    v_t.plot(pos=(t, v / 100))
    E_t.plot(pos=(t, E))
