from vpython import *
import numpy as np

#Given params. Unit: m
R = 0.12
r = 0.06
h = 0.1

#divide into m "ring" areas
m = 1000
#divide into "ds"
n = 720 #0.5 deg per segment
dtheta = 2*pi/n

print("Please wait patiently!")
#=====================method 1: divide the smaller ring=======================
dx = r / m
flux = 0
#For every dA on r
for i in range(m):
    pos_p = vec(dx * i , 0, h)
    #Sample point: any point on the ring. For the ring pair is symmetric wrt x=y=0.
    B = vec(0, 0, 0)

    #For every ds on R
    for j in range(n):
        theta = dtheta * j
        ds = R * dtheta * vec(-sin(theta), cos(theta), 0) #ds = R*dtheta, and the direction is parallel to d(pos_s)
        pos_s = R * vec(cos(theta), sin(theta), 0)
        dB = 1E-7 * cross(ds, norm(pos_p - pos_s)) / mag2(pos_p - pos_s) #Biot-Savart Law. i=1 ignored
        B += dB

    flux += 2 * pi * (dx * i)  * dx * dot(B, vec(0, 0, 1)) #2pi*x*dx = dA, x = dx*i ,vec(0,0,1) = normal vector
    if (i+1)%200 == 0:
        print("Yes, it's still running. Do not close!")
    
Flux1 = flux #M=flux



#=====================method 2: divide the larger ring=======================
dx = R / m
flux = 0
#For every dA on R
for i in range(m):
    pos_p = vec(dx * i, 0, 0)
    #Sample point: any point on the ring. For the ring pair is symmetric wrt x=y=0.
    B = vec(0, 0, 0)

    #For every ds on r
    for j in range(n):
        theta = dtheta * j
        ds = r * dtheta * vec(-sin(theta), cos(theta), 0) #ds = R*dtheta, and the direction is parallel to d(pos_s)
        pos_s = vec(r*cos(theta), r*sin(theta), h)
        dB = 1E-7 * cross(ds, norm(pos_p - pos_s)) / mag2(pos_p - pos_s) #Biot-Savart Law. i=1 ignored
        B += dB

    flux += 2 * pi * (dx * i) * dx * dot(B, vec(0, 0, 1)) #2pi*x*dx = dA, x = dx*i ,vec(0,0,1) = normal vector
    if (i+1)%200 == 0:
        print("Yes, it's still running. Do not close!")

#result printing
print("=======================================")
print("M by method 1:", Flux1) 
print("M by method 2:", flux) #M=flux
