import numpy as np
import scipy as sp
from scipy.constants import gravitational_constant as G
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import orbitalTransforms as orb
orb.startup_plotting()

# Extracting elements from TLE
e,mean_motion,mean_anom,inc, arg_peri, right_asc = orb.parameters("Users/kelly/OneDrive/Desktop/A1/twoelement.txt")
# The directory was not working, I tried everything I could, so I had to unfortunately hard code it in
# However if the directory did work than this would be the correct code:
# e,mean_motion,mean_anom,inc, arg_peri, right_asc = orb.parameters("twoelement.txt")


# Cosntants
me = 5.9722e24 # mass of earth
mu = G*me # gravitation parameter
re = 6371e3 # radius of the earth

### PART 1 ###

## Keplerian elements from TLE ##

a = (mu/mean_motion**2)**(1/3)
b = a*np.sqrt(1-e**2)
h = np.sqrt(a*mu*(1-e**2)) 


rp = ((h**2)/mu)* (1/(1+e))
ra = ((h**2)/mu)* (1/(1-e))

peri_alt = rp-re # m
apo_alt = ra-re # m

period = (2*np.pi)/mean_motion

# Computing anomalies
t_sec = np.array(range(0,int(period)))
ta,ea,ma_vec = orb.compute_anomalies(t_sec,mean_anom,e,mean_motion)
print(mean_motion)

## Graphing Tangential Velocity ##
v_r, v_n = orb.compute_orbital_velocity(a, e, ta, mu)
print(max(v_n) - min(v_n))

plt.plot(t_sec, v_n)
plt.xlabel("Time (s)")
plt.ylabel("Tangential Velocity (m/s)")
plt.show()


# Converting to Perifocal and ECI frame
p, q, w, dp, dq, dw = orb.elements_to_perifocal(ta, a, e, mu, h)
x, y, z, dx, dy, dz = orb.perifocal_to_eci(p, q, w, dp, dq, dw, inc, right_asc, arg_peri)


## Plotting the position over time
fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(t_sec,p, label = "p")
ax.plot(t_sec, q, label = "q")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Position (m)")
ax.legend()
plt.show()

# Finding magnitude of radius from centre of the earth to the satellite
rmag_list = orb.magnitude_of_radius(p,q,w,t_sec)

# Finding index of perigee
perigee = min(rmag_list)
apogee = max(rmag_list)
peri_index = rmag_list.index(perigee)
apo_index = rmag_list.index(apogee)


## Plotting orbit in Perifocal frame ##
fig = plt.figure()
ax = fig.add_subplot(111)

circle1 = plt.Circle((0, 0), re, color='g', alpha = 0.7, label = "Earth")
ax.add_patch(circle1)

ax.plot(p[peri_index],q[peri_index], color="blue",marker="o", label = "Perigee")
ax.plot(p[apo_index],q[apo_index], color = "red", marker="o", label = "Apogee")
ax.plot(p,q, color = "black", label = "ISS Orbit")
ax.set_xlabel( "P (m)")
ax.set_ylabel("Q (m)")
ax.set_aspect('equal')
ax.legend(loc='center right', bbox_to_anchor=(1.5, 0.5))
plt.show()



## Plotting orbit in ECI frame ##
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


theta = np.linspace(0, 2 * np.pi, 100)
phi = np.linspace(0, np.pi, 100)
theta, phi = np.meshgrid(theta, phi)

x_earth = re * np.sin(phi) * np.cos(theta)
y_earth = re * np.sin(phi) * np.sin(theta)
z_earth = re * np.cos(phi)

# Plotting Earth
ax.plot_surface(x_earth/1000, y_earth/1000, z_earth/1000, alpha=0.5, label = "Earth")

# Plotting orbit
ax.plot(x/1000,y/1000,z/1000, color = 'red', linewidth = 2, label = "ISS Orbit")

ax.set_xlabel('X (km) ', labelpad = 15)
ax.set_ylabel('Y (km)', labelpad = 20)
ax.set_zlabel('Z (km)', labelpad = 15)


#Set an equal aspect ratio
ax.set_aspect('equal')

plt.show()

## PART 2 ##

# Iterating through different radii to find the corresponding index for the specific angle
angle = 0
i1 = peri_index
i2 = peri_index
while angle < (170):
    i1 = i1 - 1
    i2 = i2 + 1

    r1 = np.array([p[i1], q[i1], w[i1]]) - np.array([re, 0 , 0])
    r2 = np.array([p[i2], q[i2], w[i2]]) - np.array([re, 0 , 0])
    r1dotr2 = r1[0] *r2[0] + r1[1] *r2[1] +r1[2] *r2[2]
    r1mag = np.sqrt(r1[0]**2 +r1[1]**2 + r1[2]**2)
    r2mag = np.sqrt(r2[0]**2 +r2[1]**2 + r2[2]**2)
    angle = np.arccos(r1dotr2/(r1mag* r2mag))*180/np.pi 

# Removing the angle greater than 170 degrees generated in the last interation of the while loop
i1 = i1 + 1
i2 = i2 - 1
time = t_sec[i2] - t_sec[i1]


## Plotting visible orbit of satellite ## 

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

theta = np.linspace(0, 2 * np.pi, 100)
phi = np.linspace(0, np.pi, 100)
theta, phi = np.meshgrid(theta, phi)

x_earth = re * np.sin(phi) * np.cos(theta)
y_earth = re * np.sin(phi) * np.sin(theta)
z_earth = re * np.cos(phi)


# Plot the surface
ax.plot_surface(x_earth/1000, y_earth/1000, z_earth/1000, alpha = 0.7, label = "Earth")

# Plotting Orbits
ax.plot(x/1000,y/1000,z/1000, color = 'g', linewidth = 2,label='Satellite Orbit')
ax.plot(x[i1:i2+1]/1000,y[i1:i2+1]/1000,z[i1:i2+1]/1000, color = 'b', linewidth = 2,label='Visible Orbit')

ax.set_xlabel('X (km)', labelpad = 15)
ax.set_ylabel('Y (km)', labelpad = 15)
ax.set_zlabel('Z (km)', labelpad = 15)

#Set an equal aspect ratio
ax.set_aspect('equal')
plt.show()

### PART 3 ###

# Knowns
target_apo = 20000e3 + re
target_peri = 2500e3 + re

# Eccentricity of transfer orbit
e_transfer = (target_apo - perigee)/(target_apo + perigee)

# Calculating angular momentum for each orbit
h_intial = np.sqrt(2*mu)* np.sqrt((perigee*apogee)/(perigee +apogee))
h_trans = np.sqrt(2*mu)* np.sqrt((perigee*target_apo)/(perigee+target_apo))
h_target = np.sqrt(2*mu)*np.sqrt((target_apo * target_peri)/(target_apo + target_peri))

## Calculating velocities ##
va = abs(h_trans/perigee - h_intial/perigee) # tangential velocity at perigee kick
vb = abs(h_target/target_apo - h_trans/target_apo) # tangential velocity at apogee kick

## Calculating time for transfer ##
a_tran = (h_trans**2/mu) *(1/(1-e_transfer**2))
period_tran = ((2*np.pi)/np.sqrt(mu))*a_tran**(3/2)
transfer_time = period_tran/2



# Calculating values for transfer orbit

n_tran = np.pi/transfer_time # Mean motion of transfer orbit
t_tran =  np.array(range(0,int(period_tran))) # time array for transer orbit

ta_tran,ea_tran,ma_vec_tran = orb.compute_anomalies(t_tran,0,e_transfer,n_tran)
p_tran, q_tran, w_tran = orb.elements_to_perifocal(ta_tran, a_tran, e_transfer, mu, h_trans)[0:3]


## Graphing transfer orbit ##
fig = plt.figure()
ax = fig.add_subplot(111)

circle1 = plt.Circle((0, 0), re, color='g', alpha = 0.7, label = "Earth")
ax.add_patch(circle1)

ax.plot(p_tran,q_tran, color = 'black', linewidth = 2, label = "Transfer Orbit")

ax.set_xlabel('P (m)')
ax.set_ylabel('Q (m)')
ax.legend(loc='center right', bbox_to_anchor=(0.75, 1.05))


#Set an equal aspect ratio
ax.set_aspect('equal')
plt.show()

# Calculating values for target orbit
e_target = (target_apo - target_peri)/(target_apo + target_peri)
a_target = (h_target**2/mu) *(1/(1-e_target**2))
period_target = ((2*np.pi)/np.sqrt(mu))*a_target**(3/2)
target_time = period_target/2

n_target = np.pi/target_time
t_target =  np.array(range(0,int(period_target)))

ta_target,ea_target,ma_vec_target = orb.compute_anomalies(t_target,0,e_target,n_target)
p_target, q_target, w_target= orb.elements_to_perifocal(ta_target, a_target, e_target, mu, h_target)[0:3]

## Plotting the transfer orbit ##

# Midpoint of transfer orbit time, in order to graph half the transfer orbit
mid = int((len(t_tran))/2)

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(p[0:peri_index+1],q[0:peri_index+1], color = 'r', linewidth = 2, label = "ISS Orbit")
ax.plot(p_tran[0:mid],q_tran[0:mid], color = 'g', linewidth = 2, label = "Transfer Orbit")
ax.plot(p_target,q_target, color = 'b', linewidth = 2, label = "Target Orbit")

ax.set_xlabel('P (m)')
ax.set_ylabel('Q (m)')
ax.legend(loc='center right', bbox_to_anchor=(1.15, 1.05))

#Set an equal aspect ratio
ax.set_aspect('equal')
plt.show()