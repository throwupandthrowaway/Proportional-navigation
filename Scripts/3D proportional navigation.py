#%% Basic stuff
import numpy as np
import random
import time


# Everything is in m/s, even you are in m/s...
M=1234.8/3.6
dt=0.001 #Time step
T=1000. #If the missile hasn't reached target by then, then it knows not where it is at this time.
N=3 #Navigation constant (2<=N<=5 ideally)
a=random.uniform(2,2.5)
Vc=a*M #Closing velocity
xmpn=np.array([0.,0.,0.]) #Missile position (propnav)
xmpc=np.array([0.,0.,0.]) #Missile position (pure chase)
xt=np.array([random.uniform(-35400., 35400.),random.uniform(-35400., 35400.),random.uniform(-35400., 35400.)])
while np.linalg.norm(xt)<2000. or np.linalg.norm(xt)>35400: #To make sure that it's nor too close, nor too far from the missile
    xt=np.array([random.uniform(-35400., 35400.),random.uniform(-35400., 35400.),random.uniform(-35400., 35400.)])
    if np.linalg.norm(xt)>2000.0 and np.linalg.norm(xt)<35400:
        break
b=random.uniform(1.1,2.35)
Vt=b*M #Target velocity

# Random initial target direction in 3D
u0=np.random.randn(3)
u0=u0/np.linalg.norm(u0)
vt=Vt*u0

# Initial missile velocity: aimed at target
vmpn=Vc*(xt-xmpn)/np.linalg.norm(xt-xmpn)
vmpc=Vc*(xt-xmpc)/np.linalg.norm(xt-xmpc)

#%% Big brain time
start_time=time.time()
traj_mpn=[xmpn.copy()]
traj_t=[xt.copy()]
traj_mpc=[xmpc.copy()]
text1=""

seq=np.arange(0,T,dt)
for t in seq:
    vt_u=vt/np.linalg.norm(vt) #Unitary speed vector (target)
    rnd=np.random.randn(3) #Random diretion for acceleration
    rnd=rnd-np.dot(rnd,vt_u)*vt_u
    rnd_norm=np.linalg.norm(rnd)
    if rnd_norm>1e-12:
        rnd_u=rnd/rnd_norm
    else:
        rnd_u=np.zeros(3)
    a_t_mag=abs(np.random.normal(0.,0.5*np.pi))*Vt*dt
    at=a_t_mag*rnd_u
    dtheta=np.random.normal(0.,0.5*np.pi)*dt
    vt=vt+at*dt
    vt=Vt*vt/np.linalg.norm(vt)
    xt+=vt*dt
    
    lospn=xt-xmpn #Line Of Sight vector (propnav)
    lospn_norm=np.linalg.norm(lospn)
    lospn_u=lospn/lospn_norm #lospn unitary vector direction
    
    lospc=xt-xmpc #lospn vectpr (pure chase)
    lospc_norm=np.linalg.norm(lospc)
    lospc_u=lospc/lospc_norm
    
    vmpc=Vc*lospc_u #Pure chase missile speed vector
    xmpc+=vmpc*dt
    
    v_rel=vt-vmpn #Relative speed between missile and target
    omega_los=np.cross(lospn,v_rel)/(lospn_norm**2)
    vmpn_u=vmpn/np.linalg.norm(vmpn)
    a_n=N*Vc*np.cross(omega_los,vmpn_u)
    vmpn+=a_n*dt
    xmpn+=vmpn*dt
    lambda_p=(lospn[0]*v_rel[1] - lospn[1]*v_rel[0])/lospn_norm**2 #lospn rate
        
    for traj, var in zip([traj_mpn,traj_t,traj_mpc],[xmpn, xt, xmpc]):
        traj.append(var.copy())
    
    if lospn_norm<20. and t<=T:
        text1=f"The missile knows where it is: Interception at {t:.2f}s"
        break
    
if text1=="":
   print("The missile doesn't know where it is :(")
   text1="Timeout ("+str(T)+"s), the missile doesn't know where it is :("
   
ex_time=time.time()-start_time
print("Loop execution time: %s seconds" % "{:.4f}".format(ex_time))
if ex_time>10:
    print("That's a long thinking time mister computer")
#%% Plotting
import matplotlib.pyplot as plt
traj_mpn, traj_t, traj_mpc = map(np.array, [traj_mpn, traj_t, traj_mpc])

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

ax.plot(
    [traj_mpn[0, 0]], [traj_mpn[0, 1]], [traj_mpn[0, 2]],
    "b.", label="Starting point (missiles)"
)
ax.plot(
    [traj_t[0, 0]], [traj_t[0, 1]], [traj_t[0, 2]],
    "g.", label="Starting point (target)"
)

ax.plot(
    traj_mpn[:, 0], traj_mpn[:, 1], traj_mpn[:, 2],
    label="Missile (propnav)"
)
ax.plot(
    traj_mpc[:, 0], traj_mpc[:, 1], traj_mpc[:, 2],
    label="Missile (pure chase)"
)
ax.plot(
    traj_t[:, 0], traj_t[:, 1], traj_t[:, 2],
    label="Target (semi-random trajectory)"
)

ax.plot(
    [traj_t[-1, 0]], [traj_t[-1, 1]], [traj_t[-1, 2]],
    "rD", label="Final target position"
)

ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_zlabel("z (m)")
ax.grid(True)

plt.suptitle(
    "Interception: missiles at Mach "
    + str("{:.2f}".format(a))
    + " and target at Mach "
    + str("{:.2f}".format(b))
)

ax.legend()

#%% Result analysis
distance_gain = np.linalg.norm(traj_mpn[-1] - traj_mpc[-1])

# Final distances to target
distance_mpn_t = np.linalg.norm(traj_mpn[-1] - traj_t[-1])
distance_mpc_t = np.linalg.norm(traj_mpc[-1] - traj_t[-1])

print("Minimum distance gain to be more efficient than pure chase: "
      + str("{:.2f}".format(2 * (Vc - Vt))) + " m")

if distance_mpc_t > 1e-12:
    ratio = distance_mpn_t / distance_mpc_t   # good if ratio < 1
    ratio_text = "{:.3f}".format(ratio)
else:
    ratio = np.inf
    ratio_text = "undefined"

if distance_gain > 2 * (Vc - Vt) and distance_mpn_t < distance_mpc_t:
    print("Distance gain sufficient.")
elif distance_gain < 2 * (Vc - Vt) and distance_mpn_t < distance_mpc_t:
    print("Propnav is better, but the distance gain is small.")
elif distance_mpn_t > distance_mpc_t:
    print("Pure chase missile is closer to target than propnav missile "
          + "(ratio distance_mpn_t/distance_mpc_t: " + ratio_text + ")")
else:
    print("Comparable end-state between propnav and pure chase.")

text2 = str("{:.2f}".format(distance_gain)) + " m"
text3 = str("{:.2f}".format(np.linalg.norm(traj_t[0] - traj_mpn[0]))) + " m"

plt.title(
    text1 + "\n"
    + "Distance gain between propnav and pure chase: " + text2 + "\n"
    + "Initial distance between missiles and target: " + text3
)

plt.legend()
plt.show()

# Path lengths
diffs_pn = traj_mpn[1:] - traj_mpn[:-1]
length_pn = np.sum(np.linalg.norm(diffs_pn, axis=1))

diffs_pc = traj_mpc[1:] - traj_mpc[:-1]
length_pc = np.sum(np.linalg.norm(diffs_pc, axis=1))

path_difference = length_pc - length_pn

print("Path length (propnav): " + str("{:.2f}".format(length_pn)) + " m")
print("Path length (pure chase): " + str("{:.2f}".format(length_pc)) + " m")
print("Path difference between pure chase and propnav: "
      + str("{:.2f}".format(abs(path_difference))) + " m")

if path_difference > 0:
    print("Longer path for pure chase")
elif path_difference < 0:
    print("Longer path for propnav")
else:
    print("Same path length")

# More meaningful performance metric
target_distance_gain = distance_mpc_t - distance_mpn_t
print("Gain relative to target: " + str("{:.2f}".format(target_distance_gain)) + " m")

plt.show()
