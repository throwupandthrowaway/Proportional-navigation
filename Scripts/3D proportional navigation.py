#%% Basic stuff
import numpy as np
import random
import time


# Everything is in m/s
M=1234.8/3.6
g=9.81
dt=0.001 #Time step
T=1000. #If the missile hasn't reached target by then, then it knows not where it is at this time.
N=3 #Navigation constant (2<=N<=5 from what I read online)
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

# Initial missile velocity: collision triangle (resolution of ||\vec{r}+\vec{v_c}t||=||\vec{v_t}t||)
r=xt-xmpn
A=Vc**2-Vt**2
B=-2*np.dot(r,vt)
C=-np.linalg.norm(r)**2

delta=B**2-4*A*C

if delta>=0 and abs(A)>1e-6:
    tau1=(-B-np.sqrt(delta))/(2.*A)
    tau2=(-B+np.sqrt(delta))/(2.*A)
    taus=[t for t in [tau1, tau2] if t>0.]
    tau=min(taus) if taus else np.linalg.norm(r)/Vc
else:
    tau=np.linalg.norm(r)/Vc

hyp=xt+vt*tau #Hypothetical interception point
vmpn=Vc*(hyp-xmpn)/np.linalg.norm(hyp-xmpn)

#Initial missile velocity: aimed at target
vmpc=Vc*(xt-xmpc)/np.linalg.norm(xt-xmpc)

#%% Big brain time
start_time=time.time()
traj_mpn=[xmpn.copy()]
traj_t=[xt.copy()]
traj_mpc=[xmpc.copy()]
text1=""

maneuver_dir=np.random.randn(3)
maneuver_dir-=np.dot(maneuver_dir,u0)*u0 #Keep normal to initial velocity
maneuver_dir/=np.linalg.norm(maneuver_dir)
maneuver_mag=np.random.uniform(0.,12.*g)
maneuver_timer=0.
maneuver_duration = np.random.uniform(1.,4.) #Hold maneuver 1–4 s

seq=np.arange(0,T,dt)
for t in seq:
    maneuver_timer+=dt
    if maneuver_timer>=maneuver_duration:
        maneuver_timer=0.
        maneuver_duration=np.random.uniform(1.5,4.)
        vt_u=vt/np.linalg.norm(vt)
        maneuver_dir=np.random.randn(3)
        maneuver_dir-=np.dot(maneuver_dir,vt_u)*vt_u  #Strictly normal to velocity
        n=np.linalg.norm(maneuver_dir)
        maneuver_dir=maneuver_dir/n if n > 1.e-12 else np.zeros(3)
        maneuver_mag=np.random.uniform(3.*g,9.*g)

    at=maneuver_mag*maneuver_dir    
    vt+=at*dt
    vt=Vt*vt/np.linalg.norm(vt)
    xt+=vt*dt
    
    lospn=xt-xmpn #Line Of Sight vector (propnav)
    lospn_norm=np.linalg.norm(lospn)
    lospn_u=lospn/lospn_norm #lospn unitary vector direction
    
    lospc=xt-xmpc #los vector (pure chase)
    lospc_norm=np.linalg.norm(lospc)
    lospc_u=lospc/lospc_norm
    
    vmpc=Vc*lospc_u #Pure chase missile speed vector
    xmpc+=vmpc*dt
    
    v_rel=vt-vmpn #Relative speed between missile and target
    omega_los=np.cross(lospn,v_rel)/(lospn_norm**2)
    vmpn_u=vmpn/np.linalg.norm(vmpn)
    at_normal_to_los=at-np.dot(at,lospn_u)*lospn_u  #Project out LOS component
    nt=at_normal_to_los/np.linalg.norm(at_normal_to_los) if np.linalg.norm(at_normal_to_los) > 1e-12 else np.zeros(3) #Target acceleration normal to LOS
    a_n=N*Vc*np.cross(omega_los,vmpn_u)+N*nt/2
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
    print("That's a long thinking time")
#%% Plotting
import matplotlib.pyplot as plt
traj_mpn,traj_t,traj_mpc=map(np.array,[traj_mpn,traj_t,traj_mpc])

fig=plt.figure(figsize=(10,8))
ax=fig.add_subplot(111,projection="3d")

ax.plot([traj_mpn[0, 0]],[traj_mpn[0,1]],[traj_mpn[0,2]],
        "b.",label="Starting point (missiles)")
ax.plot([traj_t[0, 0]],[traj_t[0,1]],[traj_t[0,2]],
        "g.",label="Starting point (target)")

ax.plot(traj_mpn[:,0],traj_mpn[:,1],traj_mpn[:,2],
        label="Missile (propnav)")
ax.plot(traj_mpc[:,0],traj_mpc[:,1],traj_mpc[:,2],
        label="Missile (pure chase)")
ax.plot(traj_t[:,0],traj_t[:,1],traj_t[:,2],
        label="Target (semi-random trajectory)")

ax.plot([traj_t[-1,0]],[traj_t[-1,1]],[traj_t[-1,2]],
        "rD",label="Interception point")

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
distance_gain=np.linalg.norm(traj_mpn[-1]-traj_mpc[-1])
distance_mpn_t=np.linalg.norm(traj_mpn[-1]-traj_t[-1]) #Distance between propnav missile and target
distance_mpc_t=np.linalg.norm(traj_mpc[-1]-traj_t[-1]) #Distance between pure chase missile and target

print("Minimum distance gain to be more efficient than pure chase: "+str("{:.2f}".format(2*Vc))+"m")
from fractions import Fraction
ratio=""
if int(distance_mpc_t)!=0:
    ratio=Fraction(int(distance_mpn_t),int(distance_mpc_t)) #Good if ratio<1
    ratio=str(ratio.numerator)+"/"+str(ratio.denominator)

if distance_gain>2*Vc and distance_mpn_t<distance_mpc_t:
    print("Distance gain sufficient.")
elif distance_gain<2*Vc: #It would mean that at this distance, pure chase can catch up in less than 2 seconds
    print("Insufficient distance gain, pure chase is almost as effective.")
elif distance_mpn_t>distance_mpc_t:
    print("Pure chase missile closer to target than propnav missile (ratio distance_mpn_t/distance_mpc_t: "+ratio+")")
elif distance_gain<2*Vc and distance_mpn_t>distance_mpc_t:
    print("Proportional navigation absolutely useless wtf")
text2=str("{:.2f}".format(distance_gain))+"m"
text3=str("{:.2f}".format(np.linalg.norm(traj_t[0])))+"m"
plt.title(
    text1 + "\n"
    + "Distance gain between propnav and pure chase: " + text2 + "\n"
    + "Initial distance between missiles and target: " + text3
)

plt.legend()
plt.show()

length_pn=0
diffs_pn=traj_mpn[1:]-traj_mpn[:-1]
length_pn=np.sum(np.linalg.norm(diffs_pn,axis=1))
length_pc=0
diffs_pc=traj_mpc[1:]-traj_mpc[:-1]
length_pc=np.sum(np.linalg.norm(diffs_pc,axis=1))
path_difference=abs(length_pc-length_pn)
print("Path difference between propnav and pure chase: "+str("{:.2f}".format(path_difference))+"m")
if path_difference>0:
    print("Longer path for pure chase")
else:
    print("Longer path for propnav")
    
