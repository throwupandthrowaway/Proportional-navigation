#%% Basic stuff
import numpy as np
import matplotlib.pyplot as plt
import random
import time


# Everything is in m/s, even you are in m/s...
M=1234.8/3.6
dt=0.001 #Time step
T=1000. #If the missile hasn't reached target by then, then it knows not where it is at this time.
N=3 #Navigation constant (2<=N<=5 ideally)
a=random.uniform(2,2.5)
Vc=a*M #Closing velocity
xmpn=np.array([0.,0.]) #Missile position (propnav)
xmpc=np.array([0.,0.]) #Missile position (pure chase)
xt=np.array([random.uniform(-35400., 35400.),random.uniform(-35400., 35400.)])
while np.linalg.norm(xt)<2000. or np.linalg.norm(xt)>35400: #To make sure that it's nor too close, nor too far from the missile
    xt=np.array([random.uniform(-35400., 35400.),random.uniform(-35400., 35400.)])
    if np.linalg.norm(xt)>2000.0 and np.linalg.norm(xt)<35400:
        break
vmpn=Vc*(xt-xmpn)/np.linalg.norm(xt-xmpn) #Missile (propnav) speed vector
vmpc=Vc*(xt-xmpc)/np.linalg.norm(xt-xmpc) #Missile (pure chase) speed vector
b=random.uniform(1.1,2.35)
Vt=b*M #Target velocity
theta=random.uniform(0,2*np.pi)


#%% Big brain time
start_time=time.time()
traj_mpn=[xmpn.copy()]
traj_t=[xt.copy()]
traj_mpc=[xmpc.copy()]
text1=""

seq=np.arange(0,T,dt)
for t in seq:
    dtheta=np.random.normal(0.,0.5*np.pi)*dt
    theta=theta+dtheta
    vt=Vt*np.array([np.cos(theta),np.sin(theta)]) #Target speed vector
    xt+=vt*dt
    at=vt-Vt*np.array([np.cos(theta-dtheta),np.sin(theta-dtheta)]) #Target acceleration vector
    at_u=at/np.linalg.norm(at) #Unitary vector
    if np.linalg.norm(at) > 0:
        at_u=at/np.linalg.norm(at)
    else:
        at_u=np.array([0.,0.])
    
    lospn=xt-xmpn #Line Of Sight vector (propnav)
    lospn_norm=np.linalg.norm(lospn)
    lospn_u=lospn/lospn_norm #lospn unitary vector direction
    
    lospc=xt-xmpc #lospn vectpr (pure chase)
    lospc_norm=np.linalg.norm(lospc)
    lospc_u=lospc/lospc_norm
    
    vmpc=Vc*lospc_u #Pure chase missile speed vector
    xmpc+=vmpc*dt
    
    v_rel=vt-vmpn #Relative speed between missile and target
    lambda_p=(lospn[0]*v_rel[1] - lospn[1]*v_rel[0])/lospn_norm**2 #lospn rate
    
    v_u=vmpn/np.linalg.norm(vmpn)
    nm=np.array([-v_u[1],v_u[0]]) #Normal vector to missile speed
    nt=np.array([-at_u[1],at_u[0]]) #Normal vector to target acceleration
    a_n=(N*Vc*lambda_p)*nm+N*nt/2 #Normal acceleration vector
    
    vmpn+=a_n*dt #Propnav missile speed vector
    xmpn+=vmpn*dt
    
    for traj, var in zip([traj_mpn,traj_t,traj_mpc],[xmpn, xt, xmpc]):
        traj.append(var.copy())
    
    if lospn_norm<20. and t<=T:
        text1=f"The missile knows where it is: Interception at {t:.2f}s"
        break
    
if text1=="":
   print("The missile doesn't know where it is :(")
   text1="Timeout ("+str(T)+"s), the missile doesn't know where it is :("
#%% Plotting
traj_mpn,traj_t,traj_mpc=map(np.array,[traj_mpn,traj_t,traj_mpc])

fig=plt.figure(figsize=(10,8))
plt.plot(traj_mpn[0,0],traj_mpn[0,1],"b.",label="Starting point (missiles)")
plt.plot(traj_t[0,0],traj_t[0,1],"g.",label="Starting point (target)")
plt.plot(traj_mpn[:,0],traj_mpn[:,1],label="Missile (propnav)")
plt.plot(traj_mpc[:,0],traj_mpc[:,1],label="Missile (pure chase)")
plt.plot(traj_t[:,0],traj_t[:,1],label="Target (semi-random trajectory)")
plt.plot(traj_t[-1,0],traj_t[-1,1],"rD",label="Interception point")
plt.xlabel("x(m)")
plt.ylabel("y(m)")
plt.grid(True)
plt.axis("equal")

plt.suptitle("Interception: missiles at Mach "+str("{:.2f}".format(a))+" and target at Mach "+str("{:.2f}".format(b)))

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
    print("Pure chase missile closer to target than propnav missile (ratio "+ratio+")")
elif distance_gain<2*Vc and distance_mpn_t>distance_mpc_t:
    print("Proportional navigation absolutely useless wtf")
text2=str("{:.2f}".format(distance_gain))+"m"
text3=str("{:.2f}".format(np.linalg.norm(traj_t[0])))+"m"
plt.title(text1+"\n"+
          "Distance gain between propnav and pure chase: "+text2+"\n"+
          "Original distance between missiles and target: "+text3)
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
    
ex_time=time.time()-start_time
print("Execution time: %s seconds" % "{:.4f}".format(ex_time))
if ex_time>10:
    print("That's a long thinking time mister computer")
