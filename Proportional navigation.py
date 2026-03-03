#%% Basic stuff
import numpy as np
import matplotlib.pyplot as plt
import random as random

# Everything is in m/s, even you are in m/s...
M=1234.8/3.6
dt=0.01 #Time step
T=1000. #If the missile hasn't reached target by then, then it knows not where it is at this time.
N=4 #Navigation constant (2<=N<=5 ideally)
a=random.uniform(2,2.5)
Vc=a*M #Closing velocity
xmpn=np.array([0.,0.]) #Missile position (propnav)
xmpc=np.array([0.,0.]) #Missile position (pure chase)
xt=np.array([random.uniform(-35400., 35400.),random.uniform(-35400., 35400.)])
while np.linalg.norm(xt-xmpn)<1000. and np.linalg.norm(xt-xmpn)>35400:
    xt=np.array([random.uniform(-35400., 35400.),random.uniform(-35400., 35400.)])
    if np.linalg.norm(xt-xmpn)>1000.0 and np.linalg.norm(xt-xmpn)>35400:
        break
vmpn=Vc*(xt-xmpn)/np.linalg.norm(xt-xmpn) #Missile (propnav) speed vector
vmpc=Vc*(xt-xmpc)/np.linalg.norm(xt-xmpc) #Missile (pure chase) speed vector
b=random.uniform(1.1,2.35)
Vt=b*M #Target velocity
theta=random.uniform(0,2*np.pi)


#%% The stuff that works
traj_mpn=[xmpn.copy()]
traj_t=[xt.copy()]
traj_mpc=[xmpc.copy()]
text1=""

for t in np.arange(0,T,dt):
    theta+=np.random.normal(0.,0.5*np.pi)*dt
    vt=Vt*np.array([np.cos(theta),np.sin(theta)]) #Target speed vector
    xt+=vt*dt
    
    lospn=xt-xmpn #lospn vector (propnav)
    lospn_norm=np.linalg.norm(lospn)
    lospn_u=lospn/lospn_norm #lospn unitary vector direction
    
    lospc=xt-xmpc #lospn vectpr (pure chase)
    lospc_norm=np.linalg.norm(lospc)
    lospc_u=lospc/lospc_norm
    
    vmpc=Vc*lospc_u
    xmpc+=vmpc*dt
    
    v_rel=vt-vmpn #Relative speed between missile and target
    lambda_p=np.cross(lospn,v_rel)/lospn_norm**2 #lospn rate
    
    v_u=vmpn/np.linalg.norm(vmpn)
    n=np.array([-v_u[1],v_u[0]]) #Normal vector
    a_n=(N*Vc*lambda_p)*n #Normal acceleration vector
    
    vmpn+=a_n*dt
    xmpn+=vmpn*dt
    
    traj_mpn.append(xmpn.copy())
    traj_t.append(xt.copy())
    traj_mpc.append(xmpc.copy())
    
    if lospn_norm<20. and t<=T:
        print(f"The missile knows where it is: Interception at {t:.2f}s")
        text1=f"The missile knows where it is: Interception at {t:.2f}s"
        break
    
if text1=="":
   print("The missile doesn't know where it is :(")
   text1="Timeout ("+str(T)+"s), the missile doesn't know where it is :("
#%% Plotting
traj_mpn=np.array(traj_mpn)
traj_t=np.array(traj_t)
traj_mpc=np.array(traj_mpc)

fig=plt.figure(figsize=(10,8))
plt.plot(traj_mpn[:,0],traj_mpn[:,1],label="Missile (propnav)")
plt.plot(traj_mpc[:,0],traj_mpc[:,1],label="Missile (pure chase)")
plt.plot(traj_t[:,0],traj_t[:,1],label="Target (semi-random trajectory)")
plt.plot(traj_t[-1,0],traj_t[-1,1],"rx",label="Target acquired")
plt.xlabel("x(m)")
plt.ylabel("y(m)")
plt.grid(True)
plt.axis("equal")

def result(x):
    return str("{:.2F}".format(x))

plt.suptitle("Interception: missiles at Mach "+result(a)+" and target at Mach "+result(b))
text2=result(np.linalg.norm(traj_mpn[-1]-traj_mpc[-1]))+"m (Hypothetical time saved: "+result(np.linalg.norm(traj_mpn[-1]-traj_mpc[-1])/Vc)+"s)"
# A terribly gross approximation, but it gets the job done
text3=result(np.linalg.norm(traj_mpn[0]-traj_t[0]))+"m"
plt.title(text1+"\n"+
          "Distance gain between propnav and pure chase: "+text2+"\n"+
          "Original distance between missiles and target: "+text3)
plt.legend()
plt.show()
