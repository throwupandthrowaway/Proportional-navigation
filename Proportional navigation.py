#%% Basic stuff
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random as random

dt=0.05 #Time step
T=500. #If the missile hasn't reached target by then, then it knows not where it is at this time.
N=3 #Navigation constant (2<=N<=5 ideally)
Vc=40 #Clospning velocity
xmpn=np.array([0.,0.]) #Missile position (propnav)
xmpc=np.array([0.,0.]) #Missile position (pure chase)
xt=np.array([random.uniform(-1000.,1000.),random.uniform(-1000.,1000.)]) #Target position
vmpn=Vc*(xt-xmpn)/np.linalg.norm(xt-xmpn) #Missile (propnav) speed vector
vmpc=Vc*(xt-xmpc)/np.linalg.norm(xt-xmpc) #Missile (pure chase) speed vector
Vt=30 #Target velocity
theta=random.uniform(0,2*np.pi)


#%% The stuff that works
traj_mpn=[xmpn.copy()]
traj_t=[xt.copy()]
traj_mpc=[xmpc.copy()]

for t in np.arange(0,T,dt):
    theta+=np.random.normal(0.,0.5)*dt
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
    
    if lospn_norm<5.:
        print(f"The missile knows where it is: Interception at {t:.2f}s")
        break
    elif t>=T:
        print("The missile doesn't know where it is :(")
    
#%% Plotting
traj_mpn=np.array(traj_mpn)
traj_t=np.array(traj_t)
traj_mpc=np.array(traj_mpc)

fig=plt.figure(figsize=(10,6))
plt.plot(traj_mpn[:,0],traj_mpn[:,1],label="Missile (propnav)")
plt.plot(traj_mpc[:,0],traj_mpc[:,1],label="Missile (pure chase)")
plt.plot(traj_t[:,0],traj_t[:,1],label="Target (random trajectory)")
plt.plot(traj_t[-1,0],traj_t[-1,1],"rx",label="Target acquired")
plt.xlabel("x(m)")
plt.ylabel("y(m)")
plt.grid(True)
plt.axis("equal")
plt.suptitle("Proportional navigation interception")
plt.title(f"The missile knows where it is: Interception at {t:.2f}s")
plt.legend()
plt.show()
