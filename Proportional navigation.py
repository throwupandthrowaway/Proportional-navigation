#%% Basic stuff
import numpy as np
import matplotlib.pyplot as plt
import random as random

dt=0.05 #Time step
T=1000. #If the missile hasn't reached target by then, then it knows not where it is at this time.
N=5 #Navigation constant (2<=N<=5 ideally)
Vc=40. #Closing velocity
xm=np.array([0.,0.]) #Missile position
xt=np.array([random.uniform(-500.,500.),random.uniform(-500.,500.)]) #Target position
vm=Vc*(xt-xm)/np.linalg.norm(xt-xm) #Missile speed vector
Vt=30. #Target velocity
theta=random.uniform(0,2*np.pi)


#%% The stuff that works
traj_m=[xm.copy()]
traj_t=[xt.copy()]

for t in np.arange(0,T,dt):
    theta+=np.random.normal(0.,2)*dt
    vt=Vt*np.array([np.cos(theta),np.sin(theta)]) #Target speed vector
    xt+=vt*dt
    
    los=xt-xm #LOS vector
    los_norm=np.linalg.norm(los)
    los_u=los/los_norm #LOS unitary vector direction
    
    v_rel=vt-vm #Relative speed between missile and target
    lambda_p=np.cross(los,v_rel)/los_norm**2 #LOS rate
    
    v_u=vm/np.linalg.norm(vm)
    n=np.array([-v_u[1],v_u[0]]) #Normal vector
    a_n=(N*Vc*lambda_p)*n #Normal acceleration vector
    
    vm+=a_n*dt
    xm+=vm*dt
    
    traj_m.append(xm.copy())
    traj_t.append(xt.copy())
    
    if los_norm<5.:
        print(f"Interception at {t:.2f}s")
        break
    
#%% Plotting
traj_m=np.array(traj_m)
traj_t=np.array(traj_t)

plt.figure(figsize=(10,6))
plt.plot(traj_m[:,0],traj_m[:,1],label="Missile")
plt.plot(traj_t[:,0],traj_t[:,1],label="Target")
plt.plot(traj_t[-1,0],traj_t[-1,1],"rx",label="Target acquired")
plt.xlabel("x(m)")
plt.ylabel("y(m)")
plt.grid(True)
plt.axis("equal")
plt.title("Proportional navigation interception")
plt.legend()
plt.show()