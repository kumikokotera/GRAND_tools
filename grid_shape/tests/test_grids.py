import numpy as np
import matplotlib.pyplot as plt

N = np.arange(0,100)
theta0 = 2*np.pi/((1+np.sqrt(5))/2)**2
r = 1.01
a = r**N
phi = N*theta0

x = a*np.cos(phi)
y = a*np.sin(phi)

# r = 0.99
# a = r**N
# phi = N*theta0
# x2 = a*np.cos(-phi)
# y2 = a*np.sin(-phi)

fig, axs = plt.subplots(1,1) 
axs.plot(x,y, 'k.')  
# axs.plot(x2,y2, 'r.')  
axs.axis('equal')  
plt.show()
