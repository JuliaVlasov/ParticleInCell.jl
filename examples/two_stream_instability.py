import os, time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import animation
from scipy import sparse
from scipy.sparse import linalg
from tqdm import tqdm


# Set plotting parameters

params = {'axes.labelsize': 'large',
          'xtick.labelsize': 'medium', 
          'ytick.labelsize': 'medium', 
          'font.size': 15,
          'font.family': 'sans-serif', 
          'text.usetex': False, 
          'mathtext.fontset': 'stixsans',}

plt.rcParams.update(params) 
plt.rc('animation', html='html5')

# Simulation parameters
L = 20*np.pi # 20*np.pi # Domain size 
DT = 0.005 # Time step
NT = 10000 # Number of time steps
NG = 320 # Number of grid cells
N = NG * 20 # Number of particles
WP = 1.  # Plasma frequency
QM = -1.  # Charge/mass ratio
V0 = 0.9 # Stream velocity
VT = 0.0000001 # Thermal speed

# perturbation
XP1 = 1.0 
mode = 1
Q = WP**2 / (QM*N/L) # rho0*L/N: charge carried by a single particle? 
rho_back = -Q*N/L # Background charge density?
dx = L / NG # Grid step
# Auxilliary vectors
p = np.concatenate([np.arange(N), np.arange(N)]) # Some indices up to N 
Poisson = sparse.spdiags(([1, -2, 1] * np.ones((1, NG-1), dtype=int).T).T, \
[-1, 0, 1], NG-1, NG-1) 
Poisson = Poisson.tocsc()
# Cell center coordinates
xg = np.linspace(0, L-dx, NG) + dx/2
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(0, L), ylim=(-7, 7)) 
particles1, = ax.plot([],[], 'rd', ms=0.5)
particles2, = ax.plot([],[], 'bd', ms=0.5)

def init():
    global xp, vp
    xp = np.linspace(0, L-L/N, N).T # Particle positions
    vp = VT * (1 - VT**2)**(-0.5) * np.random.randn(N) # Particle momentum, initially M 
    pm = np.arange(N)
    pm = 1 - 2 * np.mod(pm+1, 2)
    vp += pm * (V0 * (1 - V0**2)**(-0.5)) # Momentum + stream velocity
    # Add electron perturbation to excite the desired mode
    xp += XP1 * (L/N) * np.sin(2 * np.pi * xp / L * mode) 
    xp[np.where(xp < 0)] += L
    xp[np.where(xp >= L)] -= L
    global histEnergy, histPotE, histKinE, histMomentum 
    histEnergy, histPotE, histKinE, histMomentum = [], [], [], [] 
    particles1.set_data(xp[0::2], vp[0::2])
    particles2.set_data(xp[1::2], vp[1::2])
    return particles1, particles2

def update_position(xp, vp, DT): 
    # update particle position xp
    xp += vp * DT
    # Periodic boundary condition 
    xp[np.where(xp < 0)] += L 
    xp[np.where(xp >= L)] -= L 
    return xp

def compute_rho(xp):
    # Project particles->grid
    g1 = np.floor(xp/dx - 0.5)
    g = np.concatenate((g1, g1+1))
    fraz1 = 1 - np.abs(xp/dx - g1 - 0.5) 
    fraz = np.concatenate((fraz1, 1-fraz1)) 
    g[np.where(g < 0)] += NG
    g[np.where(g > NG-1)] -= NG
    mat = sparse.csc_matrix((fraz, (p, g)), shape=(N, NG)) 
    return mat, Q / dx * mat.toarray().sum(axis=0) + rho_back

def solve_poisson( rho ):
    # Compute electric field potential
    Phi = linalg.spsolve(Poisson, -dx**2 * rho[0:NG-1]) 
    Phi = np.concatenate((Phi,[0]))
    # Electric field on the grid
    return (np.roll(Phi, 1) - np.roll(Phi, -1)) / (2*dx)

init()

for it in tqdm(range(NT+1)):
    xp = update_position(xp, vp, DT) 
    mat, rho = compute_rho(xp)
    Eg = solve_poisson(rho)
    vp += mat * QM * Eg * DT
    Etot = 0.5 * (Eg**2).sum() * dx
    histEnergy.append(Etot)
    histPotE.append(0.5 * (Eg**2).sum() * dx) 
    histKinE.append(0.5 * Q/QM * (vp**2).sum()) 
    histMomentum.append(Q/QM * vp.sum())

def animate(frame): 
    global xp, vp
    global Eg
    global histEnergy, histPotE, histKinE, histMomentum
    nsteps = 500 # output period for istep in range(nsteps):
    xp = update_position(xp, vp, DT) 
    mat, rho = compute_rho(xp)
    Eg = solve_poisson(rho)
    vp += mat * QM * Eg * DT
    Etot = 0.5 * (Eg**2).sum() * dx
    histEnergy.append(Etot)
    histPotE.append(0.5 * (Eg**2).sum() * dx) 
    histKinE.append(0.5 * Q/QM * (vp**2).sum()) 
    histMomentum.append(Q/QM * vp.sum())
    #particles1.set_data(xp[0:-1:2], vp[0:-1:2]) 
    #particles2.set_data(xp[1:-1:2], vp[1:-1:2]) 
    return particles1, particles2

ani = animation.FuncAnimation(fig, animate, frames=NT+1, repeat=False, blit=True, init_func=init)

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=200, metadata=dict(artist='Me'), bitrate=1800)
ani.save('tsi.mp4', writer=writer)

t = np.arange(0.0, (NT+1)*DT, DT)

plt.rcParams['figure.figsize'] = (12,8)

fig = plt.figure()
ax1 = fig.add_subplot(221, autoscale_on=False, xlim=(0, L), ylim=(-7,7))
ax1.scatter(xp[0:-1:2], vp[0:-1:2], s=0.5, marker='.', color='blue') 
ax1.scatter(xp[1:-1:2], vp[1:-1:2], s=0.5, marker='.', color='red') 
ax1.set_xlabel('x')
ax1.set_ylabel('P')
ax1.legend((mpatches.Patch(color='w'), ), (r'$\omega_{pe}$', ), loc=1, frameon=False)
# Electric field
ax2 = fig.add_subplot(222, autoscale_on=True, xlim=(0, L))
ax2.set_xlabel('x')
ax2.plot(xg, Eg, label='E', linewidth=2) 
ax2.legend(loc=1)
# Energies
ax3 = fig.add_subplot(223, autoscale_on=True, xlim=(0, NT*DT))
ax3.set_xlabel('time')
ax3.set_yscale('log')
ax3.plot(t, histEnergy, label='Total Energy', linewidth=2) 
ax3.plot(t, histPotE, label='Potential', linewidth=2) 
ax3.plot(t, histKinE, label='Kinetic', linewidth=2) 
ax3.legend(loc=4)
# Momentum
ax4 = fig.add_subplot(224, autoscale_on=True, xlim=(0, NT*DT))
ax4.set_xlabel('time')
ax4.plot(t, histMomentum, label='Momentum', linewidth=2) 
ax4.legend(loc=1);
