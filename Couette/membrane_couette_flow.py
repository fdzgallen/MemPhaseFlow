# Import libraries
import numpy as np 
import os.path
import time
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from shutil import copyfile
from numpy import *


file=os.path.join('./inputs.txt')
with open(file,'r') as f:  
    simulation = str(f.readline().partition('#')[0]).strip() # "simulation" name of the folder where the data will be saved 
    wallspd = float(f.readline().partition('#')[0]) 	# "wallspd" Poiseuille speed in dx/dt units
    R1 = float(f.readline().partition('#')[0])   	# "R1" long ellipse radius for the initial conditions for the vesicle
    R2 = float(f.readline().partition('#')[0])	# "R2" long ellipse radius for the initial conditions for the vesicle. For Red blood cells taken as 0.23*R1
    Tend = int(f.readline().partition('#')[0])   # "Tend" total number of steps for the simulation to end
    t0 = int(f.readline().partition('#')[0])      # "t0" number of time steps to evolve the interface before starting running the fluid flow
    tdump = int(f.readline().partition('#')[0])    # "tdump" data dump period
    Nx = int(f.readline().partition('#')[0])          # "Nx" lattice size X
    Ny = int(f.readline().partition('#')[0])          # "Ny" lattice size Y
    Nx0 = int(f.readline().partition('#')[0])          # "Nx0" starting position of the cell in X axis
    Ny0 = int(f.readline().partition('#')[0])         # "Ny0" starting position of the cell in Y axis
    eps = float(f.readline().partition('#')[0])         # "eps" interfacial width, to increase eps its necessary to decrease "dt" so that the simulation converges, reducing eps gives lower resolution for stream function and voritcity fields
    C0 = float(f.readline().partition('#')[0])         # "C0" Spontaneous curvature
    M = float(f.readline().partition('#')[0])            # "M" mobility
    dt = float(f.readline().partition('#')[0])      # "dt" temporal resolution
    dx = float(f.readline().partition('#')[0])           # "dx" spatial resolution
    dy = float(f.readline().partition('#')[0])          # "dy" spatial resolution
    alpha = float(f.readline().partition('#')[0])      # "alpha" Strength of the lagrange mutiplier effect for AREA
    beta = float(f.readline().partition('#')[0])      # "beta" Strength of the lagrange mutiplier effect for volume
    nit = int(f.readline().partition('#')[0])          # "nit" number of iterations for the iterative Poisson solver. As the state of the simulation is always close to the next solution we dont need many iterations
    kappa = float(f.readline().partition('#')[0]) 	# "kappa" bending modulus of the membrane
    viscosityliq = float(f.readline().partition('#')[0])  	# "viscosityliq" viscosity value of the liquid
    viscositycell = float(f.readline().partition('#')[0])	# "viscositycell" viscosity value of the liquid inside the cell
    
viscositymemb = 0.5*(viscositycell+viscosityliq) 	# "viscositymemb" is usually taken as (viscositycell+viscosityliq)/2 if there is no interest in defining a different membrane viscosity
L=float(Nx-1)
H=float(Ny-1)
idx=1/dx
idy=1/dy
kappaee=kappa/eps/eps
iH=1/(float(H))
iL=1/(float(L))
idx=1/dx
idy=1/dy 
 
print(simulation)

########################  Library of Functions  ##################################### 

def lap(u): #laplatian of a given field U
    return idx*idx*(np.roll(u,1,0)+np.roll(u,-1,0)+np.roll(u,1,1)+np.roll(u,-1,1)-4*u)

def periodic_grad(u): #gradient of a given field U with perdiodic Boundary conditions
    return [-(np.roll(u,1,0)-np.roll(u,-1,0))*0.5*idx, -(np.roll(u,1,1)-np.roll(u,-1,1))*0.5*idy]

def periodic_gradN(u): #gradient of a given field U with Neumann Boundary Conditions
    gu = [-(np.roll(u,1,0)-np.roll(u,-1,0))*0.5*idx, -(np.roll(u,1,1)-np.roll(u,-1,1))*0.5*idy]
    gu[1][:,0]=0.0
    gu[1][:,-1]=0.0
    return gu

def compute_psi(phi0,lap_phi): # Computes the parameter Psi, which is used to facilitate formulation
    psi_sc=phi0*(-1.0+phi0*phi0)-eps*eps*lap_phi-eps*C0*(1-phi0*phi0)
    lap_psi= lap(psi_sc)
    return psi_sc,lap_psi

def compute_mu(phi0,psi_sc,lap_psi): #This function computes the functional derivative of the free energy
    mu=psi_sc*(-1.0+3.0*phi0*phi0-2*eps*C0*phi0)-eps*eps*lap_psi
    lap_mu = lap(mu)
    return mu, lap_mu

def get_area(phi): #Integrating the gradient of phi is equivalent to compute the area of the membrane because 
    # the bulk in the system has gradient 0 therefore only the membrane contributes
    aux1 = np.gradient(phi)
    return sum(aux1[0]*aux1[0]+aux1[1]*aux1[1])

def get_area2(phi): #Integrating the gradient of phi is equivalent to compute the area of the membrane because 
    # the bulk in the system has gradient 0 therefore only the membrane contributes
    return (np.absolute(phi) < 0.9).sum() # phi > 0 is a boolean array, so `sum` counts the number of True elements
    
def get_area3(phi): #Integrating the gradient of phi is equivalent to compute the area of the membrane because
    # the bulk in the system has gradient 0 therefore only the membrane contributes
    return ((phi**2-1)**2).sum() # phi > 0 is a boolean array, so `sum` counts the number of True elements

def get_volume(phi): #Summing all values of phi>0 to get an idea of the vesicle volume which should remain const
    return phi.sum() # phi > 0 is a boolean array, so `sum` counts the number of True elements

def compute_multiplier_global(A0,Area): #computation of the lagrange multiplier for area
    return  alpha * (Area - A0)
    
def compute_multiplier_volume(A0,Area): #computation of the lagrange multiplier for area
    return  beta * (Area - A0)

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    return fig.colorbar(mappable, cax=cax)

#function to plot a 2x2 Figure of the cell (phi field), the stream function, vorticity and the stream function fluctuations caused by the cell
def plot_all(phi,w,curr,t):
    f, axarr = plt.subplots(2,2)

    phiplot=axarr[0,0].imshow(phi.T,cmap='seismic',origin='lower')
    f.colorbar(phiplot, ax=axarr[0,0], ticks=[-0.95,0,0.95],fraction=0.035, orientation='horizontal')
     
    mmax=np.amax(w)
    mmin=np.amin(w)
    wplot=axarr[0,1].imshow(w.T,origin='lower', vmin = mmin, vmax = mmax)
    f.colorbar(wplot, ax=axarr[0,1], ticks=[mmin*2/3,(mmax+mmin)/2,mmax*2/3],fraction=0.035,  orientation='horizontal',format=FormatStrFormatter('%.3f'))
    
    xx = np.linspace(0, Nx-1, Nx)
    yy = np.linspace(0, Ny-1, Ny)
    X, Y = np.meshgrid(xx, yy)
     
    mmax=np.amax(curr)
    mmin=np.amin(curr)
    currplot=axarr[1,0].imshow(curr.T,cmap='plasma',origin='lower', vmin = mmin, vmax = mmax)
    axarr[1,0].quiver(X[::15,::15],Y[::15,::15],periodic_gradN(curr)[1][::15,::15].T,-periodic_gradN(curr)[0][::15,::15].T);
    f.colorbar(currplot, ax=axarr[1,0], ticks=[mmin*2/3,(mmax+mmin)/2,mmax*2/3],fraction=0.035, orientation='horizontal',format=FormatStrFormatter('%.1f'))
    
    absol = np.absolute(curr-currin)
    minmax=np.amax(absol)
    
    currplot2=axarr[1,1].imshow((curr-currin).T,cmap='seismic',origin='lower', vmin = -minmax, vmax = minmax)
    f.colorbar(currplot2, ax=axarr[1,1], ticks=[-minmax*2/3,0,minmax*2/3],fraction=0.035, orientation='horizontal',format=FormatStrFormatter('%.3f'))
    axarr[1,1].autoscale(False)
    if ( (curr-currin).sum()*(curr-currin).sum() > 0.000000):
        axarr[1,1].quiver(X[::5,::5],Y[::5,::5],periodic_gradN(curr-currin)[1][::5,::5].T,-periodic_gradN(curr-currin)[0][::5,::5].T);
        

    axarr[0,1].set_title('Vorticity ω')
    axarr[1,0].set_title('Stream function ξ')
    axarr[1,1].set_title('ξ - ξ0')

    figure = plt.gcf()
    figure.set_size_inches(6, 5.0,forward=True)
    plt.tight_layout()
    file3=os.path.join('./'+simulation+'/t='+str(t)+'.png')
    figure.savefig(file3, dpi = 200)
    plt.close('all')
 

#function to be used in case you want to start with a tilter ellipsoid instead of aligned with the x or y axis
def tilted_ellipse(a,R1,R2,x,y,x0,y0):

    return ((x-x0)*np.cos(a)-(y-y0)*np.sin(a))**2/(R1**2) + ((x-x0)*np.sin(a)+(y-y0)*np.cos(a))**2/(R2**2)
    

#function to write the data in txt files and to plot the energy, volue and area state of the simulation
def dump(t,start):
#time that took between dumps
    print('Time t='+str(float(int(t/Tend*1000))/10)+'% and took '+str(1/100*float(int(100*(time.perf_counter() - start)/1)))+'s')
#Output of the system state at the given time t
    file1=os.path.join('./'+simulation+'/phi/phi_t='+str(t)+'.txt')
    file2=os.path.join('./'+simulation+'/vorticity/vorticity_t='+str(t)+'.txt')
    file3=os.path.join('./'+simulation+'/stream/stream_t='+str(t)+'.txt')
    with open(file1,'w+') as f1:
        with open(file2,'w+') as f2:
            with open(file3,'w+') as f3:
                for i in range(0,Nx):
                    for j in range(0,Ny): 
                        f1.write('{0} {1} {2}\n'.format(i,j,phi[i,j]))
                        f2.write('{0} {1} {2}\n'.format(i,j,w[i,j]))
                        f3.write('{0} {1} {2}\n'.format(i,j,curr[i,j]))
                        
 #calling the function to plot the state of the simulation
    plot_all(phi,w,curr,t)

#we compute area volume energy etc to write in a txt to easily access the values over time
    Area=get_area(phi)
    Area2=get_area2(phi)
    energy.append(np.sum(mu**2))
    vol.append((phi>0).sum())
    areat.append(Area)
    areat2.append(Area2)
    file=os.path.join('./'+simulation+'/SimulationData.txt')
    with open(file,'w+') as f:
        f.write('t energy[t] vol[t] area[t]\n')
        for t in range(0,len(energy)):  
            f.write('{0} {1} {2} {3}\n'.format(t*tdump,energy[t],vol[t],areat[t]))
        f.close()

#now we plot the evolution of the energy volume and the two surface area parameters 
    plt.subplot(2,2,1)
    plt.plot(energy[1:])
    plt.xlabel(str(int(tdump*dt))+r'$\,\Delta t$',fontsize=12)
    plt.ylabel(r'${energy\,(t)}$',fontsize=13)
    #plt.yscale('log')

    plt.subplot(2,2,2)
    plt.plot(vol[1:])
    plt.xlabel(str(int(tdump*dt))+r'$\,\Delta t$',fontsize=12)
    plt.ylabel(r'${Volume\,(t)}$',fontsize=13)
    #plt.yscale('log')

    plt.subplot(2,2,3)
    plt.plot(areat[1:])
    plt.xlabel(str(int(tdump*dt))+r'$\,\Delta t$',fontsize=12)
    plt.ylabel(r'${Area1\,(t)}$',fontsize=13)
    #plt.yscale('log')

    plt.subplot(2,2,4)
    plt.plot(areat2[1:])
    plt.xlabel(str(int(tdump*dt))+r'$\,\Delta t$',fontsize=12)
    plt.ylabel(r'${Area2\,(t)}$',fontsize=13)
    #plt.yscale('log')

    figure = plt.gcf()  
    figure.set_size_inches(9, 9,forward=True)
    plt.tight_layout()
    file3=os.path.join('./'+simulation+'/SimulationData.png')
    figure.savefig(file3, dpi = 200)

#Code to solve the Poisson equation in 2D by an iterative solver, adapted from 2017 Lorena A. Barba, Gilbert F. Forsyth https://github.com/barbagroup/CFDPython
def vorticity_periodic(p, b):
    pn = np.empty_like(p)
    for q in range(nit):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2 +
                          (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2) /
                         (2 * (dx**2 + dy**2)) -
                         dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, 1:-1])
 

                # Periodic BC Pressure @ x = 2
        p[-1,1:-1] = (((pn[0,1:-1] + pn[-2,1:-1])* dy**2 +
                        (pn[-1,2:] + pn[-1,0:-2]) * dx**2) /
                       (2 * (dx**2 + dy**2)) -
                       dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[-1,1:-1])

        # Periodic BC Pressure @ x = 0
        p[0,1:-1] = (((pn[1,1:-1] + pn[-1,1:-1])* dy**2 +
                       (pn[0,2:] + pn[0,0:-2]) * dx**2) /
                      (2 * (dx**2 + dy**2)) -
                      dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[0,1:-1])
        p[:,-1] = -wallspd*iH  #   at y = Ny
        p[:,0]  = -wallspd*iH  #   at y = 0
    return p

#Code to solve the Poisson equation in 2D by an iterative solver, adapted from 2017 Lorena A. Barba, Gilbert F. Forsyth https://github.com/barbagroup/CFDPython
def streamfunction_periodic(p, b):
    pn = np.empty_like(p)
    for q in range(nit):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2 +
                          (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2) /
                         (2 * (dx**2 + dy**2)) -
                         dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, 1:-1])
        # Periodic BC  @ x = Nx
        p[-1, 1:-1] = (((pn[0, 1:-1] + pn[-2, 1:-1])* dy**2 +
                        (pn[-1, 2:] + pn[-1, 0:-2]) * dx**2) /
                        (2 * (dx**2 + dy**2)) -
                        dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[-1, 1:-1])
        # Periodic BC  @ x = 0
        p[0, 1:-1] = (((pn[1, 1:-1] + pn[-1, 1:-1])* dy**2 +
                        (pn[0, 2:] + pn[0, 0:-2]) * dx**2) /
                        (2 * (dx**2 + dy**2)) -
                        dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[0,1:-1])
        # Wall boundary conditions, stream function
        p[:,-1] = wallspd*H*0.5 #  at y = Ny
        p[:,0] = 0  #  at y = 0
    return p

def compute_multiplier_global2_withoutflow(gradphi,grad_lapmu,grad_laparea):
    sigma_N2 = M*np.sum(gradphi[0]*grad_lapmu[0]+gradphi[1]*grad_lapmu[1])
    sigma_D = M*np.sum(gradphi[0]*grad_laparea[0]+gradphi[1]*grad_laparea[1])
    sigma = (sigma_N2)/sigma_D 
    return sigma

def compute_multiplier_global2(gradphi,grad_lapmu,grad_laplapphi,gradcurr):
    gradiii = periodic_gradN(-gradcurr[1]*gradphi[0]+gradcurr[0]*gradphi[1])
    sigma_N1 = np.sum(-gradphi[0]*gradiii[0]-gradphi[1]*gradiii[1])
    sigma_N2 = M*np.sum(gradphi[0]*grad_lapmu[0]+gradphi[1]*grad_lapmu[1])
    sigma_D = M*np.sum(gradphi[0]*grad_laplapphi[0]+gradphi[1]*grad_laplapphi[1])
    sigma = (sigma_N1+sigma_N2)/sigma_D 
    return sigma

def compute_visco(phi):
    M=viscosityliq*0.5*(1-phi)+viscositycell*0.5*(1+phi)+(viscositymemb-0.5*viscosityliq-0.5*viscositycell)*(1-phi*phi)

    return M
 
###########################################################################################################


####################################  MAIN CODE  ##########################################################

## Initial conditions
phi=np.zeros((Nx,Ny)) #Initialization of our system matrix
w=np.zeros((Nx,Ny)) #Initialization of vorticity matrix
curr=np.zeros((Nx,Ny)) #Initialization of stream funct matrix
currin=np.zeros((Nx,Ny)) #Initialization of stream funct matrix

# Initial configuration for the system. Phi=+1 is the cell volume while phi=-1 is defined as the fluid that contains the cell.
# The membrane will form during the firsts temporal steps
for i in range(0,Nx):
    for j in range(0,Ny):
#IF phase selection for each point (phase field)
        if ((i-Nx0)*(i-Nx0)/(R2**2)+(j-Ny0)*(j-Ny0)/(R1**2)<=1.0): # ellipsoidal initial conditions1
#        if(i<float(Nx)*0.66 and i>float(Nx)*0.33 ):
#        if ( tilted_ellipse(2.2,R1,R2,i,j,Nx0,Ny0) <= 1.0 ):
            phi[i,j]= 1.0;
        else:
            phi[i,j]=-1.0;
        #stream function and vorticity for a poiseuille flow
        curr[i,j]=-cte*float(j)**2*(2*float(j)-3*(Ny-1))/6  #stream function 
        currin[i,j]=curr[i,j]                               #initial stream function 
        w[i,j] = -cte*((Ny-1)-2*float(j))                   #voricity


#folders to store data, in case of non-existing folder it will make a new one
print(os.getcwd())
if not os.path.exists('./'+simulation+'/'):
    os.makedirs('./'+simulation+'/')
    print('new folder '+simulation+'')
if not os.path.exists('./'+simulation+'/phi/'):
    os.makedirs('./'+simulation+'/phi/')
    print('new folder phi')
if not os.path.exists('./'+simulation+'/stream/'):
    os.makedirs('./'+simulation+'/stream/')
    print('new folder stream')
if not os.path.exists('./'+simulation+'/vorticity/'):
    os.makedirs('./'+simulation+'/vorticity/')
    print('new folder vorticity')

#We copy the code and the parameters used to run the simulation in the data simulation folder in case we need to check in the future.
dodo=os.path.join('./'+simulation+'/'+os.path.basename(__file__))
copyfile(__file__,dodo)    
file=os.path.join('./'+simulation+'/simulation_parameters.txt')
with open(file,'w+') as f: 
    f.write('M = {0}\n'.format(M))
    f.write('dt = {0}\n'.format(dt))
    f.write('eps = {0}\n'.format(eps))
    f.write('alpha = {0}\n'.format(alpha))
    f.write('nit = {0}\n'.format(nit))
    f.write('deltaPdl = {0}\n'.format(deltaPdl))
    f.write('kappa = {0}\n'.format(kappa)) 
    f.write('dx = {0}\n'.format(dx))
    f.write('dy = {0}\n'.format(dy))
    f.write('R1 = {0}\n'.format(R1))
    f.write('R2 = {0}\n'.format(R2))
    f.write('viscosityliq = {0}\n'.format(viscosityliq))
    f.write('viscositycell = {0}\n'.format(viscositycell))
    f.write('viscositymemb = {0}\n'.format(viscositymemb))
    f.write('Poiseuille average speed = {0}\n'.format(1/H*deltaPdl*(H)**3/(6*viscosityliq)))
    f.write('Tend = {0}\n'.format(Tend))
    f.write('t0 = {0}\n'.format(t0))
    f.write('tdump = {0}\n'.format(tdump))
    f.write('Nx = {0}\n'.format(Nx))
    f.write('Ny = {0}\n'.format(Ny))
    f.write('Nx0 = {0}\n'.format(Nx0))
    f.write('Ny0 = {0}\n'.format(Ny0)) 
    f.write('C0 = {0}\n'.format(C0))
    f.write('wallspd = {0}\n'.format(wallspd)) 
    f.close()
 
############################################ Temporal evolution of the system ####################################################

# First the relaxation of the interfase to a diffuse interfase then it goes into the normal temporal evolution with the addition of the lagrange multiplier which conserves area 
Sigma=0.0 #Lagrange multiplier for the area starts at 0 while the diffuse interfase is generated
Sigma2=0.0
Sigma3=0.0
#timer to track how long takes between each dump of new data
start = time.perf_counter()
#parameters used to track if the temporal evolution is going well
energy=[]
vt=[]
vol=[]
areat=[]
areat2=[] 
 
for t in range(0,t0):
#Compute parameters of the system and their derivatives
    lap_phi = lap(phi)
    laparea2=lap(phi*(phi**2-1))
    lap_lap_phi=lap(lap_phi)
    psi_sc,lap_psi = compute_psi(phi,lap_phi)
    mu,lap_mu = compute_mu(phi,psi_sc,lap_psi)
    if(t < 100000/2): 
        A0=get_area3(phi)
        V0=get_volume(phi)
    v=get_volume(phi)
    area2=get_area3(phi)
    if(t > 100000/2):
        gradphi = periodic_gradN(phi)
        gradmu = periodic_gradN(mu)
        grad_lapmu = periodic_gradN(lap_mu)
        grad_laplapphi = periodic_gradN(lap_lap_phi) 
        Sigma = compute_multiplier_global2_withoutflow(gradphi,grad_lapmu,grad_laplapphi)
        Sigma2 = compute_multiplier_global(A0,area2)
        Sigma3 = compute_multiplier_volume(V0,v)
    phi[:,1:-1] += dt*M*(kappa*lap_mu[:,1:-1]-Sigma*lap_lap_phi[:,1:-1]+Sigma2*laparea2[:,1:-1]+Sigma3)                                 #temporal evolution of the order parameter
    
    if (t % tdump == 0): #Dumping the data
        dump(t,start)
        start = time.perf_counter()        
         


for t in range(t0,Tend):
#Compute parameters of the system and their derivatives
    lap_phi = lap(phi)
    lap_lap_phi=lap(lap_phi)
    psi_sc,lap_psi = compute_psi(phi,lap_phi)
    mu,lap_mu = compute_mu(phi,psi_sc,lap_psi)
    viscosity = compute_visco(phi) 
    area2=get_area3(phi)
    v=get_volume(phi)
    
    gradphi = periodic_gradN(phi)
    gradmu = periodic_gradN(mu)
    grad_lapmu = periodic_gradN(lap_mu)
    grad_laplapphi = periodic_gradN(lap_lap_phi)
    gradcurr = periodic_gradN(curr)
    
    Sigma = compute_multiplier_global2(gradphi,grad_lapmu,grad_laplapphi,gradcurr) #computing of the lagrange multiplier by explicit computation 
    Sigma2 = compute_multiplier_global(A0,area2)
    Sigma3 = compute_multiplier_volume(V0,v) 
    
    w = vorticity_periodic(w, (gradmu[1]*gradphi[0] - gradmu[0]*gradphi[1])/viscosity) #cuidado signos de las variables
    curr = streamfunction_periodic(curr, - w) #cuidado signos de las variablesç
    laparea2=lap(phi*(phi**2-1))
    phi[:,1:-1] += dt*M*(kappa*lap_mu[:,1:-1]-Sigma*lap_lap_phi[:,1:-1] + Sigma2*laparea2[:,1:-1]+Sigma3)   #basic temporal evolution of the order parameter
    phi[:,1:-1] += -dt*(gradcurr[1][:,1:-1] *gradphi[0][:,1:-1] -gradcurr[0][:,1:-1] *gradphi[1][:,1:-1] ) #advection term
 

    if (t % tdump == 0): #Dumping the data 
        dump(t,start)
        start = time.perf_counter()
         