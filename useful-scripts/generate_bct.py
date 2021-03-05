import numpy as np
import matplotlib.pyplot as plt

t1 = np.zeros(128)
t2 = np.zeros(128)
t3 = np.zeros(128)
t4 = np.zeros(128)
t5 = np.zeros(128)
t6 = np.zeros(128)

for i in range(1,64):
    t1 = t1 + np.loadtxt('./dat/%d.dat' %i, skiprows=1, usecols=(5))
    t2 = t2 + np.loadtxt('./dat/%d.dat' %i, skiprows=1, usecols=(6))
    t3 = t3 + np.loadtxt('./dat/%d.dat' %i, skiprows=1, usecols=(7))
    t4 = t4 + np.loadtxt('./dat/%d.dat' %i, skiprows=1, usecols=(8))
    t5 = t5 + np.loadtxt('./dat/%d.dat' %i, skiprows=1, usecols=(9))
    t6 = t6 + np.loadtxt('./dat/%d.dat' %i, skiprows=1, usecols=(10))

# load data into arrays

T11 = t1
T12 = t2
T13 = t3
T22 = t4
T23 = t5
T33 = t6

# limiting states
x1c    = np.matrix( [ 0.0 , 0.0 ] )
x2c    = np.matrix( [ 1.0 , 0.0 ] )
x3c    = np.matrix( [ 0.5 , np.sqrt(3.0)/2.0 ] )

Xcoord = np.zeros(len(t1))
Ycoord = np.zeros(len(t1))

for _ in range(len(t1)):
    reynolds_stress_tensor = np.array([[ T11[_] , T12[_] , T13[_] ],
                                        [ T12[_] , T22[_] , T23[_] ],
                                        [ T13[_] , T23[_] , T33[_] ]])
    trace = T11[_] + T22[_] + T33[_] + 1.0e-20
    #traceTot = KE1[_] + KE2[_] + KE3[_] + T11[_] + T22[_] + T33[_] + 1.0e-20

    anisotropy_tensor = (1.0/trace)*reynolds_stress_tensor - (1.0/3.0)*np.array( [[ 1.0 , 0.0 , 0.0 ],
                                                                                   [ 0.0 , 1.0 , 0.0 ],
                                                                                   [ 0.0 , 0.0 , 1.0 ]])
    #anisotropy_tensor = (1.0/traceTot)*(reynolds_stress_tensor - (1.0/3.0)*trace*np.array( [[ 1.0 , 0.0 , 0.0 ],[ 0.0 , 1.0 , 0.0 ],[ 0.0 , 0.0 , 1.0 ]]))

    eigenvalues_anisotropy_tensor = np.linalg.eigvalsh(anisotropy_tensor)
    sorted_eigenvalues_anisotropy_tensor = sorted(eigenvalues_anisotropy_tensor, reverse=True)
    xbc = x1c*( sorted_eigenvalues_anisotropy_tensor[0] - sorted_eigenvalues_anisotropy_tensor[1] ) + x2c*( 2.0*sorted_eigenvalues_anisotropy_tensor[1] - 2.0*sorted_eigenvalues_anisotropy_tensor[2] ) + x3c*( 3.0*sorted_eigenvalues_anisotropy_tensor[2] + 1 )

    Xcoord[_] = xbc[0,0]
    Ycoord[_] = xbc[0,1]


fig = plt.figure(figsize=[1.6*6.4, 1.2*4.8])
ax = fig.add_subplot(111)
ax.axis('equal')
cm = plt.cm.get_cmap('viridis')
sc = ax.scatter(Xcoord,Ycoord,cmap=cm,marker="o") # color by y-distance
#plt.plot(xBCT,yBCT,'k--')
#sc1 = ax.scatter(xBCTp1,yBCTp1,c=Yp1,cmap=cm,marker="_") # comment out if you don't need overlaid plots
#sc2 = ax.scatter(xBCTp2,yBCTp2,c=Yp2,cmap=cm,marker="d") #
#sc3 = ax.scatter(xBCTp3,yBCTp3,c=Yp3,cmap=cm,marker="o") #
ax.plot([1, 0.5], [0, 0.866], 'k-')
ax.plot([0, 0.5], [0, 0.866], 'k-')
ax.plot([0, 1], [0, 0], 'k-')
ax.plot([0.5, 1.01, 0.1], [-0.01, 0.5, 0.88], 'w.')
ax.text(0-.075,0-.05,'1C')
ax.text(1,0-.05,'2C')
ax.text(0.5-.04,0.866+.015,'3C')
ax.patch.set_visible(False)
#ax.set_title('x/h = 7.0', y=-0.1) ##########################
ax.axis('off')
#cbar = plt.colorbar(sc)
#bar.set_label('y/h')
#plt.savefig("%s.png" %xLoc, dpi=450,bbox_inches='tight')
plt.show()
