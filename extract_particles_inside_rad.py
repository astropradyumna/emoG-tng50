import numpy as np
import illustris_python as il
import matplotlib.pyplot as plt
import math
import pdb
import sys
import os
from illustris_python.snapshot import getSnapOffsets, snapPath, getNumPart
from illustris_python.util import partTypeNum
import h5py


### uses python 3. if using python 2, add 'from future import print' ###


def fold_pos(x, lbox):
        '''
        This is to account for the preiodic box condition
        '''
        aux = x > lbox / 2.
        x[aux] = x[aux] - lbox
        return x



### TNG50 output found here ###
basePath = '/rhome/psadh003/bigdata/L35n2160TNG_fixed/output'

snpz0 = 99
g = 4.3e-6
h_small = 0.6744 #0.6744

# We will be taking input from the terminal for the FoF group number. 
fof_no = int(sys.argv[1]) #This would be the FoF number that we would be getting the files for
haloID = fof_no
fof_str = 'fof' + str(fof_no)


this_fof = il.groupcat.loadSingle(basePath, 99, haloID = fof_no)
central_sfid_99 = this_fof['GroupFirstSub']
group_pos = this_fof['GroupPos'] / h_small
group_vel = this_fof['GroupVel']
r200 = this_fof['Group_R_Crit200'] / h_small


### output will be stored into the following folder
outpath = '/rhome/psadh003/bigdata/tng50/fof_partdata_3rvir/' + fof_str + '_partdata/' 

if not os.path.exists(outpath): #If the directory does not exist, then just create it!
    os.makedirs(outpath)


subset = getSnapOffsets(basePath, 99, 0, "Subhalo")

with h5py.File(snapPath(basePath, 99), 'r') as f:
        header = dict(f['Header'].attrs.items())
        nPart = getNumPart(header)

startpoint=0
chunk_len=500000
partType = 'dm'

dm_pos_ar = np.array(([np.inf, np.inf, np.inf], [np.inf, np.inf, np.inf]))
dm_vel_ar = np.array(([np.inf, np.inf, np.inf], [np.inf, np.inf, np.inf]))
dm_id_ar = np.zeros(0)
while(startpoint < nPart[partTypeNum(partType)]):
# while(startpoint < 500000 * 10):
    len_range=np.min([chunk_len, nPart[partTypeNum(partType)]-startpoint])
    subset['offsetType'][partTypeNum(partType)] = startpoint
    subset['lenType'][partTypeNum(partType)]    = len_range
    dm = il.snapshot.loadSubset(basePath, 99, 'dm', fields=['Velocities','Coordinates', 'ParticleIDs'], subset=subset)
    startpoint=startpoint+chunk_len #updating the startpoint here
    pos = dm['Coordinates']/h_small
    dist = np.sqrt((pos[:,0] - group_pos[0])**2 + (pos[:,1] - group_pos[1])**2 + (pos[:,2] - group_pos[2])**2)
    w = np.where(dist < 3*r200)[0]
    if len(w) > 0:
        pos_wrt_fof = pos[w] - group_pos
        dm_pos_ar = np.vstack((dm_pos_ar, pos_wrt_fof))
        vel = dm['Velocities']
        vel_wrt_fof = vel[w] - group_vel
        dm_vel_ar = np.vstack((dm_vel_ar, vel_wrt_fof))
        dm_id_ar = np.append(dm_id_ar, dm['ParticleIDs'][w])


dm_pos_ar = dm_pos_ar[2:, :]
dm_vel_ar = dm_vel_ar[2:, :]

# Let us now save these as .npy files
np.save(outpath + 'dm_pos.npy', dm_pos_ar)
np.save(outpath + 'dm_vel.npy', dm_vel_ar)
np.save(outpath + 'dm_ids.npy', dm_id_ar)



# Let us now look at the stars within the FoF group

startpoint=0
chunk_len=500000
partType = 'stars'

star_pos_ar = np.array(([np.inf, np.inf, np.inf], [np.inf, np.inf, np.inf]))
star_vel_ar = np.array(([np.inf, np.inf, np.inf], [np.inf, np.inf, np.inf]))
star_id_ar = np.zeros(0)
star_mass_ar = np.zeros(0)

while(startpoint < nPart[partTypeNum(partType)]):
# while(startpoint < 500000 * 10):
    len_range=np.min([chunk_len, nPart[partTypeNum(partType)]-startpoint])
    subset['offsetType'][partTypeNum(partType)] = startpoint
    subset['lenType'][partTypeNum(partType)]    = len_range
    star = il.snapshot.loadSubset(basePath, 99, 'stars', fields=['Velocities','Coordinates', 'Masses', 'ParticleIDs'], subset=subset)
    startpoint=startpoint+chunk_len #updating the startpoint here
    pos = star['Coordinates']/h_small
    dist = np.sqrt((pos[:,0] - group_pos[0])**2 + (pos[:,1] - group_pos[1])**2 + (pos[:,2] - group_pos[2])**2)
    w = np.where(dist < 3*r200)[0]
    if len(w) > 0:
        pos_wrt_fof = pos[w] - group_pos
        star_pos_ar = np.vstack((star_pos_ar, pos_wrt_fof))
        vel = star['Velocities']
        vel_wrt_fof = vel[w] - group_vel
        star_vel_ar = np.vstack((star_vel_ar, vel_wrt_fof))
        star_id_ar = np.append(star_id_ar, star['ParticleIDs'][w])
        star_mass_ar = np.append(star_mass_ar, star['Masses'][w] * 1e10/h_small)
        

star_pos_ar = star_pos_ar[2:, :]
star_vel_ar = star_vel_ar[2:, :]

# Let us now save these as .npy files

np.save(outpath + 'star_pos.npy', star_pos_ar)
np.save(outpath + 'star_vel.npy', star_vel_ar)
np.save(outpath + 'star_ids.npy', star_id_ar)
np.save(outpath + 'star_mass.npy', star_mass_ar)














# #My machine would need 350GB of memory (95 avail) to read in the entire set of HDF particles files at once. Instead, we read 1/10th of them at a time (thus the length). 
# part_len = 15625000000.
# loops = 10.0
# bin_size = part_len/loops

# #We are looking at a 5x5x5Mpc box around each location. We need to be real careful of when "h" lives in the units, else we get this wrong.
# #For Illustris, they are pretty good at keeping everything in inverse h distance units.
# vol_size = 5000.0

# for i in range(0,int(loops)):
#  start = int(i*bin_size)
#  stop = int((i+1)*bin_size - 1)
#  # print i, start, stop
#  #we use the simple access mechaism to access the data (h5py). Usage is described here: http://www.tng-project.org/data/docs/specifications/
#  #Note that particles for any given halo can live in any of the HDF snapshot files. In other words, halo #0 might have particles that live
#  #in an HDF file that is not used for the first 156250000 particles and so we would only catch it on the second (or third...tenth) loop.
#  #For now, we therefore save all of the particles for EACH HALO at A SPECIFIC LOOP over the particles. This means that we will need to
#  #do a second round of data cleaning where we combine the particles from each loop for any specific halo back into one single particle file
#  #for  that halo. Below, str(clusters[0][j]) puts the halo id in the filename while str(i) denotes the loop #.
#  with h5py.File('simulation.hdf5','r') as f:
#     dm_positions = f['/Snapshots/99/PartType1/Coordinates'][start:stop,:]
#     for j in range(0,len(centers[0])):
#         print j
#         w = np.where((np.abs(dm_positions[:,0] - clusters[8][j][0]) < vol_size) & (np.abs(dm_positions[:,1] - clusters[8][j][1]) < vol_size) & (np.abs(dm_positions[:,2] - clusters[8][j][2]) < vol_size))[0]
#         if (len(w) > 0):
#           filename = front + '/particles/halo_' + str(clusters[0][j]) + '_positions_' + str(i) + '.fits'
#           t = Table([dm_positions[w,0]/1000.0,dm_positions[w,1]/1000.0, dm_positions[w,2]/1000.0],  names=('px','py','pz'))
#           t.write(filename, format='fits',overwrite=True)
#  #We cannot read in both the velocities and the particles at the same time or else we run out of memory. We could just shrink the size of 
#  #each chunk (i.e., increase the number of loops). Instead, i choose to write out the list of member tracers. I will use these later on the
#  #velocities.
#           filename = front + '/particles/halo_' + str(clusters[0][j]) + '_members_' + str(i) + '.fits'
#           t = Table([w],names=('w'))
#           t.write(filename, format='fits',overwrite=True)
#  #I have to delete the variable before I even try to read in the next set, else we get memory issues.
#     del dm_positions
# #To get the velocities, I use the vector of saved members, which goes a lot faster.
