import sys
import math
import h5py
import numpy as np
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation
from optparse import OptionParser

# set input options
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="particles hdf5 FILE", metavar="FILE")

(options, args) = parser.parse_args()

# load the file
h5file = h5py.File( options.filename, "r" )
num_frame = len(h5file.keys())

# setup plot
fig = plt.figure()

# bounds
xbnds = [0.0,0.1];
ybnds = [0.0,0.1];
zbnds = [0.0,0.1];

def setupFrame():
    fig.clf()
    global ax1
    ax1 = fig.add_subplot(2,2,1, projection='3d')
    ax1.set_xlabel('X Axis')
    ax1.set_ylabel('Y Axis')
    ax1.set_zlabel('Z Axis')
    ax1.set_xlim3d(xbnds)
    ax1.set_ylim3d(ybnds)
    ax1.set_zlim3d(zbnds)
    global ax2
    ax2 = fig.add_subplot(2,2,2, projection='3d')
    ax2.view_init(0, 270)
    ax2.set_xlabel('X Axis')
    ax2.set_ylabel('Y Axis')
    ax2.set_zlabel('Z Axis')
    ax2.set_xlim3d(xbnds)
    ax2.set_ylim3d(ybnds)
    ax2.set_zlim3d(zbnds)
    global ax3
    ax3 = fig.add_subplot(2,2,3, projection='3d')
    ax3.set_xlabel('X Axis')
    ax3.set_ylabel('Y Axis')
    ax3.set_zlabel('Z Axis')
    ax3.set_xlim3d(xbnds)
    ax3.set_ylim3d(ybnds)
    ax3.set_zlim3d(zbnds)
    global ax4
    ax4 = fig.add_subplot(2,2,4, projection='3d')
    ax4.view_init(0, 270)
    ax4.set_xlabel('X Axis')
    ax4.set_ylabel('Y Axis')
    ax4.set_zlabel('Z Axis')
    ax4.set_xlim3d(xbnds)
    ax4.set_ylim3d(ybnds)
    ax4.set_zlim3d(zbnds)


# animation function.  This is called sequentially
def animate(i):
    print "Animating frame", i+1, "of", num_frame
    setupFrame()
    step = "TIME_STEP_" + str(i)
    group = h5file.get( step )
    time_dset = group.get("TIME");
    time = time_dset[0];
    x0 = group.get("POS_0")
    x1 = group.get("POS_1")
    x2 = group.get("POS_2")
    color = group.get("COLOR")
    v0 = group.get("VEL_0")
    v1 = group.get("VEL_1")
    v2 = group.get("VEL_2")
    s0 = group.get("STRESS_0_0")
    s1 = group.get("STRESS_1_1")
    s2 = group.get("STRESS_2_2")
    mass = group.get("MASS")
    volume = group.get("VOLUME")
    pressure = np.zeros( len(v1) )
    v = np.zeros( len(v1) )
    for p in xrange(len(x1)):
        pressure[p] = -(s0[p] + s1[p] + s2[p])/3.0
        v[p] = math.sqrt( v0[p]**2 + v1[p]**2 + v2[p]**2 )
    min_p = min(pressure)
    max_p = max(pressure)
    print min_p, max_p
    # ax1.scatter(x0[:], x1[:], x2[:], c=pressure[:], cmap='rainbow', vmin=min_p, vmax=max_p)
    # ax2.scatter(x0[:], x1[:], x2[:], c=pressure[:], cmap='rainbow', vmin=min_p, vmax=max_p)
    ax1.scatter(x0[:], x1[:], x2[:], c=color[:] )
    ax2.scatter(x0[:], x1[:], x2[:], c=color[:] )
    ax3.scatter(x0[:], x1[:], x2[:], c=v[:], cmap='rainbow', vmin=0.0, vmax=1.0)
    ax4.scatter(x0[:], x1[:], x2[:], c=v[:], cmap='rainbow', vmin=0.0, vmax=1.0)
    # timestring = 'Time ' + str(time) + ' (s)'
    # ax1.text(3, 8, timestring , style='italic',
    #         bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, frames=num_frame)
anim.save("collision.mp4", fps=10, dpi=300)

print "Animation complete."
