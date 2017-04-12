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

def setupFrame():
    fig.clf()
    global ax1
    ax1 = fig.add_subplot(2,2,1, projection='3d')
    ax1.set_xlabel('X Axis')
    ax1.set_ylabel('Y Axis')
    ax1.set_zlabel('Z Axis')
    ax1.set_xlim3d([0.0, 1.0])
    ax1.set_ylim3d([0.0, 1.0])
    ax1.set_zlim3d([0.0, 1.0])
    global ax2
    ax2 = fig.add_subplot(2,2,2, projection='3d')
    ax2.view_init(90, 0)
    ax2.set_xlabel('X Axis')
    ax2.set_ylabel('Y Axis')
    ax2.set_zlabel('Z Axis')
    ax2.set_xlim3d([0.0, 1.0])
    ax2.set_ylim3d([0.0, 1.0])
    ax2.set_zlim3d([0.0, 1.0])
    global ax3
    ax3 = fig.add_subplot(2,2,3, projection='3d')
    ax3.set_xlabel('X Axis')
    ax3.set_ylabel('Y Axis')
    ax3.set_zlabel('Z Axis')
    ax3.set_xlim3d([0.0, 1.0])
    ax3.set_ylim3d([0.0, 1.0])
    ax3.set_zlim3d([0.0, 1.0])
    global ax4
    ax4 = fig.add_subplot(2,2,4, projection='3d')
    ax4.view_init(90, 0)
    ax4.set_xlabel('X Axis')
    ax4.set_ylabel('Y Axis')
    ax4.set_zlabel('Z Axis')
    ax4.set_xlim3d([0.0, 1.0])
    ax4.set_ylim3d([0.0, 1.0])
    ax4.set_zlim3d([0.0, 1.0])


# animation function.  This is called sequentially
def animate(i):
    setupFrame()
    step = "TIME_STEP_" + str(i)
    group = h5file.get( step )
    x0 = group.get("POS_0")
    x1 = group.get("POS_1")
    x2 = group.get("POS_2")
    v0 = group.get("VEL_0")
    v1 = group.get("VEL_1")
    v2 = group.get("VEL_2")
    v = np.zeros( len(v1) )
    for p in xrange(len(x1)):
        v[p] = math.sqrt( v0[p]**2 + v1[p]**2 + v2[p]**2 )
    matid = group.get("MATID")
    ax1.scatter(x0[:], x1[:], x2[:], c=matid[:])
    ax2.scatter(x0[:], x1[:], x2[:], c=matid[:])
    ax3.scatter(x0[:], x1[:], x2[:], c=v[:], cmap='rainbow', vmin=0, vmax=0.1)
    ax4.scatter(x0[:], x1[:], x2[:], c=v[:], cmap='rainbow', vmin=0, vmax=0.1)

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, frames=num_frame)
anim.save("collision.mp4", fps=30)

plt.show()
