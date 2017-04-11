import sys
import math
import h5py
import numpy as np
from matplotlib import pyplot as plt
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
ax = plt.axes(xlim=(0, 0.04), ylim=(0, 0.04))
particles = ax.scatter([], [], c=[],s=25, cmap="hsv", vmin=0, vmax=1)

# animation function.  This is called sequentially
def animate(i):
    step = "TIME_STEP_" + str(i)
    group = h5file.get( step )
    x0 = group.get("POS_0")
    x1 = group.get("POS_1")
    v0 = group.get("VEL_0")
    v1 = group.get("VEL_1")
    s11 = group.get("STRAIN_1_1")
    v = np.zeros( len(v1) )
    x = np.zeros( (len(x0), 2) )
    for p in xrange(len(x1)):
        v[p] = math.sqrt( v0[p]**2 + v1[p]**2 )
        x[p][0] = x0[p]
        x[p][1] = x1[p]
    particles.set_offsets( x )
    particles.set_array( v )
    return particles,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, frames=num_frame, interval=10)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#anim.save('particles.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()
