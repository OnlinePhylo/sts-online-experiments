source venv/bin/activate

module use /home/matsengrp/modules
module load bpp/2.2.0
module load beagle/2.1
module load bali-phy/2.1.1
module load nlopt/2.4.2

export R_LIBS=$PWD/venv/lib/R

# the environment modules aren't updating LD_LIBRARY_PATH properly, so hack in a fix
export LD_LIBRARY_PATH="$LD_RUN_PATH"

