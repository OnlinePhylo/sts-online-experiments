source venv/bin/activate

module use /home/matsengrp/modules
module load bpp/2.1.0
module load beagle/2.1
module load bali-phy/2.1.1

export R_LIBS=$PWD/venv/lib/R

rehash
