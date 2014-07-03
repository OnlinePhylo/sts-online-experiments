# SOURCE ME
_MODULE_INIT_PATH='/app/Modules/3.2.10/init/zsh'
if [ -f $_MODULE_INIT_PATH ]; then
  source $_MODULE_INIT_PATH
fi
unset $_MODULE_INIT_PATH

module load gcc/4.8.2
module use /home/matsengrp/modules
module load bpp/2.1.0
module load miniconda
