import os 
os.environ['LD_LIBRARY_PATH']='/usr/local/lib'

!echo $LD_LIBRARY_PATH
!sudo ln -s /usr/local/lib/libmkl_intel_lp64.so /usr/local/lib/libmkl_intel_lp64.so.1
!sudo ln -s /usr/local/lib/libmkl_intel_thread.so /usr/local/lib/libmkl_intel_thread.so.1
!sudo ln -s /usr/local/lib/libmkl_core.so /usr/local/lib/libmkl_core.so.1

!ldconfig
!ldd /usr/local/lib/python3.7/dist-packages/torch/lib/libtorch.so
