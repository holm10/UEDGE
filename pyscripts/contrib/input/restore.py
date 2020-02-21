# Restore package for UEDGE, restore.py
# Separated from 160299 master input by holm10
# Changelog
# 200213 - file created
from uedge import bbb,com,grd,flx
from uedge.hdf5 import hdf5_restore

def label(restart):
    """====================================================================================================
    RESTORE SAVE FILE
    ===================================================================================================="""
    '''
    Function to restore UEDGE save file if argument is True
    '''

    if restart==True:
        bbb.restart=1       # Flag to restart from previous case (=1)
    if bbb.restart==1:      # Restart from previous case
        hdf5_restore("../solutions/"+bbb.label[0].decode('UTF-8')+".hdf5")    # Append hdf5 restart name and relative pat and relative pathh



