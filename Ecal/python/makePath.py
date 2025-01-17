"""Helpful python configuration functions for getting the path to installed data files. 

This file was configured by cmake for the installation of ldmx-sw at
   @CMAKE_INSTALL_PREFIX@
"""

import os, sys

def makeBDTPath( BDTname ) :
    """Get the full path to the installed BDT files

    Exits entire python script if file does not exist.

    Parameters
    ----------
    BDTname : str
        Name of BDT file to make path for (no extension)

    Returns
    -------
    str
        full path to installed data file

    Examples
    --------
        ecalVeto.bdt_file = makeBDTPath( 'gabrielle' )
    """

    fullPath = '@CMAKE_INSTALL_PREFIX@/data/Ecal/' + BDTname + '.onnx'
    if not os.path.isfile( fullPath ) :
        print('ERROR: ONNX model file \'%s\' does not exist.' % ( fullPath ))
        sys.exit(1)

    return fullPath

def makeCellXYPath() :
    """Get the full path to the installed cell xy text file

    Returns
    -------
    str
        full path to installed data file

    Warnings
    --------
    - The need for the cellxy.txt file will be remove in upcoming ldmx-sw versions.

    Examples
    --------
        ecalVeto.cellxy_file = makeCellXYPath()
    """

    fullPath = '@CMAKE_INSTALL_PREFIX@/data/Ecal/cellxy.txt'
    if not os.path.isfile( fullPath ) :
        print('ERROR: Cell xy text file \'%s\' does not exist.' % ( fullPath ))
        sys.exit(1)

    return fullPath

def makeRoCPath( RoCname ) :
    """Get the full path to the RoC csv file

    Exits entire python script if file does not exist.

    Parameters
    ----------
    RoCname : str
        Name of RoC file to make path for (no extension)

    Returns
    -------
    str
        full path to installed data file

    Examples
    --------
        ecalVeto.roc_file = makeRoCPath( 'RoC_v14_8gev' )
    """

    fullPath = '@CMAKE_INSTALL_PREFIX@/data/Ecal/' + RoCname + '.csv'
    if not os.path.isfile( fullPath ) :
        print('ERROR: RoC csv file \'%s\' does not exist.' % ( fullPath ))
        sys.exit(1)

    return fullPath
