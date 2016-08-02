#############################################################################
# Motrology Mask Library (demo - test)
# v 0.11
# Date: 2014.03.25
# Authors:
# - revision
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
from numpy import loadtxt
import math
#****************************************************************************
#****************************************************************************
# Setup some transmission-type optical elements
#****************************************************************************
#****************************************************************************
def srwl_opt_setup_mask(_delta, _atten_len, _thick, _hole_sh, _hole_dim1, _hole_dim2, _pitch_x, _pitch_y, _hole_nx, _hole_ny, _mask_Nx, _mask_Ny, _hole_tilt=0,_angle=0):
    """
    Setup Transmission type Optical Element which simulates a Pinhole Mask Array (PMA) for at-wavelength metrology
    :param _delta: refractive index decrement (can be one number of array vs photon energy)
    :param _atten_len: attenuation length [m] (can be one number of array vs photon energy)
    :param _thick: thickness of mask [m]
    :param _hole_sh: hole shape. 1: Circular Holes case. 2: Rectangular Holes case. 3:square
    :param _hole_dim1: hole dimension 1, radius for circular holes or width for rectangular holes.
    :param _hole_dim2: hole dimension 2, height for rectangular holes
    :param _pitch_x: mask pitch in x-direction [m]
    :param _pitch_y: mask pitch in y-direction [m]
    :param _hole_nx: number of holes in x-direction
    :param _hole_ny: number of holes in y-direction
    :param _mask_Nx: number of pixels in x-direction  #?
    :param _mask_Ny: number of pixels in y-direction
    :param _hole_tilt: tilt angle of the mask(or the hole(?)) [rad or degree(?)]
    :return: transmission (SRWLOptT) type optical element which simulates the PMA
    """

    # Set the range of mask.
    MaskRx = _hole_nx*_pitch_x # mask range vs x [m].
    MaskRy = _hole_ny*_pitch_y # mask range vs y [m].

    # Generate OpT.
    opT = SRWLOptT(_nx=_mask_Nx, _ny=_mask_Ny, _rx=MaskRx, _ry=MaskRy, _arTr=None, _extTr=0, _x=0, _y=0)

    # distance intervial in x [m].
    hx = 0 # [m]
##    if(_mask_Nx > 1): hx = MaskRx/(_mask_Nx-1) # [m]
    if(_mask_Nx > 1): hx = MaskRx/(_mask_Nx) # [m]
    
    # distance intervial in y [m].
    hy = 0 # [m]
##    if(_mask_Ny > 1): hy = MaskRy/(_mask_Ny-1) # [m]
    if(_mask_Ny > 1): hy = MaskRy/(_mask_Ny) # [m]

    # Set the starting poision (not the center) of the top left hole.
    xStartHole = 0
    yStartHole = 0
   
    #****************************************************************************
    # Same data alignment as for wavefront: outmost loop vs y, inmost loop vs x.
    #
    # !!!! Need to Think about Rotation and other shapes.
    # Now is rectangluar with no rotation only!
    #****************************************************************************

    ofst = 0 # pointer for array opT.arTr
    
    y = -0.5*MaskRy # Mask is always centered on the grid, however grid can be shifted.
    for iy in range(_mask_Ny):
        
        # Calculate the relative position in y.
        # NOTE!!! Use round to solve the precision issue!
        nyPitchesBefore = floor(round((y - yStartHole)/_pitch_y,9)) 
        yRel = y - (yStartHole + nyPitchesBefore*_pitch_y)
        if(_hole_sh==1):
        # Make a rough judgement (J1) if this row is inside the hole according to yRel.
        #insideHoleY = True
        # NOTE!!! Use round to solve the precision issue!
        #if (round(yRel-_hole_dim2,9) >= 0): insideHoleY = False
            x = -0.5*MaskRx # Mask is always centered on the grid, however grid can be shifted.
            for ix in range(_mask_Nx):

            # Calculate the relative position in x.
            # NOTE!!! Use round to solve the precision issue!
                nxPitchesBefore = floor(round((x - xStartHole)/_pitch_x,9))
                xRel = x - (xStartHole + nxPitchesBefore*_pitch_x)

            # Make a rough judgement (J2) if this column is inside the hole according to xRel.
            #insideHoleX = True
            # NOTE!!! Use round to solve the precision issue!
            #if (round(xRel-_hole_dim1,9) >= 0): insideHoleX = False
                insideHole = True
                if ((xRel-_pitch_x/2)**2+(yRel-_pitch_y/2)**2>=_hole_dim1**2): insideHole=False
                # Give values to OpT.arTr
                # Make final judgement, based on J1, J2, J3.
                #if ((insideHoleX and insideHoleY) or ((not insideHoleX) and (not insideHoleY))):
                if(insideHole):
                    opT.arTr[ofst] = 1      # amplitude transmission.
                    opT.arTr[ofst + 1] = 0  # optical path difference.
                else:
                    #opT.arTr[ofst] = exp(-0.5*_thick/_atten_len)    # amplitude transmission.
                    opT.arTr[ofst] =0
                    opT.arTr[ofst + 1] = 0            # optical path difference.

                # Shift the pointer by 2.
                ofst += 2
                # Step x by hx.
                x += hx
        if(_hole_sh==3):
        # Make a rough judgement (J1) if this row is inside the hole according to yRel.
        #insideHoleY = True
        # NOTE!!! Use round to solve the precision issue!
        #if (round(yRel-_hole_dim2,9) >= 0): insideHoleY = False
            x = -0.5*MaskRx # Mask is always centered on the grid, however grid can be shifted.
            for ix in range(_mask_Nx):

            # Calculate the relative position in x.
            # NOTE!!! Use round to solve the precision issue!
                nxPitchesBefore = floor(round((x - xStartHole)/_pitch_x,9))
                xRel = x - (xStartHole + nxPitchesBefore*_pitch_x)

            # Make a rough judgement (J2) if this column is inside the hole according to xRel.
            #insideHoleX = True
            # NOTE!!! Use round to solve the precision issue!
            #if (round(xRel-_hole_dim1,9) >= 0): insideHoleX = False
                insideHole = False
                #if (abs(xRel-_pitch_x/2)<(_hole_dim1/(2**0.5)) and abs(yRel-_pitch_y/2)<(_hole_dim1/(2**0.5))):
                xCross1=_pitch_x/2-_hole_dim1/(2**0.5)*math.cos(_angle)
                yCross1=_pitch_y/2-_hole_dim1/(2**0.5)*math.sin(_angle)
                xCross2=_pitch_x/2+_hole_dim1/(2**0.5)*math.cos(_angle)
                yCross2=_pitch_y/2+_hole_dim1/(2**0.5)*math.sin(_angle)
                k1=math.tan(math.pi/4+_angle)
                k2=-math.tan(math.pi/4-_angle)
                k4=math.tan(math.pi/4+_angle)
                k3=-math.tan(math.pi/4-_angle)
                #print("k1: ",k1," k2: ",k2," k3: ",k3," k4: ",k4)
               
                if(yRel<(k1*xRel+(yCross1-k1*xCross1)) and yRel<(k2*xRel+(yCross2-k2*xCross2)) and yRel>(k3*xRel+(yCross1-k3*xCross1)) and yRel>(k4*xRel+(yCross2-k4*xCross2))):
                    insideHole=True
                # Give values to OpT.arTr
                # Make final judgement, based on J1, J2, J3.
                #if ((insideHoleX and insideHoleY) or ((not insideHoleX) and (not insideHoleY))):
                if(insideHole):
                    opT.arTr[ofst] = 1      # amplitude transmission.
                    opT.arTr[ofst + 1] = 0  # optical path difference.
                else:
                    #opT.arTr[ofst] = exp(-0.5*_thick/_atten_len)    # amplitude transmission.
                    opT.arTr[ofst] =0
                    opT.arTr[ofst + 1] = 0            # optical path difference.

                # Shift the pointer by 2.
                ofst += 2
                # Step x by hx.
                x += hx
        # Step y by hy.
        y += hy
    
    return opT



def srwl_opt_load_mask(_FolderName,_AmpTraFileName,_OpPaDiFileName):

    # Load Mask Range and Size from files.
    f = open(_FolderName+'/'+_OpPaDiFileName, 'r')
    NumOfLines = 10
    MaxLenthOfOneLine = 100
    ValueList=[]
    for LineNum in range(NumOfLines):
        LineStr = f.readline(MaxLenthOfOneLine)
        if (LineNum>=1 and LineNum<=9):
            EndNum = LineStr.find(' #')
            ValueList.extend([float(LineStr[1:EndNum])])
    f.close()
    MaskRx = ValueList[4]-ValueList[3]
    MaskRy = ValueList[7]-ValueList[6]
    mask_Nx = int(ValueList[5])
    mask_Ny = int(ValueList[8])
    
    # Generate OpT.
    opT = SRWLOptT(_nx=mask_Nx, _ny=mask_Ny, _rx=MaskRx, _ry=MaskRy, _arTr=None, _extTr=0, _x=0, _y=0)

    # Load data from files.
    AmpTra = loadtxt(_FolderName+'/'+_AmpTraFileName, dtype="float")    # Amplitude Transmission.
    OpPaDi = loadtxt(_FolderName+'/'+_OpPaDiFileName, dtype="float")    # Optical Path Difference.
    
    ofst = 0 # pointer for array opT.arTr
    for iy in range(mask_Ny):
        for ix in range(mask_Nx):
            opT.arTr[ofst] = AmpTra[ofst/2]      # Amplitude Transmission.
            opT.arTr[ofst + 1] = OpPaDi[ofst/2]  # Optical Path Difference.
             # Shift the pointer by 2.
            ofst += 2
    
    return opT
