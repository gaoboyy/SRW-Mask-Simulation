#############################################################################
# Metrology Mask Library
# v1.0
# Date: 2016.08.02
# Authors: Bo GAO and Lei HUANG
#############################################################################

from __future__ import print_function  # Python 2.7 compatibility
from srwlib import *
from numpy import loadtxt
import math


# ***************************************************************************
# ***************************************************************************
# Setup some transmission-type optical elements
# ***************************************************************************
# ***************************************************************************

def srwl_opt_setup_mask(_delta, _atten_len, _thick,
                        _hx, _hy, _pitch_x, _pitch_y, _mask_Nx, _mask_Ny,
                        _grid_nx, _grid_ny, _grid_sh, _grid_dx, _grid_dy=0, _grid_angle=0):
    """
    Setup Transmission type Optical Element which simulates a mask array for at-wavelength metrology.
    :param _delta: refractive index decrement (can be one number of array vs photon energy)
    :param _atten_len: attenuation length [m] (can be one number of array vs photon energy)
    :param _thick: thickness of mask [m]
    :param _hx: sampling interval in x-direction [m]
    :param _hy: sampling interval in y-direction [m]
    :param _pitch_x: grid pitch in x-direction [m]
    :param _pitch_y: grid pitch in y-direction [m]
    :param _mask_Nx: number of pixels in x-direction [1]
    :param _mask_Ny: number of pixels in y-direction [1]
    :param _grid_nx: number of grids in x-direction
    :param _grid_ny: number of grids in y-direction
    :param _grid_sh: grid shape. 0: Circular grids case. 1: Rectangular grids case. 2: 2-D phase grating
(!) :param _grid_dx: grid dimension in x-direction, width for rectangular or elliptical grids [m]
(!) :param _grid_dy: grid dimension in y-direction, height for rectangular or elliptical grids [m]
    :param _grid_angle: tilt angle of the grid [rad]
    :return: transmission (SRWLOptT) type optical element which simulates the PMA
    """
    # Check if _grid_dy is set by user.
    if _grid_dy == 0:
        _grid_dy = _grid_dx  # An ellipse becomes a circle and a rectangle becomes a square.      
    if _grid_sh == 2:
        _grid_dx = _pitch_x  #Change grid_size for 2D grating, grid_size eauql to pitch
        _grid_dy = _pitch_y
    # Calculate the range of mask.
    mask_Rx = _hx * _mask_Nx  # mask range in x-direction [m].
    mask_Ry = _hy * _mask_Nx  # mask range in y-direction [m].

    # Calculate the range of grid.
    grid_Rx = _pitch_x * _grid_nx  # grid range in x-direction [m].
    grid_Ry = _pitch_y * _grid_ny  # grid range in y-direction [m].

    # Generate Transmission Optical Element.
    trans_opt = SRWLOptT(_nx=_mask_Nx, _ny=_mask_Ny, _rx=mask_Rx, _ry=mask_Ry, _arTr=None, _extTr=0, _x=0, _y=0)

    # ****************************************************************************
    # Same data alignment as for wavefront: out-most loop vs y, in-most loop vs x.
    # ****************************************************************************

    pointer = 0  # pointer for array trans_opt.arTr

    y = - mask_Ry / 2  # Mask is always centered on the grid, however grid can be shifted.

    for iy in range(_mask_Ny):

        # Calculate the relative position in y.
        # NOTE: Use round to solve the precision issue!
        pitch_num_y = floor(round(y / _pitch_y, 9))
        y_rel = y - (pitch_num_y * _pitch_y)
        if y_rel >= _pitch_y / 2: 
            y_rel -= _pitch_y
            
        x = - mask_Rx / 2  # Mask is always centered on the grid, however grid can be shifted.
        for ix in range(_mask_Nx):

            # Calculate the relative position in x.
            # NOTE: Use round to solve the precision issue!
            pitch_num_x = floor(round(x / _pitch_x, 9))
            x_rel = x - (pitch_num_x * _pitch_x)

            if x_rel >= _pitch_x / 2:
                x_rel -= _pitch_x

            # Initialize the bool parameter.
            inside_hole = False

            # Hartmann hole in an elliptical shape.
            if _grid_sh == 0:

                if (x_rel / _grid_dx) ** 2 + (y_rel/ _grid_dy) ** 2 < 1 \
                        and not (round(x_rel - x, 9) == 0 and round(y_rel - y, 9) == 0) \
                        and abs(x) < grid_Rx / 2 and abs(y) < grid_Ry / 2:
                    inside_hole = True

            # Hartmann hole in a rectangular shape.
            elif _grid_sh == 1:

                # Calculate the equations for edges of rectangle.
                xCross1 = - _grid_dx / (2 ** 0.5) * math.cos(_grid_angle)
                yCross1 = - _grid_dx / (2 ** 0.5) * math.sin(_grid_angle)
                xCross2 = + _grid_dx / (2 ** 0.5) * math.cos(_grid_angle)
                yCross2 = + _grid_dx / (2 ** 0.5) * math.sin(_grid_angle)
                k1 = math.tan(math.pi / 4 + _grid_angle)
                k2 = -math.tan(math.pi / 4 - _grid_angle)
                k4 = math.tan(math.pi / 4 + _grid_angle)
                k3 = -math.tan(math.pi / 4 - _grid_angle)

                if (k2 * x_rel + (yCross2 - k2 * xCross2)) > y_rel > (k3 * x_rel + (yCross1 - k3 * xCross1)) \
                        and (k1 * x_rel + (yCross1 - k1 * xCross1)) > y_rel > (k4 * x_rel + (yCross2 - k4 * xCross2)) \
                        and not (abs(x) < _pitch_x / 2 and abs(y) < _pitch_y / 2) \
                        and abs(x) < grid_Rx / 2 and abs(y) < grid_Ry / 2:
                    inside_hole = True

            # Grating shearing interferometry in a 2D phase grating.
            elif _grid_sh == 2:               
                phase_shift = False
                if ((x_rel>=0 and y_rel<0) or (x_rel<0 and y_rel>=0) ):
                    phase_shift = True
                

            else:
                print('''Unknown shape code.''')  # (!)

            # Give values to trans_opt.arTr.
            if inside_hole  and not(_grid_sh == 2):
                trans_opt.arTr[pointer] = 1      # amplitude transmission.  (!) not in physics yet
                trans_opt.arTr[pointer + 1] = 0  # optical path difference. (!) not in physics yet
            else:
                trans_opt.arTr[pointer] = 0      # amplitude transmission.  (!) not in physics yet
                trans_opt.arTr[pointer + 1] = 0  # optical path difference. (!) not in physics yet
                
            if _grid_sh == 2:
                # Give values to OpT.arTr
                # Make final judgement.
                if (phase_shift):                
                    trans_opt.arTr[pointer] = exp(-0.5*_thick/_atten_len)    # amplitude transmission.                    
                    trans_opt.arTr[pointer + 1] = -_delta*_thick             # optical path difference.
                else:                   
                    trans_opt.arTr[pointer] = 1      # amplitude transmission.
                    trans_opt.arTr[pointer + 1] = 0  # optical path difference.
                if not(abs(x) < grid_Rx / 2 and abs(y) < grid_Ry / 2) : # check if it is in grid area
                    trans_opt.arTr[pointer] = 0
                    trans_opt.arTr[pointer+1] = 0   
                
            # Shift the pointer by 2.
            pointer += 2

            # Step x by _hx.
            x += _hx

        # Step y by _hy.
        y += _hy

    return trans_opt



# Do not consider this function so far.
def srwl_opt_load_mask(_FolderName, _AmpTraFileName, _OpPaDiFileName):
    # Load Mask Range and Size from files.
    f = open(_FolderName + '/' + _OpPaDiFileName, 'r')
    NumOfLines = 10
    MaxLenthOfOneLine = 100
    ValueList = []
    for LineNum in range(NumOfLines):
        LineStr = f.readline(MaxLenthOfOneLine)
        if (LineNum >= 1 and LineNum <= 9):
            EndNum = LineStr.find(' #')
            ValueList.extend([float(LineStr[1:EndNum])])
    f.close()
    mask_Rx = ValueList[4] - ValueList[3]
    mask_Ry = ValueList[7] - ValueList[6]
    mask_Nx = int(ValueList[5])
    mask_Ny = int(ValueList[8])

    # Generate OpT.
    trans_opt = SRWLOptT(_nx=mask_Nx, _ny=mask_Ny, _rx=mask_Rx, _ry=mask_Ry, _arTr=None, _extTr=0, _x=0, _y=0)

    # Load data from files.
    AmpTra = loadtxt(_FolderName + '/' + _AmpTraFileName, dtype="float")  # Amplitude Transmission.
    OpPaDi = loadtxt(_FolderName + '/' + _OpPaDiFileName, dtype="float")  # Optical Path Difference.

    pointer = 0  # pointer for array trans_opt.arTr
    for iy in range(mask_Ny):
        for ix in range(mask_Nx):
            trans_opt.arTr[pointer] = AmpTra[pointer / 2]  # Amplitude Transmission.
            trans_opt.arTr[pointer + 1] = OpPaDi[pointer / 2]  # Optical Path Difference.
            # Shift the pointer by 2.
            pointer += 2

    return trans_opt
