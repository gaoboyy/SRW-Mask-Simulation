###################################################################################
# SRWLIB ExampleMetro: 
# v 0.11
# Date: 2014.03.25
# Authors:
# - Revision
###################################################################################

# ********************** Import classes: **********************
from __future__ import print_function #Python 2.7 compatibility
from metro_mask import *
import os
import random
import copy
import math

print('SRWLIB Python Example on Metrology:')
# ********************** Setting input and output folders and files: **********************
print('Setting input and output folders and files ...', end='')

## Set data sub-folder name.
strOutDataFolderName = 'data_example_metro_output'
strInpDataFolderName = 'data_example_metro_input'

## Set data files.
# Aberration
strAmpTraAberrationFileName = 'aberration_Amp.dat'      # Set file name for aberration amplitude transmission.
strOpPathAberrationFileName = 'aberration_OPD.dat'      # Set file name for aberration optical path difference.
# mask
strAmpTraOutFileName = 'ex_mask_amp_transmission.dat'   # Set file name for mask amplitude transmission.
strOpPathOutFileName = 'ex_mask_optical_path_dif.dat'   # Set file name for mask optical path difference.
# position_00
strIntenOutFileName_00 = 'ex_mask_res_inten_00.dat'     # Set file name for output intensity data at position 00.
strPhaseOutFileName_00 = 'ex_mask_res_phase_00.dat'     # Set file name for output phase data     at position 00.
# position_01
strIntenOutFileName_01 = 'ex_mask_res_inten_01.dat'     # Set file name for output intensity data at position 01.
strPhaseOutFileName_01 = 'ex_mask_res_phase_01.dat'     # Set file name for output phase data     at position 01.


print('done.')

# ********************** Input parameters and structures: **********************
print('Setting parameters and structures ...', end='')

# Number of pixels in mask and beam matrix.
intNx = 512 # At least 8 times higher than number of period to get the double-frequency fringe.
intNy = 512 # At least 8 times higher than number of period to get the double-frequency fringe.


GsnBm = SRWLGsnBm() #Gaussian Beam structure (just parameters)
GsnBm.x = 0 #Transverse Coordinates of Gaussian Beam Center at Waist [m]
GsnBm.y = 0
GsnBm.z = 0 #Longitudinal Coordinate of Waist [m]
GsnBm.xp = 0 #Average Angles of Gaussian Beam at Waist [rad]
GsnBm.yp = 0
GsnBm.avgPhotEn = 12400 #5000 #Photon Energy [eV]
GsnBm.pulseEn = 0.001 #Energy per Pulse [J] - to be corrected
GsnBm.repRate = 1 #Rep. Rate [Hz] - to be corrected
GsnBm.polar = 1 #1- linear hoirizontal
GsnBm.sigX = 2e-09/2.35 #23e-06/2.35 #Horiz. RMS size at Waist [m]
GsnBm.sigY = GsnBm.sigX #Vert. RMS size at Waist [m]
GsnBm.sigT = 10e-15 #Pulse duration [fs] (not used?)
GsnBm.mx = 0 #Transverse Gauss-Hermite Mode Orders
GsnBm.my = 0


wfr = SRWLWfr() #Initial Electric Field Wavefront
wfr.allocate(1, intNx, intNy) #Numbers of points vs Photon Energy (1), Horizontal and Vertical Positions (dummy)
wfr.mesh.zStart = 1e6 #300 #Longitudinal Position [m] at which Electric Field has to be calculated, i.e. the position of the first optical element
wfr.mesh.eStart = GsnBm.avgPhotEn #Initial Photon Energy [eV]
wfr.mesh.eFin = GsnBm.avgPhotEn #Final Photon Energy [eV]
firstHorAp = 2*1.024e-03 #First Aperture [m]
firstVertAp = 2*1.024e-03 #[m] 
wfr.mesh.xStart = -0.5*firstHorAp #Initial Horizontal Position [m]
wfr.mesh.xFin = 0.5*firstHorAp #Final Horizontal Position [m]
wfr.mesh.yStart = -0.5*firstVertAp #Initial Vertical Position [m]
wfr.mesh.yFin = 0.5*firstVertAp #Final Vertical Position [m]

wfr.partBeam.partStatMom1.x = GsnBm.x #Some information about the source in the Wavefront structure
wfr.partBeam.partStatMom1.y = GsnBm.y
wfr.partBeam.partStatMom1.z = GsnBm.z
wfr.partBeam.partStatMom1.xp = GsnBm.xp
wfr.partBeam.partStatMom1.yp = GsnBm.yp

sampFactNxNyForProp = 0 #sampling factor for adjusting nx, ny (effective if > 0)
arPrecPar = [sampFactNxNyForProp]

print('done.')

# ********************** Calculating Initial Wavefront and extracting Intensity: **********************

srwl.CalcElecFieldGaussian(wfr, GsnBm, arPrecPar)
arI0 = array('f', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI0, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0) #extracts intensity
srwl_uti_save_intens_ascii(arI0, wfr.mesh, os.path.join(os.getcwd(), strOutDataFolderName, strIntenOutFileName_00), 0)

arI0x = array('f', [0]*wfr.mesh.nx) #array to take 1D intensity data
srwl.CalcIntFromElecField(arI0x, wfr, 6, 0, 1, wfr.mesh.eStart, 0, 0) #extracts intensity

arP0 = array('d', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" array to take 2D phase data (note it should be 'd')
srwl.CalcIntFromElecField(arP0, wfr, 0, 4, 3, wfr.mesh.eStart, 0, 0) #extracts radiation phase
srwl_uti_save_intens_ascii(arP0, wfr.mesh, os.path.join(os.getcwd(), strOutDataFolderName, strPhaseOutFileName_00), 0, ['', 'Horizontal Position', 'Vertical Position', 'Phase'], _arUnits=['', 'm', 'm', 'rad'])

mesh0 = deepcopy(wfr.mesh)

# ********************** Optical Elements: **********************

# load 2D optical aberration matrices
print('Loading the aberration matrices ...', end='')
opAberr = srwl_opt_load_mask(strInpDataFolderName,strAmpTraAberrationFileName,strOpPathAberrationFileName)
print('done')


print('Setting-up the Mask ...')
delta = 2.21555115e-06  #2.21555115e-06 #Be @ 12400 eV
                        #1.3657359E-05 #Be @ 5000 eV

attenLen = 15453.7e-06  #[m]
                        #15453.7e-06
                        #1268.65e-06
                        #1.0e-07

wavelength = (1239.852/GsnBm.avgPhotEn*1e-09)
print('The wavelength =', wavelength, 'm')

thickness = wavelength/(2*delta)
print('The thickness =', thickness, 'm')

holeSize = (32.e-06)/2 #[m]
pitch_x = holeSize*4
pitch_y = holeSize*4
maskNx = intNx
maskNy = intNy
numHoles = 3
holeAngle=math.pi/4
# Generate a 2D Mask.
opMask = srwl_opt_setup_mask(_delta=delta, _atten_len=attenLen, _thick=thickness, _hole_sh=1, _hole_dim1=holeSize, _hole_dim2=holeSize, _pitch_x=pitch_x, _pitch_y=pitch_y, _hole_nx=numHoles, _hole_ny=numHoles,_mask_Nx=maskNx,_mask_Ny=maskNy,_angle=holeAngle)
print('Done!')

# Extracting transmission data characteristic for subsequent plotting.
OpticalPathDifMask = opMask.get_data(3, 3)
AmpTransmitionMask = opMask.get_data(2, 3)
# Saving transmission data to file.
srwl_uti_save_intens_ascii(OpticalPathDifMask, opMask.mesh, os.path.join(os.getcwd(), strOutDataFolderName, strOpPathOutFileName), 0, ['', 'Horizontal Position', 'Vertical Position', 'Opt. Path Diff.'], _arUnits=['', 'm', 'm', 'm'])
srwl_uti_save_intens_ascii(AmpTransmitionMask, opMask.mesh, os.path.join(os.getcwd(), strOutDataFolderName, strAmpTraOutFileName), 0, ['', 'Horizontal Position', 'Vertical Position', 'Opt. Path Diff.'], _arUnits=['', 'm', 'm', 'm'])
opPathDifMaskx = opMask.get_data(3, 1, _y=0)
opPathDifMasky = opMask.get_data(3, 2, _x=0)
opMeshMask = opMask.mesh

# Generate aperture at Mask.
diamMask_x = (numHoles+1)*pitch_x # Mask diameter in x [m]
diamMask_y = (numHoles+1)*pitch_y # Mask diameter in y [m]
opApMask = SRWLOptA('r', 'a', diamMask_x, diamMask_y)

# Calculate fractional Talbot distance for a "plane wave configuration".
TD = 2*pitch_x*pitch_x/wavelength
D = TD/16
print('Fractional TD for a plane wavefront, D =', D, 'm')

# Calculate fractional Talbot distance for a spherical wavefront.
R = wfr.mesh.zStart # radius of the coming spherical wavwavefronte.
d = R*D/(R+D)
print('Fractional TD for a spherical wavefront, d =', d, 'm')

# Drift space from Mask to a plane where an image could be observed
opDrMask = SRWLOptD(d)

#*********** Wavefront Propagation Parameters:
#[0]: Auto-Resize (1) or not (0) Before propagation
#[1]: Auto-Resize (1) or not (0) After propagation
#[2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
#[3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
#[4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
#[5]: Horizontal Range modification factor at Resizing (1. means no modification)
#[6]: Horizontal Resolution modification factor at Resizing
#[7]: Vertical Range modification factor at Resizing
#[8]: Vertical Resolution modification factor at Resizing
#[9]: Type of wavefront Shift before Resizing (not yet implemented)
#[10]: New Horizontal wavefront Center position after Shift (not yet implemented)
#[11]: New Vertical wavefront Center position after Shift (not yet implemented)
#prParRes1 = [0, 0, 1., 0, 0, 8.0, 18.0, 8.0, 18.0, 0, 0, 0]
prParRes1 = [0, 0, 1., 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]

prParRes0 = [0, 0, 1., 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]

#"Beamline" - Container of Optical Elements (together with the corresponding wavefront propagation instructions)

##optBL = SRWLOptC([opApMask,  opMask],
##                 [prParRes1, prParRes0]
##                )

##optBL = SRWLOptC([opApMask,  opMask,    opDrMask],
##                 [prParRes1, prParRes0, prParRes0]
##                )

##optBL = SRWLOptC([opApMask,  opAberr],
##                 [prParRes1, prParRes0]
##                )

optBL = SRWLOptC([opApMask,  opAberr,   opMask,     opDrMask],
                 [prParRes1, prParRes0, prParRes0,  prParRes0]
                )

#********************** Wavefront Propagation.
print('Propagating Wavefront (through Mask and Drift)...', end='')
srwl.PropagElecField(wfr, optBL)
print('done')
print('Saving resulting data to files...', end='')
mesh1 = deepcopy(wfr.mesh)

arI1 = array('f', [0]*mesh1.nx*mesh1.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI1, wfr, 6, 0, 3, mesh1.eStart, 0, 0) #extracts intensity
srwl_uti_save_intens_ascii(arI1, mesh1, os.path.join(os.getcwd(), strOutDataFolderName, strIntenOutFileName_01), 0)

yCut = 0.5*pitch_y
arI1x = array('f', [0]*mesh1.nx) #array to take 1D intensity data
srwl.CalcIntFromElecField(arI1x, wfr, 6, 0, 1, mesh1.eStart, 0, yCut) #extracts intensity

xCut = 0.5*pitch_x
arI1y = array('f', [0]*mesh1.ny) #array to take 1D intensity data
srwl.CalcIntFromElecField(arI1y, wfr, 6, 0, 2, mesh1.eStart, xCut, 0) #extracts intensity

arP1 = array('d', [0]*mesh1.nx*mesh1.ny) #"flat" array to take 2D phase data (note it should be 'd')
srwl.CalcIntFromElecField(arP1, wfr, 0, 4, 3, mesh1.eStart, 0, 0) #extracts radiation phase
srwl_uti_save_intens_ascii(arP1, mesh1, os.path.join(os.getcwd(), strOutDataFolderName, strPhaseOutFileName_01), 0)

print('done')

#********************** Plotting results (requires 3rd party graphics package)
#print('Plotting the results (blocks script execution; close any graph windows to proceed) ... ', end='')

# Input Wavefront
plotMesh0x = [1000*mesh0.xStart, 1000*mesh0.xFin, mesh0.nx]
plotMesh0y = [1000*mesh0.yStart, 1000*mesh0.yFin, mesh0.ny]
uti_plot2d(arI0, plotMesh0x, plotMesh0y, ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity Before Propagation'])
uti_plot1d(arI0x, plotMesh0x, ['Horizontal Position [mm]', 'Intensity [a.u.]', 'Intensity Before Propagation\n(cut vs horizontal position at y = 0)'])
uti_plot2d(arP0, plotMesh0x, plotMesh0y, ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Phase Before Propagation'])

# Mask
maskPlotMeshX = [1000*opMeshMask.xStart, 1000*opMeshMask.xFin, opMeshMask.nx]
maskPlotMeshY = [1000*opMeshMask.yStart, 1000*opMeshMask.yFin, opMeshMask.ny]
uti_plot2d(AmpTransmitionMask, maskPlotMeshX, maskPlotMeshY, ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Amplitude Transmission in Mask'])
uti_plot2d(OpticalPathDifMask, maskPlotMeshX, maskPlotMeshY, ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Optical Path Difference in Mask'])

# Propagated Wavefront
plotMesh1x = [1000*mesh1.xStart, 1000*mesh1.xFin, mesh1.nx]
plotMesh1y = [1000*mesh1.yStart, 1000*mesh1.yFin, mesh1.ny]
uti_plot2d(arI1, plotMesh1x, plotMesh1y, ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity After Propagation'])
uti_plot1d(arI1x, plotMesh1x, ['Horizontal Position [mm]', 'Intensity [a.u.]', 'Intensity After Propagation\n(cut vs horizontal position at y = 0.5*pitch_y)'])
uti_plot1d(arI1y, plotMesh1y, ['Vertical Position [mm]', 'Intensity [a.u.]', 'Intensity After Propagation\n(cut vs vertical position at x = 0.5*pitch_x)'])
uti_plot2d(arP1, plotMesh1x, plotMesh1y, ['Horizontal Position [microns]', 'Vertical Position [microns]', 'Phase After Propagation'])
uti_plot_show() #show all graphs (blocks script execution; close all graph windows to proceed)
##print('done')

# sys.exit(0)



