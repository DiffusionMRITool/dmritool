#!/bin/bash

mkdir -p temp
cd temp


# data generation {{{1
# generate data in 3 shells
b=1000,2000,3000
DWISimulator ../dwi_crossing.txt --outdwi dwi.nii.gz --outodf odfTrue.nii.gz --outeap eapTrue_r0.015.nii.gz --outrto rtoTrue.nii.gz --outmsd msdTrue.nii.gz --qorientations ../Elec060.txt --bvalues ${b} --rorientations ../directions_t4.txt --rvalues 0.015 --noisesigma 0.0 --outb0 dwi_diagonal_b0.nii.gz --outputdwitype EACHSHELL

ImageInfo dwi_b1000.nii.gz 
MeshFromSphericalFunctionTessellatedSamples eapTrue_r0.015.nii.gz eapTrue_r0.015_vis.vtk ../directions_t4.txt  --scale 8e-6
vtkviewer eapTrue_r0.015_vis.vtk &

MeshFromSphericalFunctionTessellatedSamples odfTrue.nii.gz odfTrue_vis.vtk ../directions_t4.txt  --scale 1.5
vtkviewer odfTrue_vis.vtk &

# SPFI {{{1

# estimate mean diffusivity
MeanDiffusivityEstimator dwi.txt D_sh4_ra1.nii.gz --sh 4 --ra 1

# estimate SPF coefficients 
SphericalPolarFourierImaging dwi.txt --sh 8 --ra 4 --lambdaSH 0 --lambdaRA 0 --signal signalSPF.nii.gz --radius 0.015 --estimation L1_DL --lambdaL1 1e-7  --mdImage D_sh4_ra1.nii.gz 

# estimate eap profiles
SPFToProfile signalSPF.nii.gz eap_r0.015.nii.gz --sh 8 --ra 4  --fourier --mdImage D_sh4_ra1.nii.gz
MeshFromSHCoefficients eap_r0.015.nii.gz eap_r0.015_vis.vtk --tessorder 4 --scale 8e-6
vtkviewer eap_r0.015_vis.vtk &

# estimate odfs
SPFToODF signalSPF.nii.gz odf.nii.gz --sh 8 --ra 4 --mdImage D_sh4_ra1.nii.gz
MeshFromSHCoefficients odf.nii.gz odf_vis.vtk --tessorder 4 --scale 1.5
vtkviewer odf_vis.vtk &

# estimate scalar maps
SPFToScalarMap  signalSPF.nii.gz rto.nii.gz  --mapType RTO --sh 8 --ra 4  --mdImage D_sh4_ra1.nii.gz
SPFToScalarMap  signalSPF.nii.gz msd.nii.gz  --mapType MSD --sh 8 --ra 4  --mdImage D_sh4_ra1.nii.gz
SPFToScalarMap  signalSPF.nii.gz pfa.nii.gz  --mapType PFA --sh 8 --ra 4  --mdImage D_sh4_ra1.nii.gz

cd ..

