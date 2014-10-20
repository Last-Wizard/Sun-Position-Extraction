# -*- coding: utf-8 -*-
# 2014.10.16

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from math import *

def get_sun(Hs=0.0,As=0.0):
	'''
	@Hs: elevation angle of the sun (0,90)
	@As: azimuth angle of the sun (-180,180)
	'''
	# maximum degree of polarization in Rayleigh model
	deltamax = 1
	# normalized radius of celestial sphere
	R = 1
	# transform the azimuth angle from 180 to 360
	# (-pi,-pi/2)->(pi/2,pi),(-pi/2,0)->(pi,pi*3/2)
	# (0,pi/2)->(pi*3/2,2*pi),(pi/2,pi)->(0,pi/2)
	# if As in (pi/2,pi), phis = As-pi/2
	# else phis = As+3/2*pi

	# thetas: zenith angle of the sun (90-Hs)
	thetas = pi/2 - Hs*pi/180.0
	As = As*pi/180
	# phis: azimuth angle of the sun (0,360)
	# relationship between sun and the sphere coordinate
	if pi/2 <= As <= pi:
		phis = As - pi/2
	else:
		phis = As + 1.5*pi

	# dop: degree of polarization
	# aop: angle of polarization	 
	dop = np.zeros((100,400)) # np.zeros(), dtype=float
	aop = np.zeros((100,400))
	x = np.zeros((100,400))
	y = np.zeros((100,400))

	# elevation angle 'el'
	# azimuth angle 'az'
	# polarization vector 'e_vector'
	el = np.zeros(3)
	az = np.zeros(3)
	el_az = np.zeros(2)
	e_vector_3d = np.zeros((100,400,3))
	e_vector_2d = np.zeros((100,400,2))

	Px, Py = 0, 0
	Nx, Ny = 0, 0
	D = 0

	# for loop	 
	i,j = 0,0

	# Ho: elevation angle of the observer
	# Ao: azimuth angle of the observer
	for Ho in np.linspace(0,pi/2,100,endpoint=True):
		for Ao in np.linspace(0,2*pi,400,endpoint=True):
			thetao = pi/2 - Ho
			phio = Ao # for simplification
			# if pi/2 <= Ao <= pi:
			# 	phio = Ao - pi/2
			# else:
			# 	phio = Ao + 1.5*pi
			# gamma: angular difference between sky point and the sun
			gamma = acos(sin(thetas)*sin(thetao)*cos(phio-phis)+cos(thetas)*cos(thetao))

			# dop: degree of polarization
			dop[i,j] = (deltamax*(sin(gamma))**2)/(1+cos(gamma)**2)

			# aop: angle of polarization
			if sin(phio-phis)*sin(thetas) == 0:
				aop[i,j] = pi/2
			else:
				aop[i,j] = atan((sin(thetao)*cos(thetas)-cos(thetao)*cos(phio-phis)*sin(thetas))/(sin(phio-phis)*sin(thetas)))

			el[0]  = -cos(thetao)*cos(phio)
			el[1]  = -cos(thetao)*sin(phio)
			el[2]  = sin(thetao)

			az[0]  = -sin(phio)
			az[1]  = cos(phio)
			az[2]  = 0

			el_az[0] = sin(thetao)*cos(phio)
			el_az[1] = sin(thetao)*sin(phio)

			e_vector_3d[i,j,:] = cos(aop[i,j])*el + sin(aop[i,j])*az

			e_vector_2d[i,j,0] = e_vector_3d[i,j,0]
			e_vector_2d[i,j,1] = e_vector_3d[i,j,1]

			Px += np.dot(e_vector_2d[i,j],el_az)*e_vector_2d[i,j,0]
			Py += np.dot(e_vector_2d[i,j],el_az)*e_vector_2d[i,j,1]

			x[i,j] = R*sin(thetao)*cos(phio)
			y[i,j] = R*sin(thetao)*sin(phio)

			j += 1

		i += 1
		j = 0

	# can consider the dop
	Nx = np.sum(e_vector_2d[:,:,1]**2)*Px - np.sum(e_vector_2d[:,:,0]*e_vector_2d[:,:,1])*Py
	Ny = np.sum(e_vector_2d[:,:,0]**2)*Py - np.sum(e_vector_2d[:,:,0]*e_vector_2d[:,:,1])*Px
	D  = np.sum(e_vector_2d[:,:,0]**2)*np.sum(e_vector_2d[:,:,1]**2) - np.sum(e_vector_2d[:,:,0]*e_vector_2d[:,:,1])**2

	# azimuth angle of the sun
	phi0 = atan2(Ny,Nx)*180/pi

	if phi0 < 0:
		phi0 += 360
	else:
		pass

	print '-----------------------------------------------------------'
	print 'the original azimuth angle of the sun (0,360):',(phis*180/pi)
	print '-----------------------------------------------------------'
	print 'the processed azimuth angle of the sun (0,360):',phi0

# main function
if __name__=='__main__':
	print 'please input the elevation angle of the sun (0,90):',
	Hs = float(raw_input())
	print 'please input the azimuth angle of the sun (-180,180):',
	As = float(raw_input())
	get_sun(Hs,As)