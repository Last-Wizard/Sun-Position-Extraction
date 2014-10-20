Sun-Position-Extraction
=======================
>version: 1.0

>last_update_time: 2014-10-20

>reference:

>1. Stürzl W, Carey N. A fisheye camera system for polarisation detection on UAVs[C]//Computer Vision–ECCV 2012. Workshops and Demonstrations. Springer Berlin Heidelberg, 2012: 431-440.
>2. Pomozi I, Gál J, Horváth G, et al. Fine structure of the celestial polarization pattern and its temporal change during the total solar eclipse of 11 August 1999\[J\]. Remote sensing of Environment, 2001, 76(2): 181-201.

###利用全天域大气偏振模式的太阳位置获取方法
####algorithms to extract the sun position using skylight polarization pattern

>二维平面中太阳的位置获取算法 (2-dimensional,get_sun_2d.py)
>>仅获取太阳的方位角信息 (only get the azimuth angle of the sun)

>三维空间中太阳的位置获取算法 (3-dimensional,get_sun_3d.py)
>>获取太阳的方位角, 高度角信息 (can get the azimuthal angle & elevation angle of the sun)

###Environment
>python 2.7.8, numpy 1.8.1, matplotlib 1.3.1, spyder 2.3.0
