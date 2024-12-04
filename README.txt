*****Crustal magnetic field model******

Description
The equivalent source dipoles method was used to model the remanent field.  To evaluate the crustal magnetic field, B_JM and B_IN were removed from the observation data. By employing the conjugate gradient technique , the optimum m can be numerically calculated. An example for 1271 dipoles with  a buried depth of 60 km are shown as follow.


Input
data_R_fun.mat : Data set after removal of B_JM and B_IN in the longitude range 270-330°E（in planetocentric coordinate, see Li et al., 2020, https://10.26464/epp2020045）
Each column correspond to: Longitude, Latitude, Height, Bx (By in CphiO), By(-Bx in CphiO), Bz(Bz in CphiO) in  planetocentric coordinate, spacecraft position in the x-direction,  y-direction and z-direction

data_R_fun2.mat : Data set after removal of B_JM and B_IN in the range of all longitudes for all three flybys


Output
Crustal magnetic field three components and total field in planetocentric coordinate.


Running 
The basic equation for this model can be expressed simply as: b=Gm+v 
First, running the "main_matrix_R.m" to calculate the G value
Then, running the "main_gongetidu_3D_2.m" to calculate the m value
Finally, running the "magnetic_components_for_Valhalla.m" to calculate the b value

The "magnetic_components_for_Valhalla2.m" calculates the inclination and declination
The "Induced_field_Callisto.m" calculates the Induced magnetic field (B_IN) in CphiO cordinate 


Please refer for any questions  to:
Jianhong Luo, School of Earth Sciences, China University of Geosciences, Wuhan, Hubei, China, jianhong@cug.edu.cn
Zhaojin Rong, Key Laboratory of Earth and Planetary Physics, Institute of Geology and Geophysics, Chinese Academy of Sciences, Beijing, China, rongzhaojin@mail.iggcas.ac.cn