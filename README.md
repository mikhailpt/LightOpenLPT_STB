# ShakeTheBox

NOTICE: Thanks for visiting STB repository. The lastest version of OpenLPT has been uploaded, in addition to a first version of tutorial and the sample data "SD00125". The first version of tutorial introduces how to compile the code on both Windows and Linux system and how to test the code with the sample data set. For later version of tutorial which may be uploaded by the middle of Feburary, we will show how to obtain a better calibration and how to adjust the parameters in the configure file to achieve the best performance of the code.    

ShakeTheBox provides a C++ code for processing images in order to obtain tracks of particles seeded in the flow. Usually at least three cameras are needed to reconstruct 3D tracks. Anyone who uses this code and develop this code for their research should cite the following publication in their work：

Tan, S., Salibindla, A., Masuk, A.U.M. and Ni, R. Introducing OpenLPT: new method of removing ghost particles and high-concentration particle shadow tracking. Exp Fluids 61, 47 (2020). https://doi.org/10.1007/s00348-019-2875-2

Tan, S., Salibindla, A., Masuk, A.U.M. and Ni, R., 2019. An open-source shake-the-box method and its performance evaluation. In 13th International Symposium on Particle Image Velocimetry.

We really welcome any research group that is willing to improve the code with us. Please send us email to join as contributor to the code. 

Instructions on how to use this code are listed as follows:

(1) Preprocessing images, to make particles more bright and images less noise. A sample code of processing images can be found in  ./Data_analysis_process/PreprocessImage.m. Since different camera gives different images, it is not suggested to use this sample code directly, but to use it as a reference.

(2) Setting up a project folder, and put the images under that folder. It is suggested to set up folder like cam1, cam2, cam3 and so on under the project folder, and place the corresponding processed images beneath each camera folder.

(3) Creat a camXImageNames.txt under the project folder to list the relative path to the images. For cam1, cam1ImageNames.txt should be like:

cam1/cam1frame00001.tif

cam1/cam1frame00002.tif

cam1/cam1frame00003.tif

cam1/cam1frame00004.tif
