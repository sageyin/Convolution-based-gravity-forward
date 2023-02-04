# Convolution-based-gravity-forward
The code for fast computation of gravity field based on circular convolution. The code includes our proposed method and the traditional spatial domain method that can be used to calculate the gravity of a density body and its tensor.

1、 Our methond Fast computation of gravity field based on circular convolution. Main function: Cal_ModelGravityinFourie.m Calling sub-functions:
Cal_tranGraf.m : Analytic formula method for calculating the gravity of a cube
GraconvelP.m : Construct the circular kernel matrix and calculate the gravity field using FFT algorithm

2、In space domain Main function: Cal_ModelGravityinFourie.m Calling sub-functions:
Cal_tranGraf.m : Analytic formula method for calculating the gravity of a cube

Note: See code comments for detailed parameters.

The article link: https://authors.elsevier.com/c/1gXR8MMTPor16
Xianzhe Yin, Changli Yao, Yuanman Zheng, Wenqiang Xu, Guangxi Chen, Xiaoyu Yuan, 2023. A fast 3D gravity forward algorithm based on circular convolution. Computers & Geosciences 172, 105309. https://doi.org/10.1016/j.cageo.2023.105309
