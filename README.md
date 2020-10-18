# MCS-on-2D-Ising-Model
This is a C++ program to calculate some physical quantities in ferromagnetic phase transition of 2D ising model and some related materials .
The size of the system is given by the parameter you pass from the command line.
For example you can input ./a.out 20 on the command line,
then you create a configuration containing 20*20 spins.

The program contains a for loop for temperature in order to do calculation at various temperature,
but for some quantities such as correlation function,
a few temperature will be enough,
so you can modify the loop scope appropriately when calculating them.

The codes were developed by Yuliang Jin's research group (http://home.itp.ac.cn/~yuliangjin/) 
at Institute of Theoretical Physics, Chinese Academy of Sciences. 
Please contact the developers,Xinyang Li (lixinyang@mail.itp.ac.cn) 
or Dr. Yuliang Jin (yuliangjin@mail.itp.ac.cn) for questions.