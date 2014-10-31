2pg_cartesian
=============

is a framework of optimization algorithms for structural biology. 
The goal of this framework is allow to build optimization algorithms for 
structural biology area. 

INSTALATION
=============
2pg_cartesian uses cmake 2.8 or later to install. Therefore, you are 
able to execute a simple cmake process as described below:

mkdir build
cd build
cmake ..
make 
make install

If you want to assign an install directory, you may use 
DCMAKE_INSTALL_PREFIX option.
For example, you want to install 2pg_cartesian at 
/home/faccioli/Programs/ directory, you have to use:

mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/home/faccioli/Programs/2pg_cartesian/
make 
make install

After make install command, will be created bin and scripts directories 
where was specified. 

