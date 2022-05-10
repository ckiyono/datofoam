# DATOFoam - Discrete Adjiont Topology Optimization

Topology optimization software with the discrete adjoint method from DAFoam.

## Instalation
 - Clone DATOFOam_base:
 ```console
 git clone https://github.com/ckiyono/datofoam_base.git
 ```

 - cd into folder:
 ```console
 cd datofoam-base
 ```

 - Generate docker image of DATOFoam_base:
 ```console
 sudo docker build -t datofoam:base .
 ```

 - Leave folder:
 ```console
 cd ..
 ```

 - Clone this repository:
 ```console
 git clone https://github.com/ckiyono/datofoam.git
 ```

 - cd into folder:
 ```console
 cd datofoam
 ```

 - Generate docker image of DATOFoam:
 ```console
 sudo docker build -t datofoam .
 ```

## Running DATOFoam
 - cd into an example folder:
 ```console
 cd /path/to/example/folder
 ```

 - Run a container with datofoam image:
 ```console
 docker run -it --rm -u dafoamuser --mount "type=bind,src=$(pwd),target=/home/dafoamuser/mount" -w /home/dafoamuser/mount datofoam bash
 ```

 Inside the container, scripts may not have persmission to run. In this case:
 - Leave mount folder:
 ```console
 cd ..
 ```

 - Change permissions:
 ```console
 chmod -R 777 mount
 ```

 - cd to mount folder
 ```console
 cd mount
 ```

 - Run the preProcessing.sh script to generate the mesh:
 ```console
 ./preProcessing.sh
 ```

 - Then, use the following command to run the optimization with 4 CPU cores:
 ```console
 mpirun -np 4 python runScript_v2.py 2>&1 | tee logOpt.txt
 ```

 The option '--task=' can be used after 'runScript_v2.py' to choose the task:
  - "opt" (default)
  - "runPrimal"
  - "runAdjoint"
  - "calcPsi"
