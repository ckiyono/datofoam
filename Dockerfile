FROM datofoam:base

USER root

# DAFoam root path
ENV DAFOAM_ROOT_PATH=/home/dafoamuser/dafoam

# copy dafoam code
COPY dafoam $DAFOAM_ROOT_PATH/repos/dafoam/

# copy rhoSimpleFoamTopOpt
COPY rhoSimpleFoamTopOpt $DAFOAM_ROOT_PATH/OpenFOAM/datofoam/applications/solvers/rhoSimpleFoamTopOpt/

RUN chmod -R 777 $DAFOAM_ROOT_PATH/OpenFOAM/datofoam/

WORKDIR $DAFOAM_ROOT_PATH/repos/dafoam

# Swith to dafoamuser
USER dafoamuser

# install datofoam
RUN pip install .

# install rhoSimpleFoamTopOpt
WORKDIR $DAFOAM_ROOT_PATH/OpenFOAM/datofoam/applications/solvers/rhoSimpleFoamTopOpt/
RUN . $DAFOAM_ROOT_PATH/OpenFOAM/OpenFOAM-v1812/etc/bashrc && \
    wmake

#install rhoSimpleFoamPen
WORKDIR $DAFOAM_ROOT_PATH/OpenFOAM/datofoam/applications/solvers/rhoSimpleFoamTopOpt/rhoSimpleFoamPen/
RUN . $DAFOAM_ROOT_PATH/OpenFOAM/OpenFOAM-v1812/etc/bashrc && \
    wmake

#install petscSolver
WORKDIR $DAFOAM_ROOT_PATH/OpenFOAM/datofoam/applications/solvers/rhoSimpleFoamTopOpt/petscSolver/
RUN . $DAFOAM_ROOT_PATH/OpenFOAM/OpenFOAM-v1812/etc/bashrc && \
    wmake