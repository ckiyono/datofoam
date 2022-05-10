FROM datofoam:base

# DAFoam root path
ENV DAFOAM_ROOT_PATH=/home/dafoamuser/dafoam

# copy dafoam code
COPY dafoam $DAFOAM_ROOT_PATH/repos/dafoam/

WORKDIR $DAFOAM_ROOT_PATH/repos/dafoam

# Swith to dafoamuser
USER dafoamuser

RUN pip install .
