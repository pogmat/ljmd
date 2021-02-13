


OPENMPI_MAJOR_VER="4.1"
OPENMPI_MINOR_VER="4.1.0"

OPENMPI_DIR="/home/${USER}/openmpi/${OPENMPI_MINOR_VER}"


#download and install openmpi

VER="v${OPENMPI_MINOR_VER}"
TARFILE="openmpi-${OPENMPI_MINOR_VER}.tar.gz"

URl=" https://download.open-mpi.org/release/open-mpi/v${OPENMPI_MAJOR_VER}/${TARFILE}"


mkdir -p ${OPENMPI_DIR}
cd ${OPENMPI_DIR}

mkdir -p "install"


wget ${URl}

#untar openmpi
tar -xvf ./${TARFILE}

#configure openmpi
#./${OPENMPI_DIR}/install/configure --prefix="${OPENMPI_DIR}/install"

#install openmpi
#make -j 2 all
#make install
