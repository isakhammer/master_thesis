# Setting the flag
## -DUSE_SUPERBUILD=OFF \
# activates MPI no matter how USE_MPI was set
#
# Setting the flag
## -DUSE_SUPERBUILD=ON \
## -DUSE_MUMPS=ON \
# gives an parmetis not found error and also activates MPI
# no matter how USE_MPI was set

mkdir -p $BUILD_DIR
mkdir -p $INSTALL_DIR
cd $BUILD_DIR

cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
      -DUSE_JPEG=ON \
      -DUSE_MPEG=ON \
      -DUSE_OCC=ON \
      -DUSE_GUI=ON \
      -DUSE_MPI=${USE_MPI} \
      -DUSE_MUMPS=OFF \
      -DUSE_HYPRE=OFF \
      -DUSE_SUPERBUILD=ON \
      ${BASEDIR}/src

#exit 0
make -j8 VERBOSE=1

make -j8 install

PYTHONPATH_TMP=`python3 -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(1,0,''))"`

###### Uncomment if you want to generate a module file
MODULE_FILE=ngsolve-${BUILD_VERSION}
echo "#%Module 1.0"  > ./$MODULE_FILE
echo "set MODULE_DIR \"${INSTALL_DIR}\"" >> ./$MODULE_FILE
echo "setenv NETGENDIR \"\$MODULE_DIR/bin\"" >> ./$MODULE_FILE
echo "prepend-path PATH \"\$MODULE_DIR/bin\""  >> ./$MODULE_FILE
echo "prepend-path PYTHONPATH \"\$MODULE_DIR/${PYTHONPATH_TMP}\""  >> ./$MODULE_FILE

#source ~/.bashrc
#cd ${BASEDIR}/${BUILD_VERSION}ngsolve-install/share/ngsolve/py_tutorials/intro
#netgen navierstokes.py

