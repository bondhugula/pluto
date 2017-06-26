#!/bin/bash

BASE_DIR=`pwd`
PLUTO_DIR=${BASE_DIR}/pluto
PLUTO_INSTALL_DIR=${BASE_DIR}/install
LLVM_GIT=${BASE_DIR}/llvm_git
LLVM_BUILD=${BASE_DIR}/llvm_build
LLVM_INSTALL=${BASE_DIR}/llvm_install
POLLY_PATCH_FILE=${BASE_DIR}/polly.patch

LLVM_COMMIT=c40271cb7545ab19216ea964c6b3607ed825a59e
CLANG_COMMIT=b38543d2c36e7e905dbdb3ec3ace2167421b562d
POLLY_COMMIT=5556f29626b25cbba0059431c11355064ce14ffe

echo "Installing Pluto..."
# Install pluto (libpluto)
git clone https://github.com/bondhugula/pluto.git -b libpluto
cd ${PLUTO_DIR}
git submodule init
git submodule update
./apply_patches.sh
./autogen.sh
./configure --enable-debug --prefix=${PLUTO_INSTALL_DIR}
make -j4 install;
	
echo "Installing Polly with pluto support."
# Clone LLVM,Polly, clang
git clone http://llvm.org/git/llvm.git ${LLVM_GIT}
git clone http://llvm.org/git/polly.git ${LLVM_GIT}/tools/polly
git clone http://llvm.org/git/clang.git ${LLVM_GIT}/tools/clang

# Checking out to specific versions of LLVM, CLANG, POLLY
cd ${LLVM_GIT} && git checkout ${LLVM_COMMIT} && cd ${BASE_DIR}
cd ${LLVM_GIT}/tools/polly && git checkout ${POLLY_COMMIT} && cd ${BASE_DIR}
cd ${LLVM_GIT}/tools/clang && git checkout ${CLANG_COMMIT} && cd ${BASE_DIR}

# Applying the patch
cd ${LLVM_GIT}/tools/polly
echo "Getting polly.patch..."
wget pluto-compiler.sourceforge.net/polly.patch
git am polly.patch

# Build LLVM/Polly with Pluto support
mkdir ${LLVM_BUILD} ${LLVM_INSTALL}
cd ${LLVM_BUILD}
cmake -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=${LLVM_INSTALL} -DCMAKE_BUILD_TYPE=Release -DLLVM_TARGETS_TO_BUILD="X86" ${LLVM_GIT} -DPLUTO_INSTALL_DIR=${PLUTO_INSTALL_DIR}
make -j4 install
