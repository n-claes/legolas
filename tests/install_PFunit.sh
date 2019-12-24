git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git
cd pFUnit || exit
mkdir build
cd build || exit
cmake ..
make tests
make install

GREEN='\033[0;32m'
NC='\033[0m'
echo "==================================================="
echo "Installation completed."
echo "${GREEN}export PFUNIT_DIR=$PWD/installed${NC}"
echo "still needs to be added to your bash/zsh profile."
echo "==================================================="
echo ""