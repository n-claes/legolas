rm -rf .pytest_cache
rm -f regression_tests/*.dat
rm -f regression_tests/*.log
rm -rf regression_tests/.pytest_cache
echo ">> regression tests cleaned"
rm -f pylbo_tests/*.dat
rm -f pylbo_tests/*.log
rm -rf pylbo_tests/.pytest_cache
echo ">> pylbo tests cleaned"
cd core_tests || exit
make clean
echo ">> core tests cleaned"


