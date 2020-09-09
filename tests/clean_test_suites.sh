rm -rf .pytest_cache
rm -rf regression_tests/results
rm -rf regression_tests/.pytest_cache
echo ">> regression tests cleaned"
rm -f pylbo_tests/*.dat
rm -f pylbo_tests/*.log
rm -rf pylbo_tests/.pytest_cache
echo ">> pylbo tests cleaned"
cd core_tests || exit
rm -rf build
rm -f test_legolas
echo ">> core tests cleaned"
