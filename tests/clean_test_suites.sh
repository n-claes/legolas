rm -rf .pytest_cache
rm -rf regression_tests/results
rm -rf regression_tests/.pytest_cache
echo ">> regression tests cleaned"
rm -rf pylbo_tests/results
rm -rf pylbo_tests/tmp
rm -rf pylbo_tests/.pytest_cache
echo ">> pylbo tests cleaned"
cd unit_tests || exit
rm -rf build
rm -f test_legolas
echo ">> unit tests cleaned"
