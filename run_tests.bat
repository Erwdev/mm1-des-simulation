#!/bin/bash
# filepath: e:\SEMESTER 5\SCIENCE SIMULATIONO\MM1_DES_SIMULATION\mm1-des-simulation\run_tests.sh

echo "========================================"
echo "  Building and Running Unit Tests"
echo "========================================"
echo ""

# Compile test file
echo "[1/2] Compiling test_des.cpp..."
g++ -std=c++17 -O2 -Wall -o test_des test_des.cpp

if [ $? -ne 0 ]; then
    echo "ERROR: Compilation failed!"
    exit 1
fi

echo "[2/2] Running tests..."
echo ""
./test_des

if [ $? -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "  All tests completed successfully!"
    echo "========================================"
else
    echo ""
    echo "========================================"
    echo "  Some tests failed - check output"
    echo "========================================"
    exit 1
fi