import os
import subprocess
import sys
import numpy as np
import matplotlib.pyplot as plt

def install_packages():
    try:
        import numpy
        import matplotlib
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "numpy", "matplotlib"])

def compile_fortran():
    try:
        subprocess.run(["gfortran", "cycDeformationModelling.f", "-o", "aaa"], check=True)
    except subprocess.CalledProcessError:
        print("Error compiling Fortran code")
        sys.exit(1)

def run_fortran_executable():
    try:
        subprocess.run(["./aaa"], check=True)
    except subprocess.CalledProcessError:
        print("Error running Fortran executable")
        sys.exit(1)

def plot_data():
    try:
        data = np.loadtxt('se3.dat')
        x = data[:, 1]
        y = data[:, 4]

        plt.plot(x, y, marker='o')
        plt.xlabel('Column 2')
        plt.ylabel('Column 5')
        plt.title('Cyclic Stress-Strain Curves')
        plt.grid(True)
        plt.show()
    except Exception as e:
        print(f"Error plotting data: {e}")
        sys.exit(1)

def main():
    install_packages()
    compile_fortran()
    run_fortran_executable()
    plot_data()

if __name__ == "__main__":
    main()
