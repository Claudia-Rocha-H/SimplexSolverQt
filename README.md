# SimplexSolverQt

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![GitHub last commit](https://img.shields.io/github/last-commit/Claudia-Rocha-H/SimplexSolverQt)](https://github.com/Claudia-Rocha-H/SimplexSolverQt/commits/main)

## Project Overview

SimplexSolverQt is a  desktop application designed to solve Linear Programming (LP) problems efficiently. Built with C++ and the Qt framework, it provides a graphical interface for defining and solving optimization problems. A key feature of this application is its  visualization capabilities, offering a 3D graphical representation of the feasible region for problems with three variables, alongside comprehensive sensitivity analysis.

## Features

* **Linear Programming Solver:** Solves maximization and minimization problems using the Simplex Method.
* **Intuitive User Interface:** Powered by Qt, providing a clean and easy-to-use experience for inputting objective functions and constraints.
* **3D Feasible Region Visualization:** For problems with three decision variables ($x_1, x_2, x_3$), the application dynamically generates and displays a 3D plot of the feasible region, allowing users to visually grasp the solution space. This visualization clearly depicts the planes defined by each constraint (labeled as S1, S2, etc., corresponding to slack variables) and the vertices of the feasible polyhedron.
* **Optimal Solution Display:** Clearly presents the optimal value of the objective function and the values of decision variables at the optimum.
* **Sensitivity Analysis:** Provides detailed insights into how changes in the problem parameters (e.g., right-hand side of constraints, objective function coefficients) affect the optimal solution, offering a comprehensive understanding of the problem's.
* **Cross-Platform Compatibility:** Developed with CMake and Qt, ensuring potential deployment across Windows, macOS, and Linux.

## Technologies Used

* **C++:** Core application logic and UI development.
* **Qt 6:** Modern cross-platform framework for the graphical user interface.
* **Eigen:** C++ template library for linear algebra, used for efficient matrix operations in calculating feasible vertices.
* **Python 3:** Utilized for 3D plotting capabilities, invoked as a subprocess.
    * **Matplotlib:** Python plotting library for generating 3D graphs.
    * **NumPy:** Essential library for numerical operations in Python.
    * **SciPy:** Used for spatial algorithms, specifically `ConvexHull` for generating the faces of the 3D feasible region.
* **CMake:** Build system for managing the project compilation and packaging.
* **CPack:** Used for creating distributable packages (e.g., `.msi` for Windows).

## Getting Started

### Prerequisites

To build and run this project from source, you will need:

* **CMake** (version 3.16 or higher)
* **Qt 6** (with MSVC or MinGW compiler on Windows, depending on your setup)
* **Python 3.x**
* **Python dependencies:** `numpy`, `matplotlib`, `scipy`. It is recommended to install these in a virtual environment.

    ```bash
    python -m venv _python_env
    # On Windows:
    _python_env\Scripts\activate 
    # On Linux/macOS:
    # source _python_env/bin/activate 
    pip install numpy matplotlib scipy
    ```

### Building from Source

1.  **Clone the repository:**

    ```bash
    git clone [https://github.com/Claudia-Rocha-H/SimplexSolverQt.git](https://github.com/Claudia-Rocha-H/SimplexSolverQt.git)
    cd SimplexSolverQt
    ```

2.  **Create a build directory:**

    ```bash
    mkdir build
    cd build
    ```

3.  **Configure CMake:**
    Replace `C:/Qt/6.x.x/mingw_64` with your actual Qt 6 installation path and appropriate compiler (e.g., `msvc2019_64` for Visual Studio).

    ```bash
    cmake .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="C:/Qt/6.x.x/mingw_64"
    ```
    * *Note for Linux/macOS:* You might just need `cmake ..` or specify your Qt installation path differently.

4.  **Build the project:**

    ```bash
    cmake --build . --config Release
    ```

### Running the Application

After building, the executable `SimplexSolver.exe` (on Windows) will be located in the `build/Release` (or `build`) directory.

Before running, ensure that the `_python_env` directory (containing your Python installation and required libraries) and the `grafico3D.py` script are correctly placed relative to your executable, as per the `CMakeLists.txt` installation rules. For a development setup, you can often run directly from `build/Release`.

### Packaging the Application (for Distribution)

The project includes CPack configurations for creating installers.

1.  **Run `windeployqt` (Windows only):** This step copies necessary Qt DLLs next to your executable.
    Navigate to your build directory (e.g., `build/Release` or `build`).
    ```bash
    cd path/to/your/SimplexSolverQt/build/Release # or build
    windeployqt SimplexSolver.exe
    ```

2.  **Generate the package using CPack:**
    From your `build` directory:

    ```bash
    cpack -G WIX # For Windows MSI installer (requires WiX Toolset)
    # cpack -G DragNDrop # For macOS .dmg
    # cpack -G DEB # For Debian/Ubuntu .deb package
    # cpack -G TGZ # For a .tar.gz archive
    ```
    The generated installer will be in your `build` directory.

## Usage

1.  Launch the `SimplexSolver` application.
2.  Input the objective function coefficients and the number of variables/restrictions.
3.  Define each constraint by entering coefficients and the right-hand side value.
4.  Click "Calculate" to solve the Linear Programming problem.
5.  View the solution tableau, optimal values, and sensitivity analysis.
6.  The program visualize the feasible region for 2 and 3 variables problem.

## Contributing

Contributions are welcome! I'm still working on this program and it has some bug with the deployments, f you have suggestions for improvements or bug fixes, please open an issue or submit a pull request.

## License

This project is licensed under the MIT License.

## Contact

Claudia Rocha  https://github.com/Claudia-Rocha-H/ email: rochaclaudia177@gmail.com
