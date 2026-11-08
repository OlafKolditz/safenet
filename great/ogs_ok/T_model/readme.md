readme.md great

06.11.2026

(1) meshing
We want to build an OGS model for a T process. First step is to build a circular finite element mesh with radius 0.1 m. We also need an inner circle at 0.066 m for calculation of volumetric strain later on. Please build the triangle mesh with a central point at (0,0,0).

(2) T process
Now we build the OGS model for the T process. The initial temperature is 293 K and the temperature at the outer boundary is 323 K. The material properties are for granite.

(3) verification
Please build a Jupyter Notebook for verification of the numerical solution including the analytical solution.

(4) T boundary conditions
Please modify the temperature at the outer boundary as follows: linear increase of temperature from 293 K to 363 K from t = 0 to 325 seconds and than maintain the temperature at 363 K until t = 4000 seconds.

(5) recover previous Notebook
Please recover the previous notebook and create two separate notebooks.

(6) TH process
Please extend the T process to a TH process with water injection in the central point from 4000 to 4050 seconds. Provide a new OGS project file for the TH process.

(7) great project
Please move our today's session to the great project, the project folder is 94-CLD/great.