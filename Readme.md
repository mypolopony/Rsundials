 ## Rsundials Version 1.6  
###### Selwyn-Lloyd McPherson

### Contents
1. Summary
2. IDA Example
3. CVODES Example

### 1. Summary
Rsundials implements the differential algebraic equation (DAE) and ordinary differential equation (ODE) solvers in the SUNDIALS suite (version 2.3.0 – http://www.llnl.gov/CASC/sundials/). The IDA module of SUNDIALS handles DAEs and the CVODES module solves ODE systems. Both modules utilize a dense linear solver and both require a user-defined (hard-coded) residual / right hand side function. Installation of the SUNDIALS libraries is not a prerequisite for this package.

### 2. IDA Example
*The Problem* 
This example, due to Robertson [1], is a model of a three-species chemical kinetics system written in DAE form. Differential equations are given for species y1 and y2 while an algebraic equation determines y3. The equations for the system concentrations yi(t) are:

<img src="https://render.githubusercontent.com/render/math?math=\dot{y}_{1}=0.4y_{1} + 10^4y_{2}y_{3}">
<img src="https://render.githubusercontent.com/render/math?math=\dot{y}_{2}=0.5y_{1} - 10^4y_{2}y_{3} - 30\cdot10^7y_{2}^2">
<img src="https://render.githubusercontent.com/render/math?math=0 = y_{1} + y_{2} + y_{3} - 1">

The initial values are taken as y1 = 1, y2 = 0, and y3 = 0. This example computes the three
concentration components on the interval from t = 0 through t = 4×10^10


*The Residual Function*
The first step is to create a compilable file (preferably in c/++) that defines the residual (right hand side) function as called by the solver on each time step. A template for such a function can be described as:

```c++
#include "include/nvector_serial.h" #include "include/sundials_dense.h"

int resrob(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *rdata)
{
	realtype *yval, *ypval, *rval;
	yval = NV_DATA_S(yy); ypval = NV_DATA_S(yp); rval = NV_DATA_S(rr);

	/* Set values of rval[0]...rval[n] here */

	return(0); 
}
```

The function takes these arguments:

| Command | Description |
| --- | --- |
| __tres__ | the current value of the independent variable |
| __yy__ | the current value of the dependent variable vector, y(t) yp the current value of y'(t) |
| __rr__ | the output residual vector F(t, y, y') |
| __rdata__ | a pointer to user data |

The two `include`s are used to bring in the data types used by the SUNDIALS solver, namely `N_Vector`, a simple vector, and `realtype`, a primitive data type. By default, this is a double-precision floating-point numeric data type (double C-type). The header files are found in the usrfcns directory.

In order to write or read from the vectors, __NV_DATA_S__ must be used on a realtype pointer. This sets the pointer to the first value of the vector. In general, nothing must be done with __tres__ or __rdata__.

IDA works by attempting to minimize the residual of the equations presented. As such, it is important to formulate the system properly by bringing everything over to one side. The output residual vector __rr__ can be written to by writing to the elements of rval. For the above system, this can be done in the following way:

```c++
rval[0] = -0.04*yval[0] + 1.0e4*yval[1]*yval[2] - ypval[0]; 
rval[1] = 0.04*yval[0] - 1.0e4*yval[1]*yval[2] - 3.0e7*yval[1]*yval[1] - ypval[1]
rval[2] = yval[0] + yval[1] + yval[2] - 1;
```

If the complete function __resrob__ is placed in the file _idafcns.c_, it can be compiled into a
shared library useable by R using:

`
     R CMD SHLIB ...path.../idafcns.c
`	

The header files _nvector_serial.h_ and _ida_dense.h_ should be present in the same directory. This will compile the file and allow the functions within to be usable in R.

*Solving the Problem in R*
Once R is started, the defined functions must be loaded before running the IDA solver. To load dynamic libraries into R, use the command

`
dyn.load("...path.../idafcns.so") or dyn.load("...path.../idafcns.dll")
`

(depending on your system type). To execute the IDA solver in R, the following command might be used:

```R
     yvals <- ida(c(1,0,0), c(-0.04,0.04,0),
          gseq(.4,4e10,10), "idafcns", "resrob", jacfunc =
          "jacrob", rtol = 1E-4, atol = c(1E-8,1E-14,1E-6),
          verbose = TRUE)
```

The following output (or something similar) should be obtained:

```
SUNDIALS IDADENSE Linear Solver Number of Equations: 3 Integration Limits: 0.4 to 4e+10 Solver Memory Allocated
Max Number of Steps: 500
Max step size: 0
Relative Tolerance: 0.0001
Absolute Tolerances: 1e-08 1e-14 1e-06 Jacobian Function Initialized
IDADense Linear Solver Initialized
```

| t | y0 | y1 | y2 | nst | k | h |
| --- | --- | --- | --- | --- | --- | --- |
| 4.0000e+00 | 9.1172e-01 | 2.3133e-05 | 8.8253e-02 | 96 | 4 | 4.1520e-01|
| 4.0000e+01 | 7.1686e-01 | 9.2258e-06 | 2.8313e-01 | 132| 4 | 1.9421e+00|
| 4.0000e+02 | 4.5066e-01 | 3.2247e-06 | 5.4933e-01 | 177| 4 | 3.2576e+01|
| 4.0000e+03 | 1.8321e-01 | 8.9427e-07 | 8.1679e-01 | 220| 4 | 2.0527e+02|
| 4.0000e+04 | 3.8985e-02 | 1.6218e-07 | 9.6102e-01 | 266| 3 | 1.2965e+03|
| 4.0000e+05 | 4.9391e-03 | 1.9853e-08 | 9.9506e-01 | 307| 5 | 1.6802e+04|
| 4.0000e+06 | 5.1672e-04 | 2.0679e-09 | 9.9948e-01 | 366| 3 | 1.5150e+05|
| 4.0000e+07 | 5.2026e-05 | 2.0811e-10 | 9.9995e-01 | 416| 4 | 2.0066e+06|
| 4.0000e+08 | 5.2132e-06 | 2.0853e-11 | 9.9999e-01 | 452| 4 | 2.7282e+07|
| 4.0000e+09 | 5.2178e-07 | 2.0871e-12 | 1.0000e+00 | 482| 4 | 2.6574e+08|
| 4.0000e+10 | 5.1049e-08 | 2.0420e-13 | 1.0000e+00 | 503| 2 | 6.8880e+09|

```
Final Run Statistics:
Number of steps = 503
Number of residual evaluations = 703
Number of Jacobian evaluations = 78
Number of nonlinear iterations = 703
Number of error test failures = 16
Number of nonlinear conv. failures = 0
Number of root fn. evaluations = 0
```

### 3. CVODES Example

*The Problem*
This example from a presentation by Borrelli and Coleman [2] is a three-compartment model for lead in the human body. Lead is input to the system at a constant rate _L_. Three state variables, _x1_, _x2_, and _x3_ describe the concentration of lead in the blood, tissue, and bones respectively. There exist transfer rates between the compartments as well as to the external environment via urine from the blood and via hair, nails, and sweat from the tissues.

![Compartment Model](img.compartment_model.png "Compartment Model")

For _i = 1, 2, 3_, we let <img src="https://render.githubusercontent.com/render/math?math=x_{i}t"> be the amount of lead in compartment i at time t and we assume that the rate of transfer from compartment i to j is proportional to xi(t) with a proportionality constant of aji. The units for amounts of lead are micrograms and the time t is measured in days.
The vector differential equation of this problem is in the form
where A is the matrix
x! = A x + b
and b = (L,0,0)T .
a21 = 0.011 a31 = 0.0039 a01 = 0.021
Units: days-1
"!(a01 +a21 +a31) a12 a13 % $ a21 !(a02+a12) 0'
$a 0!a' # 31 13&
In a paper published by Rabinowitz and colleagues [3], measurements of the concentration of lead in these compartments in a male subject living in Los Angeles allowed for the calculation of the rates of transfer. Relatively speaking, lead is somewhat slow to enter the bones and very slow to leave them.
Lead Transfer Coefficients (Rabinowitz, et al.)
a12 = 0.012 a13 = 0.000035 a02 = 0.016
from blood to tissue and back from blood to bone and back excretion from blood and tissue
The study also showed that the average rate of ingestion of lead (L) in Los Angeles over the period studied was 49.3 micrograms per day.
The Right Hand Side Function
The first step is to create a compilable file (preferably in c) that defines the right hand side function as called by the solver on each time step. A template for such a function can be described as:
#include "include/nvector_serial.h" #include "include/sundials_dense.h"
int rhs(realtype t, N_Vector y, N_Vector ydot, void *f_data) {
realtype y1, y2, y3;
y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
/* Change values of ydot here using Ith(ydot,i) */
return(0); }
The function takes the arguments:
t the current value of the independent variable
y the current value of the dependent variable vector, y(t) ydot the current value of y'(t)
fdata a pointer to user data
The two includes are used to obtain the data types used by the SUNDIALS solver, namely N_Vector, a simple vector, and realtype, a primitive data type. By default, this is a double-precision floating-point numeric data type (double C-type). The header files are found in the usrfcns directory.
The realtypes y1, y2 and y3 are first set to the corresponding value in the y vector using the macro Ith with indices from 1 to n (not from 0 to n-1). These can then be modified to reflect the proper value. For the above problem, the values of ydot can be updated using the following code:
double L = 49.3;
double a21 = 0.011; double a12 = 0.012; double a31 = 0.0039; double a13 = 0.000035; double a01 = 0.021; double a02 = 0.016;
Ith(ydot, 1) = -(a01 + a21 + a31)*y1 + a12*y2 + a13*y3 + L; Ith(ydot, 2) = a21*y1 - (a02 + a12)*y2;
Ith(ydot, 3) = a31*y1 - a13*y3;
If the complete function rhs is placed in the file cvodesfcns.c, it can be compiled into a shared library useable by R by typing:
     R CMD SHLIB ...path.../cvodesfcns.c
in the command line as long as the header files nvector_serial.h and cvodes_dense.h are present in the same directory. This will compile the file and allow the functions within to be usable in R.
The user may supply a method for computing the Jacobian; this can be placed in the same file as the RHS method.
Solving the Problem in R
Using the model defined above, let's take a look at an individual who moves to Los Angles with no lead in her body. We can use Rsundials to determine the levels in her system after a set number of days, say 400.
Once R is started, the defined functions must be loaded before running the CVODES solver. To load dynamic libraries into R, use the command
dyn.load("...path.../cvodesfcns.so") or dyn.load("...path.../cvodesfcns.dll")
depending on your system type. To execute the CVODES solver in R, the following command might be used:
     yvals <-
          cvodes(c(0.0,0.0,0.0),seq(0,400,20),"cvodesfcns",
At t = 2.0000e+01 At t = 4.0000e+01 At t = 6.0000e+01 At t = 8.0000e+01 At t = 1.0000e+02 At t = 1.2000e+02 At t = 1.4000e+02 At t = 1.6000e+02 At t = 1.8000e+02 At t = 2.0000e+02 At t = 2.2000e+02 At t = 2.4000e+02 At t = 2.6000e+02 At t = 2.8000e+02 At t = 3.0000e+02 At t = 3.2000e+02 At t = 3.4000e+02 At t = 3.6000e+02 At t = 3.8000e+02 At t = 4.0000e+02
Final Run Statistics:
7.088272e+02 1.073435e+03 1.272270e+03 1.386870e+03 1.456042e+03 1.499400e+03 1.527143e+03 1.545141e+03 1.557178e+03 1.565419e+03 1.570973e+03 1.574665e+03 1.577182e+03 1.578955e+03 1.580190e+03 1.581042e+03 1.581648e+03 1.582106e+03 1.582467e+03 1.582768e+03
7.198569e+01 1.972720e+02 3.134134e+02 4.047548e+02 4.716999e+02 5.189813e+02 5.518129e+02 5.743737e+02 5.896541e+02 5.999263e+02 6.069031e+02 6.116752e+02 6.149088e+02 6.170854e+02 6.185755e+02 6.196108e+02 6.203337e+02 6.208375e+02 6.211930e+02 6.214457e+02
3.079803e+01
1.017541e+02
1.938759e+02
2.978071e+02
4.086610e+02
5.237255e+02
6.414478e+02
7.608681e+02
8.813493e+02
1.002498e+03
1.124084e+03
1.245947e+03
1.367965e+03
1.490058e+03
1.612180e+03
1.734300e+03
1.856392e+03
1.978440e+03
2.100434e+03
2.222365e+03
"rhs", rtol = 1E-4,atol = c(1E-8,1E-14,1E-
6),verbose = TRUE)
The following output should be obtained:
SUNDIALS CVODES Linear Solver
Number of Equations: 3
Integration Limits: 0 to 400 InitialValues: y0=0 y1=0 y2=0 Solver Memory Allocated
Relative Tolerance: 0.0001
Absolute Tolerances: 1e-08
CVDENSE Solver Initiated
Max number of steps: 500
Max step size: 0
Requesting data for all time points.
Number of steps
Number of RHS evaluations = 89 Number of linear solver setups = 31 Number of nonlinear iterations = 85 Number of error test failures = 1 Number of nonlinear conv. failures = 0
Number of root fn. evaluations
= 0 [,4]
0.00000 30.79803 101.75408 193.87588 297.80714 408.66104 523.72549 641.44779 760.86810 881.34927 1002.49768 1124.08405
      [,1]
 [1,]    0
 [2,]   20
 [3,]   40
 [4,]   60
 [5,]   80
 [6,]  100
 [7,]  120
 [8,]  140
 [9,]  160
[10,]  180
[11,]  200
[12,]  220
[,2] [,3] 0.0000 0.00000 708.8272 71.98569 1073.4351 197.27204 1272.2696 313.41340 1386.8700 404.75483 1456.0422 471.69985 1499.4004 518.98131 1527.1425 551.81292 1545.1407 574.37374 1557.1780 589.65408 1565.4186 599.92632 1570.9726 606.90308
1e-14 1e-06
= 73
[13,] 240 1574.6646 611.67520 1245.94737 [14,] 260 1577.1821 614.90883 1367.96547 [15,] 280 1578.9551 617.08536 1490.05834 [16,] 300 1580.1895 618.57552 1612.18045 [17,] 320 1581.0423 619.61078 1734.29997 [18,] 340 1581.6483 620.33370 1856.39243 [19,] 360 1582.1063 620.83745 1978.44042 [20,] 380 1582.4672 621.19300 2100.43377 [21,] 400 1582.7678 621.44568 2222.36541
 1 H. H. Robertson. The solution of a set of reaction rate equations. In J. Walsh, editor, Numerical analysis: an introduction, pages 178–182. Academ. Press, 1966.
2 Differential Equations: A Modeling Approach, by R. Borrelli and C. Coleman, Prentice-Hall, 1987.
3 Rabinowitz MB, Wetherill GW, Kopple JD. Lead metabolism in the normal human: stable isotope studies. Science. 1973 Nov 16;182(113):725–727.