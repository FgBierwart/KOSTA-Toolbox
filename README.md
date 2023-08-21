# STAKOOP (Stability Koopman)-Toolbox
The toolbox presented is based on the paper [] and estimates the basin of attraction of general vector field (polynomial and non-polynomial) using the Koopman operator framework. Hereafter, we detailed the main functions needed for the construction of the Lyapunov functions and its validation. 

# Construction of The Lyapunov function. 
The construction of the Lyapunov function is based on functions : main.m and Eigenfunction.m. In function main.m: 
- The $n$ dimension vector field is initialized as well as the domain of interest.\
  Example:  
  ```ruby
  K = 0.2; F = [K*sin(x(1)-x(2))-sin(x(1)) K*sin(x(2)-x(1))-sin(x(2))];
  w = 3.5;
  ```
- The parameter approx.flag has to be set as 0 for polynomial vector field and 1 for non-polynomial vector field. For the latter, a polynomial approximation $P(x) = \left[P_1(x),...,P_n(x)\right]$ (of fixed order $d$) is given using the Taylor series of a min-max approximation. Note that for the Taylor series, constants $c_i$ are given such that $|F_i(x)-P_i(x)|<c_i\|x\|^{d}$ for some d. 
  ```ruby
  choice = 'minimax'; order_rem = 12; 
  choice = 'Taylor'; order_tayl = 5; c = [0.7;0.7]
  ```
  Note that 

