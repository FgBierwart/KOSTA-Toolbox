# STAKOOP (Stability Koopman)-Toolbox
The toolbox presented is based on the paper [] and estimates the basin of attraction of general vector field (polynomial and non-polynomial) using the Koopman operator framework. Hereafter, we detailed the main functions needed for the construction of the Lyapunov functions and its validation. 

# Construction of The Lyapunov function. 
The construction of the Lyapunov function is based on functions : main.m and Eigenfunction.m. In function main.m: 
- The vector field is initialized as well as the domain of interest.\
  Example:  
  ```ruby
  F = [x(2) -2*x(1)+(1/3)*x(1)^3-x(2)];
  w = 4;
  ```
- The parameter approx.flag has to be set as 0 for polynomial vector field and 1 for non-polynomial vector field. For the latter, a polynomial approximation is given using the Taylor series of a min-max approximation.
- choice = 'minimax';
- choice = 'minimax';
  ```ruby
choice = 'minimax';
choice = 'Taylor';
```

