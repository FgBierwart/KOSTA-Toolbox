# STAKOOP (Stability Koopman)-Toolbox
The toolbox presented is based on the paper [] and estimates the basin of attraction of general vector field (polynomial and non-polynomial) using the Koopman operator framework. Hereafter, we detailed the main functions needed for the construction of the Lyapunov functions and its validation.  

# Construction of the Lyapunov function
The construction of the Lyapunov function is based on functions : main.m and Eigenfunction.m. In function main.m: 
- The $n$ dimension vector field is initialized as well as the domain of interest.


  Example:  
```ruby
  K = 0.2; F = [K*sin(x(1)-x(2))-sin(x(1)) K*sin(x(2)-x(1))-sin(x(2))];
  w = 3.5;
```
  
- The parameter *approx.flag* has to be set as 0 for polynomial vector field and 1 for non-polynomial vector field. For the latter, a polynomial approximation $P(x) = [P_1(x),...,P_n(x)]$ (of fixed order $d$) is given using the Taylor series of a min-max approximation. Note that for the Taylor series, constants $c_i$ are given such that $|F_i(x)-P_i(x)|< c_i\lVertx\rVert^d$ for some odd $d$ and $i=1,...,n$. Choising the polynomial approximation is done by the parameter *choice*:
   
```ruby
  choice = 'minimax'; order_rem = 12; 
  choice = 'Taylor'; order_tayl = 5; c = [0.7;0.7]
```
- According to parameter *basis* and the structure *s*, the different paremeter of the basis functions are given. We precise here that up to know, only monomials and gaussian have been implemented.   

- *trunction* is a parameter which indicates which projection we are using. This is a summarize by the following flowchart where dashed arrow indicate polynomial approximation computed using either Taylor or a min-max approximation.

<img src="https://github.com/FgBierwart/STAK-Toolbox/assets/142835014/f6c583be-ada8-4391-a5ea-8c652e92d738" width="700" height="230">


- Finally, *Eigenfunction.m* computes the eigenfunctions of the approximated Koopman operator. More precisely, 
  *  *Vec* contains the coefficient of the approximate eigenfunctions for the basis considered,
  *  *indx* contains the column of the eigenvalues closest associated to eigenvalues closest to those of the jacobian,
  Those two inputs allows to compute a Lyapunov candidate and will be used for the validation step described hereafter.
 
# Validation of the Lyapunov function 
to do
