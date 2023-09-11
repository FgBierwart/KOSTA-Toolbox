# STAKOOP (Stability with Koopman)  -  Toolbox
The STAKOOP toolbox is built upon the results presented in [1]. Its aims at estimating the region of attraction of equilibria for general vector fields (polynomial and non-polynomial) using the Koopman operator framework (see [2] for further information). Hereafter, we detail the main functions needed to construct candidate Lyapunov functions, validate them, and estimate the region of attraction.  

### Initialization of the dynamical system. 
To do 

### Construction of the Lyapunov function
The construction of the Lyapunov function is based on the functions *main.m* and *Eigenfunction.m*. In *main.m*:  
- The $n$-dimensional vector field is initialized as well as the domain of interest $\mathbb{X} = [-w,w]^n$. 


  Example ( $n = 2$ ):  
```ruby
  K = 0.2; F = [K*sin(x(1)-x(2))-sin(x(1)) K*sin(x(2)-x(1))-sin(x(2))];
  w = 3.5;
```
  
- The parameter *approx.flag* has to be set to 0 for polynomial vector fields and 1 for non-polynomial vector fields. In the latter case, a polynomial approximation $P(x) = [ P_1(x),...,P_n(x) ]$ (of fixed order $d$) is given using Taylor series or min-max approximation. Note that for Taylor series, constants $c_i$ are given such that $|F_i(x)-P_i(x)|< c_i\lVert x\rVert^{d+1}$ for some odd $d$ and $i=1,...,n$. The polynomial approximation is chosen through the parameter *choice*:
   
```ruby
  choice = 'minimax'; order_rem = 12; 
  choice = 'Taylor'; order_tayl = 5; c = [0.7;0.7]
```
- The different paremeters of the basis functions are provided via the parameter *basis* and the structure *s*. In the current version, only monomials and Gaussian radial basis functions are implemented.

- The projection used to approximate the Koopman operator is given with the parameter *trunction*.

The following flowchart shows 4 different ways to construct a Lyapunov candidate, depending on the chosen basis functions and projection operator. The dashed arrows depict the polynomial approximation. 

<img src="https://github.com/FgBierwart/STAK-Toolbox/assets/142835014/f6c583be-ada8-4391-a5ea-8c652e92d738" width="700" height="230">

&nbsp;

- Finally, the Lyapunov candidate is computed with the function *Eigenfunction.m*. The function gives as an output the matrix of approximated eigenvectors and indicates those associated to the eigenvalues closest to the eigenvalues of the Jacobian matrix. See [] for more details and documentation of the related function.     
 
### Validation of the Lyapunov function 

In this section, we present the two main functions used for validation. They are based on (i) Sum-Of-Squares (SOS) programming and (ii) a “worst case” approach combined with an adaptive grid.

- SOS-based validation

  The region of attraction is estimated with the function *Lyap_certificate*. This function returns the value of $\gamma_1$ and $\gamma_2$ such that the largest set $`\{x\in\mathbb{X}~|~\gamma_1\leq V(x) \leq \gamma_2\}`$ is in the validity region $\mathcal{S} =$ $`\{x\in\mathbb{X}~|~\dot{V} < 0\}`$.

**Remark**. SOS-based validation techniques requires the SOSTOOLS toolbox from [3] 
  
- Grid validation

  An estimation of the region of attraction using an adaptive grid is performed in two steps. The first step provides an approximation of the validity region $\mathcal{S}$. The associated function *AGM.m* returns a matrix containing the coordinates and size of grid cells lying within the validity set. Finally, the region of attraction is computed through a bisection method with *bisection_grid.m*.

  **Please note that, in the current version, *AGM.m* and *bisection_grid.m* are implemented only for 2D polynomial systems.**

# References 
[1] Bierwart and Mauroy. “A numerical Koopman-based framework to estimate regions of attraction for general vector fields.” Submitted to Communications in Nonlinear Science and Numerical Simulation

[2] A. Mauroy, Y. Susuki, and I. Mezi´c, The Koopman operator in systems and control, Springer, 2020.

[3] S. Prajna, A. Papachristodoulou, and P. A. Parrilo, Introducing sostools: A general purpose sum of squares programming solver, in Proceedings of the 41st IEEE Conference on Decision and Control,
2002., vol. 1, IEEE, 2002, pp. 741–746

