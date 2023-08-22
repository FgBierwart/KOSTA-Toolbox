# STAKOOP (Stability Koopman)-Toolbox
The toolbox presented here is based on the pre-print [1] and estimates the basin of attraction of general vector fields (polynomial and non-polynomial) using a Koopman operator framework (see [2] for further information). Hereafter, we detailed the main functions needed for the construction of the Lyapunov functions and its validation.  

# Construction of the Lyapunov function
The construction of the Lyapunov function is based on functions : $\texttt{main.m}$ and  $\texttt{Eigenfunction.m}$. In function main.m: 
- The $n$ dimension vector field is initialized as well as the domain of interest $\mathbb{X} = [-w,w]^n$. 


  Example ( $n = 2$ ):  
```ruby
  K = 0.2; F = [K*sin(x(1)-x(2))-sin(x(1)) K*sin(x(2)-x(1))-sin(x(2))];
  w = 3.5;
```
  
- The parameter $\color{royalblue}\texttt{approx.flag}$ has to be set as 0 for polynomial vector field and 1 for non-polynomial vector field. For the latter, a polynomial approximation $P(x) = [ P_1(x),...,P_n(x) ]$ ( of fixed order $d$ ) is given using the Taylor series of a min-max approximation. Note that for the Taylor series, constants $c_i$ are given such that $|F_i(x)-P_i(x)|< c_i\lVert x\rVert^{d+1}$ for some odd $d$ and $i=1,...,n$. Choising the polynomial approximation is done by the parameter $\color{royalblue}\textit{choice}$:
   
```ruby
  choice = 'minimax'; order_rem = 12; 
  choice = 'Taylor'; order_tayl = 5; c = [0.7;0.7]
```
- According to parameter *basis* and the structure *s*, the different paremeter of the basis functions are given. We precise here that up to know, only monomials and gaussian have been implemented.   

- *trunction* is a parameter which indicates which projection we are using. This is a summarize by the following flowchart where dashed arrow indicate polynomial approximation computed using either Taylor or a min-max approximation. So, there is 4 different way to construct a Lyapunov candidate. 

<img src="https://github.com/FgBierwart/STAK-Toolbox/assets/142835014/f6c583be-ada8-4391-a5ea-8c652e92d738" width="700" height="230">

&nbsp;

- Finally, according to the function *Eigenfunction.m*, this lyapunov candidate is estimated. More precisely, This function gives as an output the matrix of approximated eigenvectors and points out those associated to the eigenvalues closest the one of the Jacobian matrix. See [] for more details and documentations of the related function.     
 
# Validation of the Lyapunov function 

In this section, we present the two main functions developped for two validation techniques based on (i) Sum-Of-Square (SOS) programming and (ii) a “worst case” approach combined with an
adaptive grid

- SOS validation

  An estimation of the basin of attraction is given by the function *Lyap_certificate*. This function return the value of $\gamma_1$ and $\gamma_2$ such that the $`\{x\in\mathbb{X}~|~\gamma_1\leq V(x) \leq \gamma_2\}`$ is in the validity region delimited by the set $\mathcal{S} =$ $`\{x\in\mathbb{X}~|~\dot{V} < 0\}`$. As illustrated on the flowchart, when a Lyapunov candidate is computed using other basis functions than monomials, a twlever order (default) min-max polynomial approximation of it is constructed for SOS use.      
  
- Grid validation

  An estimation of the basin of attraction using the adaptive grid is done in two step. The first one is for the construction of an approximation of the validity region $\mathcal{S}$. This is done by constructing an adaptive grid where each cells are well in the validity region. The associated function is given by *AGM.m* and returns a matrix containing coordinates and associated stepsize of valid cells. Fianlly, the ROA is computed according to *bisection_grid.m* among the valid cells. 

  **Please note that we only tested *AGM.m* and implemented *bisection_grid.m* in 2D only for polynomial systems. We leave their application of non-polynomial vector fields for future reasearch.**

# References 
[1] 

[2] A. Mauroy, Y. Susuki, and I. Mezi ́c, The Koopman operator in systems and control, Springer, 2020.

