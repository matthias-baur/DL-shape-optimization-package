----------------------------
|          README          |
----------------------------
Matlab scripts for Numerical Optimization of consecutive sums of eigenvalues 
of the Dirichlet Laplacian

by Matthias Baur, July 2021


This package contains MATLAB files for the numerical optimization of 2D shapes with respect to consecutive sums of eigenvalues of the Dirichlet Laplacian (MATLAB version R2020a). The shape optimization algorithms and methods implemented are all described in

Pedro R. S. Antunes and Pedro Freitas. "Numerical Optimization of Low Eigenvalues of the Dirichlet and Neumann Laplacians". In: Journal of Optimization Theory and Applications 154 (2012), pp. 235-257. doi: 10.1007/s10957-011-9983-3.

The package includes three main optimization scripts:

1. genetic_optimization.m 
    -> the genetic algorithm that is used to find a good initial domain for the gradient descent method
2. gradient_descent_optimization.m
    -> main part of the optimization scheme
3. random_optimization.m 
    -> fine tuning of domain

Most other scripts are functions that are used in those optimization scripts. A detailed descriptions of what the files do is given below.


---------------------------
|          Usage          |
---------------------------

---    Note   ---
Domains are always described by a pair of vectors, a and b. The vector a contains the Fourier coefficients to cosine terms and b contains Fourier coefficients to sine terms of the radial function r(t). The splitting of the full coefficient vector c=(a,b) into two coefficient vectors of half the length is purely a matter of taste.


-------------------------------
---    1. Set parameters    ---
-------------------------------
First thing to do before using any of the optimization scripts is to run the set_parameters.m file. Here, the sum that will be minimized is specified via the parameters n and sum_length. Minimization will then be performed on the objective function

J = lambda_{n-sum_length+1} + ... + lambda_n

Setting sum_length = n corresponds to minimizing

J = lambda_1 + ... + lambda_n

i.e. the sum of all eigenvalues up to n. Setting sum_length = 1 corresponds to minimizing

J =  lambda_n

i.e. just the n-th eigenvalue.

The next important parameter is len_ab which sets the number of Fourier coefficients used to describe the domains (number of sine coefficients = number of cosine coefficients = len_ab). 

The parameter interval_end has to be chosen carefully. interval_end specifies an upper cutoff for eigenfrequencies, i.e. when computing eigenvalue approximations, no eigenfrequencies above interval_end can be found. In this regard, interval_end should not be chosen too low to avoid the program getting stuck because it cannot find the n-th eigenfrequency of a domain. interval_end needs to be chosen >sqrt(lambda_n) to be able to find all n eigenvalues up to lambda_n. At the same time, choosing interval_end too big should not be done either, because this leads to unnecessary long computations when unneeded eigenfrequencies are computed. Another problem that can arise is that some eigenfrequencies can be missed when two or more of them lie close together and interval_end>sqrt(lambda_{n+1}). For a fixed domain with simple eigenvalues, interval_end should be chosen such that

sqrt(lambda_n) < interval_end < sqrt(lambda_{n+1})

to be sure to find all eigenvalues up to lambda_n.

For gradient_descent_optimization.m and random_optimization.m it is recommended to set interval_end = sqrt(lambda_n of initial domain) + some small constant (0.3-0.5). For genetic_optimization.m set interval_end = sqrt(lambda_n of circle) + some large constant (3 to 5).

Another parameter to emphasize here is p_beta which sets the displacement length of the source points of MFS eigenfunctions with respect to the collocation points. It is recommended to use p_beta = 0.03; for gradient_descent_optimization.m and random_optimization.m and to use p_beta = 0.01; for genetic_optimization.m. While p_beta = 0.03 generally leads to more accurate eigenvalue approximations for relatively close to circular domains, the first eigenvalue for some wild shapes in the genetic algorithm can be garbage when using p_beta = 0.03. Here, p_beta = 0.01 seems to prevent non-sensical eigenvalues, but decreases the overall accuracy of the eigenvalues slightly. For the phase of exploration, this reduced accuracy is still sufficient.

Other parameters listed in set_parameters.m can be left to their default values.

-------------------------------------------
---    2. Run genetic_optimization.m    ---
-------------------------------------------
With all parameters set, running this script will yield a final population of close to optimal domains (w.r.t the objective function that was set). The Fourier coefficients are row-wise stored in A_pop and B_pop, where A_pop stores the cosine coefficients and B_pop stores the sine coefficients for each domain. The coefficients are sorted such that the uppermost row of A_pop and B_pop describe the best domain in the population. The associated value of the objective function to each domain is stored in sum_pop.

----------------------------------------------------
---    3. Run gradient_descent_optimization.m    ---
----------------------------------------------------
After finishing genetic optimization, one can proceed with gradient descent. The initial coefficients a, b are set in the beginning. The default setting is

a = A_pop(1,:);  
b = B_pop(1,:);

which just takes over the coefficients of the best domain found in the genetic algorithm. Before running this script, run set_parameters.m with p_beta = 0.03; to ensure more accurate eigenvalue approximations.

Coefficient vectors of intermediate steps and their associated value for the objective function are stored in a_iter, b_iter and sum_n_iter.

------------------------------------------
---    4. Run random_optimization.m    ---
------------------------------------------
For fine tuning. After gradient descent has finished its loop, it is possible to fine tune the coefficients by small random deviations. The default setting takes over the final coefficients of the gradient descent loop and randomizes these. This can be run as long as one likes, but the gains become very small after a few hundred iterations.


---------------------------------------
|          Individual Files           |
---------------------------------------

- define_domain.m 
Given Fourier coefficients a, b, this function returns all important variables that describe a domain. This includes the functions $r(t)$ and $r'(t)$, collocation points on the boundary of the set, outer normals at the collocation points, source points for the Method of Fundamental solutions and the volume and center of mass of the domain.

- define_r.m 
Defines function r(t) and its derivative r'(t) for the domain boundary with given Fourier coefficients.

- direct_problem.m 
Used to solve the eigenvalue problem for given Fourier coefficients. Gives back the value of the sum of eigenvalue in question, multiplicity of the largest eigenvalue (stored in mult(1)), all eigenfrequencies up to index n (stored in mult(2:n+1)) and one coefficient vector $\alpha*$ that defines an eigenfunction in the MFS scheme for each eigenvalue found (stored in V, sorted columnwise, largest to smallest eigenvalue).

- dn_eigf.m 
Computes the normal derivative of an approx. eigenfunction in points x,y on the boundary, needs coefficients alpha, eigenfrequency approximation kappa and source points p_source. This method uses the exact derivative of the first Hankel function.

- dn_eigf_FD.m 
Computes the normal derivative of an approx. eigenfunction in points x, y on the boundary, needs an eigenfunction eigf, normal vectors normals and a small step width dx. This method uses finite difference with step width dx to compute normal derivatives. Currently unused, dn_eigf is preferred.

- eigf.m 
Used to evaluate an approximate eigenfunction on a meshgrid X, Y. Needs eigenfrequency approximation kappa, coefficients alpha and source points p_source.

- eigf_norm.m 
Computes approximate L2-norm of an approximate eigenfunction by evaluating it on a meshgrid X, Y with meshwidth dh and doing a basic quadrature. Needs eigenfrequency approximation kappa, coefficients alpha and source points p_source and meshwidth dh.

- eigval_gradient.m 
Computes the gradient of eigenvalue lambda = k_crit^2 wrt Fourier coefficients for given Fourier coefficients a,b. Needs parameter struct param, eigenfrequency approximation k_crit and coefficients alpha that define the approximate eigenfunction.

- genetic_optimization.m 
Performs genetic optimization. At the beginning of the script one can set additional parameters: population size (N_pop), number of fittest domains selected for reproduction (N_best) and number of iterations (N_iter). The Fourier coefficients are stored in A_pop and B_pop where A_pop row-wise stores the cosine coefficients and B_pop stores the sine coefficients for each domain. The coefficients are sorted such that the upmost row of A_pop and B_pop describe the best domain in the population. The associated value of the objective function to each domain is stored in sum_pop

- gradient_descent_optimization.m
Performs gradient descent. Set initial coefficients a, b first. The default setting is takes over the coefficients of the best domain found in the genetic algorithm. Automatically adjusts objective function if the next lower eigenvalues becomes closer than mult_tol (set in set_parameters.m ) to the lowest eigenvalue in the sum in question. If sum starts with lambda_1, this is irrelevant. Coefficient vectors of intermediate steps and their associated value for the objective function are stored in a_iter, b_iter and sum_n_iter. The final coefficient vectors are found in the last row of a_iter resp. b_iter.

- gr_max_search.m 
Searches for a local maximizer of a real function f:R->R in the interval interval by the golden ratio method. The accuracy is set by tol.

- gr_min_search.m 
Searches for a local minimizer of a real function f:R->R in the interval interval by the golden ratio method. The accuracy is set by tol.

- gr_min_line_search.m 
Searches for the step width beta_opt in beta_interval that minimizes the objective function along the half line that is defined by the negative gradient d. The minimization is done with the golden ratio method. 

- parameter_struct.m 
Puts all important parameters set in set_parameters.m in a struct.

- random_optimization.m 
Fine tunes Fourier coefficients by iteratively applying random deviations and keeping coefficients when they turn out to be better than before. The default setting takes over the final coefficients of the gradient descent loop and randomizes these. The maximal size of relative deviations is set width dev_max. dev_max = 0.025 means that each randomized coefficient does not deviate more than 2.5% from the resp. coefficient from the current best domain.

- rescale_domain.m 
Rescales Fourier coefficients a, b such that the area of the associated domain is normalized to one.

- search_eigval_interval.m 
Searches for local minima of f = log(|det(A(kappa))|) in interval. Nk gives the number of equidistant grid points on which f is evaluated first. Then the function searches for minima between the local maxima on the coarse grid and interval endpoints (using golden ratio search). The higher Nk, the unlikelier it becomes that local minima very close to each other will be missed. There is no guarantee that the function finds all local minima in the interval when two minima are very close and the choice of Nk is not large enough to clearly resolve those two minima.

- set_parameters.m (see Usage above for more details)
Used to set all basic parameters needed for optimization procedures. Most importantly here one defines the objective function J(\lambda_1,...,\lambda_n) = lambda_{n-sum_length+1} + ... + lambda_n that is used and the number of Fourier coefficients for domain boundary parametrization. Sets also parameters for gradient descent loop and parameters that define the collocation and source points for the Method of Fundamental solutions.
