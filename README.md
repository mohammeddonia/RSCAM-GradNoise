# RSCAM-GradNoise
# Intro about the project by course organizers
When Langevin schemes are used in sampling the Bayesian posterior distribution defined
by a large (huge) data set, the cost of each step can be very high. This has led to a lot of
work to try to reduce the cost by simplifying the gradient. The most popular method of this
type is stochastic gradient Langevin dynamics. The idea behind these methods is that the
error in computing the gradient using a subsampling of the data set can be viewed as a
white noise perturbation (like the one in Langevin). The complication is that this gradient
noise may not actually be white (iid from step to step) in practice, and, even if nearly white,
its properties (e.g. variance) are not known. Several papers have tried to estimate the
gradient noise adaptively.

[4.1] Bayesian Learning via Stochastic Gradient Langevin Dynamics, Max Welling and Yee
Whye Teh, Proceedings of the 28th International Conference on Machine Learning, 2011.
https://www.ics.uci.edu/~welling/publications/papers/stoclangevin_v6.pdf
[4.2] Exploration of the (Non-)Asymptotic Bias and Variance of Stochastic Gradient
Langevin Dynamics, Sebastian J. Vollmer, Konstantinos C. Zygalakis, and Yee Whye Teh,
J. Mach. Learn. Res., 17, 1−48, 2016
http://jmlr.org/papers/volume17/15-494/15-494.pdf
[4.3] Stochastic Gradient Hamiltonian Monte Carlo, Tianqi Chen, Emily B. Fox and Carlos
Guestrin, Proceedings of the 31st International Conference on Machine Learning, Beijing,
China, 2014.
https://arxiv.org/pdf/1402.4102.pdf
