# Parametrization of the density functions

------------------------------------------------------------------------

## TL; DR

This document states the parametrization used for the positive density
functions used in the package.

## Log-normal

Let $Y$ be a random variable following a Log-normal distribution with
parameters $\eta$ and $\tau$. We denote $Y \sim LN(\eta,\tau)$. The
density function associated with such a variable is
$$g(y) = (2\pi\tau)^{-1/2}y^{-1}\exp\left\{ -\frac{\left( \log y - \eta \right)^{2}}{2\tau} \right\}.$$
The expected value (theoretical mean) of such a variable is
$$E\lbrack Y\rbrack = \exp\{\eta + \tau/2\},$$ while the variance is:
$$V\lbrack Y\rbrack = \left\lbrack \exp\{\tau\} - 1 \right\rbrack\exp\{ 2\eta + \tau\}.$$

Since we are usually interested in parametrization regression models
through the mean of a distribution, we work with an alternative
parametrization of this distribution. In particular, to parametrize the
distribution in terms of the mean (denoted $\mu$) on the original scale.
The simplest way to achieve this is as follows by solving the following
system of equations for $\eta$ and $\tau$$$\begin{aligned}
\mu & {= \exp\{\eta + \tau/2\}} \\
\phi & {= 2\tau,}
\end{aligned}$$ which yields $$\begin{array}{r}
{\eta = \log\mu - \phi} \\
{\tau = \frac{\phi}{2}.}
\end{array}$$

Alternatively, one may set
$$\phi = \left\lbrack \exp\{\tau\} - 1 \right\rbrack\exp\{ 2\eta + \tau\}.$$
Which yields: $$\begin{array}{r}
{\eta = \log\left( \frac{\mu^{2}}{\sqrt{\phi + \mu^{2}}} \right)} \\
{\tau = \log\left( \frac{\phi + \mu^{2}}{\mu^{2}} \right).}
\end{array}$$

For simplicity, we work with the first option, which yields the
following density:
$$g(y) = (\pi\phi)^{-1/2}y^{-1}\exp\left\{ -\frac{1}{\phi}\left( \log y - \log\mu + \phi \right)^{2} \right\}.$$

## Gamma

Similarly, let~$Y \sim Gamma(\alpha,\beta)$, with density given by:
$$g(y) = \frac{\beta^{\alpha}}{\Gamma(\alpha)}y^{\alpha - 1}\exp\{-\beta y\},$$
while its mean and variance are $$E\lbrack Y\rbrack = \alpha/\beta$$ and
$$V\lbrack Y\rbrack = \alpha/\beta^{2},$$ respectively.

Using a similar strategy as before, we get a parametrization in terms of
the mean as follows: $$\begin{aligned}
\mu & {= \alpha/\beta} \\
\phi & {= \alpha.}
\end{aligned}$$ Solving for $\alpha$ and $\beta$ gives:
$$\begin{aligned}
\alpha & {= \phi} \\
\beta & {= \phi/\mu.}
\end{aligned}$$ which yields the following density (parametrized in
terms of the mean)
$$g(y) = \frac{1}{\Gamma(\phi)}\frac{\phi^{\phi}}{\mu^{\phi}}y^{\phi - 1}\exp\left\{ -\frac{\phi}{\mu}y \right\},$$

## References
