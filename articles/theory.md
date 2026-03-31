# Theoretical Background

------------------------------------------------------------------------

## TL; DR

The table below details each parameter in the DRM along with their code
representation and meaning.

|  Parameter   |   Code    |                Meaning                 |
|:------------:|:---------:|:--------------------------------------:|
| $\beta_{t}$  | `beta_t`  |    Reg. coef. for proba. of absence    |
| $\beta_{r}$  | `beta_r`  |       Reg. coef. for recruitment       |
|    $\phi$    |   `phi`   |          Dispersion parameter          |
| $\beta_{s}$  | `beta_s`  |        Reg. coef. for survival         |
|    $\xi$     |   `xi`    | ${logit}(\rho) = \cdots\xi{\log}(\mu)$ |
|   $\zeta$    |  `zeta`   |  Prob. of remaining in the same patch  |
|   $\alpha$   |  `alpha`  |       AR(1) temporal correlation       |
| $\sigma_{t}$ | `sigma_t` |          AR(1) conditional SD          |
|   $z_{t}$    |   `z_t`   |          AR(1) random effects          |
| $\sigma_{i}$ | `sigma_i` |    Patch level iid random effect SD    |
|   $z_{i}$    |   `z_i`   |     Patch level iid random effects     |
| $\sigma_{s}$ | `sigma_s` |   Patch level ICAR random effect SD    |
|   $z_{s}$    |   `z_s`   |    Patch level ICAR random effects     |

## Notation

Before diving into coding, let us establish some common notation so we
understand what we are estimating.

Ecological data on species density often exhibits zeros. This occurs
because species may be patchily distributed or because sampling methods
may fail to detect individuals at low densities. Since densities are
continuous non-negative quantities, we cannot use discrete probability
distributions (e.g., Poisson or Binomial) to model this type of
phenomena. However, standard continuous probability distributions assign
a zero probability to specific values (that is, ${Pr}(Y_{t,i} = c) = 0$,
for any constant $c$ and continuous probability distribution). To
address this, we use zero-augmented (Hurdle/Delta) models, which
explicitly account for the excess zeros by combining a discrete and a
continuous probability distribution through two components:

1.  **A probability of observing a zero (or probability of absence)**:
    This accounts for the “true” absence or non-detection of the
    species.

2.  **A continuous distribution for positive values**: This models the
    density of the species when it is present.

Formally, let $Y_{t,i}$ be a random variable representing the density
(individuals per unit of area) of a focal species at time $t$ and
patch/site $i$. We denote realizations (i.e., what we observe) of
$Y_{t,i}$ by $y_{t,i}$. We define a *zero-augmented* probability density
function (pdf) as follows:
$$f(y_{t,i} \mid \mu_{t,i},\phi,\rho_{t,i}) = \begin{cases}
{\rho_{t,i},} & {{\mspace{6mu}\text{if}\mspace{6mu}}y_{t,i} = 0,} \\
{(1 - \rho_{t,i})g(y_{t,i} \mid \mu_{t,i},\phi),} & {{\mspace{6mu}\text{if}\mspace{6mu}}y_{t,i} > 0,}
\end{cases}\qquad(1)$$ where

- $\rho_{t,i} = {Pr}(Y_{t,i} = 0)$ is the probability of observing a
  zero density.

- $g(\cdot \mid \mu_{t,i},\phi)$ is the pdf of a continuous probability
  distribution with expected value (or theoretical mean) $\mu_{t,i}$ and
  additional parameter $\phi$. This distribution represents the
  densities of the species when they are present.

In [Equation 1](#eq-za_dens), we can choose different probability
distributions to specify the pdf $g(\cdot \mid \mu_{t,i},\phi)$. Popular
choices include the Log-Normal, Log-Logistic, and Gamma distributions,
which can be parametrized in terms of their mean facilitating modeling
([Ye et al. 2021](#ref-ye2021comparisons)).

## Species Distribution Models

A typical SDM based on [Equation 1](#eq-za_dens) relates the parameters
of the zero-augmented distribution to environmental covariates. Let
$t = 1,\ldots,T$ represent time points and $i = 1,\ldots,P$ represents
patches/sites. The first assumptions of SDMs is that
$Y_{1,1},\ldots,Y_{T,P}$ are independent random variables when
conditioned on the environmental factors. From here on, only SDMs based
on (generalized) linear models will be discussed.

We define as $\mathbf{x}_{t,i}^{(1)}$ a vector of $p_{1}$ covariates
related to the probability of observing a zero density at time $t$ and
site $i$, and $\mathbf{x}_{t,i}^{(2)}$ as a vector of $p_{2}$ covariates
related to the (theoretical) mean density when the species is present.
Note that, $\mathbf{x}_{t,i}^{(1)}$ and $\mathbf{x}_{t,i}^{(2)}$ may
depend on a different set of environmental factors.

We define linear-predictors as follows: $$\begin{aligned}
\eta_{t,i}^{(1)} & {= {\mathbf{β}}^{(1)}\mathbf{x}_{t,i}^{(1)}} \\
\eta_{t,i}^{(1)} & {= {\mathbf{β}}^{(2)}\mathbf{x}_{t,i}^{(2)},}
\end{aligned}\qquad(2)$$ where ${\mathbf{β}}^{(1)}$ and
${\mathbf{β}}^{(2)}$ are vectors of regression coefficients. Finally,
through the introduction of suitable link functions, $h_{1}$ and
$h_{2}$, we relate the linear predictors (and consequently environmental
factors) to the probability of absence $\rho_{t,i}$ and mean density
$\mu_{t,i}$ as follows: $$\begin{aligned}
{h_{1}(\rho_{t,i})} & {= \eta_{t,i}^{(1)}} \\
{h_{2}(\mu_{t,i})} & {= \eta_{t,i}^{(2)}.}
\end{aligned}\qquad(3)$$

The link functions will ensure that the estimated absence probabilities
and mean densities are on the appropriate scale. For instance, $h_{1}$
should be a function that take a number between 0 and 1 as input and
outputs a real number. Common choices for $h_{1}$ are the $logit$ and
the $cloglog$ functions. The latter is asymmetric and preferred when the
proportion of zeros is very high or low. For $h_{2}$, the $\log$
function is often used as it ensures that the predicted mean densities
are positive.

The assumption of independence can be relaxed by introducing random
effects into the model.

## Dynamic Range Models

Similar to SDMs, a DRM based [Equation 1](#eq-za_dens) also relates the
parameters of the zero-augmented distribution to environmental
covariates. However, the DRMs offer a more comprehensive approach by
incorporating explicit representations of the underlying biological
processes that govern species distributions. This advancement enables
the integration of expert knowledge and empirical data regarding species
demography, dispersal patterns, and mechanistic interactions with
environmental factors ([Pagel and Schurr
2012](#ref-pagel2012forecasting)).

The assumptions of our DRM are placed on the mean density $\mu_{t,i}$,
with the probability of absence remaining the same as in the SDM. In
fact, those assumptions can be interpreted as a non-linear function of
the environmental variables. In particular, we make assumptions about
the (unobserved[¹](#fn1)) age-structure behind the densities at a given
time point and patch.

Define $Y_{a,t,i}$ as a random variable representing the *unobserved*
density (individuals per unit of area) for individuals of age $a$, where
$a \in \{ 1,\ldots,A\}$. Similarly, we denote the expected density (or
theoretical mean for a given age, time, and patch) as:
$$\lambda_{a,t,i} = {\mathbb{E}}\lbrack Y_{a,t,i}\rbrack.$$ We do not
observe $Y_{a,i,t}$ and assume the previously defined (overall) density,
which we actually observe, to be the sum of the densities across all age
groups. In other words, $Y_{t,i} = \sum_{a}Y_{a,t,i}$. Consequently,
$$\mu_{t,i} = \sum\limits_{a}\lambda_{a,t,i}.\qquad(4)$$

Biological processes and additional assumptions are encoded through the
expected “age densities” $\lambda_{a,t,i}$. Some key process we consider
are: recruitment, survival, and movement. They are described below.

### Recruitment

We begin by defining the expected density for recruits (fish at age
$a = 1$) at time $t$ and patch $i$ as:
$${\mathbb{E}}\lbrack Y_{1,i,t}\rbrack = \lambda_{1,i,t} = {\exp}\{\psi\},\qquad(5)$$

where ${\exp}\{\psi\}$ represents the overall mean recruitment per unit
of area. The assumption in [Equation 5](#eq-recruits) is quite
restrictive. To make the model more realistic, we can replace the
constant $\psi$ with

$$\psi_{t,i} = {\mathbf{β}}_{r}\mathbf{x}_{t,i}^{(r)} + z_{t,i}^{(r)},\qquad(6)$$

where~$\mathbf{x}_{t,i}^{(r)}$ is the vector of environmental drivers
associated to recruitment at time~$t$ and patch~$i$, ${\mathbf{β}}_{r}$
is a vector of corresponding regression coefficients. This depicts a
log-linear relationship between recruitment and the environment.

To account for the influence of time on recruitment (or SDM densities),
we introduce a random effect in Equation [Equation 6](#eq-lmrec). In the
current version of the model, this random effect is the same across all
groups for a given time point, meaning $z_{t,i}^{(r)} = z_{t}^{(r)}$ for
every $i$. This creates a temporal dependence, where the recruitment or
density at one point in time is related to its value at the previous
time point. We model this temporal dependence using an autoregressive
process of order 1 (AR(1)). Formally, this is defined as:

$$z_{t}^{(r)} = \alpha z_{t - 1}^{(r)} + \varepsilon_{t},\qquad(7)$$

where $\varepsilon_{1},\ldots,\varepsilon_{T}$ are independent and
identically distributed residual terms from a zero-mean Normal
distribution with variance $\tau^{2}$, and $\alpha$ is the temporal
autocorrelation parameter. In this model, the recruitment (or density
from [Equation 2](#eq-linpreds)) at time $t$ is a function of its value
at time $t - 1$, plus a random error term.

### Survival

In the proposed model, the expected density at age $a$, time $t$, and
patch $i$ evolves over time according to survival rates denoted
$s_{a,t,i}$. Formally, we have

$${\mathbb{E}}\lbrack Y_{a,t,i}\rbrack = \lambda_{a,t,i} = \lambda_{a - 1,t - 1,i}s_{a - 1,t - 1},$$

As with recruitment, the survival rates may vary by patch/site $i$ and
time $t$:

$$s_{a,t,i} = {\exp}\{{\mathbf{β}}_{m}^{\top}\mathbf{x}_{t,i}^{(m)}\},\qquad(8)$$

where $\mathbf{x}_{t,i}^{(m)}$ represents environmental drivers
potentially different from those affecting recruitment, and
${\mathbf{β}}_{m}$ are associated regression coefficients. Note that,
despite the age-specific index~$a$, the survival rates are constant
across age-groups.

In [Equation 8](#eq-surv), we may use external information to allow for
different survival rates across age-groups. For instance, in fisheries
research, stock assessments often provide estimates of age-specific
fishing mortality instantaneous rates for several species. It is
reasonable, if not logical, to assume that the survival rates for an
age-group with a higher fishing mortality is lower. In that case,
Equation~ may be updated to:

$$s_{a,t,i} = {\exp}\{{\mathbf{β}}_{m}^{\top}\mathbf{x}_{t,i}^{(m)} - f_{a,t}\},$$

where~$f_{a,t}$ is the instantaneous rate of fishing mortality for the
age-group~$a$ at time~$t$.

In our current implementation, we do not allow for a random effect on
mortality/survival yet.

### Movement

Movement is an important demographic process affecting species’
densities. In our package, we implement a simplistic movement routine
described as follows. Denote $P$ the total number of patches in the
present study. For each age-group, we define a $P \times P$ movement
matrix denoted $\mathbf{M}_{a}$. The element $\{ i,j\}$ in
$\mathbf{M}_{a}$ represents the probability that individuals of age $a$
move from patch $i$ at time $t - 1$ to patch $j$ at time $t$. The
movement matrix $\mathbf{M}_{a}$ is derived from adjacency matrix
$\mathbf{A}$, which encodes connections between patches. In
$\mathbf{A}$, the element $\{ i,j\}$ is given by

$$\mathbf{A}_{ij} = \begin{cases}
{1/N(i),} & {{\mspace{6mu}\text{if}\mspace{6mu}}i \sim j,} \\
{0,} & {{\mspace{6mu}\text{otherwise}},}
\end{cases}$$

where $N(i)$ is the number of neighbors of patch $i$, and $i \sim j$
indicates that $i$ and $j$ are neighbors (i.e., share borders).

In addition to the adjacency matrix $\mathbf{A}$, the user must also
indicate the age-groups \`\`allowed’’ to move. If an age group is not
assumed to move, the movement matrix becomes the identity matrix. Given
those inputs, we define $\zeta$ as the probability of individuals
remaining in the same patch between two time points. The movement matrix
is then constructed as:

$$\mathbf{M}_{a} = \zeta\mathbf{I}_{P} + (1 - \zeta)\mathbf{A},\qquad(9)$$

where $\mathbf{I}_{P}$ is a $P \times P$ identity matrix. In other
words, while the probability of remaining on the same patch if fixed,
the probability of moving~(i.e., $(1 - \zeta)$) is evenly distributed
across neighbors. We could extend this setting in, at least, three ways:
(1) allowing $\zeta$ to vary with time and patch; (2) having different
movement probabilities for different age-groups, and (3) make the
movement probabilities to vary according to environmental conditions~.
Any of these extensions, however, would represent significant increase
in computation. Therefore, for now, we opted to stick with the simplest
option.

Once a movement matrix is defined, we need to adjust the age-specific
expected densities accordingly. Let ${\mathbf{λ}}_{a,t,\cdot}$ be a
vector of length $P$ (number of patches) representing the expected
densities of age $a$ at time $t$ for every patch in the given study.
After movement, the expected density is:

$${\widetilde{\mathbf{λ}}}_{a,t,\cdot} = \mathbf{M}_{a}{\mathbf{λ}}_{a - 1,t - 1,\cdot},\qquad(10)$$

where $\mathbf{M}_{a}$ is the aforementioned movement matrix.

### Selectivity at age

Our main data source, trawl surveys, often exhibit age-based
selectivity, where younger fish have a lower probability of being caught
due to their smaller size. In particular, the probability of capturing
an animal of a given length given it encounters the gear is denoted
*selectivity at length*. On the other hand, *selectivity at age* is the
probability of capturing a fish of a certain age, provided it
“encountered the gear”. We represent selectivity at age by the vector
$\mathbf{v} = \lbrack v_{1},\ldots,v_{A}\rbrack^{\top}$, where $A$ is
the total number of age groups. Each element $v_{a}$ of this vector is a
value between 0 and 1, inclusive, representing the selectivity for age
group $a$. Our package does not estimate $\mathbf{v}$ and, unless
provided by the user, assumes the selectivity is 1 for every age group.

When a vector $\mathbf{v}$ is provided, our model assumes:

$$Y_{t,i} = \sum\limits_{a}v_{a}Y_{a,t,i}$$

and, consequently,

$$\mu_{t,i} = \sum\limits_{a}v_{a}\lambda_{a,t,i}.$$

Note that, our model cannot estimate selectivity. Instead, the model
assumes it is provided as known without error. Of course, this is an
oversimplification, however, we seldomly have the necessary data to
estimate selectivity readily available. cvf

## Bayesian Inference and Forecasting

To complete the Bayesian model specification, we need to specify prior
distributions for all model parameters. This involves defining our prior
beliefs about the likely values of these parameters before observing any
data.

First, we define the priors for the fixed effects regression
coefficients, typically denoted by $\mathbf{β}$ with superscrip. The
default priors for these coefficients are uncorrelated zero-mean Normal
distributions with marginal variances of 1. The parameter $\phi$, which
controls the scale of the non-zero part of the pdfs in
[Equation 1](#eq-za_dens) (e.g., Log-Normal, Gamma, or Log-Logistic),
has a Gamma prior with shape and rate parameters equal to 2 and 1,
respectively. In general, the scale of this parameter varies with the
pdf chosen.

Next, we specify the priors for the parameters associated with the AR(1)
process in Equation [Equation 7](#eq-ar1). We take a hierarchical
approach, placing a Normal prior on the log of the square root of the
conditional variance (i.e., the log of the conditional standard
deviation). This prior has a default mean of $-2$ and and a standard
deviation of 0.25, favoring a model with no AR(1) term. For the
autocorrelation term, we employ a beta prior on $\alpha$. The default
parameters on this prior are $0.5$ and $0.5$. Lastly, we the same Beta
prior on $\zeta$. The priors and default hyperparameters are highlighted
in the table below:

|  Parameter   |   Code    | Priors                           |   Default hyperparameters    |
|:------------:|:---------:|----------------------------------|:----------------------------:|
| $\beta_{t}$  | `beta_t`  | $N(m_{t},s_{t})$                 |    $m_{t} = 0,s_{t} = 1$     |
| $\beta_{r}$  | `beta_r`  | $N(m_{r},s_{r})$                 |    $m_{r} = 0,s_{r} = 1$     |
|    $\phi$    |   `phi`   | ${Gamma}(a_{p},b_{p})$           |    $a_{p} = 2,b_{p} = 1$     |
| $\beta_{s}$  | `beta_s`  | $N(m_{s},s_{s})$                 |    $m_{s} = 0,s_{s} = 1$     |
|    $\xi$     |   `xi`    | \${- \rm LN}(m_lxi, s_lxi)\$     |  $m_{lxi} = 0,s_{lxi} = 1$   |
|   $\zeta$    |  `zeta`   | ${Beta}(a_{z},b_{z})$            |   $a_{z} = .5,b_{z} = .5$    |
|   $\alpha$   |  `alpha`  | ${Beta}(a_{a},b_{a})$            |   $a_{a} = .5,b_{a} = .5$    |
| $\sigma_{t}$ | `sigma_t` | ${LN}(m_{l}st,s_{l}st)$          | $m_{lst} = -2,s_{lst} = .25$ |
|   $z_{t}$    |   `z_t`   | $z_{t} \sim AR(1)$               |             N/A              |
| $\sigma_{i}$ | `sigma_i` | ${LN}(m_{l}si,s_{l}si)$          | $m_{lsi} = -2,s_{lsi} = .25$ |
|   $z_{i}$    |   `z_i`   | $z_{i} \sim N(0,\sigma_{i}^{2})$ |             N/A              |
| $\sigma_{s}$ | `sigma_s` | ${LN}(m_{l}ss,s_{l}ss)$          | $m_{lss} = -2,s_{lss} = .25$ |
|   $z_{s}$    |   `z_s`   | $z_{s} \sim N(0,\sigma_{i}^{2})$ |             N/A              |

Denote by $\mathbf{θ}$ a set containing all the parameters associated
with the model we have chosen. The posterior distribution of
$\mathbf{θ}$ (and $\mathbf{z}$) given the observed data $\mathbf{y}$ and
environmental factors $\mathbf{X}$ is

$$\pi({\mathbf{θ}},\mathbf{z} \mid \mathbf{y}) \propto p(\mathbf{y} \mid \mathbf{X},\mathbf{z},{\mathbf{β}},{\mathbf{γ}})\pi(\mathbf{z} \mid {\mathbf{σ}},{\mathbf{δ}})\pi({\mathbf{θ}}),\qquad(11)$$

where $\pi(\cdot)$ denotes prior and posterior distributions of its
arguments. Importantly, the distribution placed on the random effects
$\mathbf{z}$ also functions as a prior, as we draw samples from these
random effects distribution conditioned on the observed data.

### Parameter estimation

The estimation of the model parameters is obtained from the posterior
distribution defined in [Equation 11](#eq-posterior). As the posterior
distribution is intractable, we use the No-U-Turn ([Hoffman and Gelman
2014](#ref-hoffman2014no)) Markov Chain Monte Carlo~(MCMC) sampler
available in `Stan` to draw samples from it. The No-U-Turn algorithm
samples all the model parameters jointly and efficiently exploits the
parameter space using automatic differentiation. Additionally, it
eliminates the need for hand-tuning, making it a highly convenient
sampler. We initialize the parameters with samples from their respective
prior distributions. The latent random effects are initialized from a
standard Gaussian distribution. The number of samples and the warm-up
period of the MCMC algorithm are application dependent. We assess the
convergence of the chains using the split-$\widehat{R}$ diagnostic
([Vehtari et al. 2021](#ref-vehtari2021rank)). Finally, the parameters’
point estimates and 95% credible intervals~(CI) are obtained as,
respectively, the median and percentiles ($2.5$ and $97.5$, unless
otherwise stated) of the marginal MCMC samples.

### Model comparison

The assessment of goodness-of-fit GoF is carried out using the
leave-one-out information criterion \[Vehtari et al.
([2017](#ref-vehtari2017practical)); LOOIC\]. Lower LOOIC values
indicate a better fit.

### Predictions/Forecasting

Given the posterior MCMC samples, we can obtain Monte Carlo samples from
the posterior predictive distribution under different environmental
conditions and make forecasts. Suppose we wish to compute forecasts $k$
time points ahead. We define an $k$-dimensional variable
$\mathbf{Y}^{*}$ representing the densities at these time points.
Assuming environmental at these time points are also available, we can
generate predictions accounting for uncertainty using the posterior
predictive distribution of $\mathbf{Y}^{*}$:

$$p(\mathbf{y}^{*} \mid \mathbf{y}) = \int p(\mathbf{y}^{*} \mid \mathbf{z}^{*},{\mathbf{θ}})p(\mathbf{z}^{*} \mid \mathbf{z},{\mathbf{θ}})\pi({\mathbf{θ}} \mid \mathbf{y})d{\mathbf{θ}},\qquad(12)$$

where

$\mathbf{z}^{*} = {(z_{T + 1}),\ldots,z_{T + k}))}^{\top}$ is a
realization of the AR(1) process at the new time points, and
$\mathbf{θ}$ represents the model parameters. Moreover, the pdf
$p(\mathbf{y}^{*} \mid \mathbf{z}^{*},{\mathbf{θ}})$ can be obtained as
in [Section 3](#sec-drm) or [Section 2](#sec-sdm).

## Pontential issues with estimation & Closed formulas for special cases

Define

$$s_{a,t,i} = {\exp}\{ m_{i,t} - f_{a,t}\},$$

as the environment dependent survival for age $a$ at time $t$ and site
$i$.

Then, in the scenario without movement, the expected density for age $a$
at time $t$ and site $i$ can be written as:

$$\begin{aligned}
\lambda_{a,t,i} & {= {\exp}\{\psi_{t,i} + z_{t}^{(r)}\prod\limits_{k = 1}^{a - 1}s_{k,t,i}\}} \\
 & {= {\exp}\{\psi_{t,i} + z_{t}^{(r)}\}{\exp}\left\{ (a - 1)m_{t,i} - \sum\limits_{k = 1}^{a - 1}f_{a,t} \right\}.}
\end{aligned}$$

Therefore, a closed form for the expected density $\mu_{i,t}$ is:
$$\begin{aligned}
\mu_{t,i} & {= \sum\limits_{a = 1}^{A}\lambda_{a,t,i}} \\
 & {= {\exp}\{\psi_{t,i} + z_{t}^{(r)}\}\left( 1 + \sum\limits_{a = 2}^{A}\exp\left\{ (a - 1)m_{t,i} - \sum\limits_{k = 1}^{a - 1}f_{k,t} \right\} \right).}
\end{aligned}$$

Or, with selectivity,

$$\begin{aligned}
\mu_{t,i} & {= \sum\limits_{a = 1}^{A}v_{a}\lambda_{a,t,i}} \\
 & {= {\exp}\{\psi_{t,i} + z_{t}^{(r)}\}\left( v_{1} + \sum\limits_{a = 2}^{A}v_{a}\exp\left\{ (a - 1)m_{t,i} - \sum\limits_{k = 1}^{a - 1}f_{k,t} \right\} \right).}
\end{aligned}$$

## References

Hoffman, Matthew D., and Andrew Gelman. 2014. “The No-U-turn Sampler:
Adaptively Setting Path Lengths in Hamiltonian Monte Carlo.” *Journal of
Machine Learning Research* 15 (1): 1593–623.

Pagel, Jörn, and Frank M. Schurr. 2012. “Forecasting Species Ranges by
Statistical Estimation of Ecological Niches and Spatial Population
Dynamics.” *Global Ecology and Biogeography* 21 (2): 293–304.
https://doi.org/<https://doi.org/10.1111/j.1466-8238.2011.00663.x>.

Vehtari, Aki, Andrew Gelman, and Jonah Gabry. 2017. “Practical Bayesian
Model Evaluation Using Leave-One-Out Cross-Validation and WAIC.”
*Statistics and Computing* 27: 1413–32.

Vehtari, Aki, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
Paul-Christian Bürkner. 2021. “Rank-Normalization, Folding, and
Localization: An Improved $\widehat{R}$ for Assessing Convergence of
MCMC (with Discussion).” *Bayesian Analysis* 16 (2): 667–718.

Ye, Tairan, Victor H Lachos, Xiaojing Wang, and Dipak K Dey. 2021.
“Comparisons of Zero-Augmented Continuous Regression Models from a
Bayesian Perspective.” *Statistics in Medicine* 40 (5): 1073–100.

------------------------------------------------------------------------

1.  By “unobserved”, I mean something we do not actually measure or have
    records of (although it can be inferred to some extent from the
    length composition)
