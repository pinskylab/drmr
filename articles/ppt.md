# Parameters, priors, and toggles

The `fit_drm` function is designed to be flexible, allowing different
model features to be enabled or disabled through a set of toggles. This
vignette provides a comprehensive guide to the parameters used in the
DRM, their meaning, the priors they are assigned, and the specific
toggles required to activate them.

### DRM model parameters and priors

The table below details each parameter in the DRM. The “Toggle” column
indicates the condition needed to include a parameter in the model. If a
parameter is not activated by its toggle, it is excluded from the model,
and its prior is not used.

|  Parameter   |   Code    |               Meaning                |         Toggle          | Priors                                        |   Default hyperparameters    |
|:------------:|:---------:|:------------------------------------:|:-----------------------:|-----------------------------------------------|:----------------------------:|
| $\beta_{t}$  | `beta_t`  |   Reg. coef. for proba. of absence   |           N/A           | $N\left( m_{t},s_{t} \right)$                 |    $m_{t} = 0,s_{t} = 1$     |
| $\beta_{r}$  | `beta_r`  |      Reg. coef. for recruitment      |           N/A           | $N\left( m_{t},s_{t} \right)$                 |    $m_{r} = 0,s_{r} = 1$     |
|    $\phi$    |   `phi`   |         Dispersion parameter         |           N/A           | ${Gamma}\left( a_{p},b_{p} \right)$           |    $a_{p} = 2,b_{p} = 1$     |
| $\beta_{s}$  | `beta_s`  |       Reg. coef. for survival        |       `est_surv`        | $N\left( m_{s},s_{s} \right)$                 |    $m_{s} = 0,s_{s} = 1$     |
|    $\xi$     |   `xi`    | ${logit}(\rho) = \cdots\xi\log(\mu)$ |        `rho_mu`         | \${- \rm LN}(m_lxi, s_lxi)\$                  |  $m_{lxi} = 0,s_{lxi} = 1$   |
|   $\zeta$    |  `zeta`   | Prob. of remaining in the same patch |       `movement`        | ${Beta}\left( a_{z},b_{z} \right)$            |   $a_{z} = .5,b_{z} = .5$    |
|   $\alpha$   |  `alpha`  |      AR(1) temporal correlation      | `ar_re` $\neq$`"none"`  | ${Beta}\left( a_{a},b_{a} \right)$            |   $a_{a} = .5,b_{a} = .5$    |
| $\sigma_{t}$ | `sigma_t` |         AR(1) conditional SD         | `ar_re` $\neq$`"none"`  | ${LN}\left( m_{l}st,s_{l}st \right)$          | $m_{lst} = -2,s_{lst} = .25$ |
|   $z_{t}$    |   `z_t`   |         AR(1) random effects         | `ar_re` $\neq$`"none"`  | $z_{t} \sim AR(1)$                            |             N/A              |
| $\sigma_{i}$ | `sigma_i` |   Patch level iid random effect SD   | `iid_re` $\neq$`"none"` | ${LN}\left( m_{l}si,s_{l}si \right)$          | $m_{lsi} = -2,s_{lsi} = .25$ |
|   $z_{i}$    |   `z_i`   |    Patch level iid random effects    | `iid_re` $\neq$`"none"` | $z_{i} \sim N\left( 0,\sigma_{i}^{2} \right)$ |             N/A              |
| $\sigma_{s}$ | `sigma_s` |  Patch level ICAR random effect SD   | `sp_re` $\neq$`"none"`  | ${LN}\left( m_{l}ss,s_{l}ss \right)$          | $m_{lss} = -2,s_{lss} = .25$ |
|   $z_{s}$    |   `z_s`   |   Patch level ICAR random effects    | `sp_re` $\neq$`"none"`  | $z_{s} \sim N\left( 0,\sigma_{i}^{2} \right)$ |             N/A              |

### SDM

TBD

## References
