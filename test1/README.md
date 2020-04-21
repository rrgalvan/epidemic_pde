# Test 1. 1D SIR + reaction + diffusion

Numerical resolution of a SIR + reaction + difusión problem:
$$
   S_t - alpha \Delta S = mu N - mu S - beta S I
   I_t - alpha \Delta I = -(mu + nu) I  + beta S I
   + homogeneous b.c.
$$

Here:
  * S, I: susceptible and infectious individuals,
  * N: total population (fixed positive constant)
From this model, one can compute R (recovered individuals) as
  * $R_t - alpha I = nu I + mu R$

Test based on:
> [1] Numerical modelling of an SIR epidemic model with diffusion
> Settapat Chinviriyasit, Wirawan Chinviriyasit
> Applied Mathematics and Computation 216 (2010) 395–409
