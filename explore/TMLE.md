# TMLE Implementation for `seine`

## Overview

This document develops a TMLE for the global ecological estimands in `seine`.
In addition to targeting the usual efficient influence function (EIF), the
estimator should respect the aggregate accounting identity

```
px' B = py,
```

where `B[j, k]` is the estimated outcome-`k` mean for predictor group `j`,
`px` is the population-weighted predictor margin, and `py` is the
population-weighted outcome margin.

This identity should be part of the model and targeting equations, rather than
enforced by projecting the final matrix of estimates. The key implementation
change is to maintain one collection of counterfactual group curves throughout
targeting and to construct the factual regression as their mixture. The current
implementation instead updates the factual and counterfactual regressions
separately. Since its update links are nonlinear, it need not preserve their
mixture relationship.

The derivation below assumes the conditional-aggregate-representativeness
condition needed for a population-weighted target. Equivalently, the nuisance
regression is understood to be the appropriate `N`-weighted conditional mean.
Including `N` among the covariates can make this assumption more plausible and
can reduce the accounting discrepancy of plug-in and one-step estimators. It
does not fix the present TMLE because the separate targeting updates can
reintroduce incoherence.

## Notation and coherent nuisance model

Let there be `J` predictor groups and `K` outcomes. Define

```
W       = N / E[N]
p_j     = E[W X_j]
q_k     = E[W Y_k]
b_jk(z) = Q_k(e_j, z)
Q_k(x,z) = sum_j x_j b_jk(z)
B_jk    = E[W X_j b_jk(Z)] / p_j.
```

Here `p = px`, `q = py`, and `E[W] = 1`. Because the components of `X`
partition the population, `sum_j X_j = 1` and `sum_j p_j = 1`.

The factual regression must always be constructed from the group curves:

```
Q_k(X_i,Z_i) = sum_j X_ij b_jk(Z_i).
```

It then follows mechanically that

```
sum_j p_j B_jk
    = E[W sum_j X_j b_jk(Z)]
    = E[W Q_k(X,Z)].
```

Thus, the accounting identity follows from nuisance coherence plus the weighted
margin equation

```
E[W {Y_k - Q_k(X,Z)}] = 0.
```

In a sample, a coherent targeted regression satisfying

```
P_n[W {Y_k - Q*_k(X,Z)}] = 0
```

produces plug-in estimates satisfying `p_hat' B_hat = q_hat` exactly, up to
optimization tolerance.

## Unrestricted EIF and the accounting restriction

Because `p_j` is itself random, the ratio-form EIF for `B_jk` is

```
D_B,jk(O) =
    W X_j / p_j {b_jk(Z) - B_jk}
    + alpha_j(O) {Y_k - Q_k(X,Z)}.
```

The margin EIFs are

```
D_p,j(O) = W (X_j - p_j),
D_q,k(O) = W (Y_k - q_k).
```

Let

```
r_k(O) = Y_k - Q_k(X,Z),
s(O)   = sum_j p_j alpha_j(O).
```

Differentiating `sum_j p_j B_jk - q_k = 0` and substituting the EIFs gives

```
sum_j p_j D_B,jk
    + sum_j B_jk D_p,j
    - D_q,k
  = {s(O) - W} r_k(O).
```

Consequently, compatible population Riesz representers obey

```
sum_j p_j alpha_j(O) = W
```

in the relevant tangent space. This is the population identity behind exact
accounting for an unregularized DML estimator. Estimated, penalized Riesz
representers need only satisfy it approximately.

### Independent coordinates and restricted canonical gradient

For simplex outcomes, use the nonredundant joint parameter

```
theta = (
    vec(B[, 1:(K-1)]),
    p[1:(J-1)],
    q[1:(K-1)]
).
```

The `K-1` independent restrictions are

```
g_k(theta) = sum_j p_j B_jk - q_k = 0,
k = 1, ..., K-1.
```

Writing `p_J = 1 - sum_(a<J) p_a`, the Jacobian `A = d g / d theta`
has entries

```
A[k, B_jl] = p_j 1(k = l),
A[k, p_a]  = B_ak - B_Jk,
A[k, q_l]  = -1(k = l).
```

The `-I` block makes `A` full row rank. If `D` is the canonical gradient in a
larger unrestricted model and `Sigma = E[D D']`, the canonical gradient after
imposing the smooth equality restriction is

```
D_c = {I - Sigma A' (A Sigma A')^(-1) A} D,
```

with covariance

```
Var(D_c) =
    Sigma - Sigma A' (A Sigma A')^(-1) A Sigma.
```

A generalized inverse may be used with redundant coordinates. This projection
is useful for deriving and checking the constrained EIF. It is not, by itself,
a substitution estimator: the targeted nuisance functions must also remain in
the restricted model.

For componentwise outcomes without a sum-to-one restriction, all `K`
accounting equations are independent. Use

```
theta = (vec(B), p[1:(J-1)], q)
```

and the same Jacobian formulas for `k,l = 1, ..., K`. In that case `A` has
`K` rows and rank `K`.

## Coherent fluctuation in general form

Only the group curves are fluctuated. For a support-preserving inverse link
`G`, write

```
b_jk^epsilon(z) =
    G_k( eta_jk^0(z) + sum_a epsilon_a h_ajk(z) ),
```

and always define

```
Q_k^epsilon(x,z) = sum_j x_j b_jk^epsilon(z).
```

For a direction `a`, let

```
dot_b,ajk = d b_jk^epsilon / d epsilon_a at epsilon = 0,
dot_Q,ak  = sum_j X_j dot_b,ajk.
```

Under a loss whose score has the form

```
S_a(O) = sum_k {Y_k - Q_k} dot_Q,ak / V_k(Q),
```

a desired residual score

```
R_a(O) = sum_k c_ak(O) {Y_k - Q_k}
```

is obtained by choosing coherent group-curve directions satisfying

```
sum_j X_j dot_b,ajk = V_k(Q) c_ak.
```

The margin directions use `c_ak = W 1(a = k)`. The remaining directions use
the residual coefficients from an independent set of constrained EIFs. There
are generally many group-curve directions producing the same factual score;
a minimum-norm solution, with regularization where necessary, gives a
well-defined local fluctuation.

These displays take the targeting loss to have uniform observation weights.
If the loss instead uses observation weight `omega`, replace a desired
coefficient `c_a` in the fluctuation equations by `c_a / omega`. In
particular, the weighted-margin directions then use `W / omega`, so that the
resulting score remains `W(Y-Q)`. The loss weights used for numerical
efficiency and the target-population weights `W` need not be the same.

This construction is naturally iterative because `V(Q)`, the group-link
derivatives, and the restricted EIF coefficients change after each update.

## Sum-to-one simplex outcomes

Assume

```
b_j(z) in Delta^(K-1),
Q(x,z) in Delta^(K-1).
```

Use group-specific multinomial logits:

```
b_jk^epsilon(z) =
    softmax_k(
        log b_j^0(z) + sum_a epsilon_a h_aj(z)
    ).
```

The factual prediction is the mixture of the updated group probabilities:

```
Q_k^epsilon(x,z) = sum_j x_j b_jk^epsilon(z).
```

The observed-data loss is the multinomial mixture loss

```
L(epsilon) =
    -sum_i omega_i sum_k Y_ik log Q_ik^epsilon.
```

It is not a multinomial loss obtained by separately applying softmax to the
factual regression. For a logit direction `h_ajk`,

```
dot_b,ajk =
    b_jk {h_ajk - sum_l b_jl h_ajl},
dot_Q,ak = sum_j X_j dot_b,ajk.
```

The score is

```
S_a(O) =
    sum_k {Y_k - Q_k} dot_Q,ak / Q_k.
```

Only `K-1` outcome directions are independent. Include margin directions
whose induced factual derivatives satisfy, modulo an outcome-common term,

```
dot_Q,ak / Q_k = W 1(a = k),  a = 1, ..., K-1.
```

Solving their empirical scores gives

```
P_n[W(Y_k - Q*_k)] = 0,  k = 1, ..., K-1.
```

The final equation follows because both `Y` and `Q*` sum to one. Combined with
mixture coherence, this yields exact accounting and row-simplex constraints for
`B_hat`.

## Unconstrained componentwise outcomes

For outcomes supported on the real line, use an identity-link fluctuation:

```
b_jk^epsilon(z) =
    b_jk^0(z) + sum_a epsilon_a h_ajk(z),
Q_k^epsilon(x,z) = sum_j x_j b_jk^epsilon(z).
```

With squared-error loss

```
L(epsilon) =
    (1/2) sum_i omega_i sum_k
    {Y_ik - Q_ik^epsilon}^2,
```

the score is, up to sign,

```
S_a(O) =
    sum_k {Y_k - Q_k} sum_j X_j h_ajk(Z).
```

The margin directions can be especially simple: take

```
h_ajk(Z,N) = W 1(a = k)
```

for every group `j`. Since `sum_j X_j = 1`, their induced factual direction is
`dot_Q,ak = W 1(a = k)`. The `K` fitted score equations then impose

```
P_n[W(Y_k - Q*_k)] = 0
```

for every outcome. No accounting equation is redundant unless the outcomes
have an additional linear restriction.

## Two-sided componentwise bounds

Suppose every component lies in `[L, U]`, without a sum-to-one restriction.
Let `S = U - L` and use group-specific logistic curves:

```
b_jk^epsilon(z) =
    L + S expit(
        logit{(b_jk^0(z)-L)/S}
        + sum_a epsilon_a h_ajk(z)
    ).
```

Their mixture remains in `[L,U]`:

```
Q_k^epsilon(x,z) = sum_j x_j b_jk^epsilon(z).
```

Define

```
Ybar_k = (Y_k-L)/S,
Qbar_k = (Q_k-L)/S.
```

Use the componentwise quasi-binomial mixture loss

```
L(epsilon) =
    -sum_i omega_i sum_k [
        Ybar_ik log Qbar_ik^epsilon
        + (1-Ybar_ik) log(1-Qbar_ik^epsilon)
    ].
```

For a group-logit direction,

```
dot_b,ajk =
    {(b_jk-L)(U-b_jk)/S} h_ajk,
dot_Q,ak = sum_j X_j dot_b,ajk.
```

Up to a constant scale, the score is

```
S_a(O) =
    sum_k (Y_k-Q_k)
        dot_Q,ak / {(Q_k-L)(U-Q_k)}.
```

Thus a desired residual coefficient `c_ak` is obtained by solving

```
sum_j X_j dot_b,ajk =
    (Q_k-L)(U-Q_k) c_ak
```

up to the irrelevant common scale. In particular, the `K` margin directions
use `c_ak = W 1(a = k)`. Their score equations impose the weighted outcome
margins while the group-logistic parameterization preserves the bounds.

## One-sided componentwise bounds

### Finite lower bound

For support `[L, Inf)`, use a log link for the distance above the bound:

```
b_jk^epsilon(z) =
    L + (b_jk^0(z)-L)
        exp{sum_a epsilon_a h_ajk(z)}.
```

Then

```
Q_k^epsilon(x,z) =
    L + sum_j x_j {b_jk^epsilon(z)-L}
```

remains above `L`. With

```
Ytilde_k = Y_k-L,
Qtilde_k = Q_k-L,
```

use the Poisson quasi-loss

```
L(epsilon) =
    sum_i omega_i sum_k [
        Qtilde_ik^epsilon
        - Ytilde_ik log Qtilde_ik^epsilon
    ].
```

Here

```
dot_b,ajk = (b_jk-L) h_ajk,
dot_Q,ak = sum_j X_j dot_b,ajk,
```

and the score is, up to sign,

```
S_a(O) =
    sum_k (Y_k-Q_k) dot_Q,ak / (Q_k-L).
```

Choose directions satisfying

```
sum_j X_j dot_b,ajk = (Q_k-L) c_ak.
```

Setting `c_ak = W 1(a = k)` supplies the `K` accounting-margin equations.

### Finite upper bound

For support `(-Inf, U]`, apply the same construction to the distance below the
upper bound:

```
b_jk^epsilon(z) =
    U - (U-b_jk^0(z))
        exp{sum_a epsilon_a h_ajk(z)}.
```

Use the Poisson quasi-loss for

```
Ytilde_k = U-Y_k,
Qtilde_k = U-Q_k.
```

Now

```
dot_b,ajk = -(U-b_jk) h_ajk,
dot_Q,ak = sum_j X_j dot_b,ajk.
```

Accounting directions may be obtained by choosing

```
sum_j X_j dot_b,ajk = (U-Q_k) c_ak
```

with the sign of `h` chosen consistently with the reflected Poisson score.
Again, `c_ak = W 1(a = k)` targets the weighted factual residual margins.

## Direct constrained optimization versus margin directions

For every support case, an equivalent way to ensure accounting during
targeting is to minimize the appropriate coherent mixture loss subject to

```
P_n[W {Y_k-Q_k^epsilon}] = 0.
```

There are `K-1` constraints in the simplex case and `K` otherwise. A
sequential quadratic programming or augmented-Lagrangian solver could impose
these equations jointly with the restricted-EIF targeting directions.

This is not post-estimation calibration: the constraints are imposed on the
coherent nuisance model inside the likelihood optimization. Explicit margin
directions are preferable when they can be made least favorable, while direct
constraints may provide a simpler and more stable first implementation.

## Implementation sketch

1. **Construct a coherent initial nuisance representation.**
   Store only the counterfactual group curves `b0[[j]]`. Reconstruct the
   factual regression as `Q0 = sum_j X_j b0[[j]]`; do not independently use or
   update a stored factual fit. Check the initial mixture discrepancy and warn
   if the supplied regression model is not coherent.

2. **Select the support-specific link and mixture loss.**
   Use softmax/multinomial loss for simplex outcomes, identity/squared loss for
   unconstrained outcomes, logistic/quasi-binomial loss for two-sided bounds,
   and log/Poisson quasi-loss for one-sided bounds.

3. **Work in independent target coordinates.**
   Use `K-1` outcomes for the simplex case and all `K` outcomes otherwise.
   Parameterize `J-1` independent group contrasts plus the outcome margins, or
   form the joint EIF and its accounting-restricted projection.

4. **Make the estimated Riesz directions compatible.**
   Check the empirical analogue of

   ```
   sum_j p_hat_j alpha_hat_j = W.
   ```

   Estimate the representers jointly in the restricted tangent space, or
   project their residual-score coefficients using the loss-induced inner
   product. Do not project the final `B` matrix.

5. **Add explicit weighted-margin directions.**
   Include `K-1` directions for simplex outcomes and `K` otherwise. As an
   initial implementation, impose their score equations as equality
   constraints in the coherent mixture-loss optimization.

6. **Target iteratively.**
   At each iteration, compute the current restricted-EIF residual
   coefficients, solve for support-preserving group-curve directions, take a
   small loss-minimizing fluctuation step, reconstruct the factual mixture,
   and repeat until all empirical EIF and margin scores are below tolerance.

7. **Form the substitution estimate.**
   Compute

   ```
   B_hat_jk =
       P_n[W X_j b*_jk(Z)] / P_n[W X_j].
   ```

   This estimate inherits its support from the group curves and satisfies
   `p_hat' B_hat = q_hat` through nuisance coherence and the targeted margin
   equations.

8. **Compute covariance in independent coordinates.**
   Evaluate the restricted EIF, verify the differentiated accounting identity
   observation by observation, compute its empirical covariance, and map back
   to the full `J` by `K` matrix. The full covariance is singular in the
   simplex case because of both row-sum and accounting identities.

9. **Add invariant-based tests.**
   Test, to numerical tolerance,

   ```
   Q*_i = sum_j X_ij b*_j(Z_i),
   P_n[W(Y_i-Q*_i)] = 0,
   p_hat' B_hat = q_hat,
   A D_c(O_i) = 0.
   ```

   Repeat these tests for simplex, unconstrained, two-sided, lower-bounded,
   and upper-bounded outcomes, including subsets and nonuniform totals.
