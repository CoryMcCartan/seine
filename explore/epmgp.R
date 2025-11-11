# rank deficient
n = 3
mu = c(0.95, 0.8, 0.7)
# Sigma = matrix(0.01*c(2, -1, -1, 2), 2)
Sigma = 0.09 * (tcrossprod(c(2, 1, -2)) + tcrossprod(c(0, -1, 1)))
# L = chol_pivot(Sigma)[c(3, 1, 2), ]
L = t(chol(Sigma + diag(n)*1e-12))[, 1:(n-1), drop=FALSE]

N = 1e5
# B = R_ess_tmvn(N + 10, mu, L, init=mu + 0.3*L[, 1])[10 + 1:N, ]
B = R_ess_tmvn(N + 10, mu, L, init=mu)[10 + 1:N, ]
Bu = t(mu + L %*% matrix(rnorm((n-1)*N), ncol=N)) # untrunc

plot(B, xlim=c(-0.8, 1.0), ylim=0:1, pch=16, cex=1, col="#00000020"); abline(h=0:1, v=0:1, lty="dashed")
plot(Bu, xlim=c(-0.8, 1.0), ylim=0:1, pch=16, cex=1, col="#00000020"); abline(h=0:1, v=0:1, lty="dashed")
mean(rowSums(Bu < 0 | Bu > 1) == 0)

x = lm.fit(B, rep(1, N)) |> coef() |> prop.table() # race %
y = sum(x * B[1, ])
S = diag(n)
S[n, ] = x
Sinv = solve(S)
C = t(Sinv[, -n, drop=FALSE])

L0 = L
L = (S %*% L0)[-n, , drop=FALSE]
Kinv = solve(L %*% t(L))

tau = numeric(n)
nu = numeric(n)
cav_loc = numeric(n)
cav_var = numeric(n)

m1 = mu[-n]
prev_m1 = m1
m1_fac = numeric(n)
m2 = L %*% t(L)
log_Z = 0.0

for (i in 1:50) {
    m2_diag = diag(t(C) %*% m2 %*% C)
    cav_var = 1 / (1 / m2_diag - tau)
    cav_loc = cav_var * (c(t(C) %*% m1) / m2_diag - nu)
    log_Z = 0.0

    for (j in 1:n) {
        if (j < n) {
            ms = R_utn_moments(cav_loc[j], cav_var[j])
        } else {
            ms = R_utn_moments(cav_loc[j] + y/x[n], cav_var[j])
            ms[2] = ms[2] - y/x[n]
        }
        log_Z = log_Z + ms[1]
        tau[j] = 1 / ms[3] - 1 / cav_var[j]
        nu[j] = ms[2] / ms[3] - cav_loc[j] / cav_var[j]
    }

    # m2 = solve(Kinv + diag(tau))
    m2 = solve(Kinv + C %*% diag(tau) %*% t(C))
    m1 = c(m2 %*% (Kinv %*% mu[-n] + C %*% nu))

    if (max(abs(m1 - prev_m1)) < 1e-8) break
    prev_m1 = m1
}

log_Z = log_Z + 0.5 * (
    # log(det(m2)) - 2*sum(log(diag(L))) +
    log(det(m2)) + log(det(Kinv)) +
        (t(m1) %*% solve(m2) %*% m1) - (t(mu[-n]) %*% Kinv %*% mu[-n]) +
        sum(log1p(tau*cav_var)) +
        sum((cav_loc^2*tau - 2*cav_loc*nu - nu^2*cav_var) / (1 + tau*cav_var))
)

m_cpp = R_ep_moments(mu[-n], L, C[, n], y/x[n], tol=1e-7)

c(log_Z)
m_cpp[[1]]
log(mean(rowSums(Bu < 0 | Bu > 1) == 0))

c(m1)
m_cpp[[2]]
colMeans(B)

m2
# m_cpp[[3]] %*% t(m_cpp[[3]])
m_cpp[[3]]
cov(B)
