x = with(elec_1968, cbind(vap_white, vap_black, vap_other))
y = elec_1968$pres_rep_nix
# z = model.matrix(~ 0 + (state + pop_city + pop_urban + farm + educ_elem +
#                 educ_hsch + inc_00_03k + inc_03_08k + inc_08_25k), elec_1968)
# z = scale(z)
z = with(elec_1968, bases::b_bart(
    model.matrix(~ 0 + state), pop_city, pop_urban, pop_rural, farm, educ_elem,
    educ_hsch, educ_coll, inc_00_03k, inc_03_08k, inc_08_25k, inc_25_99k, log(pop), pres_turn,
    trees = 250
))
colnames(z) = paste0("z_", seq_len(ncol(z)))
w = sqrt(elec_1968$pres_total) / mean(sqrt(elec_1968$pres_total))

n = nrow(x)
p = ncol(z)
n_x = ncol(x)
m0 = ei_ridge_impl(x, y, z)

xz = row_kronecker(x, z, int_scale = 1e4)
udv = svd(xz)
dvec = crossprod(xz, y)
Dmat = crossprod(xz) + diag(ncol(xz)) * m0$penalty * 1#0.8
Amat = matrix(0, nrow = nrow(Dmat), ncol = n * n_x * 2)
for (j in seq_len(k)) {
    use = n_x + p*(j-1) + seq_len(p)
    idx = (j - 1) * n + seq_len(n)
    Amat[j, idx] = 1e4
    Amat[j, idx + n*k] = -1e4
    Amat[use, idx] = t(z)
    Amat[use, idx + n*k] = -t(z)
}
bvec = c(rep(0, n*k), rep(-1, n*k))

n_nz = sum(Amat[, 1] != 0)
Amat2 = apply(Amat, 2, \(x) x[x != 0])
Aind = matrix(0, nrow = n_nz + 1, ncol = n * n_x * 2)
Aind[1, ] = n_nz
Aind[1 + seq_len(n_nz), ] = matrix(which(Amat != 0, arr.ind=TRUE)[, 1], nrow=n_nz)

# Dmat2 = t(scale_cols(udv$v, 1 / sqrt(udv$d^2 + m0$penalty)))
# Dmat2 = solve(chol(Dmat))
# bench::mark(
# res = quadprog::solve.QP(Dmat, dvec, Amat, bvec),
# res = quadprog::solve.QP.compact(Dmat, dvec, Amat2, Aind, bvec),
# chol = {
Dmat2 = solve(chol(Dmat))
# Dmat2 = backsolve(qr.R(qr(Dmat)), diag(nrow(Dmat)))
res = quadprog::solve.QP.compact(Dmat2, dvec, Amat2, Aind, bvec, factorized = TRUE)
# }
# )
all.equal(c(m0$coef), res$unconstrained.solution, check.attributes=FALSE)

matplot(cbind(m0$coef, res$solution), cex=0.5)
plot(fitted(m0), xz %*% res$solution, pch=1, cex=0.4*w,
    xlab="Unconstrained", ylab="Constrained", main=expr(hat(y)))
plot(xz %*% res$solution, y)
r2 = print(c(
    unconstrained = cor(fitted(m0), y)^2,
    constrained = cor(xz %*% res$solution, y)^2
))

preds0 = matrix(crossprod(Amat[, 1:(n*n_x)], m0$coef), ncol=n_x)
preds_qp = matrix(crossprod(Amat[, 1:(n*n_x)], res$solution), ncol=n_x)
plot(preds0[, 1], preds_qp[, 1], pch=1, cex=0.4*w,
    xlab="Unconstrained", ylab="Constrained", main="E[Nixon | White]"); abline(a=0, b=1, col="red")
plot(preds0[, 2], preds_qp[, 2], pch=1, cex=0.4*w,
    xlab="Unconstrained", ylab="Constrained", main="E[Nixon | Black]"); abline(a=0, b=1, col="red")
plot(preds0[, 3], preds_qp[, 3], pch=1, cex=0.4*w,
    xlab="Unconstrained", ylab="Constrained", main="E[Nixon | Other]"); abline(a=0, b=1, col="red")

est = cbind(
    unconstrained = colSums(preds0 * x * elec_1968$pres_total) / colSums(x * elec_1968$pres_total),
    constrained = colSums(preds_qp * x * elec_1968$pres_total) / colSums(x * elec_1968$pres_total)
)
print(est)

p_r = colSums(x * elec_1968$pres_total) / sum(elec_1968$pres_total)
p_r %*% est
weighted.mean(xz %*% res$solution, elec_1968$pres_total)
weighted.mean(fitted(m0), elec_1968$pres_total)
weighted.mean(y, elec_1968$pres_total)

