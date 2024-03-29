# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

dAMH <- function(u, rho, logf) {
    .Call(`_clusterid_dAMH`, u, rho, logf)
}

vdAMH <- function(u1, u2, rho, logf) {
    .Call(`_clusterid_vdAMH`, u1, u2, rho, logf)
}

dClayton <- function(u, rho, logf) {
    .Call(`_clusterid_dClayton`, u, rho, logf)
}

dFrank <- function(u, rho, logf) {
    .Call(`_clusterid_dFrank`, u, rho, logf)
}

vdFrank <- function(u1, u2, rho, logf) {
    .Call(`_clusterid_vdFrank`, u1, u2, rho, logf)
}

dcopf <- function(u, rho, logf, copula) {
    .Call(`_clusterid_dcopf`, u, rho, logf, copula)
}

vdcopf <- function(u1, u2, rho, logf, copula) {
    .Call(`_clusterid_vdcopf`, u1, u2, rho, logf, copula)
}

dccopf <- function(u, u0, rho, logf, copula) {
    .Call(`_clusterid_dccopf`, u, u0, rho, logf, copula)
}

vdccopf <- function(u1, u2, u0, rho, logf, copula) {
    .Call(`_clusterid_vdccopf`, u1, u2, u0, rho, logf, copula)
}

pClayton <- function(u, rho) {
    .Call(`_clusterid_pClayton`, u, rho)
}

pFrank <- function(u, rho) {
    .Call(`_clusterid_pFrank`, u, rho)
}

pAMH <- function(u, rho) {
    .Call(`_clusterid_pAMH`, u, rho)
}

pccopf <- function(u, u0, rho, copula) {
    .Call(`_clusterid_pccopf`, u, u0, rho, copula)
}

pcopf <- function(u, rho, copula) {
    .Call(`_clusterid_pcopf`, u, rho, copula)
}

hf <- function(u, del, rho, copula) {
    .Call(`_clusterid_hf`, u, del, rho, copula)
}

hcf <- function(u, del, u0, rho, copula) {
    .Call(`_clusterid_hcf`, u, del, u0, rho, copula)
}

hf_Clayton <- function(u, del, rho) {
    .Call(`_clusterid_hf_Clayton`, u, del, rho)
}

hf_Frank <- function(u, del, theta) {
    .Call(`_clusterid_hf_Frank`, u, del, theta)
}

h1 <- function(u1, u2, rho, copula) {
    .Call(`_clusterid_h1`, u1, u2, rho, copula)
}

hc1 <- function(u1, u2, u0, rho, copula) {
    .Call(`_clusterid_hc1`, u1, u2, u0, rho, copula)
}

vh1 <- function(u1, u2, rho, copula) {
    .Call(`_clusterid_vh1`, u1, u2, rho, copula)
}

vhc1 <- function(u1, u2, u0, rho, copula) {
    .Call(`_clusterid_vhc1`, u1, u2, u0, rho, copula)
}

hpc <- function(x, levels, cuts, logf) {
    .Call(`_clusterid_hpc`, x, levels, cuts, logf)
}

Hpc <- function(x, levels, cuts, logf) {
    .Call(`_clusterid_Hpc`, x, levels, cuts, logf)
}

ppc <- function(q, levels, cuts, lower, logf) {
    .Call(`_clusterid_ppc`, q, levels, cuts, lower, logf)
}

vppc <- function(q, levels, cuts, lower, logf) {
    .Call(`_clusterid_vppc`, q, levels, cuts, lower, logf)
}

dpc <- function(x, levels, cuts, logf) {
    .Call(`_clusterid_dpc`, x, levels, cuts, logf)
}

vdpc <- function(x, levels, cuts, logf) {
    .Call(`_clusterid_vdpc`, x, levels, cuts, logf)
}

order_cpp <- function(x) {
    .Call(`_clusterid_order_cpp`, x)
}

pG0 <- function(r_id, G, p) {
    .Call(`_clusterid_pG0`, r_id, G, p)
}

pG <- function(r_id, G, p) {
    .Call(`_clusterid_pG`, r_id, G, p)
}

loglikFD2_pch <- function(par, theta, Y_F, X_F, Y_proband, X_proband, Age, Cal, cut_F, lam03, fgau, combn, copula) {
    .Call(`_clusterid_loglikFD2_pch`, par, theta, Y_F, X_F, Y_proband, X_proband, Age, Cal, cut_F, lam03, fgau, combn, copula)
}

loglikFD2_pch_gene <- function(par, theta, Y_F, X_F, Y_proband, X_proband, Age, Cal, cut_F, lam03, fgau, combn, copula) {
    .Call(`_clusterid_loglikFD2_pch_gene`, par, theta, Y_F, X_F, Y_proband, X_proband, Age, Cal, cut_F, lam03, fgau, combn, copula)
}

loglikR_pch <- function(par, cut_F, Y_R, X_R, LAM03R, LAM12R, cutR, fgau) {
    .Call(`_clusterid_loglikR_pch`, par, cut_F, Y_R, X_R, LAM03R, LAM12R, cutR, fgau)
}

loglikR_pch_gene <- function(par, cut_F, Y_R, X_R, LAM03R, LAM12R, cutR, fgau) {
    .Call(`_clusterid_loglikR_pch_gene`, par, cut_F, Y_R, X_R, LAM03R, LAM12R, cutR, fgau)
}

loglikS_pch <- function(par, cut_F, Y_S, LAM03S, LAM12S, cutS, fgau) {
    .Call(`_clusterid_loglikS_pch`, par, cut_F, Y_S, LAM03S, LAM12S, cutS, fgau)
}

loglikS_pch_gene <- function(par, cut_F, Y_S, LAM03S, LAM12S, cutS, fgau) {
    .Call(`_clusterid_loglikS_pch_gene`, par, cut_F, Y_S, LAM03S, LAM12S, cutS, fgau)
}

