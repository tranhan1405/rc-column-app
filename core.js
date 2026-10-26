/* RC column seismic model: computations (units N, mm, MPa unless stated; forces reported in kN).
   Model coefficients come from model.js (window.APP_MODEL), written by validation/app_export.py. */
(function (root) {
  "use strict";
  const LS_KEYS = ["LS1", "LS80"];

  function clip(x, lo, hi) { return Math.min(Math.max(x, lo), hi); }

  // standard normal CDF via erfc (Numerical Recipes erfcc, fractional error < 1.2e-7)
  function erfc(x) {
    const z = Math.abs(x), t = 1 / (1 + 0.5 * z);
    const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
      t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 +
      t * (-0.82215223 + t * 0.17087277)))))))));
    return x >= 0 ? r : 2 - r;
  }
  function phi(x) { return 0.5 * erfc(-x / Math.SQRT2); }

  // flexural reference shear V_p and ASCE 41-17 Eq. (10-3) shear strength ratio (k_nl = 1, lambda = 1), as in the study
  function shearRatio(c) {
    const n = c.nu * c.fc * c.b * c.h;
    const a = Math.min(n / (0.85 * c.fc * c.b), c.h);
    const vp = (n * (c.h - a) / 2 + 0.4 * c.rho_l * c.b * c.h * c.h * c.fyl) / c.ls;
    const d = 0.8 * c.h;
    const mvd = clip(c.ls / d, 2, 4);
    const av = c.rho_sh * c.b * c.s;
    const vs = av * c.fyt * d / c.s;
    const vc = 0.5 * Math.sqrt(c.fc) / mvd * Math.sqrt(1 + n / (0.5 * Math.sqrt(c.fc) * c.b * c.h)) * 0.8 * c.b * c.h;
    return { vp: vp / 1e3, vo: (vs + vc) / 1e3, vp_vo: vp / (vs + vc) };
  }

  function features(c, vp_vo) {
    return { ln_nu: Math.log(c.nu + 0.05), ln_fc: Math.log(c.fc), ln_rho_l: Math.log(c.rho_l),
      ln_rho_sh: Math.log(c.rho_sh), ln_ls_h: Math.log(c.ls / c.h), ln_fyl: Math.log(c.fyl),
      ln_fyt: Math.log(c.fyt), ln_vp_vo: Math.log(vp_vo) };
  }

  function ll(m, x, feats) {
    let y = m.intercept;
    for (const f of feats) y += m.coef[f] * (Number.isFinite(x[f]) ? x[f] : m.impute[f]);
    return y;
  }

  function ec(fc) { return 4700 * Math.sqrt(fc); }

  // ASCE/SEI 41-17 Table 10-8 (columns not controlled by inadequate development or splicing) and Table 10-5
  function asce41(c, vp, vp_vo) {
    const nu = c.nu, rt = Math.min(c.rho_sh, 0.0175), ratio = Math.max(vp_vo, 0.2);
    const a = Math.max(0.042 - 0.043 * nu + 0.63 * rt - 0.023 * ratio, 0);
    const bOf = (n) => 0.5 / (5 + n / 0.8 * (1 / rt) * c.fc / c.fyt) - 0.01;
    const nb = Math.max(nu, 0.1);
    let b = nb <= 0.5 ? bOf(nb) : bOf(0.5) * Math.max(0.7 - nb, 0) / 0.2;
    b = Math.max(b, a);
    const cc = Math.max(0.24 - 0.4 * nu, 0);
    const k = clip(0.3 + 0.4 * (nu - 0.1) / 0.4, 0.3, 0.7);
    const thy = vp * 1e3 * c.ls * c.ls / (3 * k * ec(c.fc) * c.b * Math.pow(c.h, 3) / 12);
    const valid = c.rho_sh >= 0.0005;
    return { theta_y: thy, a: a, b: b, c: cc, theta_C: valid ? thy + a : NaN, theta_E: valid ? thy + b : NaN, valid: valid };
  }

  // EN 1998-3 Annex A, Eq. (A.1) with gamma_el = 1, symmetric reinforcement of equal bars, no diagonal bars;
  // omega includes the web (intermediate side-face) bars; Eq. (A.2) for alpha
  function ec8(c) {
    const bo = c.b - 2 * c.cover - c.d_hoop, ho = c.h - 2 * c.cover - c.d_hoop;
    const gb = c.n_perp + 1, gh = c.n_par + 1;
    const sb2 = 2 * bo * bo / gb + 2 * ho * ho / gh;
    const alpha = Math.max((1 - c.s / (2 * bo)) * (1 - c.s / (2 * ho)) * (1 - sb2 / (6 * bo * ho)), 0);
    const nFace = c.n_perp + 2, nWeb = Math.max(c.n_par, 0), nTot = 2 * nFace + 2 * nWeb;
    const wc = c.rho_l * nFace / nTot * c.fyl / c.fc, wt = c.rho_l * (nFace + 2 * nWeb) / nTot * c.fyl / c.fc;
    const th = 0.016 * Math.pow(0.3, c.nu) * Math.pow(Math.max(0.01, wc) / Math.max(0.01, wt) * c.fc, 0.225) *
      Math.pow(Math.min(9, c.ls / c.h), 0.35) * Math.pow(25, alpha * c.rho_sh * c.fyt / c.fc);
    return { theta_um: th, theta_um_primary: th / 1.5, alpha: alpha };
  }

  // Phase 16: shift of the log-median by k tests on columns of the same programme (shrinkage with the between- and
  // within-programme variances), scaled bounds and reduction factors phi_q for q = 5 and 10 %
  function anchored(mod, phiT, median, tests) {
    const v = (tests || []).filter((x) => Number.isFinite(x) && x > 0).slice(0, 2);
    const k = v.length, t2 = mod.tau2, s2 = mod.sigma2;
    const u = k ? v.reduce((a, x) => a + Math.log(x / median), 0) * t2 / (s2 + k * t2) : 0;
    const sk = Math.sqrt(s2 + t2 * s2 / (s2 + k * t2)) / Math.sqrt(s2 + t2) * mod.sigma_ln;
    const m = median * Math.exp(u);
    return { median: m, lo: m * Math.exp(-sk), hi: m * Math.exp(sk), sigma: sk, k: k,
      d05: phiT["0.05"][k] * m, d10: phiT["0.10"][k] * m };
  }

  function shearClass(vp_vo) { return vp_vo <= 0.6 ? "Flexure" : vp_vo <= 1.0 ? "Flexure-Shear" : "Shear"; }

  function backbone(ty, tp, t80, v, ls) {
    const t1 = 0.75 * ty;
    const p = Math.max(tp, 1.2 * t1);
    const q = Math.max(t80, 1.1 * p);
    const t4 = q + (q - p);
    return [[t1 * ls, 0.75 * v], [p * ls, v], [q * ls, 0.8 * v], [t4 * ls, 0.6 * v]];
  }

  // linear interpolation in T on the model grid; outside the grid the nearest period is used
  function atPeriod(M, T, fn) {
    const P = M.periods;
    const t = clip(T, P[0], P[P.length - 1]);
    let i = 0;
    while (i < P.length - 2 && t > P[i + 1]) i++;
    const w = (t - P[i]) / (P[i + 1] - P[i]);
    return (1 - w) * fn(P[i]) + w * fn(P[i + 1]);
  }

  // tests: optional {theta_y: [..], theta_p: [..], theta_80: [..], v_peak: [..]} measured on up to two columns of the
  // same design from the same test programme (drifts as ratios, force in kN)
  function compute(c, M, demand, tests) {
    const sr = shearRatio(c);
    const x = features(c, sr.vp_vo);
    const F = M.features, mod = M.models;
    const T = tests || {};
    const out = { vp: sr.vp, vo: sr.vo, vp_vo: sr.vp_vo, cls: shearClass(sr.vp_vo) };
    const med = (k) => Math.exp(ll(mod[k], x, F));
    out.strength = anchored(mod.strength, M.phi.strength, sr.vp * med("strength"), T.v_peak);
    out.v_peak = out.strength.median;
    out.v_sigma = out.strength.sigma;
    out.deform = {};
    for (const k of ["theta_y", "theta_p", "theta_80"]) out.deform[k] = anchored(mod[k], M.phi[k], med(k), T[k]);
    out.backbone = backbone(out.deform.theta_y.median, out.deform.theta_p.median, out.deform.theta_80.median,
      out.v_peak, c.ls);
    out.cyclic = M.cyclic_class_medians[out.cls];
    out.asce41 = asce41(c, sr.vp, sr.vp_vo);
    out.ec8 = ec8(c);
    if (demand) {
      out.fragility = {};
      for (const ls of LS_KEYS) {
        const lnTheta = atPeriod(M, demand.T, (T) => ll(mod["R_" + ls + "_T" + T], x, F));
        const sR = atPeriod(M, demand.T, (T) => mod["R_" + ls + "_T" + T].sigma_ln);
        const bR = atPeriod(M, demand.T, (T) => M.beta_rtr[ls + "_T" + T]);
        const beta = Math.sqrt(bR * bR + sR * sR + out.v_sigma * out.v_sigma);
        const thetaR = Math.exp(lnTheta);
        const sa50 = thetaR * out.v_peak / (demand.mass_t * 9.80665);       // g
        const R = demand.mass_t * demand.sa_g * 9.80665 / out.v_peak;
        out.fragility[ls] = { thetaR: thetaR, beta: beta, beta_rtr: bR, sigma_pred: sR, sa50: sa50, R: R,
          p: phi(Math.log(R / thetaR) / beta) };
      }
    }
    return out;
  }

  function warnings(c, r, M, demand) {
    const w = [];
    const rg = M.ranges;
    const vals = { b: c.b, h: c.h, ls: c.ls, ls_h: c.ls / c.h, fc: c.fc, fyl: c.fyl, fyt: c.fyt, rho_l: c.rho_l,
      rho_sh: c.rho_sh, nu: c.nu, s: c.s, vp_vo: r.vp_vo, cover: c.cover };
    for (const k in vals) {
      if (rg[k] && (vals[k] < rg[k][0] || vals[k] > rg[k][1])) w.push({ level: "range", key: k, value: vals[k], range: rg[k] });
    }
    if (c.nu > 0.7) w.push({ level: "limit", key: "nu", text: "Axial load ratio above 0.7." });
    if (c.rho_sh < 0.0005) w.push({ level: "limit", key: "rho_sh", text: "Transverse ratio below 0.0005: ASCE 41-17 Table 10-8 does not apply." });
    if (c.ls / c.h < 1.5) w.push({ level: "limit", key: "ls_h", text: "Shear-span ratio below 1.5." });
    if (demand && (demand.T < M.periods[0] || demand.T > M.periods[M.periods.length - 1]))
      w.push({ level: "limit", key: "T", text: "Period outside 0.3-1.0 s: fragility taken at the nearest analysed period." });
    const hi80 = r.deform.theta_80.hi;
    if (Number.isFinite(r.asce41.theta_C) && r.asce41.theta_C > hi80)
      w.push({ level: "code", key: "asce41", text: "ASCE 41-17 drift at point C exceeds the 84 % bound of the predicted drift at 20 % strength loss." });
    if (r.ec8.theta_um > hi80)
      w.push({ level: "code", key: "ec8", text: "EN 1998-3 ultimate chord rotation exceeds the 84 % bound of the predicted drift at 20 % strength loss." });
    return w;
  }

  function pinching4Args(tag, bb, cyc) {
    const env = [];
    for (const p of bb) env.push(p[1], p[0]);
    const neg = [];
    for (const p of bb) neg.push(-p[1], -p[0]);
    return [tag, ...env, ...neg, cyc.rDisp, cyc.rForce, cyc.uForce, cyc.rDisp, cyc.rForce, cyc.uForce,
      cyc.gK, 0, 1, 0, 0.95, cyc.gD, 0, 1, 0, 0.95, 0, 0, 1, 0, 0.95, 10];
  }

  root.ColumnModel = { compute, warnings, pinching4Args, phi, shearRatio, asce41, ec8, backbone, anchored };
})(typeof window !== "undefined" ? window : globalThis);
