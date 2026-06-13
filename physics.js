// Physics utilities, Nilsson model calculations, Strutinsky shell correction, and LDM
const fCache = [1, 1];

// Proton & Neutron Nilsson model parameter sets (kappa, mu) indexed by principal quantum number N (0 to 7)
export const PARAMETER_SETS = {
    universal: {
        name: "Lund Universal (Default)",
        proton_kappa: [0.05, 0.05, 0.05, 0.05, 0.05, 0.0637, 0.0637, 0.06],
        proton_mu:    [0.00, 0.00, 0.00, 0.35, 0.625, 0.600, 0.600, 0.54],
        neutron_kappa: [0.05, 0.05, 0.05, 0.05, 0.05, 0.0637, 0.0637, 0.06],
        neutron_mu:    [0.00, 0.00, 0.00, 0.25, 0.450, 0.450, 0.450, 0.40]
    },
    rare_earth: {
        name: "Rare Earth (A ≈ 150 - 180)",
        proton_kappa: [0.05, 0.05, 0.05, 0.06, 0.0637, 0.0637, 0.0637, 0.06],
        proton_mu:    [0.00, 0.00, 0.00, 0.35, 0.600, 0.600, 0.600, 0.54],
        neutron_kappa: [0.05, 0.05, 0.05, 0.06, 0.0637, 0.0637, 0.0637, 0.06],
        neutron_mu:    [0.00, 0.00, 0.00, 0.25, 0.390, 0.420, 0.440, 0.35]
    },
    actinide: {
        name: "Actinides (A ≈ 250)",
        proton_kappa: [0.05, 0.05, 0.05, 0.05, 0.0577, 0.0577, 0.0577, 0.0577],
        proton_mu:    [0.00, 0.00, 0.00, 0.35, 0.650, 0.650, 0.650, 0.650],
        neutron_kappa: [0.05, 0.05, 0.05, 0.05, 0.0635, 0.0635, 0.0635, 0.062],
        neutron_mu:    [0.00, 0.00, 0.00, 0.25, 0.390, 0.400, 0.400, 0.30]
    }
};

export function fact(n) { if (n < 0) return 0; while (fCache.length <= n) fCache.push(fCache[fCache.length - 1] * fCache.length); return fCache[n]; }
function cgDelta(j1, j2, j3) { return Math.sqrt((fact(j1 + j2 - j3) * fact(j1 - j2 + j3) * fact(-j1 + j2 + j3)) / fact(j1 + j2 + j3 + 1)); }
function CG(j1, j2, j3, m1, m2, m3) {
    if (m1 + m2 !== m3) return 0;
    if (Math.abs(m1) > j1 || Math.abs(m2) > j2 || Math.abs(m3) > j3) return 0;
    if (j3 > j1 + j2 || j3 < Math.abs(j1 - j2)) return 0;
    let term1 = Math.sqrt(2 * j3 + 1) * cgDelta(j1, j2, j3) * Math.sqrt(fact(j1 + m1) * fact(j1 - m1) * fact(j2 + m2) * fact(j2 - m2) * fact(j3 + m3) * fact(j3 - m3));
    let sum = 0, kmin = Math.max(0, j2 - j3 - m1, j1 - j3 + m2), kmax = Math.min(j1 + j2 - j3, j1 - m1, j2 + m2);
    for (let k = kmin; k <= kmax; k++) sum += ((k % 2 === 0) ? 1 : -1) / (fact(k) * fact(j1 + j2 - j3 - k) * fact(j1 - m1 - k) * fact(j2 + m2 - k) * fact(j3 - j2 + m1 + k) * fact(j3 - j1 - m2 + k));
    return term1 * sum;
}
export function get_cg(j1_d, j2_d, j3_d, m1_d, m2_d) { return CG(j1_d / 2, j2_d / 2, j3_d / 2, m1_d / 2, m2_d / 2, (m1_d + m2_d) / 2); }

export function jacobiEigvals(A, maxIter = 100) {
    let n = A.length, D = A.map(row => [...row]);
    for (let it = 0; it < maxIter; it++) {
        let max = 0, p = 0, q = 1;
        for (let i = 0; i < n - 1; i++) for (let j = i + 1; j < n; j++) if (Math.abs(D[i][j]) > max) { max = Math.abs(D[i][j]); p = i; q = j; }
        if (max < 1e-10) break;
        let theta = (D[q][q] - D[p][p]) / (2 * D[p][q]);
        let t = (theta === 0) ? 1 : Math.sign(theta) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        let c = 1 / Math.sqrt(t * t + 1), s = c * t;
        for (let i = 0; i < n; i++) {
            if (i !== p && i !== q) {
                let D_ip = D[i][p], D_iq = D[i][q];
                D[i][p] = D[p][i] = c * D_ip - s * D_iq;
                D[i][q] = D[q][i] = s * D_ip + c * D_iq;
            }
        }
        let D_pp = D[p][p], D_qq = D[q][q], D_pq = D[p][q];
        D[p][p] = c * c * D_pp - 2 * s * c * D_pq + s * s * D_qq;
        D[q][q] = s * s * D_pp + 2 * s * c * D_pq + c * c * D_qq;
        D[p][q] = D[q][p] = 0;
    }
    let evals = []; for (let i = 0; i < n; i++) evals.push(D[i][i]); return evals.sort((a, b) => a - b);
}

export function generateAsymptoticLabels(N, Omega) {
    let valid_states = [];
    for (let nz = N; nz >= 0; nz--) {
        let n_perp = N - nz;
        for (let sigma of [0.5, -0.5]) {
            let Lambda = Omega - sigma;
            if (Math.abs(Lambda) <= n_perp && (n_perp - Math.floor(Math.abs(Lambda))) % 2 === 0) valid_states.push([nz, Math.round(Lambda)]);
        }
    }
    return valid_states.sort((a, b) => b[0] - a[0]).map(s => `${Math.round(Omega * 2)}/2[${N}${s[0]}${s[1]}]`);
}
export function nilssonEnergies(beta, gamma_deg, A, nucleonType = 'proton', paramSet = 'universal') {
    const pSet = PARAMETER_SETS[paramSet] || PARAMETER_SETS.universal;
    const rkappa = (nucleonType === 'proton') ? pSet.proton_kappa : pSet.neutron_kappa;
    const rmu = (nucleonType === 'proton') ? pSet.proton_mu : pSet.neutron_mu;
    let hw0_base = 41.0 * Math.pow(A, -1.0 / 3.0);
    
    // Spheroidal volume conservation factor (approximate)
    let fdel = Math.pow(Math.pow(1.0 + (2.0 / 3.0) * beta * Math.cos(gamma_deg * Math.PI / 180.0), 2.0) * (1.0 - (4.0 / 3.0) * beta * Math.cos(gamma_deg * Math.PI / 180.0)), -1.0 / 6.0);
    let hw0_mev = hw0_base * fdel;
    let orbitals = [];
    
    const gamma_rad = gamma_deg * Math.PI / 180.0;
    const cosG = Math.cos(gamma_rad);
    const sinG = Math.sin(gamma_rad);

    // Solve for major shells N = 0 to 9 (n = 0 to 18)
    for (let n = 0; n <= 18; n += 2) {
        let N_idx = Math.floor(n / 2);
        
        // Safe lookup for higher shells parameters
        const kappa_val = rkappa[Math.min(N_idx, rkappa.length - 1)];
        const mu_val = rmu[Math.min(N_idx, rmu.length - 1)];
        let c = -2.0 * hw0_base * kappa_val, d = mu_val * c / 2.0;

        // Build Block A basis states (separating degenerate Kramer's pairs via projection parity)
        let basis = [];
        for (let l = n; l >= 0; l -= 4) {
            for (let lam = -l; lam <= l + 1; lam += 2) {
                for (let isig of [1, -1]) {
                    let iom = lam + isig;
                    // Condition for Block A: ((Omega - 1/2) / 2) is an even integer
                    if (Math.abs(((iom - 1) / 2) % 2) === 0) {
                        basis.push([l, lam, isig]);
                    }
                }
            }
        }
        
        let nbas = basis.length;
        if (nbas === 0) continue;

        let H = Array(nbas).fill(0).map(() => Array(nbas).fill(0));
        for (let i = 0; i < nbas; i++) {
            let [l, lam, isig] = basis[i];
            for (let j = i; j < nbas; j++) {
                let [l1, lam1, isig1] = basis[j], h00 = 0.0, hl2 = 0.0, hls = 0.0, hr2 = 0.0, hy20 = 0.0, hy22 = 0.0;
                
                if (i === j) { 
                    h00 = (n + 3) / 2.0 * hw0_mev; 
                    hl2 = l * (l + 2) / 4.0; 
                }
                
                // Spin-orbit coupling (l.s)
                if (l1 === l) {
                    if (lam1 === lam && isig1 === isig) {
                        hls = (lam * isig) / 4.0;
                    } else if (lam1 === lam + 2 && isig1 === isig - 2) {
                        hls = Math.sqrt((l - lam) * (l + lam + 2)) / 4.0;
                    } else if (lam1 === lam - 2 && isig1 === isig + 2) {
                        hls = Math.sqrt((l + lam) * (l - lam + 2)) / 4.0;
                    }
                }
                
                // Radial matrix elements of r^2
                if (l1 === l) {
                    hr2 = (n + 3) / 2.0;
                } else if (l1 === l - 4) {
                    hr2 = Math.sqrt((n - l + 4) * (n + l + 2)) / 2.0;
                } else if (l1 === l + 4) {
                    hr2 = Math.sqrt((n - l) * (n + l + 6)) / 2.0;
                }

                // Axial Quadrupole coupling (Y20)
                if (Math.abs(hr2) > 1e-5 && lam1 === lam && isig1 === isig) {
                    hy20 = Math.sqrt((l + 1) / (l1 + 1)) * get_cg(l, 4, l1, lam, 0) * get_cg(l, 4, l1, 0, 0);
                }
                
                // Triaxial Quadrupole coupling (Y22 and Y2,-2)
                if (Math.abs(hr2) > 1e-5 && isig1 === isig) {
                    if (lam1 === lam + 4) {
                        hy22 = Math.sqrt((l + 1) / (l1 + 1)) * get_cg(l, 4, l1, lam, 4) * get_cg(l, 4, l1, 0, 0);
                    } else if (lam1 === lam - 4) {
                        hy22 = Math.sqrt((l + 1) / (l1 + 1)) * get_cg(l, 4, l1, lam, -4) * get_cg(l, 4, l1, 0, 0);
                    }
                }

                // Quadrupole interaction term: -beta * r^2 * (cos(gamma)*Y20 + sin(gamma)/sqrt(2)*(Y22 + Y2,-2))
                let v_quad = - beta * hw0_mev * (2.0 / 3.0) * hr2 * (cosG * hy20 + (sinG / Math.sqrt(2.0)) * hy22);
                H[i][j] = H[j][i] = h00 + v_quad + c * hls + d * hl2;
            }
        }
        
        let evals_mev = jacobiEigvals(H);
        for (let i = 0; i < evals_mev.length; i++) {
            orbitals.push([evals_mev[i], 0.5, `N=${N_idx}`]);
        }
    }
    orbitals.sort((a, b) => a[0] - b[0]); return { hw0_base, orbitals };
}

// error function approximation
export function erf(x) {
    const a1 =  0.254829592, a2 = -0.284496736, a3 =  1.421413741, a4 = -1.453152027, a5 =  1.061405429, p  =  0.3275911;
    const sign = x >= 0 ? 1 : -1;
    const absx = Math.abs(x);
    const t = 1.0 / (1.0 + p * absx);
    const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-absx * absx);
    return sign * y;
}

// hermite polynomials
export function hermite(x, n) {
    if (n === 0) return 1;
    if (n === 1) return 2*x;
    if (n === 2) return 4*x*x - 2;
    if (n === 3) return 8*x*x*x - 12*x;
    if (n === 4) return 16*Math.pow(x,4) - 48*x*x + 12;
    if (n === 5) return 32*Math.pow(x,5) - 160*x*x*x + 120*x;
    if (n === 6) return 64*Math.pow(x,6) - 480*Math.pow(x,4) + 720*x*x - 120;
    return 0;
}

// smooth level density integral (Strutinsky occupation sum)
export function smoothDensityIntegral(lambda, evals, gamma, pOrder = 6) {
    let sum = 0;
    for (let i = 0; i < evals.length; i++) {
        let u = (lambda - evals[i]) / gamma;
        let hsum = 0;
        for (let m = 2; m <= pOrder; m += 2) {
            let Cm = (m === 2) ? -0.25 : (m === 4) ? 1.0/32.0 : -1.0/384.0;
            hsum += Cm * hermite(u, m - 1);
        }
        let term = 0.5 * (1.0 + erf(u)) - Math.exp(-u*u) * hsum / Math.sqrt(Math.PI);
        sum += term;
    }
    return sum;
}

// Strutinsky smooth Fermi level search
export function findSmoothLambda(evals, numParticles, gamma, pOrder = 6) {
    let low = evals[0] - 2 * gamma;
    let high = evals[evals.length - 1] + 2 * gamma;
    let iter = 0;
    let lambda = 0.5 * (low + high);
    while (high - low > 1e-6 && iter < 100) {
        lambda = 0.5 * (low + high);
        let val = smoothDensityIntegral(lambda, evals, gamma, pOrder);
        if (val > numParticles) {
            high = lambda;
        } else {
            low = lambda;
        }
        iter++;
    }
    return lambda;
}

// Strutinsky smooth energy calculation
export function calculateSmoothEnergy(lambda, evals, gamma, pOrder = 6) {
    let sum = 0;
    let Cp = (pOrder === 6) ? -1.0 / 384.0 : (pOrder === 4) ? 1.0 / 32.0 : -0.25;
    for (let i = 0; i < evals.length; i++) {
        let u = (lambda - evals[i]) / gamma;
        let hsum = 0;
        for (let m = 2; m <= pOrder; m += 2) {
            let Cm = (m === 2) ? -0.25 : (m === 4) ? 1.0/32.0 : -1.0/384.0;
            hsum += Cm * hermite(u, m - 1);
        }
        let occupation = 0.5 * (1.0 + erf(u)) - Math.exp(-u*u) * hsum / Math.sqrt(Math.PI);
        let termE = evals[i] * occupation - gamma * Math.exp(-u*u) * Cp * hermite(u, pOrder) / (2.0 * Math.sqrt(Math.PI));
        sum += termE;
    }
    return sum;
}

// Complete Strutinsky Shell Correction
export function calculateShellCorrection(orbitals, numParticles, gamma, pOrder = 6) {
    // 1. Duplicate orbitals to represent degenerate Kramers pairs (each holds 2 particles)
    let evals = [];
    for (let i = 0; i < orbitals.length; i++) {
        evals.push(orbitals[i][0]);
        evals.push(orbitals[i][0]);
    }
    evals.sort((a, b) => a - b);

    // 2. Find lambda
    let lambda = findSmoothLambda(evals, numParticles, gamma, pOrder);

    // 3. Sum microscopic energies
    let E_shell = 0;
    for (let i = 0; i < Math.min(numParticles, evals.length); i++) {
        E_shell += evals[i];
    }

    // 4. Calculate smooth energy
    let E_tilde = calculateSmoothEnergy(lambda, evals, gamma, pOrder);

    // 5. Delta
    let deltaE = E_shell - E_tilde;

    return {
        E_shell,
        E_tilde,
        deltaE,
        lambda
    };
}

// Liquid Drop Model surface & Coulomb energy
export function calculateLDM(beta, gamma_deg, Z, N) {
    let A = Z + N;
    // Standard LDM coefficients (in MeV)
    const a_surf = 17.94;
    const a_coul = 0.705;

    let gamma_rad = gamma_deg * Math.PI / 180.0;
    let cos3G = Math.cos(3.0 * gamma_rad);

    // Spheroidal shape functions with triaxial gamma parameter (3rd order expansion)
    let B_surf = 1.0 + (2.0 / 5.0) * beta * beta - (4.0 / 105.0) * beta * beta * beta * cos3G;
    let B_coul = 1.0 - (1.0 / 5.0) * beta * beta - (4.0 / 105.0) * beta * beta * beta * cos3G;

    let E_surf_0 = a_surf * Math.pow(A, 2.0 / 3.0);
    let E_coul_0 = a_coul * Z * (Z - 1) / Math.pow(A, 1.0 / 3.0);

    let E_surf = E_surf_0 * B_surf;
    let E_coul = E_coul_0 * B_coul;

    let E_ldm = E_surf + E_coul;
    let E_ldm_0 = E_surf_0 + E_coul_0;
    
    // Return relative energy compared to spherical LDM to show deformation effect
    let E_ldm_rel = E_ldm - E_ldm_0;

    return {
        E_surf,
        E_coul,
        E_ldm,
        E_ldm_rel
    };
}

// Preset Isotopes for educational demonstration
export const PRESETS = [
    { name: "Sm152 (Sm, Z=62, N=90)", Z: 62, N: 90, type: "deformed" },
    { name: "Sn120 (Sn, Z=50, N=70)", Z: 50, N: 70, type: "spherical" },
    { name: "Pb208 (Pb, Z=82, N=126)", Z: 82, N: 126, type: "double-magic" },
    { name: "Nd150 (Nd, Z=60, N=90)", Z: 60, N: 90, type: "deformed" }
];
