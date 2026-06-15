// Canvas and WebGL rendering engine using Three.js for the Strutinsky Shell Correction tool
import { nilssonEnergies, calculateShellCorrection, calculateLDM, thermalDampingFactor, fermiDiracOccupation } from './physics.js';

// Setup high-DPI 2D canvas to ensure crisp drawings
function getScaledContext(canvas) {
    const rect = canvas.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    
    canvas.width = rect.width * dpr;
    canvas.height = rect.height * dpr;
    
    const ctx = canvas.getContext('2d');
    ctx.resetTransform();
    ctx.scale(dpr, dpr);
    
    return { ctx, width: rect.width, height: rect.height };
}

// Global cache for Nilsson level calculations to keep dragging smooth
let cachedNilssonGrid = {
    proton: null,
    neutron: null,
    paramSet: '',
    gamma: -1,
    A: 0
};

// Calculate and cache the Nilsson levels grid for beta in [0, 0.6] at the current shape gamma
export function updateNilssonCache(A, paramSet, gamma) {
    const betas = [];
    const steps = 60;
    for (let i = 0; i <= steps; i++) {
        betas.push((i / steps) * 0.6);
    }
    
    const protonGrid = betas.map(b => ({
        beta: b,
        data: nilssonEnergies(b, gamma, A, 'proton', paramSet)
    }));
    
    const neutronGrid = betas.map(b => ({
        beta: b,
        data: nilssonEnergies(b, gamma, A, 'neutron', paramSet)
    }));
    
    cachedNilssonGrid = {
        proton: protonGrid,
        neutron: neutronGrid,
        paramSet,
        gamma,
        A,
        betas
    };
}

// Draw the Nilsson level splitting diagram
export function drawNilssonDiagram(canvas, state) {
    const { ctx, width, height } = getScaledContext(canvas);
    if (!ctx) return;
    
    const A = state.Z + state.N;
    
    // Regenerate cache if parameter set, nucleons, or active triaxiality changes
    if (cachedNilssonGrid.A !== A || 
        cachedNilssonGrid.paramSet !== state.paramSet || 
        cachedNilssonGrid.gamma !== state.gamma_shape) {
        updateNilssonCache(A, state.paramSet, state.gamma_shape);
    }
    
    const grid = state.nucleonType === 'proton' ? cachedNilssonGrid.proton : cachedNilssonGrid.neutron;
    const numParticles = state.nucleonType === 'proton' ? state.Z : state.N;
    const accentColor = state.nucleonType === 'proton' ? '#00f0ff' : '#10b981';
    
    ctx.clearRect(0, 0, width, height);
    
    const margin = { top: 30, right: 60, bottom: 40, left: 50 };
    const chartWidth = width - margin.left - margin.right;
    const chartHeight = height - margin.top - margin.bottom;
    
    // Find min and max energies in the grid to scale Y-axis
    let minE = Infinity, maxE = -Infinity;
    grid.forEach(point => {
        point.data.orbitals.forEach(o => {
            if (o[0] < minE) minE = o[0];
            if (o[0] > maxE) maxE = o[0];
        });
    });
    
    minE = Math.max(0, minE - 2.0);
    maxE = Math.min(100, maxE + 2.0);
    
    const getX = (b) => margin.left + (b / 0.6) * chartWidth;
    const getY = (e) => margin.top + chartHeight - ((e - minE) / (maxE - minE)) * chartHeight;
    
    // Draw Gridlines
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.05)';
    ctx.lineWidth = 1;
    
    // X Grid
    for (let b = 0.0; b <= 0.61; b += 0.1) {
        const x = getX(b);
        ctx.beginPath();
        ctx.moveTo(x, margin.top);
        ctx.lineTo(x, margin.top + chartHeight);
        ctx.stroke();
        
        ctx.fillStyle = '#6b7280';
        ctx.font = '10px Inter';
        ctx.textAlign = 'center';
        ctx.fillText(b.toFixed(1), x, margin.top + chartHeight + 15);
    }
    
    // Y Grid
    const roundedMinE = Math.ceil(minE / 5) * 5;
    for (let e = roundedMinE; e <= maxE; e += 5) {
        const y = getY(e);
        ctx.beginPath();
        ctx.moveTo(margin.left, y);
        ctx.lineTo(margin.left + chartWidth, y);
        ctx.stroke();
        
        ctx.fillStyle = '#6b7280';
        ctx.font = '10px Inter';
        ctx.textAlign = 'right';
        ctx.fillText(e.toFixed(0), margin.left - 10, y + 3);
    }
    
    ctx.fillStyle = '#9ca3af';
    ctx.font = '11px Outfit';
    ctx.textAlign = 'center';
    ctx.fillText('Quadrupole Deformation (β)', margin.left + chartWidth / 2, margin.top + chartHeight + 35);
    
    ctx.save();
    ctx.translate(15, margin.top + chartHeight / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText('Energy (MeV)', 0, 0);
    ctx.restore();
    
    // Draw Thermal Level Density Heatmap in the background (smeared quasi-continuum) if T > 0
    if (state.T > 0.01) {
        const gridX = 40;
        const gridY = 50;
        const tempCanvas = document.createElement('canvas');
        tempCanvas.width = gridX;
        tempCanvas.height = gridY;
        const tempCtx = tempCanvas.getContext('2d');
        const tempImg = tempCtx.createImageData(gridX, gridY);
        const tempImgData = tempImg.data;

        const T_eff = Math.max(0.05, state.T);

        for (let c = 0; c < gridX; c++) {
            const betaVal = (c / (gridX - 1)) * 0.6;
            const gridPointIdx = Math.min(grid.length - 1, Math.floor((betaVal / 0.6) * (grid.length - 1)));
            const pointOrbitals = grid[gridPointIdx].data.orbitals;

            for (let r = 0; r < gridY; r++) {
                const energyVal = maxE - (r / (gridY - 1)) * (maxE - minE);
                
                let density = 0.0;
                for (let o = 0; o < pointOrbitals.length; o++) {
                    const diffE = energyVal - pointOrbitals[o][0];
                    if (Math.abs(diffE) < 6.0 * T_eff) {
                        const arg = diffE / (2.0 * T_eff);
                        const ex = Math.exp(arg);
                        const cosh = 0.5 * (ex + 1.0 / ex);
                        const sech2 = 1.0 / (cosh * cosh);
                        density += sech2 / (4.0 * T_eff);
                    }
                }

                // Normalise density scaling based on temperature
                const maxDensityVal = 2.2 / Math.sqrt(T_eff);
                const intensity = Math.min(1.0, density / maxDensityVal);
                
                const idx = (r * gridX + c) * 4;
                const isProton = state.nucleonType === 'proton';
                if (isProton) {
                    tempImgData[idx] = 0;
                    tempImgData[idx + 1] = 240;
                    tempImgData[idx + 2] = 255;
                } else {
                    tempImgData[idx] = 16;
                    tempImgData[idx + 1] = 185;
                    tempImgData[idx + 2] = 129;
                }
                // Alpha glows/smeares with level density
                tempImgData[idx + 3] = Math.round(intensity * 150);
            }
        }
        tempCtx.putImageData(tempImg, 0, 0);

        ctx.save();
        ctx.imageSmoothingEnabled = true;
        ctx.imageSmoothingQuality = 'high';
        ctx.drawImage(tempCanvas, margin.left, margin.top, chartWidth, chartHeight);
        ctx.restore();
    }

    // Draw Level Curves    
    const numOrbitals = grid[0].data.orbitals.length;
    const fermiIdx = Math.ceil(numParticles / 2) - 1;

    for (let oIdx = 0; oIdx < numOrbitals; oIdx++) {
        // Draw level curve segment-by-segment to show Fermi-Dirac smearing along the deformation path
        for (let pIdx = 0; pIdx < grid.length - 1; pIdx++) {
            const p1 = grid[pIdx];
            const p2 = grid[pIdx + 1];
            
            const e1 = p1.data.orbitals[oIdx][0];
            const e2 = p2.data.orbitals[oIdx][0];
            
            // Calculate Fermi level (chemical potential) at p1 and p2 beta positions
            const p1Highest = p1.data.orbitals[fermiIdx] ? p1.data.orbitals[fermiIdx][0] : 0;
            const p1Lowest = p1.data.orbitals[fermiIdx + 1] ? p1.data.orbitals[fermiIdx + 1][0] : p1Highest;
            const pLambda1 = (p1Highest + p1Lowest) / 2.0;

            const p2Highest = p2.data.orbitals[fermiIdx] ? p2.data.orbitals[fermiIdx][0] : 0;
            const p2Lowest = p2.data.orbitals[fermiIdx + 1] ? p2.data.orbitals[fermiIdx + 1][0] : p2Highest;
            const pLambda2 = (p2Highest + p2Lowest) / 2.0;

            const n1 = fermiDiracOccupation(e1, pLambda1, state.T);
            const n2 = fermiDiracOccupation(e2, pLambda2, state.T);
            const avgN = (n1 + n2) / 2.0;
            
            ctx.beginPath();
            ctx.moveTo(getX(p1.beta), getY(e1));
            ctx.lineTo(getX(p2.beta), getY(e2));
            
            // Set styles dynamically based on local occupation factor, fading out discrete lines at high T
            const lineFade = Math.max(0.05, 1.0 - state.T / 1.6);
            ctx.globalAlpha = (0.25 + 0.6 * avgN) * lineFade;
            ctx.strokeStyle = avgN > 0.5 ? accentColor : '#4b5563';
            ctx.lineWidth = 1.0 + 1.0 * avgN;
            ctx.stroke();
        }
        
        // Draw orbital shells labels on the right edge for levels near Fermi surface
        const lastPoint = grid[grid.length - 1];
        const lastVal = lastPoint.data.orbitals[oIdx][0];
        const lastLabel = lastPoint.data.orbitals[oIdx][2];
        const isNearFermi = Math.abs(oIdx - Math.ceil(numParticles / 2)) <= 4;
        
        if (isNearFermi) {
            const lastHighest = lastPoint.data.orbitals[fermiIdx] ? lastPoint.data.orbitals[fermiIdx][0] : 0;
            const lastLowest = lastPoint.data.orbitals[fermiIdx + 1] ? lastPoint.data.orbitals[fermiIdx + 1][0] : lastHighest;
            const lastLambda = (lastHighest + lastLowest) / 2.0;
            const lastN = fermiDiracOccupation(lastVal, lastLambda, state.T);

            ctx.fillStyle = lastN > 0.5 ? accentColor : '#9ca3af';
            ctx.globalAlpha = 0.35 + 0.55 * lastN;
            ctx.font = '9px Fira Code';
            ctx.textAlign = 'left';
            ctx.fillText(lastLabel, margin.left + chartWidth + 5, getY(lastVal) + 3);
        }
    }
    ctx.globalAlpha = 1.0;
    
    // Draw Current Beta vertical cursor line
    const curX = getX(state.beta);
    ctx.strokeStyle = '#fbbf24';
    ctx.lineWidth = 1.5;
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(curX, margin.top);
    ctx.lineTo(curX, margin.top + chartHeight);
    ctx.stroke();
    ctx.setLineDash([]);
    
    // Draw intersection dots and highlight Fermi Level with temperature-dependent sizes & opacities
    const currentEnergies = nilssonEnergies(state.beta, state.gamma_shape, A, state.nucleonType, state.paramSet).orbitals;
    
    const curHighest = currentEnergies[fermiIdx] ? currentEnergies[fermiIdx][0] : 0;
    const curLowest = currentEnergies[fermiIdx + 1] ? currentEnergies[fermiIdx + 1][0] : curHighest;
    const lambdaSystem = (curHighest + curLowest) / 2.0;
    let fermiEnergy = currentEnergies[fermiIdx] ? currentEnergies[fermiIdx][0] : 0;
    
    currentEnergies.forEach((orb, oIdx) => {
        const y = getY(orb[0]);
        const n_i = fermiDiracOccupation(orb[0], lambdaSystem, state.T);
        
        if (n_i > 0.005) {
            // Draw filled dot representing orbital occupancy
            ctx.fillStyle = accentColor;
            ctx.globalAlpha = 0.15 + 0.75 * n_i; // scale opacity with occupancy
            ctx.beginPath();
            ctx.arc(curX, y, 1.5 + 3.0 * n_i, 0, 2 * Math.PI); // scale size with occupancy
            ctx.fill();
            
            if (oIdx === fermiIdx) {
                ctx.strokeStyle = '#fbbf24';
                ctx.lineWidth = 1.5;
                ctx.globalAlpha = 0.85;
                ctx.beginPath();
                ctx.arc(curX, y, 7, 0, 2 * Math.PI);
                ctx.stroke();
                
                ctx.fillStyle = '#fbbf24';
                ctx.font = '10px Fira Code';
                ctx.textAlign = state.beta > 0.3 ? 'right' : 'left';
                const xOffset = state.beta > 0.3 ? -12 : 12;
                ctx.fillText(`Fermi Level: ${fermiEnergy.toFixed(2)} MeV`, curX + xOffset, y - 6);
            }
        }
    });
    ctx.globalAlpha = 1.0;
    
    // Title inside chart (appended with temperature)
    ctx.fillStyle = '#fff';
    ctx.font = '12px Outfit';
    ctx.textAlign = 'left';
    ctx.fillText(`${state.nucleonType.toUpperCase()} LEVELS (γ=${state.gamma_shape}°, T=${state.T.toFixed(2)} MeV)`, margin.left + 10, margin.top + 20);
}

// Draw the Strutinsky Plateau Diagram (deltaE_shell vs. gamma_smear)
export function drawPlateauDiagram(canvas, state) {
    const { ctx, width, height } = getScaledContext(canvas);
    if (!ctx) return;
    
    ctx.clearRect(0, 0, width, height);
    
    const margin = { top: 30, right: 30, bottom: 40, left: 60 };
    const chartWidth = width - margin.left - margin.right;
    const chartHeight = height - margin.top - margin.bottom;
    
    const A = state.Z + state.N;
    const hw0_base = 41.0 * Math.pow(A, -1.0 / 3.0);
    
    // Compute levels at current shape (beta, gamma_shape)
    const currentProtons = nilssonEnergies(state.beta, state.gamma_shape, A, 'proton', state.paramSet).orbitals;
    const currentNeutrons = nilssonEnergies(state.beta, state.gamma_shape, A, 'neutron', state.paramSet).orbitals;
    
    // Generate data points: gamma_factor from 0.3 to 2.2
    const gammas = [];
    const shellCorrections = [];
    const steps = 40;
    
    const fDamp = thermalDampingFactor(state.T, A);
    for (let i = 0; i <= steps; i++) {
        const gFact = 0.3 + (i / steps) * 1.9;
        const gMeV = gFact * hw0_base;
        
        const shp = calculateShellCorrection(currentProtons, state.Z, gMeV, state.pOrder);
        const shn = calculateShellCorrection(currentNeutrons, state.N, gMeV, state.pOrder);
        
        gammas.push(gFact);
        shellCorrections.push((shp.deltaE + shn.deltaE) * fDamp);
    }
    
    let minCorrection = Math.min(...shellCorrections);
    let maxCorrection = Math.max(...shellCorrections);
    
    let diff = maxCorrection - minCorrection;
    if (diff < 2.0) {
        const mid = (maxCorrection + minCorrection) / 2.0;
        minCorrection = mid - 1.0;
        maxCorrection = mid + 1.0;
    } else {
        minCorrection -= diff * 0.15 + 1e-5;
        maxCorrection += diff * 0.15 + 1e-5;
    }
    
    const getX = (g) => margin.left + ((g - 0.3) / 1.9) * chartWidth;
    const getY = (e) => margin.top + chartHeight - ((e - minCorrection) / (maxCorrection - minCorrection)) * chartHeight;
    
    // Draw Gridlines
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.05)';
    ctx.lineWidth = 1;
    
    for (let g = 0.5; g <= 2.2; g += 0.5) {
        const x = getX(g);
        ctx.beginPath();
        ctx.moveTo(x, margin.top);
        ctx.lineTo(x, margin.top + chartHeight);
        ctx.stroke();
        
        ctx.fillStyle = '#6b7280';
        ctx.font = '10px Inter';
        ctx.textAlign = 'center';
        ctx.fillText(g.toFixed(1), x, margin.top + chartHeight + 15);
    }
    
    const step = (maxCorrection - minCorrection) / 5;
    for (let e = minCorrection; e <= maxCorrection; e += step) {
        const y = getY(e);
        ctx.beginPath();
        ctx.moveTo(margin.left, y);
        ctx.lineTo(margin.left + chartWidth, y);
        ctx.stroke();
        
        ctx.fillStyle = '#6b7280';
        ctx.font = '10px Inter';
        ctx.textAlign = 'right';
        ctx.fillText(e.toFixed(1) + ' MeV', margin.left - 10, y + 3);
    }
    
    ctx.fillStyle = '#9ca3af';
    ctx.font = '11px Outfit';
    ctx.textAlign = 'center';
    ctx.fillText('Smearing Parameter (γ_smear / ħω₀)', margin.left + chartWidth / 2, margin.top + chartHeight + 35);
    
    ctx.save();
    ctx.translate(15, margin.top + chartHeight / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText('Shell Correction Energy δE_shell (MeV)', 0, 0);
    ctx.restore();
    
    // Draw the curve
    ctx.strokeStyle = '#10b981';
    ctx.lineWidth = 3;
    ctx.beginPath();
    gammas.forEach((g, idx) => {
        const x = getX(g);
        const y = getY(shellCorrections[idx]);
        if (idx === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
    });
    ctx.stroke();
    
    // Highlight stable plateau region
    const platX1 = getX(1.0);
    const platX2 = getX(1.5);
    ctx.fillStyle = 'rgba(16, 185, 129, 0.06)';
    ctx.fillRect(platX1, margin.top, platX2 - platX1, chartHeight);
    
    ctx.strokeStyle = 'rgba(16, 185, 129, 0.2)';
    ctx.lineWidth = 1;
    ctx.setLineDash([3, 3]);
    ctx.beginPath();
    ctx.moveTo(platX1, margin.top); ctx.lineTo(platX1, margin.top + chartHeight);
    ctx.moveTo(platX2, margin.top); ctx.lineTo(platX2, margin.top + chartHeight);
    ctx.stroke();
    ctx.setLineDash([]);
    
    ctx.fillStyle = 'rgba(16, 185, 129, 0.7)';
    ctx.font = '9px Outfit';
    ctx.textAlign = 'center';
    ctx.fillText('STABLE PLATEAU REGION', (platX1 + platX2) / 2, margin.top + 15);
    
    // Highlight Current Smearing Width
    const curG = state.gamma; // factor
    const curGMeV = curG * hw0_base;
    const curCorrectP = calculateShellCorrection(currentProtons, state.Z, curGMeV, state.pOrder).deltaE * fDamp;
    const curCorrectN = calculateShellCorrection(currentNeutrons, state.N, curGMeV, state.pOrder).deltaE * fDamp;
    const curCorrect = curCorrectP + curCorrectN;
    
    const curX = getX(curG);
    const curY = getY(curCorrect);
    
    ctx.fillStyle = '#fbbf24';
    ctx.beginPath();
    ctx.arc(curX, curY, 5, 0, 2 * Math.PI);
    ctx.fill();
    ctx.strokeStyle = '#fff';
    ctx.lineWidth = 1.5;
    ctx.stroke();
    
    ctx.fillStyle = '#fbbf24';
    ctx.font = '10px Fira Code';
    ctx.textAlign = curG > 1.3 ? 'right' : 'left';
    const labelOffset = curG > 1.3 ? -10 : 10;
    ctx.fillText(`δE = ${curCorrect.toFixed(2)} MeV`, curX + labelOffset, curY - 6);
    
    ctx.fillStyle = '#fff';
    ctx.font = '12px Outfit';
    ctx.textAlign = 'left';
    ctx.fillText(`PLATEAU (β = ${state.beta.toFixed(2)}, γ = ${state.gamma_shape}°, T = ${state.T.toFixed(2)} MeV)`, margin.left + 10, margin.top + 20);
}

// -------------------------------------------------------------
// Interactive 3D Rendering engine using Three.js
// -------------------------------------------------------------

// WebGL Global instances
let nScene, nCamera, nRenderer, nMesh, nControls;
let pScene, pCamera, pRenderer, pMesh, pMarker, pStar, pControls;
let borderLines = [];

// Initialize 3D Viewports
export function init3D(nucleusContainer, pesContainer) {
    if (!window.THREE) {
        console.error("Three.js not loaded.");
        return;
    }

    // 1. Initialize 3D Nucleus Shape Viewer
    try {
        if (nucleusContainer && !nRenderer) {
            const width = nucleusContainer.clientWidth;
            const height = nucleusContainer.clientHeight;

            nScene = new THREE.Scene();
            nScene.background = null; // transparent background

            nCamera = new THREE.PerspectiveCamera(40, width / height, 0.1, 10);
            // Position camera and adjust dynamically for narrow viewports
            const aspect = width / height;
            const baseDistance = 3.6;
            nCamera.position.set(0, 0, aspect < 1.0 ? baseDistance / aspect : baseDistance);

            nRenderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
            nRenderer.setSize(width, height);
            nRenderer.setPixelRatio(window.devicePixelRatio);
            nucleusContainer.appendChild(nRenderer.domElement);

            if (THREE.OrbitControls) {
                nControls = new THREE.OrbitControls(nCamera, nRenderer.domElement);
                nControls.enableZoom = true;
                nControls.enablePan = false;
                nControls.autoRotate = true;
                nControls.autoRotateSpeed = 1.5;
            } else {
                console.warn("THREE.OrbitControls not found. Orbit controls disabled for nucleus viewer.");
            }

            // Lighting
            const ambientLight = new THREE.AmbientLight(0xffffff, 0.45);
            nScene.add(ambientLight);

            const dirLight1 = new THREE.DirectionalLight(0x00f0ff, 0.8); // Cyan highlight
            dirLight1.position.set(2, 3, 4);
            nScene.add(dirLight1);

            const dirLight2 = new THREE.DirectionalLight(0xf43f5e, 0.45); // Rose fill
            dirLight2.position.set(-2, -3, -4);
            nScene.add(dirLight2);

            // Deformable Sphere Mesh (reduced size to fit comfortably in responsive window)
            const sphereGeom = new THREE.SphereGeometry(0.9, 48, 48);
            sphereGeom.userData.originalPosition = sphereGeom.attributes.position.clone();

            const sphereMat = new THREE.MeshStandardMaterial({
                color: 0x00f0ff,
                roughness: 0.2,
                metalness: 0.1,
                emissive: 0x0c1328,
                flatShading: false
            });

            nMesh = new THREE.Mesh(sphereGeom, sphereMat);
            nScene.add(nMesh);
            
            // Render loop for Orbit rotation
            const animateNucleus = () => {
                requestAnimationFrame(animateNucleus);
                if (nControls) nControls.update();
                nRenderer.render(nScene, nCamera);
            };
            animateNucleus();
        }
    } catch (e) {
        console.error("Failed to initialize 3D Nucleus Viewer:", e);
    }

    // 2. Initialize 3D PES Landscape Viewer
    try {
        if (pesContainer && !pRenderer) {
            const width = pesContainer.clientWidth;
            const height = pesContainer.clientHeight;

            pScene = new THREE.Scene();
            pScene.background = new THREE.Color(0x07090e); // Dark space bg

            pCamera = new THREE.PerspectiveCamera(45, width / height, 0.1, 50);
            pCamera.position.set(4, -8, 5); // Isometric angle looking at the sector
            pCamera.up.set(0, 0, 1);

            pRenderer = new THREE.WebGLRenderer({ antialias: true });
            pRenderer.setSize(width, height);
            pRenderer.setPixelRatio(window.devicePixelRatio);
            pesContainer.appendChild(pRenderer.domElement);

            if (THREE.OrbitControls) {
                pControls = new THREE.OrbitControls(pCamera, pRenderer.domElement);
                pControls.enableZoom = true;
                pControls.maxPolarAngle = Math.PI / 2 - 0.05; // Prevent camera from going under the base plane
                pControls.minDistance = 3;
                pControls.maxDistance = 15;
            } else {
                console.warn("THREE.OrbitControls not found. Orbit controls disabled for PES viewer.");
            }

            // Lighting
            const ambientLight = new THREE.AmbientLight(0xffffff, 0.35);
            pScene.add(ambientLight);

            const dirLight = new THREE.DirectionalLight(0xffffff, 0.95);
            dirLight.position.set(5, -5, 10);
            pScene.add(dirLight);

            // Current state marker (yellow glowing ball)
            const markerGeom = new THREE.SphereGeometry(0.12, 16, 16);
            const markerMat = new THREE.MeshBasicMaterial({ color: 0xfbbf24 });
            pMarker = new THREE.Mesh(markerGeom, markerMat);
            pScene.add(pMarker);

            // Ground state stable minimum marker (glowing octahedron)
            const starGeom = new THREE.OctahedronGeometry(0.18, 0);
            const starMat = new THREE.MeshBasicMaterial({ color: 0xf43f5e });
            pStar = new THREE.Mesh(starGeom, starMat);
            pScene.add(pStar);

            // Render loop
            const animatePES = () => {
                requestAnimationFrame(animatePES);
                if (pControls) pControls.update();
                pRenderer.render(pScene, pCamera);
            };
            animatePES();
        }
    } catch (e) {
        console.error("Failed to initialize 3D PES Viewer:", e);
    }
}

// Resize Three.js viewports
export function resize3D(nucleusContainer, pesContainer) {
    if (nRenderer && nucleusContainer) {
        const w = nucleusContainer.clientWidth;
        const h = nucleusContainer.clientHeight;
        nRenderer.setSize(w, h);
        
        const aspect = w / h;
        nCamera.aspect = aspect;
        nCamera.updateProjectionMatrix();
        
        // Adjust camera Z based on aspect ratio to prevent clipping when narrow
        const baseDistance = 3.6;
        nCamera.position.z = aspect < 1.0 ? baseDistance / aspect : baseDistance;
    }
    if (pRenderer && pesContainer) {
        const w = pesContainer.clientWidth;
        const h = pesContainer.clientHeight;
        pRenderer.setSize(w, h);
        pCamera.aspect = w / h;
        pCamera.updateProjectionMatrix();
    }
}

// Update 3D Nucleus mesh shape deformation in real-time
export function updateNucleus3D(beta, gamma) {
    if (!nMesh) return;
    
    const geom = nMesh.geometry;
    const pos = geom.attributes.position;
    const orig = geom.userData.originalPosition;
    
    if (!orig) return;
    
    const gamma_rad = gamma * Math.PI / 180.0;
    const cosG = Math.cos(gamma_rad);
    const sinG = Math.sin(gamma_rad);
    
    const origV = new THREE.Vector3();
    
    for (let i = 0; i < pos.count; i++) {
        origV.fromBufferAttribute(orig, i);
        
        const r = origV.length();
        if (r < 1e-4) continue;
        
        const cosTheta = origV.z / r;
        const theta = Math.acos(Math.max(-1.0, Math.min(1.0, cosTheta)));
        const phi = Math.atan2(origV.y, origV.x);
        
        // Spheroidal quadrupole expansions
        const y20 = 0.31539 * cosG * (3.0 * cosTheta * cosTheta - 1.0);
        const y22 = 0.54627 * sinG * (1.0 - cosTheta * cosTheta) * Math.cos(2.0 * phi);
        const scale = 1.0 + beta * (y20 + y22);
        
        pos.setXYZ(i, origV.x * scale, origV.y * scale, origV.z * scale);
    }
    
    pos.needsUpdate = true;
    geom.computeVertexNormals();
}

// Custom color palette generator for the PES surface (purple -> blue -> green -> yellow -> red)
function getPESColor(E) {
    const min = -18.0;
    const max = 22.0;
    const t = Math.max(0.0, Math.min(1.0, (E - min) / (max - min)));
    
    const color = new THREE.Color();
    if (t < 0.25) {
        color.setRGB(0.5 - t * 2.0, 0.0, 0.5 + t * 2.0); // purple to blue
    } else if (t < 0.5) {
        const f = (t - 0.25) * 4.0;
        color.setRGB(0.0, f, 1.0 - f);                  // blue to green
    } else if (t < 0.75) {
        const f = (t - 0.5) * 4.0;
        color.setRGB(f, 1.0, 0.0);                      // green to yellow
    } else {
        const f = (t - 0.75) * 4.0;
        color.setRGB(1.0, 1.0 - f, 0.0);                // yellow to red
    }
    return color;
}

// Catmull-Rom cubic interpolation helper
function cubicInterpolate(p0, p1, p2, p3, t) {
    return p1 + 0.5 * t * (p2 - p0 + t * (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3 + t * (3.0 * p1 - p0 - 3.0 * p2 + p3)));
}

// 2D bicubic interpolation from the polar grid
function bicubicInterpolate(grid, x, y, nBeta, nGamma) {
    const x_idx = Math.floor(x);
    const y_idx = Math.floor(y);
    
    const u = x - x_idx;
    const v = y - y_idx;
    
    const getVal = (i, j) => {
        const bi = Math.max(0, Math.min(nBeta, i));
        const gi = Math.max(0, Math.min(nGamma, j));
        return grid[bi][gi];
    };
    
    const arr = [];
    for (let j = -1; j <= 2; j++) {
        const p0 = getVal(x_idx - 1, y_idx + j);
        const p1 = getVal(x_idx,     y_idx + j);
        const p2 = getVal(x_idx + 1, y_idx + j);
        const p3 = getVal(x_idx + 2, y_idx + j);
        arr.push(cubicInterpolate(p0, p1, p2, p3, u));
    }
    
    return cubicInterpolate(arr[0], arr[1], arr[2], arr[3], v);
}

// Generate the 3D potential energy landscape surface (relative to spherical shape energy)
export function updatePES3D(state) {
    if (!pScene) return;

    // Clean up previous mesh and boundary lines
    if (pMesh) {
        pScene.remove(pMesh);
        pMesh.geometry.dispose();
        pMesh = null;
    }
    borderLines.forEach(l => pScene.remove(l));
    borderLines = [];

    // Regenerate cache if needed
    if (cachedPESGrid.Z !== state.Z ||
        cachedPESGrid.N !== state.N ||
        cachedPESGrid.paramSet !== state.paramSet ||
        cachedPESGrid.gamma !== state.gamma ||
        cachedPESGrid.pOrder !== state.pOrder ||
        !cachedPESGrid.ldmGrid) {
        updatePESCache(state);
    }

    const A = state.Z + state.N;
    const { energyGrid, minE_rel, minB, minG, nBeta, nGamma } = getPESEnergyGrid(state.T, A);

    const vertices = [];
    const colors = [];
    const indices = [];

    // High resolution 3D surface mesh (doubled resolution relative to calculation density)
    const mBeta = 36;
    const mGamma = 24;
    const db = 0.6 / nBeta;
    const dg = 60.0 / nGamma;

    // Build vertices and colors from the cached grid using bicubic interpolation for maximum fidelity
    for (let bi = 0; bi <= mBeta; bi++) {
        const b = (bi / mBeta) * 0.6;
        for (let gi = 0; gi <= mGamma; gi++) {
            const g = (gi / mGamma) * 60.0;
            
            const x = b * Math.cos(g * Math.PI / 180.0);
            const y = b * Math.sin(g * Math.PI / 180.0);
            
            const x_val = b / db;
            const y_val = g / dg;
            const E_rel = bicubicInterpolate(energyGrid, x_val, y_val, nBeta, nGamma);
            const z = E_rel * 0.12;
            
            // Scale X and Y by 5 for visibility in 3D scene (max beta = 3.0 units)
            vertices.push(x * 5.0, y * 5.0, z);
            
            const color = getPESColor(E_rel);
            colors.push(color.r, color.g, color.b);
        }
    }

    // Build mesh faces indices
    for (let bi = 0; bi < mBeta; bi++) {
        for (let gi = 0; gi < mGamma; gi++) {
            const i00 = bi * (mGamma + 1) + gi;
            const i10 = (bi + 1) * (mGamma + 1) + gi;
            const i01 = bi * (mGamma + 1) + gi + 1;
            const i11 = (bi + 1) * (mGamma + 1) + gi + 1;

            indices.push(i00, i10, i01);
            indices.push(i10, i11, i01);
        }
    }

    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
    geometry.setIndex(indices);
    geometry.computeVertexNormals();

    const material = new THREE.MeshStandardMaterial({
        vertexColors: true,
        roughness: 0.35,
        metalness: 0.1,
        side: THREE.DoubleSide
    });

    pMesh = new THREE.Mesh(geometry, material);
    pScene.add(pMesh);

    // Add high-tech wireframe grid lines overlay to the surface
    const wireGeom = new THREE.WireframeGeometry(geometry);
    const wireMat = new THREE.LineBasicMaterial({
        color: 0xffffff,
        transparent: true,
        opacity: 0.08,
        linewidth: 1
    });
    const wireGrid = new THREE.LineSegments(wireGeom, wireMat);
    pMesh.add(wireGrid);

    // Place the Ground State star at the minimum
    const starX = minB * Math.cos(minG * Math.PI / 180.0) * 5.0;
    const starY = minB * Math.sin(minG * Math.PI / 180.0) * 5.0;
    const starZ = minE_rel * 0.12;
    pStar.position.set(starX, starY, starZ + 0.12);

    // Draw Lund Plane grid boundaries at base (z = -2.2) using a closed loop (no crossing lines!)
    const baseZ = -2.2;
    const borderGeom = new THREE.BufferGeometry();
    const borderVertices = [];

    // Origin
    borderVertices.push(0, 0, baseZ);
    
    // Wedge arc boundary from 0 to 60 deg (which connects back to the origin, drawing prolate/oblate axes naturally)
    for (let g = 0; g <= 60; g += 2) {
        const r = g * Math.PI / 180.0;
        borderVertices.push(0.6 * 5.0 * Math.cos(r), 0.6 * 5.0 * Math.sin(r), baseZ);
    }

    borderGeom.setAttribute('position', new THREE.Float32BufferAttribute(borderVertices, 3));
    const borderMat = new THREE.LineBasicMaterial({ color: 0x374151 });
    const borderLine = new THREE.LineLoop(borderGeom, borderMat);
    pScene.add(borderLine);
    borderLines.push(borderLine);

    // Position Current Marker
    updateMarker3D(state);
}

// Update the position of the yellow marker sphere on the 3D surface
export function updateMarker3D(state) {
    if (!pMarker) return;
    
    // Ensure cache is populated
    if (cachedPESGrid.Z !== state.Z ||
        cachedPESGrid.N !== state.N ||
        cachedPESGrid.paramSet !== state.paramSet ||
        cachedPESGrid.gamma !== state.gamma ||
        cachedPESGrid.pOrder !== state.pOrder ||
        !cachedPESGrid.ldmGrid) {
        updatePESCache(state);
    }
    
    const A = state.Z + state.N;
    const { energyGrid, nBeta, nGamma } = getPESEnergyGrid(state.T, A);
    const db = 0.6 / nBeta;
    const dg = 60.0 / nGamma;
    
    // Bicubic interpolation
    const x_val = state.beta / db;
    const y_val = state.gamma_shape / dg;
    const curE_rel = bicubicInterpolate(energyGrid, x_val, y_val, nBeta, nGamma);
    
    const curX = state.beta * Math.cos(state.gamma_shape * Math.PI / 180.0) * 5.0;
    const curY = state.beta * Math.sin(state.gamma_shape * Math.PI / 180.0) * 5.0;
    const curZ = curE_rel * 0.12;
    
    pMarker.position.set(curX, curY, curZ + 0.12);
}

// Global cache for PES calculation to keep slider dragging at 60 FPS
let cachedPESGrid = {
    Z: -1,
    N: -1,
    paramSet: '',
    gamma: -1,
    pOrder: -1,
    ldmGrid: null,
    shellCorrGrid: null,
    nBeta: 18,
    nGamma: 12
};

// Compute and cache the entire 2D Potential Energy Surface (PES) at T = 0
export function updatePESCache(state) {
    const A = state.Z + state.N;
    const hw0_base = 41.0 * Math.pow(A, -1.0 / 3.0);
    const gMeV = state.gamma * hw0_base;

    // Calculate spherical reference energy at T = 0
    const pSph = nilssonEnergies(0, 0, A, 'proton', state.paramSet).orbitals;
    const nSph = nilssonEnergies(0, 0, A, 'neutron', state.paramSet).orbitals;
    const shpSph = calculateShellCorrection(pSph, state.Z, gMeV, state.pOrder).deltaE;
    const shnSph = calculateShellCorrection(nSph, state.N, gMeV, state.pOrder).deltaE;
    const E_sph = shpSph + shnSph;

    const nBeta = 18;
    const nGamma = 12;
    const ldmGrid = Array(nBeta + 1).fill(0).map(() => Array(nGamma + 1).fill(0));
    const shellCorrGrid = Array(nBeta + 1).fill(0).map(() => Array(nGamma + 1).fill(0));

    for (let bi = 0; bi <= nBeta; bi++) {
        const b = (bi / nBeta) * 0.6;
        for (let gi = 0; gi <= nGamma; gi++) {
            const g = (gi / nGamma) * 60.0;

            const ldm = calculateLDM(b, g, state.Z, state.N).E_ldm_rel;
            const pO = nilssonEnergies(b, g, A, 'proton', state.paramSet).orbitals;
            const nO = nilssonEnergies(b, g, A, 'neutron', state.paramSet).orbitals;
            const sP = calculateShellCorrection(pO, state.Z, gMeV, state.pOrder).deltaE;
            const sN = calculateShellCorrection(nO, state.N, gMeV, state.pOrder).deltaE;

            ldmGrid[bi][gi] = ldm;
            shellCorrGrid[bi][gi] = sP + sN - E_sph;
        }
    }

    cachedPESGrid = {
        Z: state.Z,
        N: state.N,
        paramSet: state.paramSet,
        gamma: state.gamma,
        pOrder: state.pOrder,
        ldmGrid,
        shellCorrGrid,
        nBeta,
        nGamma
    };
}

// Compute the temperature-damped potential energy surface dynamically
export function getPESEnergyGrid(T, A) {
    const fDamp = thermalDampingFactor(T, A);
    const { ldmGrid, shellCorrGrid, nBeta, nGamma } = cachedPESGrid;

    const energyGrid = Array(nBeta + 1).fill(0).map((_, bi) => 
        Array(nGamma + 1).fill(0).map((_, gi) => 
            ldmGrid[bi][gi] + shellCorrGrid[bi][gi] * fDamp
        )
    );

    let minE_rel = Infinity;
    let minB = 0;
    let minG = 0;
    const db = 0.6 / nBeta;
    const dg = 60.0 / nGamma;

    for (let bi = 0; bi <= nBeta; bi++) {
        const b = bi * db;
        for (let gi = 0; gi <= nGamma; gi++) {
            const E = energyGrid[bi][gi];
            if (E < minE_rel) {
                minE_rel = E;
                minB = b;
                minG = gi * dg;
            }
        }
    }

    return { energyGrid, minE_rel, minB, minG, nBeta, nGamma };
}

// Draw the 2D polar heatmap projection of the Potential Energy Surface
export function drawPESHeatmap(canvas, state) {
    const { ctx, width, height } = getScaledContext(canvas);
    if (!ctx) return;
    
    ctx.clearRect(0, 0, width, height);
    
    const margin = { top: 35, right: 40, bottom: 45, left: 60 };
    const chartWidth = width - margin.left - margin.right;
    const chartHeight = height - margin.top - margin.bottom;
    
    // Bounding polar sector center and dimensions
    const centerX = margin.left + 30;
    const centerY = margin.top + chartHeight - 20;
    const maxRadius = Math.min(chartWidth - 120, chartHeight - 40);
    
    // Regenerate cache if needed
    if (cachedPESGrid.Z !== state.Z ||
        cachedPESGrid.N !== state.N ||
        cachedPESGrid.paramSet !== state.paramSet ||
        cachedPESGrid.gamma !== state.gamma ||
        cachedPESGrid.pOrder !== state.pOrder ||
        !cachedPESGrid.ldmGrid) {
        updatePESCache(state);
    }
    
    const A = state.Z + state.N;
    const { energyGrid, minE_rel, minB, minG, nBeta, nGamma } = getPESEnergyGrid(state.T, A);
    
    // Create offscreen canvas for rendering the interpolated polar wedge
    const offscreenSize = 250;
    const offCanvas = document.createElement('canvas');
    offCanvas.width = offscreenSize;
    offCanvas.height = offscreenSize;
    const offCtx = offCanvas.getContext('2d');
    
    const offCX = 0;
    const offCY = offscreenSize;
    const offR = offscreenSize;
    
    const offImg = offCtx.createImageData(offscreenSize, offscreenSize);
    const data = offImg.data;
    
    const db = 0.6 / nBeta;
    const dg = 60.0 / nGamma;
    
    // Grid to keep track of discretized level indices for contour drawing
    const levelGrid = Array(offscreenSize).fill(0).map(() => Array(offscreenSize).fill(null));
    
    // Color mapping helper with 14 discrete bands (contourf filled contour style)
    const getPESColorValues = (E) => {
        const min = -18.0;
        const max = 22.0;
        let t = Math.max(0.0, Math.min(1.0, (E - min) / (max - min)));
        
        // Discretize into 14 bands
        const levels = 14;
        t = Math.floor(t * levels) / levels;
        
        let r, g, b;
        if (t < 0.25) {
            r = Math.round((0.5 - t * 2.0) * 255);
            g = 0;
            b = Math.round((0.5 + t * 2.0) * 255);
        } else if (t < 0.5) {
            const f = (t - 0.25) * 4.0;
            r = 0;
            g = Math.round(f * 255);
            b = Math.round((1.0 - f) * 255);
        } else if (t < 0.75) {
            const f = (t - 0.5) * 4.0;
            r = Math.round(f * 255);
            g = 255;
            b = 0;
        } else {
            const f = (t - 0.75) * 4.0;
            r = 255;
            g = Math.round((1.0 - f) * 255);
            b = 0;
        }
        return { r, g, b };
    };
    
    // Evaluate color at each pixel inside the polar wedge on offscreen using bicubic interpolation for high-fidelity smooth borders
    for (let py = 0; py < offscreenSize; py++) {
        const dy = offCY - py; // screen y goes down
        for (let px = 0; px < offscreenSize; px++) {
            const dx = px - offCX;
            const r = Math.sqrt(dx * dx + dy * dy);
            
            if (r > offR) continue;
            
            const beta = (r / offR) * 0.6;
            const gamma = Math.atan2(dy, dx) * 180 / Math.PI;
            
            if (beta <= 0.6 && gamma <= 60.001) {
                // Bicubic interpolation coordinates
                const x_val = beta / db;
                const y_val = gamma / dg;
                const E = bicubicInterpolate(energyGrid, x_val, y_val, nBeta, nGamma);
                
                // Store level index for contour line border drawing
                const min = -18.0;
                const max = 22.0;
                const t = Math.max(0.0, Math.min(1.0, (E - min) / (max - min)));
                const levels = 14;
                levelGrid[px][py] = Math.floor(t * levels);
                
                const c = getPESColorValues(E);
                const idx = (py * offscreenSize + px) * 4;
                data[idx] = c.r;
                data[idx + 1] = c.g;
                data[idx + 2] = c.b;
                data[idx + 3] = 255;
            }
        }
    }
    
    // Draw thin contour boundary lines between different energy levels
    for (let py = 1; py < offscreenSize - 1; py++) {
        for (let px = 1; px < offscreenSize - 1; px++) {
            const currentLevel = levelGrid[px][py];
            if (currentLevel !== null) {
                const rightLevel = levelGrid[px + 1][py];
                const bottomLevel = levelGrid[px][py + 1];
                
                if ((rightLevel !== null && rightLevel !== currentLevel) || 
                    (bottomLevel !== null && bottomLevel !== currentLevel)) {
                    const idx = (py * offscreenSize + px) * 4;
                    // Dark navy charcoal border line
                    data[idx] = 13;
                    data[idx + 1] = 18;
                    data[idx + 2] = 28;
                    data[idx + 3] = 255;
                }
            }
        }
    }
    
    offCtx.putImageData(offImg, 0, 0);
    
    // Draw the smooth interpolated offscreen canvas onto the main canvas
    ctx.save();
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(offCanvas, 0, 0, offscreenSize, offscreenSize, centerX, centerY - maxRadius, maxRadius, maxRadius);
    ctx.restore();
    
    // 3. Draw Grid Lines (Matplotlib polar contourf grid lines style - increased visibility)
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.28)';
    ctx.lineWidth = 1.0;
    ctx.setLineDash([2, 3]);
    
    // Interior beta arcs at 0.1, 0.2, 0.3, 0.4, 0.5
    for (let b = 0.1; b < 0.6; b += 0.1) {
        const r = (b / 0.6) * maxRadius;
        ctx.beginPath();
        ctx.arc(centerX, centerY, r, 0, -Math.PI / 3, true);
        ctx.stroke();
        
        // Label values along 15° ray
        ctx.fillStyle = 'rgba(255, 255, 255, 0.7)';
        ctx.font = '8px Fira Code';
        ctx.textAlign = 'center';
        const labelRad = 15 * Math.PI / 180.0;
        ctx.fillText(b.toFixed(1), centerX + r * Math.cos(labelRad), centerY - r * Math.sin(labelRad) + 3);
    }
    
    // Interior gamma rays at 10°, 20°, 30°, 40°, 50°
    for (let g = 10; g < 60; g += 10) {
        const rad = g * Math.PI / 180.0;
        ctx.beginPath();
        ctx.moveTo(centerX, centerY);
        ctx.lineTo(centerX + maxRadius * Math.cos(rad), centerY - maxRadius * Math.sin(rad));
        ctx.stroke();
    }
    
    ctx.setLineDash([]); // Reset line dash
    
    // Draw solid outer boundaries (Wedge borders)
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.65)';
    ctx.lineWidth = 1.5;
    
    ctx.beginPath();
    // Prolate axis
    ctx.moveTo(centerX, centerY);
    ctx.lineTo(centerX + maxRadius, centerY);
    // Oblate axis
    ctx.moveTo(centerX, centerY);
    ctx.lineTo(centerX + maxRadius * Math.cos(Math.PI / 3), centerY - maxRadius * Math.sin(Math.PI / 3));
    ctx.stroke();
    
    // Outer arc
    ctx.beginPath();
    ctx.arc(centerX, centerY, maxRadius, 0, -Math.PI / 3, true);
    ctx.stroke();
    
    // Labels for gamma rays around the outer boundary
    ctx.fillStyle = '#e5e7eb';
    ctx.font = '9px Outfit';
    
    for (let g = 0; g <= 60; g += 10) {
        const rad = g * Math.PI / 180.0;
        const labelR = maxRadius + 8;
        const lx = centerX + labelR * Math.cos(rad);
        const ly = centerY - labelR * Math.sin(rad);
        
        ctx.textAlign = (g === 0) ? 'left' : (g === 60) ? 'right' : 'center';
        
        let labelText = `${g}°`;
        if (g === 0) labelText = 'Prolate (γ = 0°)';
        if (g === 60) labelText = 'Oblate (γ = 60°)';
        if (g === 30) labelText = 'Triaxial (γ = 30°)';
        
        ctx.fillText(labelText, lx, ly + 3);
    }
    
    // 4. Draw Ground State Minimum (Red Star)
    const starX = centerX + (minB / 0.6) * maxRadius * Math.cos(minG * Math.PI / 180.0);
    const starY = centerY - (minB / 0.6) * maxRadius * Math.sin(minG * Math.PI / 180.0);
    
    ctx.fillStyle = '#f43f5e';
    ctx.beginPath();
    ctx.arc(starX, starY, 6, 0, 2 * Math.PI);
    ctx.fill();
    ctx.strokeStyle = '#fff';
    ctx.lineWidth = 1;
    ctx.stroke();
    
    ctx.fillStyle = '#fbbf24';
    ctx.font = '9px Outfit';
    ctx.textAlign = 'center';
    ctx.fillText('Ground State', starX, starY - 9);
    
    // 5. Draw Current Position (Yellow Dot)
    const curX = centerX + (state.beta / 0.6) * maxRadius * Math.cos(state.gamma_shape * Math.PI / 180.0);
    const curY = centerY - (state.beta / 0.6) * maxRadius * Math.sin(state.gamma_shape * Math.PI / 180.0);
    
    ctx.fillStyle = '#fbbf24';
    ctx.beginPath();
    ctx.arc(curX, curY, 5, 0, 2 * Math.PI);
    ctx.fill();
    ctx.strokeStyle = '#fff';
    ctx.lineWidth = 1.5;
    ctx.stroke();
    
    ctx.fillStyle = '#fff';
    ctx.font = '12px Outfit';
    ctx.textAlign = 'left';
    ctx.fillText(`2D PES HEATMAP (E_sph = 0)`, margin.left + 10, margin.top + 20);
}

// Translate mouse coordinates to polar coordinates (beta, gamma_shape)
export function getLundCoordinatesFromEvent(canvas, event) {
    const rect = canvas.getBoundingClientRect();
    
    // Client event coordinates mapped to canvas bounds
    const x = event.clientX - rect.left;
    const y = event.clientY - rect.top;
    
    const margin = { top: 35, right: 40, bottom: 45, left: 60 };
    const chartWidth = rect.width - margin.left - margin.right;
    const chartHeight = rect.height - margin.top - margin.bottom;
    
    const centerX = margin.left + 30;
    const centerY = margin.top + chartHeight - 20;
    const maxRadius = Math.min(chartWidth - 120, chartHeight - 40);
    
    const dx = x - centerX;
    const dy = centerY - y; // Screen coordinate system inversion
    
    const r = Math.sqrt(dx * dx + dy * dy);
    const beta = Math.max(0.0, Math.min(0.6, (r / maxRadius) * 0.6));
    
    let gamma = Math.atan2(dy, dx) * 180 / Math.PI;
    if (gamma < -180) gamma += 360;
    if (gamma > 180) gamma -= 360;
    
    const gamma_shape = Math.round(Math.max(0.0, Math.min(60.0, gamma)));
    
    return { beta, gamma_shape };
}
