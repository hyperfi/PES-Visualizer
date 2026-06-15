// Main application orchestrator for the Strutinsky Shell Correction tool
import { PRESETS } from './physics.js';
import { 
    updateNilssonCache, 
    drawNilssonDiagram, 
    drawPlateauDiagram, 
    init3D, 
    resize3D, 
    updateNucleus3D, 
    updatePES3D, 
    updateMarker3D,
    drawPESHeatmap,
    updatePESCache,
    getLundCoordinatesFromEvent
} from './render.js';
import { nilssonEnergies, calculateShellCorrection, calculateLDM, thermalDampingFactor } from './physics.js';

// Application state
const state = {
    Z: 62,               // Protons
    N: 90,               // Neutrons
    beta: 0.20,          // Quadrupole deformation magnitude
    gamma_shape: 0,      // Triaxiality parameter (degrees)
    gamma: 1.20,         // Smearing width (factor of hbar*omega_0)
    T: 0.00,             // Nuclear temperature (MeV)
    nucleonType: 'proton', // 'proton' or 'neutron' for level diagram
    paramSet: 'universal', // Parameter set ID
    pOrder: 6,           // Hermite polynomial order
    pesViewMode: '3d'    // '3d' or '2d' view mode for PES
};

// DOM elements
const el = {
    presetSelect: document.getElementById('preset-select'),
    paramSelect: document.getElementById('param-select'),
    zSlider: document.getElementById('z-slider'),
    zVal: document.getElementById('z-val'),
    nSlider: document.getElementById('n-slider'),
    nVal: document.getElementById('n-val'),
    betaSlider: document.getElementById('beta-slider'),
    betaVal: document.getElementById('beta-val'),
    gammaShapeSlider: document.getElementById('gamma-shape-slider'),
    gammaShapeVal: document.getElementById('gamma-shape-val'),
    gammaSmearSlider: document.getElementById('gamma-smear-slider'),
    gammaSmearVal: document.getElementById('gamma-smear-val'),
    pOrderSelect: document.getElementById('porder-select'),
    
    // Toggle for nucleon type
    toggleProton: document.getElementById('toggle-proton'),
    toggleNeutron: document.getElementById('toggle-neutron'),
    
    // Toggle for PES representation (3D landscape vs 2D heatmap)
    togglePes3D: document.getElementById('toggle-pes-3d'),
    togglePes2D: document.getElementById('toggle-pes-2d'),
    
    // Canvas & 3D Containers
    nilssonCanvas: document.getElementById('nilsson-canvas'),
    plateauCanvas: document.getElementById('plateau-canvas'),
    pes2DCanvas: document.getElementById('pes-2d-canvas'),
    nucleus3DContainer: document.getElementById('nucleus-3d-container'),
    pes3DContainer: document.getElementById('pes-3d-container'),
    
    // Metrics elements
    metricZP: document.getElementById('m-zp'),
    metricNN: document.getElementById('m-nn'),
    metricA: document.getElementById('m-a'),
    metricShellP: document.getElementById('m-shell-p'),
    metricShellN: document.getElementById('m-shell-n'),
    metricShellTot: document.getElementById('m-shell-tot'),
    metricLdmRel: document.getElementById('m-ldm-rel'),
    metricEtot: document.getElementById('m-etot'),
    metricStableB: document.getElementById('m-stable-b'),
    metricStableG: document.getElementById('m-stable-g'),
    
    tempSlider: document.getElementById('temp-slider'),
    tempVal: document.getElementById('temp-val'),
    
    // Mobile Drawer elements
    drawerToggleBtn: document.getElementById('drawer-toggle-btn'),
    drawerOverlay: document.getElementById('drawer-overlay'),
    sidebar: document.getElementById('sidebar')
};

// Initialize event listeners and 3D scenes
function init() {
    console.log("PES Visualizer: init() started.");
    
    // 1. Populate presets select options
    PRESETS.forEach((preset, idx) => {
        const opt = document.createElement('option');
        opt.value = idx;
        opt.textContent = preset.name;
        el.presetSelect.appendChild(opt);
    });
    
    // Select first preset by default
    el.presetSelect.value = 0;
    
    // 2. Initialize Three.js viewport renderers
    init3D(el.nucleus3DContainer, el.pes3DContainer);
    
    // Load first preset initial state
    loadPreset(PRESETS[0]);
    
    // 3. Event listeners for inputs
    el.presetSelect.addEventListener('change', (e) => {
        const idx = parseInt(e.target.value);
        if (!isNaN(idx) && PRESETS[idx]) {
            loadPreset(PRESETS[idx]);
        }
    });
    
    el.paramSelect.addEventListener('change', (e) => {
        state.paramSet = e.target.value;
        onStructureChanged();
    });
    
    el.zSlider.addEventListener('input', (e) => {
        state.Z = parseInt(e.target.value);
        el.zVal.textContent = state.Z;
        onStructureChanged();
    });
    
    el.nSlider.addEventListener('input', (e) => {
        state.N = parseInt(e.target.value);
        el.nVal.textContent = state.N;
        onStructureChanged();
    });
    
    el.betaSlider.addEventListener('input', (e) => {
        state.beta = parseFloat(e.target.value);
        el.betaVal.textContent = state.beta.toFixed(2);
        onDeformationChanged();
    });
    
    el.gammaShapeSlider.addEventListener('input', (e) => {
        state.gamma_shape = parseInt(e.target.value);
        el.gammaShapeVal.textContent = `${state.gamma_shape}°`;
        onDeformationChanged();
    });
    
    el.gammaSmearSlider.addEventListener('input', (e) => {
        state.gamma = parseFloat(e.target.value);
        el.gammaSmearVal.textContent = state.gamma.toFixed(2);
        onSmearingChanged();
    });
    
    if (el.pOrderSelect) {
        el.pOrderSelect.addEventListener('change', (e) => {
            state.pOrder = parseInt(e.target.value);
            onSmearingChanged();
        });
    }
    
    if (el.tempSlider) {
        el.tempSlider.addEventListener('input', (e) => {
            state.T = parseFloat(e.target.value);
            console.log("PES Visualizer: Temperature slider input triggered. T =", state.T);
            el.tempVal.textContent = `${state.T.toFixed(2)} MeV`;
            onTemperatureChanged();
        });
    }
    
    console.log("PES Visualizer: init() completed. All event listeners registered.");
    
    // Nucleon level display toggle
    el.toggleProton.addEventListener('click', () => {
        state.nucleonType = 'proton';
        el.toggleProton.classList.add('active');
        el.toggleNeutron.classList.remove('active');
        drawAll();
    });
    
    el.toggleNeutron.addEventListener('click', () => {
        state.nucleonType = 'neutron';
        el.toggleNeutron.classList.add('active');
        el.toggleProton.classList.remove('active');
        drawAll();
    });
    
    // PES representation toggles
    el.togglePes3D.addEventListener('click', () => {
        state.pesViewMode = '3d';
        el.togglePes3D.classList.add('active');
        el.togglePes2D.classList.remove('active');
        el.pes3DContainer.style.display = 'block';
        el.pes2DCanvas.style.display = 'none';
        
        // Re-fit 3D viewport sizes and update
        resize3D(el.nucleus3DContainer, el.pes3DContainer);
        updatePES3D(state);
    });
    
    el.togglePes2D.addEventListener('click', () => {
        state.pesViewMode = '2d';
        el.togglePes2D.classList.add('active');
        el.togglePes3D.classList.remove('active');
        el.pes3DContainer.style.display = 'none';
        el.pes2DCanvas.style.display = 'block';
        
        drawPESHeatmap(el.pes2DCanvas, state);
    });
    
    // 2D Heatmap interactivity (mouse down/move/up and touch start/move/end)
    let isDraggingHeatmap = false;
    
    const handleHeatmapPointer = (e) => {
        const coords = getLundCoordinatesFromEvent(el.pes2DCanvas, e);
        state.beta = coords.beta;
        state.gamma_shape = coords.gamma_shape;
        
        // Update slider values and text
        el.betaSlider.value = state.beta;
        el.betaVal.textContent = state.beta.toFixed(2);
        el.gammaShapeSlider.value = state.gamma_shape;
        el.gammaShapeVal.textContent = `${state.gamma_shape}°`;
        
        onDeformationChanged();
    };
    
    el.pes2DCanvas.addEventListener('mousedown', (e) => {
        isDraggingHeatmap = true;
        handleHeatmapPointer(e);
    });
    
    window.addEventListener('mousemove', (e) => {
        if (isDraggingHeatmap) {
            handleHeatmapPointer(e);
        }
    });
    
    window.addEventListener('mouseup', () => {
        isDraggingHeatmap = false;
    });

    // Touch events for mobile responsiveness
    const handleTouchPointer = (e) => {
        if (e.touches.length > 0) {
            handleHeatmapPointer(e.touches[0]);
            e.preventDefault();
        }
    };

    el.pes2DCanvas.addEventListener('touchstart', (e) => {
        isDraggingHeatmap = true;
        handleTouchPointer(e);
    }, { passive: false });

    window.addEventListener('touchmove', (e) => {
        if (isDraggingHeatmap) {
            handleTouchPointer(e);
        }
    }, { passive: false });

    window.addEventListener('touchend', () => {
        isDraggingHeatmap = false;
    });
    
    // 4. Tab switching
    const tabButtons = document.querySelectorAll('.tab-btn');
    const viewPanels = document.querySelectorAll('.view-panel');
    
    tabButtons.forEach(btn => {
        btn.addEventListener('click', () => {
            const tabId = btn.getAttribute('data-tab');
            
            // Set button active
            tabButtons.forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            
            // Set panel active
            viewPanels.forEach(p => {
                if (p.id === `${tabId}-panel`) {
                    p.classList.add('active');
                } else {
                    p.classList.remove('active');
                }
            });
            
            // Re-fit 3D viewport sizes and update
            resize3D(el.nucleus3DContainer, el.pes3DContainer);
            drawAll();
            updatePESRepresentation();
        });
    });
    
    // 5. Mobile Drawer toggling
    if (el.drawerToggleBtn && el.sidebar && el.drawerOverlay) {
        el.drawerToggleBtn.addEventListener('click', (e) => {
            e.stopPropagation();
            el.sidebar.classList.toggle('active');
            el.drawerOverlay.classList.toggle('active');
        });
        
        el.drawerOverlay.addEventListener('click', () => {
            closeMobileDrawer();
        });
    }
    
    // 6. Initialize explainer states based on screen size
    adjustExplainersForScreenSize();
    
    // 7. Resize event
    window.addEventListener('resize', debounce(() => {
        resize3D(el.nucleus3DContainer, el.pes3DContainer);
        drawAll();
        updatePESRepresentation();
        adjustExplainersForScreenSize();
    }, 150));
}

// Close mobile parameters drawer
function closeMobileDrawer() {
    if (el.sidebar && el.drawerOverlay) {
        el.sidebar.classList.remove('active');
        el.drawerOverlay.classList.remove('active');
    }
}

// Automatically collapse or expand details based on window size
function adjustExplainersForScreenSize() {
    const isMobile = window.innerWidth <= 1024;
    document.querySelectorAll('.explainer-details').forEach(details => {
        details.open = !isMobile;
    });
}

// Load isotope preset
function loadPreset(preset) {
    state.Z = preset.Z;
    state.N = preset.N;
    
    // Default deformation depending on isotope characteristics
    if (preset.type === 'deformed') {
        state.beta = 0.25;
        state.gamma_shape = 0;
    } else if (preset.type === 'spherical') {
        state.beta = 0.05;
        state.gamma_shape = 0;
    } else {
        state.beta = 0.0;
        state.gamma_shape = 0;
    }
    
    state.T = 0.00;
    
    // Update sliders and text values
    el.zSlider.value = state.Z;
    el.zVal.textContent = state.Z;
    el.nSlider.value = state.N;
    el.nVal.textContent = state.N;
    el.betaSlider.value = state.beta;
    el.betaVal.textContent = state.beta.toFixed(2);
    el.gammaShapeSlider.value = state.gamma_shape;
    el.gammaShapeVal.textContent = `${state.gamma_shape}°`;
    if (el.tempSlider) {
        el.tempSlider.value = 0.00;
        el.tempVal.textContent = "0.00 MeV";
    }
    
    // Determine appropriate param set
    if (state.Z > 80 || state.N > 120) {
        state.paramSet = 'actinide';
    } else if (state.Z > 50 && state.N > 80) {
        state.paramSet = 'rare_earth';
    } else {
        state.paramSet = 'universal';
    }
    el.paramSelect.value = state.paramSet;
    
    closeMobileDrawer();
    onStructureChanged();
}

// Event triggered when Z, N, or paramSet changes
function onStructureChanged() {
    const A = state.Z + state.N;
    // Pre-cache Nilsson splitting grid at current gamma_shape
    updateNilssonCache(A, state.paramSet, state.gamma_shape);
    // Pre-cache the PES grid
    updatePESCache(state);
    updateMetrics();
    
    // Update energy surface visualizer
    updatePESRepresentation();
    drawAll();
}

// Event triggered when quadrupole beta or triaxiality gamma sliders change
function onDeformationChanged() {
    updateMetrics();
    
    // Nilsson levels: if gamma_shape changed, we need a cache update.
    drawNilssonDiagram(el.nilssonCanvas, state);
    drawPlateauDiagram(el.plateauCanvas, state);
    
    // Update 3D nucleus shape vertex positions
    updateNucleus3D(state.beta, state.gamma_shape);
    
    // Move yellow marker on 3D surface or redraw 2D marker
    if (state.pesViewMode === '3d') {
        updateMarker3D(state);
    } else {
        drawPESHeatmap(el.pes2DCanvas, state);
    }
}

// Event triggered when smearing width changes
function onSmearingChanged() {
    updatePESCache(state);
    updateMetrics();
    drawPlateauDiagram(el.plateauCanvas, state);
    updatePESRepresentation();
}

// Event triggered when nuclear temperature changes
function onTemperatureChanged() {
    updateMetrics();
    drawNilssonDiagram(el.nilssonCanvas, state);
    drawPlateauDiagram(el.plateauCanvas, state);
    updatePESRepresentation();
}

// Triggers active representation updates (3D surface mesh vs. 2D polar projection canvas)
function updatePESRepresentation() {
    if (state.pesViewMode === '3d') {
        updatePES3D(state);
    } else {
        drawPESHeatmap(el.pes2DCanvas, state);
    }
}

// Calculate values for current state and update dashboard metrics
function updateMetrics() {
    const A = state.Z + state.N;
    const hw0_base = 41.0 * Math.pow(A, -1.0 / 3.0);
    const gMeV = state.gamma * hw0_base;
    const fDamp = thermalDampingFactor(state.T, A);
    
    // Calculate current coordinates energy values
    const ldm = calculateLDM(state.beta, state.gamma_shape, state.Z, state.N);
    const pOrbs = nilssonEnergies(state.beta, state.gamma_shape, A, 'proton', state.paramSet).orbitals;
    const nOrbs = nilssonEnergies(state.beta, state.gamma_shape, A, 'neutron', state.paramSet).orbitals;
    
    const shp = calculateShellCorrection(pOrbs, state.Z, gMeV, state.pOrder);
    const shn = calculateShellCorrection(nOrbs, state.N, gMeV, state.pOrder).deltaE;
    
    const shpDamped = shp.deltaE * fDamp;
    const shnDamped = shn * fDamp;
    const shellCorr = shpDamped + shnDamped;
    const totalEnergy = ldm.E_ldm_rel + shellCorr;
    
    // Update numerical values
    el.metricZP.textContent = state.Z;
    el.metricNN.textContent = state.N;
    el.metricA.textContent = A;
    el.metricShellP.textContent = `${shpDamped.toFixed(2)} MeV`;
    el.metricShellN.textContent = `${shnDamped.toFixed(2)} MeV`;
    el.metricShellTot.textContent = `${shellCorr.toFixed(2)} MeV`;
    el.metricLdmRel.textContent = `${ldm.E_ldm_rel.toFixed(2)} MeV`;
    el.metricEtot.textContent = `${totalEnergy.toFixed(2)} MeV`;
    
    // Find ground state stable deformation by evaluating 2D grid in the Lund Plane
    let minEtot = Infinity;
    let stableBeta = 0;
    let stableGamma = 0;
    
    const stepsBeta = 20;
    const stepsGamma = 8;
    
    for (let bi = 0; bi <= stepsBeta; bi++) {
        const b = (bi / stepsBeta) * 0.6;
        for (let gi = 0; gi <= stepsGamma; gi++) {
            const g = (gi / stepsGamma) * 60.0;
            
            const curL = calculateLDM(b, g, state.Z, state.N).E_ldm_rel;
            const pO = nilssonEnergies(b, g, A, 'proton', state.paramSet).orbitals;
            const nO = nilssonEnergies(b, g, A, 'neutron', state.paramSet).orbitals;
            const sP = calculateShellCorrection(pO, state.Z, gMeV, state.pOrder).deltaE * fDamp;
            const sN = calculateShellCorrection(nO, state.N, gMeV, state.pOrder).deltaE * fDamp;
            const eTot = curL + sP + sN;
            
            if (eTot < minEtot) {
                minEtot = eTot;
                stableBeta = b;
                stableGamma = g;
            }
        }
    }
    
    el.metricStableB.textContent = stableBeta.toFixed(2);
    el.metricStableG.textContent = `${Math.round(stableGamma)}°`;
    
    // Visual indicators for deformed shapes
    if (stableBeta > 0.1) {
        el.metricStableB.style.color = '#fbbf24';
        el.metricStableG.style.color = '#fbbf24';
    } else {
        el.metricStableB.style.color = '#10b981';
        el.metricStableG.style.color = '#10b981';
    }
}

// Redraw all 2D charts and updates WebGL objects
function drawAll() {
    drawNilssonDiagram(el.nilssonCanvas, state);
    drawPlateauDiagram(el.plateauCanvas, state);
    updateNucleus3D(state.beta, state.gamma_shape);
}

// Debounce helper
function debounce(func, wait) {
    let timeout;
    return function(...args) {
        clearTimeout(timeout);
        timeout = setTimeout(() => func.apply(this, args), wait);
    };
}

// Load script on DOM ready (with fallback if DOM is already interactive/complete)
if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
} else {
    init();
}
