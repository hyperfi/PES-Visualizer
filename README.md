# Strutinsky Shell Correction & Potential Energy Surfaces (PES) Visualizer

An interactive, responsive 3D web application designed to demonstrate the **Microscopic-Macroscopic (Micro-Macro) approach** to nuclear deformation energy. This tool calculates Nilsson single-particle levels, performs Strutinsky shell corrections, computes Liquid Drop Model (LDM) macroscopic energies, and projects the resulting Potential Energy Surface (PES) in real-time.

---

## ⚛️ Physical Background

The total energy of a deformed nucleus is calculated as a sum of macroscopic bulk energy and microscopic quantum corrections:

$$E_{\text{total}}(\beta, \gamma) = E_{\text{LDM, rel}}(\beta, \gamma) + \delta E_{\text{shell}, p}(\beta, \gamma) + \delta E_{\text{shell}, n}(\beta, \gamma)$$

### 1. The Nilsson Model (Microscopic)
In a deformed nucleus, the spherical shell model degeneracies are broken. The Nilsson Model describes how single-particle states split under quadrupole deformation magnitude ($\beta$) and triaxiality asymmetry ($\gamma$). The Hamiltonian solved is:

$$H = H_0 - \frac{2}{3}\beta \hbar\omega_0 r^2 \left[ \cos\gamma Y_{20} + \frac{\sin\gamma}{\sqrt{2}} (Y_{22} + Y_{2,-2}) \right] + C \mathbf{l}\cdot\mathbf{s} + D \mathbf{l}^2$$

### 2. Strutinsky Shell Correction (Quantum Correction)
Extracts the quantum mechanical shell fluctuation energy by subtracting a smooth average energy level density $\tilde{E}$ from the discrete sum of occupied energies:

$$\delta E_{\text{shell}} = \sum_{i=1}^{A} e_i - \tilde{E}(\gamma_{\text{smear}}, p)$$

To verify that the calculated correction is physically sound, the application plots $\delta E_{\text{shell}}$ vs. the smearing width $\gamma_{\text{smear}}$ to visualize the **Strutinsky Plateau Region** where $\frac{\partial (\delta E)}{\partial \gamma_{\text{smear}}} \approx 0$.

### 3. Liquid Drop Model (Macroscopic)
Describes the bulk properties of the nucleus (volume energy, surface tension, and Coulomb repulsion) as a function of deformation.

---

## 🌟 Interactive Features

* **Real-time Level Splitting**: Visualize proton and neutron Nilsson orbitals crossing as you adjust deformation. A dynamic 3D WebGL nucleus shape preview displays the resulting spheroid (prolate, oblate, or triaxial).
* **Strutinsky Plateau Plot**: Interactively adjust the smearing width and correction order to verify the plateau condition.
* **Energy Surface Exploration**: Render the PES landscape in interactive **3D WebGL** or switch to a **2D Polar Heatmap** (Lund coordinate projection).
* **Ground State Search**: The application evaluates the entire PES grid in real-time to find and label the predicted stable ground-state shape ($\beta, \gamma$) marked with a gold star.
* **Summary Metrics & Tooltips**: Tap or hover over any metric card at the top to display a descriptive pop-up explaining the physical calculations behind the displayed numbers.

---

## 📱 Mobile & Touch Optimizations

The application features a responsive design layout optimized for mobile viewports (phones and tablets) and touch targets:
* **Off-Canvas controls drawer**: On screens <= 1024px, the sidebar is converted into a side-sliding parameters drawer, triggered by a `⚙️ Parameters` button in the header.
* **Smart auto-close**: Choosing an isotope preset automatically dismisses the drawer, updating the charts in full view.
* **Collapsible Theory cards**: Long explanatory text blocks are wrapped in native details accordions that collapse by default on mobile.
* **Screen-constrained tooltips**: Pop-up tooltips project downward below the metrics cards and automatically align to the inner viewport margins, preventing any text truncation.
* **Large Touch Targets**: Custom-styled range slider inputs feature expanded touch thumbs for precise fingers-on shape modifications.

---

## 🛠️ Technology Stack

* **Structure**: HTML5 Semantic markup
* **Styling**: Modern, fluid CSS3 with custom variables, glassmorphic panels, and flexbox/grid responsive layouts
* **Logic**: Vanilla ES Modules (JavaScript)
* **3D Rendering**: [Three.js](https://threejs.org/) & [OrbitControls](https://threejs.org/docs/#examples/en/controls/OrbitControls)
* **Math Rendering**: [MathJax v3](https://www.mathjax.org/) (LaTeX formatting)

---

## 🚀 How to Run Locally

Because the application uses native **ES Modules** (`import`/`export`), opening the `index.html` file directly in the browser via `file://` protocol will trigger CORS restrictions. 

You must serve the files using a local web server.

### Options:
1. **VS Code Live Server**: Right-click `index.html` and select **Open with Live Server**.
2. **Python HTTP Server**:
   ```bash
   # Python 3
   python -m http.server 8000
   ```
   Open `http://localhost:8000` in your web browser.
3. **NPM Static Server**:
   ```bash
   npx serve .
   ```
   Open the address provided in your terminal.

---

*Developed by **Dr. Abhishek** ([www.dr-abhishek.com](https://www.dr-abhishek.com/))*
