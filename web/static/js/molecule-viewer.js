/**
 * MoleculeViewer — SmilesDrawer molecular structure rendering
 * Lightweight, client-side 2D molecule display from SMILES strings
 */
class MoleculeViewer {
    constructor(options = {}) {
        this.defaultWidth = options.width || 280;
        this.defaultHeight = options.height || 220;
        this.drawer = null;
        this._initDrawer();
    }

    _initDrawer() {
        if (typeof SmilesDrawer === 'undefined') {
            console.warn('SmilesDrawer not loaded. Include CDN: https://unpkg.com/smiles-drawer@2.0.3/dist/smiles-drawer.min.js');
            return;
        }
        this.drawer = new SmilesDrawer.Drawer({
            width: this.defaultWidth,
            height: this.defaultHeight,
            bondThickness: 1.2,
            bondLength: 18,
            shortBondLength: 0.85,
            bondSpacing: 4.5,
            atomVisualization: 'default',
            isomeric: true,
            debug: false,
            terminalCarbons: false,
            explicitHydrogens: false,
            overlapSensitivity: 0.42,
            overlapResolutionIterations: 1,
            compactDrawing: true,
            fontSizeLarge: 11,
            fontSizeSmall: 4,
            padding: 12,
        });
    }

    /**
     * Render a molecule to a canvas element
     * @param {string} canvasId - Canvas element ID
     * @param {string} smiles - SMILES string
     */
    render(canvasId, smiles) {
        if (!this.drawer || !smiles) return;
        const canvas = document.getElementById(canvasId);
        if (!canvas) return;

        SmilesDrawer.parse(smiles, (tree) => {
            try {
                this.drawer.draw(tree, canvas, 'light', false);
            } catch (e) {
                this._drawFallback(canvas, smiles);
            }
        }, (err) => {
            this._drawFallback(canvas, smiles);
        });
    }

    /**
     * Render molecule to a new canvas and return as data URL
     * @param {string} smiles - SMILES string
     * @returns {string} data URL or empty string
     */
    renderToDataURL(smiles) {
        if (!this.drawer || !smiles) return '';
        const canvas = document.createElement('canvas');
        canvas.width = this.defaultWidth;
        canvas.height = this.defaultHeight;
        SmilesDrawer.parse(smiles, (tree) => {
            try { this.drawer.draw(tree, canvas, 'light', false); } catch (e) {}
        }, () => {});
        return canvas.toDataURL('image/png');
    }

    /**
     * Batch render molecules into a grid container
     * @param {string} containerId - Container element ID
     * @param {Array} molecules - [{smiles, name, ...}]
     * @param {Object} options - {columns}
     */
    renderGrid(containerId, molecules, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) return;
        const minCardWidth = options.minCardWidth || 190;
        container.innerHTML = '';
        container.classList.add('molecule-gallery');
        container.style.setProperty('--molecule-card-min', `${minCardWidth}px`);

        molecules.forEach((mol, idx) => {
            if (!mol.smiles) return;
            const card = document.createElement('div');
            card.className = 'molecule-card';
            card.innerHTML = `
                <div class="molecule-canvas-wrap"><canvas id="mol-canvas-${idx}" width="${this.defaultWidth}" height="${this.defaultHeight}"></canvas></div>
                <div class="molecule-label">${mol.name || mol.smiles.substring(0, 30)}</div>
            `;
            container.appendChild(card);
            setTimeout(() => this.render(`mol-canvas-${idx}`, mol.smiles), 10);
        });
    }

    _drawFallback(canvas, smiles) {
        const ctx = canvas.getContext('2d');
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        ctx.fillStyle = '#94a3b8';
        ctx.font = '12px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('Invalid SMILES', canvas.width / 2, canvas.height / 2);
    }
}

// Global instance
window.MoleculeViewer = MoleculeViewer;
