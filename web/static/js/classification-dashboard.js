/**
 * ClassificationDashboard — Classification metrics visualization
 * Per-class accuracy bars, confusion heatmap, confidence distribution
 */
class ClassificationDashboard {
    constructor() {
        this.classColors = {
            flavonoid: '#f59e0b', alkaloid: '#8b5cf6', terpenoid: '#ef4444',
            steroid: '#f97316', saponin: '#22c55e', phenolic_acid: '#06b6d4',
            coumarin: '#ec4899', lignan: '#6366f1', quinone: '#14b8a6',
            sugar: '#a3e635', fatty_acyl: '#fb923c', lipid: '#818cf8',
            amino_acid: '#2dd4bf', other: '#94a3b8',
        };
    }

    /**
     * Render per-class accuracy bar chart
     * @param {string} containerId
     * @param {Object} classAccuracy - {className: accuracy, ...}
     */
    renderClassAccuracy(containerId, classAccuracy) {
        const container = document.getElementById(containerId);
        if (!container) return;

        const sorted = Object.entries(classAccuracy).sort((a, b) => b[1] - a[1]);
        const names = sorted.map(([n]) => this._formatName(n));
        const values = sorted.map(([, v]) => v);
        const colors = sorted.map(([n]) => this.classColors[n] || '#6b7280');

        // Target line at 85%
        const targetLine = {
            type: 'line', xref: 'paper', x0: 0, x1: 1,
            y0: 0.85, y1: 0.85,
            line: { color: '#22c55e', width: 2, dash: 'dash' },
        };

        Plotly.newPlot(container, [{
            x: names, y: values,
            type: 'bar',
            marker: { color: colors, opacity: 0.85, line: { width: 1, color: '#fff' } },
            text: values.map(v => (v * 100).toFixed(1) + '%'),
            textposition: 'outside',
            textfont: { size: 10 },
            hovertemplate: '%{x}: %{y:.1%}<extra></extra>',
        }], {
            yaxis: { title: 'Accuracy', range: [0, 1.15], tickformat: '.0%', showgrid: true, gridcolor: '#f0f0f0' },
            xaxis: { title: '', tickangle: -45, tickfont: { size: 10 } },
            shapes: [targetLine],
            annotations: [{
                x: 1, y: 0.87, xref: 'paper', yref: 'y',
                text: 'Target 85%', showarrow: false,
                font: { size: 10, color: '#22c55e' },
            }],
            margin: { t: 20, b: 100, l: 50, r: 20 },
            height: 380,
            paper_bgcolor: 'white',
        }, { responsive: true, modeBarButtonsToAdd: ['toImage'] });
    }

    /**
     * Render confusion matrix heatmap
     * @param {string} containerId
     * @param {Array} matrix - 2D array [true][pred] normalized
     * @param {Array} classNames
     */
    renderConfusionMatrix(containerId, matrix, classNames) {
        const container = document.getElementById(containerId);
        if (!container) return;

        const names = classNames.map(n => this._formatName(n));
        const textLabels = matrix.map(row => row.map(v => (v * 100).toFixed(1) + '%'));

        Plotly.newPlot(container, [{
            z: matrix, x: names, y: names,
            type: 'heatmap',
            colorscale: [[0, '#f0f9ff'], [0.3, '#7dd3fc'], [0.6, '#0284c7'], [1, '#1e3a5f']],
            text: textLabels,
            texttemplate: '%{text}',
            textfont: { size: 8 },
            hovertemplate: 'True: %{y}<br>Pred: %{x}<br>%{z:.1%}<extra></extra>',
            xgap: 2, ygap: 2,
        }], {
            xaxis: { title: 'Predicted', side: 'bottom', tickangle: -45, tickfont: { size: 9 } },
            yaxis: { title: 'True', autorange: 'reversed', tickfont: { size: 9 } },
            margin: { t: 20, b: 100, l: 100, r: 20 },
            height: 500,
        }, { responsive: true, modeBarButtonsToAdd: ['toImage'] });
    }

    /**
     * Render confidence distribution as box plot
     * @param {string} containerId
     * @param {Object} confidences - {className: [confidence_values], ...}
     */
    renderConfidenceDistribution(containerId, confidences) {
        const container = document.getElementById(containerId);
        if (!container) return;

        const traces = Object.entries(confidences).map(([name, vals]) => ({
            y: vals,
            type: 'box',
            name: this._formatName(name),
            marker: { color: this.classColors[name] || '#6b7280' },
            boxmean: 'sd',
        }));

        Plotly.newPlot(container, traces, {
            yaxis: { title: 'Confidence', range: [0, 1], tickformat: '.0%' },
            xaxis: { tickangle: -45, tickfont: { size: 10 } },
            margin: { t: 20, b: 100, l: 50, r: 20 },
            height: 380,
            showlegend: false,
        }, { responsive: true });
    }

    /**
     * Render all three charts for a task's classification results
     * @param {string} taskResults - results array from API
     */
    renderFullDashboard(taskResults) {
        // Calculate class accuracy
        const classCorrect = {};
        const classTotal = {};
        taskResults.forEach(r => {
            const cls = r.classification || r.Classification || 'Unknown';
            if (!classTotal[cls]) { classTotal[cls] = 0; classCorrect[cls] = 0; }
            classTotal[cls]++;
            if (r.confidence && r.confidence > 0.5) classCorrect[cls]++;
        });

        const classAccuracy = {};
        Object.keys(classTotal).forEach(cls => {
            classAccuracy[cls] = classCorrect[cls] / classTotal[cls];
        });

        this.renderClassAccuracy('class-accuracy-chart', classAccuracy);
    }

    _formatName(name) {
        return name.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase());
    }
}

window.ClassificationDashboard = ClassificationDashboard;
