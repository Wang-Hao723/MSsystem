# Browser drawing dependencies

These are local copies of the exact versions previously loaded by index.html from external CDNs. No dependency version was upgraded.

| Bundle | Version | Package |
| --- | --- | --- |
| SmilesDrawer | 2.0.3 | smiles-drawer |
| Plotly | 2.35.2 | plotly.js-dist-min |
| Cytoscape | 3.28.1 | cytoscape |

VISUAL_DEPENDENCIES.json records the official npm tarball URLs, verified SHA-512 integrity values and extracted bundle SHA-256 hashes. Each package's original license is included beside its bundle. These third-party notices are separate from the project's prepared software license.

Molecule cards render their saved candidate SMILES locally. When the browser renderer is unavailable, a saved-result-index structure image is requested from the local API. Spectrum panels draw saved measured peaks; a mirror comparison requires saved reference peaks and peak-pair evidence.
