# GIZMO biochemistry substrate (v1)

`graph.json` is the 38,148-node merged biochemistry graph that defines the coordinate system used by GIZMO Paper 1.

## Composition

- **16,343** gene nodes
- **6,406** metabolite nodes
- **15,399** reaction nodes (reified as full nodes, not pass-through edges; validated by the IDH-glioma contraction test — Manuscript §1)
- ancillary nodes: pathway, disease, drug, etc.

## Source databases

- **Reactome** (15,399 reactions, 6,406 metabolites, 2,382 pathway nodes)
- **StringDB** (PPI edges at confidence ≥ 700)
- **HMDB** (metabolite identifiers and synonyms)
- **KEGG** (reaction directions and EC-class annotations)

## Format

`graph.json` is a NetworkX node-link JSON. Load via:

```python
from gizmo.export.json_export import read_json
mg = read_json("substrate/graph.json")
# mg.graph is the underlying nx.Graph
```

The biochem subgraph used for propagation (`biochem_subgraph(mg, hub_cap=200)` from `gizmo`) caps high-degree promiscuous nodes (water, ATP, NADH) at 200 incoming edges to prevent Laplacian dominance.

## License

CC-BY 4.0. When using the substrate, cite **GIZMO Paper 1** (Zenodo DOI on acceptance) AND the four source databases — see `../LICENSE` for upstream license details.
