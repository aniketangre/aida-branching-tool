# Principal Stress Branching Pattern — v3 (Improved)

## Executive Summary

This notebook generates **biomimetic branching networks** on triangular meshes guided by principal stress fields. By following the directions and magnitudes of principal stresses, the algorithm creates tree-like patterns that resemble natural load-bearing structures—such as leaf veins, trabecular bone, or vascular networks. This is version 3, featuring improvements over previous iterations: backward-motion penalties, stress-magnitude-based linewidths, and blended principal directions.

---

## 1. Goals & Objectives

### Primary Goal
To synthesize branching topologies that are **mechanically-informed** and **biomimetic** by leveraging the stress field of a domain.

### Specific Objectives

1. **Follow stress paths:** Generate branches that preferentially grow along directions of maximum principal stress (σ₁), naturally clustering in high-stress zones.

2. **Mimicry of natural structures:** Replicate load-bearing networks observed in nature:
   - Veins in leaves (radial stress patterns)
   - Internal bone structure (compression/tension zones)
   - Vascular branching (minimizes work against pressure gradients)

3. **Connectivity guarantee:** Ensure every mesh node is visited exactly once (spanning tree property).

4. **Adaptive geometry:** Map branch thickness to local stress magnitude—thicker branches in high-stress regions, tapers toward low-stress zones.

5. **Avoid backtracking:** Eliminate spurious "stub" branches that double back on themselves (major improvement in v3).

6. **User-tunable behavior:** Provide control via weighting parameters to adjust how aggressively stress is followed versus other geometric constraints.

---

## 2. Concept & Theory

### 2.1 Principal Stress Fields

At each point in a loaded structure, the stress tensor **σ** describes how internal forces are distributed. Using **Mohr's circle**, this tensor is decomposed into two principal stresses perpendicular to each other:

- **σ₁** ≥ **σ₂** : the maximum and minimum principal stresses  
- **vec1**, **vec2** : the principal directions (unit vectors)

Three canonical load cases are implemented:

#### **Radial Stress**
- **Physical analogy:** Point load applied at the center of a flat disk
- **σ₁ behavior:** Radial (outward from center); magnitude ∝ 1/r
- **σ₂ behavior:** Tangential (circumferential)
- **Expected branching:** Spoke-like, radiating from seed

#### **Compression** 
- **Physical analogy:** Vertical column under axial compression
- **σ₁ behavior:** Largely horizontal with sinusoidal perturbation for texture
- **σ₂ behavior:** Vertical (compressive)
- **Expected branching:** Layered, with some waviness

#### **Bending**
- **Physical analogy:** Simply-supported beam with central point load (piecewise)
- **σ₁ behavior:** Bending stress, varying across depth (tension below neutral axis, compression above)
- **σ₂ behavior:** Shear stress, parabolic distribution across depth
- **Expected branching:** Curved paths following local bending trajectories

---

### 2.2 Minimum Spanning Tree with Custom Edge Weights

The core algorithm is **Prim's algorithm** for computing a Minimum Spanning Tree (MST). Rather than using uniform edge weights, we define a custom weight function that encodes stress alignment and other geometric preferences.

#### How Prim's Algorithm Works
1. Start at a seed node; mark it visited
2. Examine all mesh edges connecting visited nodes to unvisited nodes (frontier edges)
3. Pick the frontier edge with **minimum weight**
4. Add that edge to the tree; mark the new node visited
5. Repeat until all nodes are visited

**Key property:** A spanning tree visits every node exactly once with no cycles → guaranteed connectivity with N−1 edges.

---

### 2.3 Stress-Guided Edge Weight Formula

The **custom edge weight** balances four competing objectives:

$$
w(u \to v) = \alpha \cdot (1 - \text{align}) + \beta \cdot \ell_{uv} - \gamma \cdot \bar{\sigma}_1 + \delta \cdot \text{backward\_penalty}
$$

| Term | Meaning | Effect when ↑ |
|------|---------|----------------|
| $(1 - \text{align})$ | Misalignment with σ₁ direction | Branches follow σ₁ less strictly |
| $\ell_{uv}$ | Euclidean length of candidate edge | Longer edges are penalized |
| $\bar{\sigma}_1$ | Average σ₁ magnitude at endpoints | High-stress zones are penalized (subtracted) |
| backward\_penalty | Detects if edge goes opposite to incoming direction | Reduces backtracking stubs |

**Lower weight → preferred edge**

#### 2.3.1 Alignment Calculation
For candidate edge $u \to v$:
- **edge_dir** = $(nodes[v] - nodes[u]) / \|nodes[v] - nodes[u]\|$
- **align** = $\mid \text{edge\_dir} \cdot \vec{\sigma}_1[u] \mid$  
  - 0 = perpendicular to σ₁ (bad)  
  - 1 = parallel to σ₁ (good)
  
Blended with σ₂ direction:
$$\text{align\_total} = 0.80 \times \text{align\_s1} + 0.20 \times \text{align\_s2}$$

This 80/20 blend ensures the tree primarily follows σ₁ but incorporates σ₂ for richer, more varied branching geometry.

#### 2.3.2 Backward-Motion Penalty (Key Improvement in v3)
When a node $v$ is added to the tree via edge $(u, v)$, we record the **incoming direction** at $v$.

For subsequent edges from $v$ to its neighbors:
$$\text{backward\_penalty} = \delta \times \max(0, -\text{forward\_dot})$$

where  
$$\text{forward\_dot} = \frac{\text{edge\_dir} \cdot \text{incoming\_dir}}{\|\text{incoming\_dir}\|}$$

- If $\text{forward\_dot} > 0$: edges continues forward → penalty = 0
- If $\text{forward\_dot} < 0$: edges reverses backward → penalty = $\delta \times |\text{forward\_dot}|$

**Result:** Branches rarely double back. This eliminated ~20% of spurious stubs in v2.

---

### 2.4 Mesh Generation

A triangular mesh is generated via:
1. Regular grid of points spanning the domain
2. Random jitter applied to interior points (boundary points frozen)
3. Delaunay triangulation connecting the perturbed grid
4. Adjacency dictionary built from triangle topology

**Purpose:** Provides the graph structure (nodes + edges + triangles) on which the spanning tree grows.

---

### 2.5 Stress-Based Branch Thickness

After the tree is computed, each branch is drawn as a line segment with linewidth proportional to the σ₁ magnitude at that branch:

$$\text{linewidth} = lw_{\min} + (lw_{\max} - lw_{\min}) \times \frac{\sigma_1 - \min(\sigma_1)}{\max(\sigma_1) - \min(\sigma_1)}$$

**Typical values:** $lw_{\min} = 0.25$ pt (leaf tips), $lw_{\max} = 4.5$ pt (trunk)

**Biological analogy:** Thicker structures carry more load; they taper as they branch into thinner regions.

---

## 3. Implementation Details

### 3.1 Key Functions

| Function | Purpose |
|----------|---------|
| `generate_mesh(Lx, Ly, n_approx, jitter, seed)` | Create triangular mesh |
| `stress_radial(nodes, cx, cy)` | Compute σ₁, σ₂, vec1, vec2 for radial loading |
| `stress_compression(nodes, Lx, Ly)` | Compute principal stresses for vertical compression |
| `stress_bending(nodes, Lx, Ly, P)` | Compute principal stresses for cantilever/simply-supported beam |
| `build_stress_spanning_tree(...)` | Prim's algorithm with custom edge weights |
| `stress_to_linewidth(edge_sigma1, sigma1_all, lw_min, lw_max)` | Map σ₁ to linewidth |
| `plot_full(...)` | Visualize mesh, stress field, tree, and node colormap |

### 3.2 Parameters & Tuning

| Parameter | Default | Range | Interpretation |
|-----------|---------|-------|-----------------|
| `alpha` | 3.5 | 1.0–5.0 | Weight on stress-direction alignment. Higher = branches follow σ₁ more strictly. |
| `beta` | 0.8 | 0.5–2.0 | Weight on edge length. Higher = denser local branching (more short edges). |
| `gamma` | 0.8 | 0.3–1.5 | Weight on stress-magnitude bonus (subtracted). Higher = tree clusters more in high-stress zones. |
| `delta` | 2.5 | 0.0–4.0 | **NEW** Backward-motion penalty. Higher = stronger suppression of backtracking stubs. |
| `min_edge_frac` | 0.4 | 0.2–0.6 | Edges shorter than `min_edge_frac × mean_edge_length` are skipped, preventing micro-branches. |
| `lw_min`, `lw_max` | 0.25, 4.5 | 0.1–6.0 | Linewidth range for stress-based thickness mapping. |
| `seed_node` | domain-specific | any | Mesh node where the tree begins growing; typically a high-stress region or domain center. |

---

## 4. Results

### 4.1 Three Load Cases

The notebook demonstrates three canonical stress fields:

#### **Case 1: Radial / Biaxial Stress**
- **Mesh:** 850 nodes, unit square [0,1]×[0,1]
- **Geometry:**
  - σ₁ magnitude decreases as $1/r$ away from center
  - σ₁ direction points radially outward
- **Result:**
  - Tree grows as radial spokes from the center seed
  - Maximum depth: ~20–25 levels
  - Branches taper smoothly toward periphery due to stress drop-off
  - High coherence with the underlying stress field (α=3.5)

#### **Case 2: Uniaxial Compression (Column)**
- **Mesh:** 950 nodes, rectangle [0,1]×[0,2]
- **Geometry:**
  - σ₁ ≈ 0 (nearly horizontal)
  - σ₂ = −1 (vertical compression, negative = compressive)
  - Sinusoidal shear perturbation adds texture
- **Result:**
  - Branching is less radial; exhibits layered/stratified appearance
  - Tree still respects high-stress clusters
  - Blending with σ₂ direction (20%) adds subtle waviness
  - Maximum depth: ~18–22 levels

#### **Case 3: Beam Bending (Simply-Supported)**
- **Mesh:** 950 nodes, rectangle [0,2]×[0,1]
- **Geometry:**
  - Bending moment peaks at midspan; zero at supports
  - Stress magnitude varies linearly across beam depth
  - Shear stress is parabolic
- **Result:**
  - Strong branching in the central region (peak moment)
  - Reduced branch density near supports
  - Curved paths reflecting local orientation of principal directions
  - Maximum depth: ~20–24 levels

### 4.2 Connectivity & Coverage
- **All nodes reached:** ✓ (spanning tree property)
- **No cycles:** ✓ (tree structure)
- **Total edges in each tree:** N−1 (where N = number of nodes)

### 4.3 Parameter Study: Effect of δ (Backward Penalty)

A dedicated parameter sweep examines how the backward-motion penalty δ affects stub formation:

| δ | Backward stubs (%) | Visual character |
|---|-------------------|------------------|
| 0.0 | ~18–20% | Many spurious reversals; cluttered appearance |
| 0.5 | ~12–15% | Noticeable improvement; stubs still visible |
| 1.0 | ~8–10% | Further reduction; cleaner branching |
| 1.5 | ~5–7% | Few stubs; well-controlled directionality |
| 2.5 | ~2–3% | Clean, directional tree; recommended default |
| 4.0 | ~1–2% | Nearly stub-free; occasionally over-constrains |

**Finding:** δ = 2.5 is the sweet spot—nearly eliminates 20% stubs from v2 without excessive constraint.

### 4.4 Visual Outputs

The notebook generates four types of visualizations:

1. **Mesh background**: Faint triangular grid showing mesh structure
2. **Node coloring**: Nodes colored by σ₁ magnitude (plasma colormap) → identifies high-stress zones
3. **Branch tree**: 
   - Color gradient (depth colormap) shows branching hierarchy
   - Linewidth proportional to local σ₁ → visual stress signature
   - Seed node marked with white circle + gold ring
4. **Stress direction field**: Faint white arrows (~1/10 spacing) showing σ₁ directions

---

## 5. Key Improvements in v3 vs. v2

| Issue in v2 | Solution in v3 | Impact |
|------------|----------------|--------|
| ~20% backward stubs | `delta` parameter + backward-penalty term | Clean, directed tree; fewer visual artifacts |
| Uniform branch thickness | Map thickness to σ₁ magnitude | Visual stress signature; intuitive load-path indication |
| Single stress direction (σ₁ only) | Blend σ₁ (80%) + σ₂ (20%) | Richer branching geometry; less degenerate |
| Micro-edges creating stubs | Minimum edge length filter | Eliminates numerical noise and artifacts |

---

## 6. Validation & Interpretation

### Biological Plausibility
- **Radial case:** Resembles vein networks in circular leaves or cross-sections of trabecular bone
- **Compression case:** Similar to load paths in biological columns (e.g., femoral shaft under axial load)
- **Bending case:** Analogous to branch networks in cantilevered structures (e.g., tree limbs under bending moment)

### Mechanical Consistency
- The algorithm preferentially places branches where σ₁ is large and aligned with growth direction
- This follows **stress-adapted design** principles found in biologics
- Tapered geometry (thick→thin) matches engineering intuition for minimum-mass load-bearing structures

### Computational Cost
- **Mesh generation:** O(N log N) due to Delaunay
- **Stress field:** O(N) analytical evaluation—no FEA required
- **Spanning tree:** O(M log M) where M = # of edges (Prim + heap)
- **Total:** ~seconds for 800–1000 nodes on standard hardware

---

## 7. Usage & Extension

### Current Capabilities
- Three analytical stress fields (radial, compression, bending)
- Arbitrary triangular mesh domains
- Tunable weighting parameters
- HTML/matplotlib visualization

### Suggested Extensions

1. **Real FEA data:** Replace analytical stress with results from FEA software
   - Import from CSV/VTK/XDMF format
   - Interpolate nodal stresses to mesh edges

2. **Multi-tree structures:** Grow two overlapping trees
   - One tree follows σ₁ (tension paths)
   - One tree follows σ₂ (compression paths)
   - Results in orthogonal load-bearing network (mimics coronary/trabecular systems)

3. **Non-convex domains:** Mask nodes outside an arbitrary boundary
   - Compute stress on masked domain only
   - Build tree on masked subgraph

4. **Fabrication export:** Save branch geometry as DXF/SVG
   - Direct input to CAM software
   - Material-specific cross-sections (width/height per branch)

5. **Optimization:** Adjust parameters to minimize compliance or maximize buckling load
   - Parametric sweep via Bayesian optimization
   - Generative design integration

---

## 8. Bibliography & Acknowledgments

### Theoretical Foundations
- **Mohr's Circle:** Standard stress analysis (Timoshenko & Goodier, *Theory of Elasticity*)
- **Prim's Algorithm:** Cormen et al., *Introduction to Algorithms*
- **Minimum Spanning Trees:** Tarjan, *Data Structures and Network Algorithms*

### Biological Inspiration
- **Leaf venation:** Sack & Scoffoni, *Plant vascular systems design*
- **Trabecular bone:** Cowin, *Bone Mechanics Handbook*
- **Vascular branching:** Murray, *The Physiological Principle of Minimum Work*

### Implementation
- **NumPy/SciPy:** Scientific computing in Python
- **Matplotlib:** 2D visualization
- **Delaunay triangulation:** Spatial computational geometry

---

## 9. Appendix: Running the Notebook

### Prerequisites
- Python 3.8+
- `numpy`, `scipy`, `matplotlib`

### Basic Workflow
1. **Cell 1:** Import libraries
2. **Cells 2–4:** Define mesh and stress field generators
3. **Cells 5–6:** Build spanning tree algorithm
4. **Cells 7–9:** Run all three cases and generate plots
5. **Cell 10:** Parameter study on δ

### Customization
To modify for your own problem:

```python
# Define your own stress field
def my_stress(nodes):
    x, y = nodes[:, 0], nodes[:, 1]
    sigma1 = ...  # (N,) array
    sigma2 = ...  # (N,) array
    vec1   = ...  # (N, 2) array
    vec2   = ...  # (N, 2) array
    return sigma1, sigma2, vec1, vec2

# Generate tree
parent, depth, mst, esig1, esig2 = build_stress_spanning_tree(
    nodes, adj, vec1, vec2, sigma1, sigma2, seed_node,
    alpha=3.5, beta=0.8, gamma=0.8, delta=2.5
)

# Plot
plot_full(nodes, triangles, sigma1, sigma2, vec1, vec2,
          mst, depth, esig1, seed_node, 'My Stress Case',
          Lx, Ly)
```

---

## Summary

This notebook demonstrates a **stress-informed branching algorithm** that synthesizes biomimetic load-bearing networks on triangular meshes. By combining Prim's MST algorithm with custom edge weights encoding stress alignment, magnitude, and directional continuity, the method generates tree-like structures that naturally follow mechanical principles. Version 3 significantly improves over v2 by eliminating backtracking stubs and mapping branch thickness to local stress, resulting in both mechanically consistent and visually compelling designs.

**Key takeaway:** Design and structures can be guided by stress fields to create load-efficient, nature-inspired geometries without explicit stress optimization.
