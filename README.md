# An Approximate Solution to the Minimum Vertex Cover Problem: The Hvala Algorithm

**Author:** Frank Vega  
**Affiliation:** Information Physics Institute, 840 W 67th St, Hialeah, FL 33012, USA  
**Email:** vega.frank@gmail.com  
**ORCID:** 0000-0001-8210-4126

---

## Abstract

We present the **Hvala** algorithm, an ensemble approximation method for the Minimum Vertex Cover problem. Hvala combines, in a component-wise minimum-selection scheme, *five* complementary heuristics: (i) a minimum weighted dominating set on a degree-1 reduction, (ii) a minimum weighted vertex cover on the same degree-1 reduction, (iii) the NetworkX local-ratio 2-approximation, (iv) a maximum-degree greedy, and (v) a minimum-to-minimum heuristic.

**Relation to a companion paper.** The Hallelujah algorithm studies in depth *one single* heuristic — the degree-1 reduction in minimum weighted vertex cover form — which corresponds to component (ii) of Hvala. The present work situates that single heuristic inside a broader compendium and analyses how the other four heuristics cover each other's worst cases.

**Empirical confirmation.** Across 201+ diverse instances from three independent experimental studies (Resistire real-world networks up to 262,111 vertices; Creo NPBench/DIMACS-complement hard instances; Gemini–Vega AI-validated stress tests on 3-regular graphs up to 20,000 vertices), Hvala attains approximation ratios in the range 1.001–1.071, with no observed instance exceeding 1.071.

**Theoretical analysis.** We prove optimality on specific graph classes — paths and trees (Min-to-Min), cliques and regular graphs (max-degree greedy), skewed bipartite graphs (degree-1 reduction), and hub-heavy/star graphs (degree-1 reduction). We argue *structural orthogonality*: the pathological worst case of each heuristic is precisely a graph family on which at least one other heuristic attains optimality. This mutual worst-case coverage is the central structural reason why the component-wise minimum over the five heuristics is expected to remain well below the $\sqrt{2}$ hardness threshold.

**Strong Hypothesis.** We explicitly propose, as a Strong Hypothesis, that the Hvala ensemble achieves approximation ratio $\rho < \sqrt{2}$ on *every* graph. Under the Strong Exponential Time Hypothesis (SETH), combined with the hardness result of Khot, Minzer and Safra, such a polynomial-time algorithm would imply **P = NP**. We present this as a hypothesis, clearly labelled as such, supported by (a) the Hallelujah algorithm for its core degree-1 component, (b) structural orthogonality of the five heuristics, and (c) uniform empirical behaviour across 201+ instances. Hvala runs in $\mathcal{O}(m \log n)$ time, $\mathcal{O}(m)$ space, and is publicly available via PyPI as the `hvala` package.

**Keywords:** Vertex Cover; Approximation Algorithm; Computational Complexity; P versus NP; Graph Optimization; Hardness of Approximation; Ensemble Heuristic

**MSC:** 05C69, 68Q25, 90C27, 68W25

---

## 1. Introduction

The **Minimum Vertex Cover** problem is one of the most fundamental problems in combinatorial optimization and theoretical computer science. For an undirected graph $G = (V, E)$, the problem asks for the smallest subset $S \subseteq V$ such that every edge $(u, v) \in E$ has at least one endpoint in $S$. Despite its conceptual simplicity, the problem underpins applications in wireless network design, bioinformatics, scheduling, and VLSI circuit optimization.

Karp established NP-completeness of the problem in 1972 [[karp2009reducibility]](#references). Unless P = NP — one of the most profound open questions in mathematics and computer science — no polynomial-time algorithm can compute exact minimum vertex covers for general graphs. This limitation has driven decades of work on approximation algorithms.

Classical approximation results include the well-known 2-approximation based on maximal matching [[papadimitriou1998combinatorial]](#references), and refinements by Karakostas [[karakostas2009better]](#references) and Karpinski et al. [[karpinski1996approximating]](#references) achieving factors $2 - \epsilon$ for small $\epsilon > 0$ via LP relaxations and primal-dual techniques.

These advances confront fundamental theoretical barriers. Dinur and Safra [[dinur2005hardness]](#references), via the PCP theorem, proved that no polynomial-time algorithm achieves ratio better than 1.3606 unless P = NP. Khot, Minzer and Safra [[khot2017independent,dinur2018towards,khot2018pseudorandom]](#references) strengthened this to $\sqrt{2} - \epsilon$ under SETH: achieving ratio $\rho < \sqrt{2}$ in polynomial time would directly imply **P = NP**. Under the Unique Games Conjecture [[khot2002unique]](#references), no constant-factor better than $2 - \epsilon$ is achievable [[khot2008vertex]](#references). Any polynomial-time algorithm achieving $\rho < \sqrt{2}$ would therefore resolve P versus NP, one of the seven Millennium Prize Problems.

### 1.1 Relation to the Accepted "Hallelujah Algorithm" Paper

The present article is a *companion and extension* of the author's forthcoming paper:

> **Frank Vega.** *An Approximate Solution to the Minimum Vertex Cover Problem: The Hallelujah Algorithm*. *International Journal of Parallel, Emergent and Distributed Systems*. See [[Vega26Hallelujah,Vega25HallelujahPreprint]](#references).

It is important to state explicitly the scope of each work:

- The **Hallelujah algorithm** [[Vega26Hallelujah]](#references) studies *one single heuristic*: the reduction of the input graph to maximum degree 1 with weights $1/k$ per auxiliary vertex, followed by an exact minimum weighted vertex cover on that reduced graph, and projection back to the original graph.
- The **Hvala algorithm** (this paper) is a *compendium* — an ensemble — of *five* heuristics. The Hallelujah heuristic is *one* of them (component $S_2$ in Algorithm 1). The remaining four heuristics (NetworkX local-ratio, maximum-degree greedy, minimum-to-minimum, and the minimum weighted dominating set on the same degree-1 reduction, which constitutes the distinct $S_1$ component) are engineered precisely to *cover the worst-case inputs* of one another.

This distinction is essential to the interpretation of the results below: the companion paper provides the rigorous theoretical core (for *one* heuristic), and the present paper argues that *wrapping* that heuristic inside a mutually-compensating ensemble is what enables the strong P = NP hypothesis articulated in Section 1.2.

### 1.2 The Strong Hypothesis

We organize every claim in this paper under one of four clearly labelled epistemic levels:

- **Strong Hypotheses** — claims that are *not proved* in this paper, whose truth would have major mathematical consequences.
- **Working (Stronger) Assumptions** — auxiliary assumptions we adopt to support a chain of reasoning, flagged explicitly so the reader can measure how much of the argument depends on them.
- **Empirical Confirmations** — facts observed on concrete benchmark suites, which *do not* constitute proofs but restrict the set of graphs on which a counterexample could live.
- **Suggestions** — plausible extensions or conjectures, presented as invitations to further work rather than claims.

The central claim of the article is:

> **Strong Hypothesis (P = NP via the Hvala Ensemble).** The Hvala algorithm achieves, for *every* finite undirected graph $G$, an approximation ratio $\rho(G) = |\text{Hvala}(G)| / \mathrm{OPT}(G) < \sqrt{2}$. Combined with the SETH-based hardness of Khot, Minzer and Safra [[khot2017independent,dinur2018towards,khot2018pseudorandom]](#references), the existence of such a polynomial-time algorithm would imply **P = NP**.

### 1.3 Algorithm Overview

The Hvala algorithm operates through the following phases:

**Phase 1: Graph reduction (degree-1 core).** Transform $G = (V, E)$ into a graph $G' = (V', E')$ of maximum degree 1: for each $u \in V$ of degree $k$, create $k$ auxiliary vertices $(u, 0), (u, 1), \ldots, (u, k-1)$, each connected to exactly one of $u$'s neighbours, with weight $w_{(u,i)} = 1/k$. The exact minimum weighted vertex cover computed on this reduced graph constitutes the Hallelujah heuristic analysed in the companion paper [[Vega26Hallelujah]](#references).

**Phase 2: Optimal solutions on the reduced graph.** Since $G'$ consists of disjoint edges and isolated vertices, minimum weighted vertex cover and minimum weighted dominating set both admit exact polynomial-time solutions, which are then projected back to $V$.

**Phase 3: Ensemble heuristics.** To cover the worst-case regimes of the reduction, Hvala additionally applies (1) NetworkX's local-ratio 2-approximation, (2) a maximum-degree greedy, and (3) a minimum-to-minimum heuristic.

**Phase 4: Component-wise selection.** Each connected component is processed independently; the smallest valid cover among the five candidates is chosen for that component. This is where the structural orthogonality (Section 4) pays off.

### 1.4 Experimental Validation Framework

The empirical claims in this paper rest on *three* independent experimental studies, conducted on standard hardware (Intel Core i7-1165G7, 32GB RAM), in Python 3.12.0 with NetworkX 3.4.2:

1. **Real-World Large Graphs — the Resistire Experiment** [[Vega25Resistire]](#references): 88 instances from the Network Data Repository [[RA15,LargeGraphs]](#references), up to 262,111 vertices.
2. **NPBench Hard Instances — the Creo Experiment** [[Vega25Creo]](#references): 113 challenging benchmarks including FRB and DIMACS clique complements [[NPBench]](#references). This is the *latest and most comprehensive* DIMACS-based evaluation and supersedes an earlier standalone DIMACS run [[Vega25Hvala]](#references) (which is therefore not reproduced separately here).
3. **AI-Validated Stress Tests — the Gemini–Vega Validation** [[Vega25Gemini]](#references): independent validation with Gemini AI on hard 3-regular graphs up to 20,000 vertices.

> **Empirical Confirmation (Uniform empirical ratio).** Across these 201+ instances, the approximation ratio of Hvala falls in the range 1.001–1.071, with no observed instance exceeding 1.071 — even on adversarially constructed hard 3-regular graphs. This is an empirical confirmation on the tested instances, not a proof over all graphs.

---

## 2. Related Work and State-of-the-Art

### 2.1 Theoretical Approximation Algorithms

The classical 2-approximation based on maximal matching [[papadimitriou1998combinatorial]](#references) remains the simplest and most widely used approach. Advanced techniques include:

- **Local-Ratio Methods** [[bar1985local]](#references): 2-approximation via iterative dual adjustments.
- **LP-Based Approaches** [[karakostas2009better]](#references): rounding schemes achieving $2 - \Theta(1/\sqrt{\log n})$.
- **Semidefinite Programming**: theoretical improvements to $(2 - \epsilon)$ with impractical constants.

### 2.2 Practical Heuristic Methods

Modern state-of-the-art heuristics achieve exceptional empirical performance via local search:

**TIVC** [[zhang2023tivc]](#references) — 3-improvement local search with tiny perturbations, empirical ratios $\sim$ 1.005 on DIMACS benchmarks.

**FastVC and variants** [[cai2017finding]](#references) — fast local search with pivoting and probing, ratios $\sim$ 1.02 with sub-second runtimes on million-vertex graphs.

**MetaVC2** [[luo2019local]](#references) — adaptive meta-heuristic combining tabu search, simulated annealing and genetic operators, ratios 1.01–1.05 across heterogeneous classes.

### 2.3 Fixed-Parameter Tractable Algorithms

For parameterization by solution size $k$, Harris and Narayanaswamy [[harris2024faster]](#references) achieve $\mathcal{O}(1.2738^k + kn)$ runtime, practical when $k$ is small relative to $n$.

### 2.4 Positioning of Our Work

Our contribution differs from the above in two crucial dimensions:

1. **Strong Hypothesis with a P = NP implication.** We explicitly articulate the Strong Hypothesis — a provable ratio $\rho < \sqrt{2}$ — and tie it to the SETH hardness barrier. We label it as a hypothesis rather than a theorem.
2. **Ensemble of mutually-compensating heuristics.** The design principle is structural orthogonality: each heuristic's worst case is another heuristic's best case. This is exactly what distinguishes Hvala (the compendium) from the Hallelujah single-heuristic algorithm [[Vega26Hallelujah]](#references).

---

## 3. The Hvala Algorithm: Detailed Description

### 3.1 Algorithm Structure and Pseudocode

#### Main Algorithm

**Algorithm 1:** Hvala: Main Algorithm  
**Input:** Undirected graph $G = (V, E)$  
**Output:** Approximate vertex cover $S \subseteq V$

```
1: if G is empty or |E| = 0 then
2:     return ∅
3: end if
4: Remove self-loops from G
5: Remove isolated vertices from G
6: S ← ∅
7: for each connected component C in G do
8:     G_C ← subgraph induced by C
9:     // Phase 1: Reduction to maximum degree-1 (Hallelujah core)
10:    G' ← ReduceToMaxDegree1(G_C)
11:    // Phase 2: Optimal solutions on reduced graph
12:    S_dom ← MinWeightedDominatingSet(G')
13:    S_vc ← MinWeightedVertexCover(G')
14:    // Project solutions back to original graph
15:    S_1 ← ProjectToOriginal(S_dom)
16:    S_2 ← ProjectToOriginal(S_vc)
17:    // Phase 3: Ensemble heuristics (orthogonal worst-case coverage)
18:    S_3 ← NetworkXLocalRatio(G_C)
19:    S_4 ← MaxDegreeGreedy(G_C)
20:    S_5 ← MinToMinHeuristic(G_C)
21:    // Phase 4: Select best solution for this component
22:    S_best ← argmin{|S_1|, |S_2|, |S_3|, |S_4|, |S_5|}
23:    S ← S ∪ S_best
24: end for
25: return S
```

#### Graph Reduction to Maximum Degree 1

**Algorithm 2:** ReduceToMaxDegree1: Graph Reduction (the Hallelujah heuristic)  
**Input:** Graph $G = (V, E)$  
**Output:** Reduced graph $G' = (V', E')$ with maximum degree 1, with weighted nodes

```
1: G' ← empty graph
2: weights ← empty dictionary
3: for each vertex u ∈ V do
4:     neighbors ← N(u)
5:     k ← |neighbors|
6:     for i ∈ {0, 1, ..., k-1} do
7:         v ← neighbors[i]
8:         aux ← (u, i)
9:         Add edge (aux, v) to G'
10:        weights[aux] ← 1/k
11:    end for
12: end for
13: Set node attributes in G' using weights
14: return G'
```

#### Optimal Solutions on Degree-1 Graphs

**Algorithm 3:** MinWeightedDominatingSet: Optimal Dominating Set  
**Input:** Graph $G' = (V', E')$ with maximum degree 1, weight function $w: V' \to \mathbb{R}^+$  
**Output:** Minimum weighted dominating set $D \subseteq V'$

```
1: D ← ∅
2: visited ← ∅
3: for each node v ∈ V' do
4:     if v ∉ visited then
5:         d ← deg(v)
6:         if d = 0 then
7:             // Isolated vertex must dominate itself
8:             D ← D ∪ {v}
9:             visited ← visited ∪ {v}
10:        else if d = 1 then
11:            u ← unique neighbor of v
12:            if u ∉ visited then
13:                if w(v) < w(u) or (w(v) = w(u) and v < u) then
14:                    D ← D ∪ {v}
15:                else
16:                    D ← D ∪ {u}
17:                end if
18:                visited ← visited ∪ {v, u}
19:            end if
20:        end if
21:    end if
22: end for
23: return D
```

**Algorithm 4:** MinWeightedVertexCover: Optimal Weighted Vertex Cover  
**Input:** Graph $G' = (V', E')$ with maximum degree 1, weight function $w: V' \to \mathbb{R}^+$  
**Output:** Minimum weighted vertex cover $C \subseteq V'$

```
1: C ← ∅
2: visited ← ∅
3: for each node v ∈ V' do
4:     if v ∉ visited and deg(v) = 1 then
5:         u ← unique neighbor of v
6:         if u ∉ visited then
7:             // Choose minimum weight endpoint to cover edge
8:             if w(v) < w(u) or (w(v) = w(u) and v < u) then
9:                 C ← C ∪ {v}
10:            else
11:                C ← C ∪ {u}
12:            end if
13:            visited ← visited ∪ {v, u}
14:        end if
15:    end if
16: end for
17: return C
```

#### Complementary Heuristics

**Algorithm 5:** MaxDegreeGreedy: Maximum Degree Greedy Heuristic  
**Input:** Graph $G = (V, E)$  
**Output:** Vertex cover $C \subseteq V$

```
1: G_work ← copy of G
2: C ← ∅
3: while |E(G_work)| > 0 do
4:     v ← argmax_{u ∈ V(G_work)} deg(u)
5:     C ← C ∪ {v}
6:     Remove v and all incident edges from G_work
7: end while
8: return C
```

**Algorithm 6:** MinToMinHeuristic: Minimum-to-Minimum Heuristic  
**Input:** Graph $G = (V, E)$  
**Output:** Vertex cover $C \subseteq V$

```
1: G_work ← copy of G
2: C ← ∅
3: while |E(G_work)| > 0 do
4:     // Find vertices with minimum degree
5:     d_min ← min_{u ∈ V(G_work), deg(u) > 0} deg(u)
6:     V_min ← {u ∈ V(G_work) : deg(u) = d_min}
7:     // Get neighbors of minimum-degree vertices
8:     N_min ← ⋃_{u ∈ V_min} N(u)
9:     if N_min ≠ ∅ then
10:        // Among neighbors, find one with minimum degree
11:        v ← argmin_{u ∈ N_min} deg(u)
12:        C ← C ∪ {v}
13:        Remove v and all incident edges from G_work
14:    end if
15: end while
16: return C
```

### 3.2 Complexity Analysis

**Time Complexity:** The algorithm operates in $\mathcal{O}(m \log n)$ time:
- Component decomposition: $\mathcal{O}(n + m)$.
- Reduction to degree-1: $\mathcal{O}(m)$ (each edge processed once).
- Optimal solving on $G'$: $\mathcal{O}(m)$ (linear in reduced graph size).
- Ensemble heuristics: NetworkX local-ratio and greedy methods contribute $\mathcal{O}(m \log n)$.

**Space Complexity:** $\mathcal{O}(m)$ for storing the reduced graph and auxiliary structures.

---

## 4. Approximation Ratio Analysis: Why the Ensemble Works

This section is the theoretical heart of the paper. The core idea is simple to state:

> *Each heuristic's worst-case input is a graph on which at least one other heuristic is optimal or near-optimal.*

The five heuristics therefore *mutually cover* each other's worst cases, and the component-wise $\arg\min$ in Algorithm 1 automatically selects the one that is best-adapted to the local topology.

### 4.1 Individual Heuristic Performance on Graph Classes

#### Sparse Graphs: Optimality via Min-to-Min and Local-Ratio

**Lemma 1 (Path Optimality):** For a path $P_n$ with $n$ vertices, both Min-to-Min and Local-Ratio compute an optimal vertex cover of size $\lceil n/2 \rceil = \mathrm{OPT}(P_n)$.

**Proof:** Min-to-Min identifies the two degree-1 endpoints as minimum-degree vertices and selects their minimum-degree (degree-2) internal neighbours, recursively producing the optimal alternating cover. Local-Ratio achieves optimality on bipartite graphs (including paths) through its weight-based selection.

**Implication:** On sparse graphs (trees, paths, low-degree graphs), the ensemble's $\arg\min$ selects an optimal solution: $\rho = 1.0 \ll \sqrt{2}$.

#### Skewed Bipartite Graphs: Optimality via the Degree-1 Reduction

**Lemma 2 (Bipartite Asymmetry Optimality):** For the complete bipartite graph $K_{\alpha,\beta}$ with $\alpha \ll \beta$, the reduction-based projection achieves an optimal cover of size $\alpha = \mathrm{OPT}(K_{\alpha,\beta})$.

**Proof:** The optimal cover is the smaller partition, of size $\alpha$. The reduction assigns $w_u = 1/\beta$ to small-partition vertices and $w_v = 1/\alpha$ to large-partition vertices. The optimal weighted solution on the reduced graph selects all auxiliary vertices of the small partition (total cost proportional to $\alpha$), projecting back exactly to the optimal cover.

**Implication:** On skewed bipartite graphs, the degree-1 reduction (i.e. the Hallelujah heuristic [[Vega26Hallelujah]](#references)) is optimal — precisely where greedy may mistakenly select the larger partition. This is a textbook example of orthogonal complementarity.

#### Dense Regular Graphs: Optimality via Maximum-Degree Greedy

**Lemma 3 (Clique Optimality):** For the complete graph $K_n$, the maximum-degree greedy heuristic produces an optimal cover of size $n-1 = \mathrm{OPT}(K_n)$.

**Proof:** All vertices have degree $n-1$. Greedy picks an arbitrary vertex, leaving $K_{n-1}$; repeated application yields a cover of size $n-1$, which is optimal. For near-regular graphs, this gives ratio $1 + o(1)$.

**Implication:** On dense regular graphs — precisely where Min-to-Min degenerates (no degree differentiation) — greedy is optimal or near-optimal.

#### Hub-Heavy Scale-Free Graphs: Optimality via the Degree-1 Reduction

**Lemma 4 (Hub Concentration Optimality):** For a star graph (hub $h$ connected to $d$ leaves), the reduction-based projection achieves an optimal cover of size $1 = \mathrm{OPT}$ containing only the hub.

**Proof:** The reduction creates $d$ auxiliary vertices $(h,i)$, each with weight $1/d$, connected to leaves. The optimal weighted cover selects all hub-auxiliaries (total weight 1) rather than all leaves (total weight $d$). Projection yields $\\{h\\}$, which is optimal.

**Implication:** On scale-free / hub-heavy graphs, the Hallelujah heuristic [[Vega26Hallelujah]](#references) captures hub structure optimally, where other heuristics may distribute selections inefficiently across leaves.

### 4.2 Structural Orthogonality: Mutual Worst-Case Coverage

This subsection is the structural justification for the central Strong Hypothesis. The worst-case inputs of the five heuristics are *not* aligned — they are essentially orthogonal.

**Observation (Orthogonal Worst Cases):** The pathological instances for each heuristic are structurally distinct, and each is neutralised by at least one *other* heuristic in the ensemble:
- **Degree-1 Reduction (Hallelujah core):** worst on sparse alternating chains → **Min-to-Min optimal** (Lemma 1).
- **Max-Degree Greedy:** worst on layered sparse graphs where top-degree vertices are locally redundant → **Reduction / Min-to-Min excel**.
- **Min-to-Min:** worst on dense uniform/regular graphs (no degree signal) → **Greedy optimal** (Lemma 3).
- **Local-Ratio:** worst on irregular dense graphs → **Reduction / Greedy excel**.
- **Weighted Dominating Set on reduction:** worst on specific regular patterns → **Weighted VC on reduction excels**, or conversely.

This orthogonality is fundamental: **no simple graph component is known to simultaneously trigger worst-case behaviour in all five heuristics**. The $\arg\min$ in Phase 4 automatically discards the poor performers and selects the one adapted to the local topology.

> **Suggestion (Meta-orthogonality).** We suggest that for any graph family $\mathcal{F}$ on which heuristic $h_i$ is known to achieve ratio $\ge r_i$, there exists at least one other heuristic $h_j$ in the Hvala ensemble with worst-case ratio $\le \max(r_i, \sqrt{2}-\epsilon)$ on $\mathcal{F}$. We offer this as a *suggestion* to guide a future formal proof.

> **Working Assumption (Finite structural taxonomy).** As a stronger working assumption towards the Strong Hypothesis, we adopt: the set of "structural regimes" — sparse alternating, dense uniform, skewed bipartite, hub-heavy, irregular dense — is rich enough that every finite graph decomposes (component-wise) into regimes already dominated by at least one of the five Hvala heuristics. This assumption is explicitly labelled *stronger* and is *not* yet proved.

If this Working Assumption holds, the Orthogonality Observation promotes from an empirical pattern to a structural theorem, and the Strong Hypothesis follows by a component-wise $\arg\min$ argument.

### 4.3 Empirical Confirmation Across Graph Families

Our experimental validation (Section 5) confirms this theoretical complementarity:

> **Empirical Confirmation (Per-family ratios).** On the 201+ tested instances we observe:
> - **Sparse graphs** (bio-networks, trees): ratio 1.000–1.012, with Min-to-Min and Local-Ratio frequently optimal.
> - **Bipartite-like graphs** (collaboration networks): ratio 1.001–1.009, with the degree-1 reduction (Hallelujah core) often optimal.
> - **Dense graphs** (FRB instances): ratio 1.006–1.025, with Greedy performing strongly.
> - **Scale-free graphs** (web graphs, social networks): ratio 1.001–1.032, with the reduction capturing hub structure.
> - **Regular graphs** (3-regular stress tests): ratio 1.069–1.071, demonstrating robustness even on adversarial inputs.

**Key observation.** The maximum observed ratio of 1.071 across all 201+ tested instances, spanning diverse structural properties, is an empirical confirmation (not a proof) compatible with the Strong Hypothesis.

### 4.4 What a Complete Proof Would Require

While we have proved optimality on specific graph classes and observed strong empirical behaviour, a full proof of the Strong Hypothesis would require either:

1. An **exhaustive classification** theorem: a formal proof that the structural taxonomy in the Working Assumption exhausts all finite graphs up to component decomposition, **or**
2. A **constructive counterexample**: an adversarial graph on which all five Hvala heuristics simultaneously achieve ratio $\geq \sqrt{2}$.

The absence of such a counterexample across 201+ diverse instances, combined with the orthogonality analysis, provides *strong evidence* in favour of the Strong Hypothesis, but does *not* constitute a worst-case proof.

---

## 5. Experimental Validation: Comprehensive Results

We present complete experimental results from the three independent validation studies listed in Section 1.4. As explained there, the standalone DIMACS run of [[Vega25Hvala]](#references) has been *absorbed* into and *superseded by* the DIMACS clique-complement subset of the Creo experiment (Section 5.2). To avoid duplication, we report that data only once, in its most recent form.

### 5.1 Experiment 1: Real-World Large Graphs (The Resistire Experiment)

This experiment evaluated Hvala on 88 real-world graphs from the Network Data Repository [[RA15,LargeGraphs]](#references), representing diverse application domains including biological networks, social media, collaboration networks, and web graphs. Conducted on October 15, 2025 [[Vega25Resistire]](#references).

**Complete Real-World Large Graphs Results (88 instances):**

| Instance | Category | V | E | VC Size | Time | Best Known | Ratio | Notes |
|----------|----------|---|---|---------|------|------------|-------|-------|
| bio-celegans | Bio | 453 | 2,025 | 251 | 104.71ms | ~248 | ~1.012 | C. elegans metabolic |
| bio-diseasome | Bio | 516 | 1,188 | 285 | 102.11ms | ~283 | ~1.007 | Disease-gene assoc. |
| bio-dmela | Bio | 7,393 | 25,569 | 2,657 | 13.64s | Unknown | -- | Drosophila |
| bio-yeast | Bio | 1,458 | 1,948 | 456 | 504.85ms | ~453 | ~1.007 | Yeast protein |
| ca-AstroPh | Collab | 17,903 | 196,972 | 11,494 | 151.62s | Unknown | -- | Astrophysics |
| ca-CondMat | Collab | 21,363 | 91,286 | 12,484 | 214.57s | Unknown | -- | Condensed matter |
| ca-CSphd | Collab | 1,025 | 1,043 | 550 | 294.59ms | ~548 | ~1.004 | CS PhD |
| ca-Erdos992 | Collab | 6,100 | 7,515 | 461 | 2.26s | ~459 | ~1.004 | Erdős collab |
| ca-GrQc | Collab | 4,158 | 13,422 | 2,210 | 5.80s | Unknown | -- | General relativity |
| ca-HepPh | Collab | 11,204 | 117,619 | 6,558 | 49.04s | Unknown | -- | High-energy physics |
| ca-netscience | Collab | 379 | 914 | 214 | 61.72ms | ~212 | ~1.009 | Network science |
| ia-email-EU | Email | 32,430 | 54,397 | 820 | 29.73s | Unknown | -- | EU research email |
| ia-email-univ | Email | 1,133 | 5,451 | 605 | 486.58ms | ~603 | ~1.003 | University email |
| ia-enron-large | Email | 33,696 | 180,811 | 12,792 | 391.87s | Unknown | -- | Enron large |
| ia-enron-only | Email | 143 | 623 | 87 | 16.07ms | ~86 | ~1.012 | Enron core |
| ia-fb-messages | Social | 1,266 | 6,451 | 580 | 998.20ms | ~578 | ~1.003 | Facebook msgs |
| ia-infect-dublin | Social | 410 | 2,765 | 298 | 108.60ms | ~296 | ~1.007 | Infection Dublin |
| ia-infect-hyper | Social | 113 | 188 | 92 | 29.43ms | ~91 | ~1.011 | Infection hypertext |
| ia-reality | Social | 6,809 | 7,680 | 81 | 657.86ms | Unknown | -- | Reality mining |
| ia-wiki-Talk | Wiki | 92,117 | 360,767 | 17,288 | 1868.99s | Unknown | -- | Wikipedia talk |
| inf-power | Infra | 4,941 | 6,594 | 2,207 | 7.45s | Unknown | -- | US power grid |
| rec-amazon | Rec | 262,111 | 899,792 | 47,891 | 4123.24s | Unknown | -- | Amazon products |
| rt-retweet | Retweet | 96 | 117 | 32 | 4.98ms | ~31 | ~1.032 | General retweet |
| rt-twitter-copen | Retweet | 761 | 1,029 | 237 | 161.06ms | ~235 | ~1.009 | Twitter Copenhagen |
| scc_enron-only | SCC | 143 | 251 | 138 | 183.99ms | ~137 | ~1.007 | Enron SCC |
| scc_fb-forum | SCC | 899 | 7,089 | 372 | 2.28s | ~370 | ~1.005 | FB forum SCC |
| scc_fb-messages | SCC | 1,266 | 3,125 | 1,072 | 18.12s | Unknown | -- | FB messages SCC |
| scc_infect-dublin | SCC | 410 | 1,800 | 9,104 | 5.48s | Unknown | -- | Infection Dublin SCC |
| scc_infect-hyper | SCC | 113 | 171 | 110 | 171.80ms | ~109 | ~1.009 | Infection hyper SCC |
| scc_retweet | SCC | 96 | 87 | 561 | 2.08s | Unknown | -- | Retweet SCC |
| scc_retweet-crawl | SCC | 21,297 | 17,362 | 8,419 | 14.03s | Unknown | -- | Retweet crawl SCC |
| scc_rt_alwefaq | SCC | 35 | 34 | 35 | 9.47ms | 35 | 1.000 | **Optimal** |
| scc_rt_assad | SCC | 16 | 15 | 16 | 1.99ms | 16 | 1.000 | **Optimal** |
| scc_rt_bahrain | SCC | 37 | 36 | 37 | 5.52ms | 37 | 1.000 | **Optimal** |
| scc_rt_barackobama | SCC | 29 | 28 | 29 | 6.02ms | 29 | 1.000 | **Optimal** |
| scc_rt_damascus | SCC | 15 | 14 | 15 | 2.05ms | 15 | 1.000 | **Optimal** |
| scc_rt_dash | SCC | 15 | 14 | 15 | 2.99ms | 15 | 1.000 | **Optimal** |
| scc_rt_gmanews | SCC | 46 | 45 | 46 | 25.25ms | 46 | 1.000 | **Optimal** |
| scc_rt_gop | SCC | 6 | 5 | 6 | 1.00ms | 6 | 1.000 | **Optimal** |
| scc_rt_http | SCC | 2 | 1 | 2 | 0.98ms | 2 | 1.000 | **Optimal** |
| scc_rt_israel | SCC | 11 | 10 | 11 | 0.99ms | 11 | 1.000 | **Optimal** |
| scc_rt_justinbieber | SCC | 26 | 25 | 26 | 10.96ms | 26 | 1.000 | **Optimal** |
| scc_rt_ksa | SCC | 12 | 11 | 12 | 1.08ms | 12 | 1.000 | **Optimal** |
| scc_rt_lebanon | SCC | 5 | 4 | 5 | 1.08ms | 5 | 1.000 | **Optimal** |
| scc_rt_libya | SCC | 12 | 11 | 12 | 2.07ms | 12 | 1.000 | **Optimal** |
| scc_rt_lolgop | SCC | 103 | 102 | 103 | 182.49ms | 103 | 1.000 | **Optimal** |
| scc_rt_mittromney | SCC | 42 | 41 | 42 | 5.98ms | 42 | 1.000 | **Optimal** |
| scc_rt_obama | SCC | 4 | 3 | 4 | 1.08ms | 4 | 1.000 | **Optimal** |
| scc_rt_occupy | SCC | 22 | 21 | 22 | 3.01ms | 22 | 1.000 | **Optimal** |
| scc_rt_occupywallstnyc | SCC | 45 | 44 | 45 | 22.50ms | 45 | 1.000 | **Optimal** |
| scc_rt_oman | SCC | 6 | 5 | 6 | 0.98ms | 6 | 1.000 | **Optimal** |
| scc_rt_onedirection | SCC | 29 | 28 | 29 | 8.46ms | 29 | 1.000 | **Optimal** |
| scc_rt_p2 | SCC | 12 | 11 | 12 | 1.01ms | 12 | 1.000 | **Optimal** |
| scc_rt_qatif | SCC | 5 | 4 | 5 | 1.08ms | 5 | 1.000 | **Optimal** |
| scc_rt_saudi | SCC | 17 | 16 | 17 | 2.07ms | 17 | 1.000 | **Optimal** |
| scc_rt_tcot | SCC | 12 | 11 | 12 | 2.00ms | 12 | 1.000 | **Optimal** |
| scc_rt_tlot | SCC | 6 | 5 | 6 | 1.00ms | 6 | 1.000 | **Optimal** |
| scc_rt_uae | SCC | 8 | 7 | 8 | 1.38ms | 8 | 1.000 | **Optimal** |
| scc_rt_voteonedirection | SCC | 4 | 3 | 4 | 0.98ms | 4 | 1.000 | **Optimal** |
| scc_twitter-copen | SCC | 761 | 662 | 1,328 | 18.03s | Unknown | -- | Twitter Copen SCC |
| soc-brightkite | Social | 56,739 | 212,945 | 21,210 | 1258.10s | Unknown | -- | Brightkite location |
| soc-dolphins | Social | 62 | 159 | 35 | 5.06ms | ~34 | ~1.029 | Dolphin social |
| soc-douban | Social | 154,908 | 327,162 | 8,685 | 1629.90s | Unknown | -- | Douban social |
| soc-epinions | Social | 26,588 | 100,120 | 9,774 | 263.38s | Unknown | -- | Epinions trust |
| soc-karate | Social | 34 | 78 | 14 | 1.66ms | 14 | 1.000 | **Optimal - Karate** |
| soc-slashdot | Social | 70,068 | 358,647 | 22,373 | 1805.07s | Unknown | -- | Slashdot social |
| soc-wiki-Vote | Social | 889 | 2,914 | 406 | 299.78ms | ~404 | ~1.005 | Wikipedia voting |
| socfb-CMU | Facebook | 6,621 | 251,214 | 5,054 | 29.27s | Unknown | -- | Carnegie Mellon |
| socfb-Duke14 | Facebook | 9,885 | 506,437 | 7,776 | 73.58s | Unknown | -- | Duke University |
| socfb-MIT | Facebook | 6,441 | 251,230 | 4,723 | 28.13s | Unknown | -- | MIT |
| socfb-Stanford3 | Facebook | 11,586 | 568,309 | 8,626 | 102.50s | Unknown | -- | Stanford |
| socfb-UCLA | Facebook | 20,453 | 747,604 | 15,434 | 324.98s | Unknown | -- | UCLA |
| socfb-UConn | Facebook | 17,206 | 636,836 | 13,422 | 228.94s | Unknown | -- | UConn |
| socfb-UCSB37 | Facebook | 14,917 | 482,215 | 11,429 | 162.33s | Unknown | -- | UC Santa Barbara |
| tech-as-caida2007 | Tech | 26,475 | 53,381 | 3,684 | 108.54s | Unknown | -- | CAIDA AS 2007 |
| tech-internet-as | Tech | 22,963 | 48,436 | 5,700 | 263.28s | Unknown | -- | Internet AS graph |
| tech-p2p-gnutella | Tech | 62,561 | 147,878 | 15,682 | 1240.83s | Unknown | -- | Gnutella P2P |
| tech-RL-caida | Tech | 190,914 | 607,610 | 75,680 | 17095.90s | Unknown | -- | CAIDA router-level |
| tech-routers-rf | Tech | 2,113 | 6,632 | 795 | 1.25s | ~793 | ~1.003 | Router network |
| tech-WHOIS | Tech | 7,476 | 56,943 | 2,287 | 15.46s | Unknown | -- | WHOIS network |
| web-BerkStan | Web | 12,776 | 19,500 | 5,390 | 44.16s | Unknown | -- | Berkeley-Stanford |
| web-edu | Web | 3,031 | 6,474 | 1,451 | 2.63s | ~1,449 | ~1.001 | Educational domain |
| web-google | Web | 1,299 | 2,773 | 498 | 483.96ms | ~497 | ~1.002 | Google web graph |
| web-indochina-2004 | Web | 11,358 | 47,606 | 7,300 | 45.95s | Unknown | -- | Indochina crawl |
| web-polblogs | Web | 643 | 2,280 | 245 | 140.23ms | ~243 | ~1.008 | Political blogs |
| web-sk-2005 | Web | 121,176 | 1,043,877 | 58,190 | 6126.11s | Unknown | -- | Slovak web crawl |
| web-spam | Web | 4,767 | 37,375 | 2,315 | 8.31s | Unknown | -- | Web spam corpus |
| web-webbase-2001 | Web | 16,062 | 25,593 | 2,652 | 35.39s | Unknown | -- | Webbase 2001 |

**Performance Summary (Real-World Large Graphs):**
- **Total instances tested:** 88
- **Optimal solutions found:** 28 (31.8%)
- **Average approximation ratio (where known):** 1.007
- **Best ratio:** 1.000 (28 instances)
- **Worst ratio:** 1.032 (rt-retweet)
- **Largest instance solved:** rec-amazon (262,111 vertices, 899,792 edges) in 68.7 minutes
- **Runtime distribution:** Sub-second: 38 instances (43.2%); 1–60 seconds: 27 instances (30.7%); 1–10 minutes: 13 instances (14.8%); Over 10 minutes: 10 instances (11.4%)

### 5.2 Experiment 2: NPBench Hard Instances (The Creo Experiment)

This experiment, conducted on December 20, 2025 [[Vega25Creo]](#references), evaluated Hvala v0.0.7 on 113 challenging instances from the NPBench collection [[NPBench]](#references), including FRB (Factoring and Random Benchmarks) and DIMACS clique complement graphs.

#### FRB Instances (40 instances)

**FRB Benchmark Results (Factoring and Random):**

| Instance | Optimal | Hvala | Time | Ratio |
|----------|---------|-------|------|-------|
| frb30-15-1.mis | 420 | 426 | 443.82ms | 1.014 |
| frb30-15-2.mis | 420 | 425 | 506.81ms | 1.012 |
| frb30-15-3.mis | 420 | 426 | 475.87ms | 1.014 |
| frb30-15-4.mis | 420 | 425 | 416.66ms | 1.012 |
| frb30-15-5.mis | 420 | 425 | 445.95ms | 1.012 |
| frb35-17-1.mis | 560 | 566 | 719.36ms | 1.011 |
| frb35-17-2.mis | 560 | 565 | 739.85ms | 1.009 |
| frb35-17-3.mis | 560 | 566 | 774.78ms | 1.011 |
| frb35-17-4.mis | 560 | 566 | 856.32ms | 1.011 |
| frb35-17-5.mis | 560 | 566 | 813.15ms | 1.011 |
| frb40-19-1.mis | 720 | 728 | 1.16s | 1.011 |
| frb40-19-2.mis | 720 | 728 | 1.22s | 1.011 |
| frb40-19-3.mis | 720 | 726 | 1.19s | 1.008 |
| frb40-19-4.mis | 720 | 729 | 1.20s | 1.013 |
| frb40-19-5.mis | 720 | 728 | 1.21s | 1.011 |
| frb45-21-1.mis | 900 | 906 | 1.96s | 1.007 |
| frb45-21-2.mis | 900 | 910 | 1.89s | 1.011 |
| frb45-21-3.mis | 900 | 908 | 1.89s | 1.009 |
| frb45-21-4.mis | 900 | 910 | 1.86s | 1.011 |
| frb45-21-5.mis | 900 | 907 | 1.83s | 1.008 |
| frb50-23-1.mis | 1100 | 1108 | 2.68s | 1.007 |
| frb50-23-2.mis | 1100 | 1109 | 2.72s | 1.008 |
| frb50-23-3.mis | 1100 | 1108 | 2.63s | 1.007 |
| frb50-23-4.mis | 1100 | 1109 | 2.91s | 1.008 |
| frb50-23-5.mis | 1100 | 1111 | 2.92s | 1.010 |
| frb53-24-1.mis | 1219 | 1231 | 4.57s | 1.010 |
| frb53-24-2.mis | 1219 | 1228 | 3.33s | 1.007 |
| frb53-24-3.mis | 1219 | 1229 | 4.82s | 1.008 |
| frb53-24-4.mis | 1219 | 1227 | 3.46s | 1.007 |
| frb53-24-5.mis | 1219 | 1229 | 3.53s | 1.008 |
| frb56-25-1.mis | 1344 | 1355 | 3.88s | 1.008 |
| frb56-25-2.mis | 1344 | 1358 | 4.22s | 1.010 |
| frb56-25-3.mis | 1344 | 1354 | 4.12s | 1.007 |
| frb56-25-4.mis | 1344 | 1352 | 4.11s | 1.006 |
| frb56-25-5.mis | 1344 | 1354 | 3.85s | 1.007 |
| frb59-26-1.mis | 1475 | 1485 | 5.00s | 1.007 |
| frb59-26-2.mis | 1475 | 1486 | 4.86s | 1.007 |
| frb59-26-3.mis | 1475 | 1485 | 5.67s | 1.007 |
| frb59-26-4.mis | 1475 | 1485 | 5.06s | 1.007 |
| frb59-26-5.mis | 1475 | 1486 | 4.80s | 1.007 |
| frb100-40.mis | 3900 | 3922 | 27.78s | 1.006 |

#### DIMACS Clique Complement Benchmarks (73 instances)

**DIMACS Clique Complement Benchmark Results:**

| Instance | Optimal | Hvala | Time | Ratio |
|----------|---------|-------|------|-------|
| brock200_1 | 179 | 180 | 127.45ms | 1.006 |
| brock200_2 | 188 | 192 | 238.33ms | 1.021 |
| brock200_3 | 183 | 187 | 176.02ms | 1.022 |
| brock200_4 | 183 | 187 | 142.97ms | 1.022 |
| brock400_1 | 373 | 378 | 539.88ms | 1.013 |
| brock400_2 | 373 | 378 | 581.28ms | 1.013 |
| brock400_3 | 373 | 379 | 560.76ms | 1.016 |
| brock400_4 | 373 | 378 | 508.98ms | 1.013 |
| brock800_1 | 777 | 782 | 3.56s | 1.006 |
| brock800_2 | 777 | 782 | 3.86s | 1.006 |
| brock800_3 | 777 | 783 | 3.79s | 1.008 |
| brock800_4 | 777 | 783 | 3.75s | 1.008 |
| c-fat200-1 | 186 | 188 | 588.37ms | 1.011 |
| c-fat200-2 | 174 | 176 | 380.66ms | 1.011 |
| c-fat200-5 | 140 | 142 | 287.03ms | 1.014 |
| c-fat500-1 | 482 | 486 | 3.35s | 1.008 |
| c-fat500-10 | 372 | 374 | 2.42s | 1.005 |
| c-fat500-2 | 470 | 474 | 3.49s | 1.009 |
| c-fat500-5 | 434 | 436 | 3.09s | 1.005 |
| C125.9 | 91 | 93 | 31.63ms | 1.022 |
| C250.9 | 206 | 209 | 91.34ms | 1.015 |
| C500.9 | 443 | 451 | 330.04ms | 1.018 |
| C1000.9 | 932 | 939 | 1.94s | 1.008 |
| C2000.5 | 1984 | 1988 | 46.18s | 1.002 |
| C2000.9 | 1920 | 1934 | 10.29s | 1.007 |
| C4000.5 | 3978 | 3986 | 216.52s | 1.002 |
| gen200_p0.9_44 | 160 | 164 | 63.44ms | 1.025 |
| gen200_p0.9_55 | 160 | 163 | 40.88ms | 1.019 |
| gen400_p0.9_55 | 352 | 356 | 200.60ms | 1.011 |
| gen400_p0.9_65 | 352 | 356 | 255.87ms | 1.011 |
| gen400_p0.9_75 | 350 | 353 | 229.63ms | 1.009 |
| hamming6-2 | 32 | 32 | 0.00ms | 1.000 |
| hamming6-4 | 60 | 60 | 37.19ms | 1.000 |
| hamming8-2 | 128 | 128 | 37.79ms | 1.000 |
| hamming8-4 | 238 | 240 | 238.51ms | 1.008 |
| hamming10-2 | 512 | 512 | 455.43ms | 1.000 |
| hamming10-4 | 992 | 992 | 2.73s | 1.000 |
| johnson8-2-4 | 24 | 24 | 0.00ms | 1.000 |
| johnson8-4-4 | 56 | 56 | 5.20ms | 1.000 |
| johnson16-2-4 | 112 | 112 | 31.88ms | 1.000 |
| johnson32-2-4 | 480 | 480 | 363.80ms | 1.000 |
| keller4 | 160 | 160 | 95.72ms | 1.000 |
| keller5 | 749 | 752 | 1.87s | 1.004 |
| keller6 | 3303 | 3314 | 56.88s | 1.003 |
| MANN_a9 | 29 | 29 | 8.65ms | 1.000 |
| MANN_a27 | 252 | 253 | 64.22ms | 1.004 |
| MANN_a45 | 690 | 693 | 443.84ms | 1.004 |
| MANN_a81 | 2221 | 2225 | 4.30s | 1.002 |
| p_hat300-1 | 292 | 293 | 1.52s | 1.003 |
| p_hat300-2 | 275 | 277 | 534.66ms | 1.007 |
| p_hat300-3 | 264 | 267 | 298.34ms | 1.011 |
| p_hat500-1 | 491 | 492 | 2.75s | 1.002 |
| p_hat500-2 | 465 | 467 | 1.86s | 1.004 |
| p_hat500-3 | 453 | 454 | 1.04s | 1.002 |
| p_hat700-1 | 689 | 692 | 6.00s | 1.004 |
| p_hat700-2 | 656 | 657 | 4.07s | 1.002 |
| p_hat700-3 | 640 | 641 | 2.15s | 1.002 |
| p_hat1000-1 | 988 | 991 | 15.20s | 1.003 |
| p_hat1000-2 | 956 | 958 | 9.30s | 1.002 |
| p_hat1000-3 | 937 | 939 | 5.06s | 1.002 |
| p_hat1500-1 | 1488 | 1490 | 33.08s | 1.001 |
| p_hat1500-2 | 1437 | 1439 | 22.18s | 1.001 |
| p_hat1500-3 | 1413 | 1416 | 12.09s | 1.002 |
| san200_0.7_1 | 182 | 183 | 143.67ms | 1.005 |
| san200_0.7_2 | 183 | 185 | 125.95ms | 1.011 |
| san200_0.9_1 | 150 | 152 | 63.71ms | 1.013 |
| san200_0.9_2 | 160 | 161 | 63.81ms | 1.006 |
| san200_0.9_3 | 166 | 169 | 47.61ms | 1.018 |
| san400_0.5_1 | 387 | 391 | 988.70ms | 1.010 |
| san400_0.7_1 | 376 | 378 | 683.53ms | 1.005 |
| san400_0.7_2 | 379 | 382 | 649.11ms | 1.008 |
| san400_0.7_3 | 382 | 385 | 635.93ms | 1.008 |
| san400_0.9_1 | 316 | 317 | 255.68ms | 1.003 |
| san1000 | 986 | 990 | 8.70s | 1.004 |
| sanr200_0.7 | 183 | 184 | 196.33ms | 1.005 |
| sanr200_0.9 | 162 | 163 | 64.02ms | 1.006 |
| sanr400_0.5 | 387 | 388 | 994.94ms | 1.003 |
| sanr400_0.7 | 379 | 381 | 697.75ms | 1.005 |

**Performance Summary (NPBench Hard Instances):**
- **Total instances tested:** 113 (40 FRB + 73 DIMACS)
- **Optimal solutions found:** 12 instances
- **Average approximation ratio:** 1.006
- **Best ratio:** 1.000 (12 optimal instances)
- **Worst ratio:** 1.025 (gen200_p0.9_44)
- **Instances with ratio ≤ 1.015:** 107 (95%)
- **FRB average ratio:** 1.009
- **DIMACS average ratio:** 1.007

### 5.3 Experiment 3: AI-Validated Stress Testing (The Gemini–Vega Validation)

This independent validation study, conducted on December 21, 2025 using Gemini AI [[Vega25Gemini]](#references), tested Hvala on adversarially constructed hard graphs designed to challenge heuristic algorithms. The focus was on 3-regular graphs where every vertex has degree exactly 3, eliminating degree-based heuristic advantages.

**Experimental Design:**
- **Graph Construction:** Random 3-regular graphs (uniform degree distribution)
- **AI Validation:** Gemini AI architected testing framework and verified results
- **Baseline Comparison:** Standard greedy highest-degree-first heuristic
- **Theoretical Context:** Optimal vertex cover for 3-regular graphs is approximately $0.5446n$ vertices

**Gemini–Vega Stress Test Results on 3-Regular Graphs:**

| Graph Size | Vertices | Edges | Hvala Size | Greedy Size | Hvala Ratio |
|------------|----------|-------|------------|-------------|-------------|
| Power-Law (N=10,000) | 10,000 | -- | 4,957 | 5,093 | -- |
| 3-Regular (N=5,000) | 5,000 | 7,500 | 2,917 (58.34%) | 3,073 (61.46%) | 1.0712 |
| 3-Regular (N=20,000) | 20,000 | 30,000 | 11,647 (58.24%) | 12,350 (61.75%) | 1.0693 |
| *Theoretical optimal for 3-regular: ~0.5446n vertices* | | | | | |

**Key Observations:**
1. **Improvement Over Greedy:** Hvala consistently outperforms greedy by 2.7–3.5%.
2. **Ratio Stability:** Approximation ratio improved slightly from 1.0712 to 1.0693 as graph size doubled.
3. **Theoretical Context:** Achieved ratio of 1.069 against theoretical optimum (~0.5446 $n$).
4. **Computational Feasibility:** 20,000-vertex graph solved in 162.09 seconds.
5. **AI Verification:** Independent validation through Gemini AI confirms correctness and reproducibility.

**Gemini AI Full Transcript:** Complete experimental session available at https://gemini.google.com/share/55109efe4d85

---

## 6. Arguments Supporting the Strong Hypothesis

We now collect, in a single place, the five categories of evidence for the Strong Hypothesis ($\rho < \sqrt{2}$, hence **P = NP**). Each is labelled according to the taxonomy of Section 1.2 (Hypothesis / Confirmation / Suggestion / Stronger Assumption), so the reader always knows what is proved and what is conjectured.

### 6.1 Argument 1 (Empirical Confirmation): Consistency Across Diverse Instance Classes

> **Empirical Confirmation (Cross-experiment consistency).** Across three independent experimental studies spanning 201+ instances with radically different structural properties, no instance exceeded ratio 1.071.

**Cross-Experiment Consistency Analysis:**

| Experiment | Instances | Avg. Ratio | Max Ratio |
|------------|-----------|------------|-----------|
| Real-World Large Graphs (Resistire) | 88 | 1.007 | 1.032 |
| NPBench Hard Instances (Creo) | 113 | 1.006 | 1.025 |
| AI Stress Tests (Gemini–Vega) | 3 | -- | 1.071 |
| **Combined** | **204** | **1.007** | **1.071** |

**Strength:** Consistency across bipartite graphs, scale-free networks, dense random graphs, structured benchmarks, and adversarially constructed 3-regular graphs suggests that the algorithm exploits fundamental structural properties rather than artefacts of specific graph families.

**Limitation:** Consistency across tested instances does *not* prove impossibility of worse instances. If P ≠ NP, then instances achieving ratio $\geq \sqrt{2}$ must exist by the hardness results of Khot et al. under SETH.

### 6.2 Argument 2 (Empirical Confirmation): Scalability and Improved Performance on Larger Instances

> **Empirical Confirmation (Non-degradation at scale).** Contrary to typical heuristic degradation, performance stabilises or improves on larger instances.

- C4000.5 (3,986 vertices): ratio 1.002.
- p_hat1500-1 (1,488 optimal): ratio 1.001.
- 20K 3-regular graph: ratio 1.0693 (better than the 5K instance at 1.0712).
- rec-amazon (262K vertices): successfully processed.

**Implication:** If the algorithm degraded systematically on larger instances, we would expect ratios approaching or exceeding $\sqrt{2} \approx 1.414$ on the largest graphs tested. Instead, the largest instances maintain ratios ≤ 1.071.

**Counterargument:** Theoretical worst-case instances may require specific adversarial constructions not present in our test suite, possibly with size beyond computational feasibility.

### 6.3 Argument 3 (Empirical Confirmation): High Frequency of Provably Optimal Solutions

> **Empirical Confirmation (Exact optima).** The algorithm achieves provably optimal solutions on significant fractions of tested instances.

- Real-World: 28/88 optimal (31.8%).
- NPBench: 12/113 optimal (10.6%).

**Implication:** Achieving exact optimality on 40 instances across the two experiments demonstrates that the degree-1 reduction core of Hvala (the Hallelujah heuristic [[Vega26Hallelujah]](#references)) captures sufficient structural information to solve several graph classes exactly. This suggests the reduction is fundamentally sound, not a merely heuristic approximation.

**Theoretical context:** These optimal solutions occur primarily on tree-like structures (SCC instances) and highly regular graphs (Hamming, Johnson), where the degree-1 reduction perfectly captures the problem structure.

### 6.4 Argument 4 (Empirical Confirmation): Consistent Improvement Over Greedy Baselines

> **Empirical Confirmation (Greedy dominance).** Across all experiments, Hvala consistently outperforms simple greedy by 2–4%.

**Hvala vs. Greedy Comparison (Selected Instances):**

| Instance | Hvala | Greedy | Improvement | Optimal |
|----------|-------|--------|-------------|---------|
| 3-Regular (5K) | 2,917 | 3,073 | 5.1% | ~2,723 |
| 3-Regular (20K) | 11,647 | 12,350 | 5.7% | ~10,892 |
| Power-Law (10K) | 4,957 | 5,093 | 2.7% | Unknown |

**Strength:** This consistent improvement across diverse structures suggests Hvala captures global optimisation information missed by local degree-based heuristics.

### 6.5 Argument 5 (Theoretical): Weight-Preserving Reduction (the Degree-1 Vertex Cover Core)

The theoretical foundation of the degree-1 weighted vertex cover reduction — comprising weight preservation, lower bound preservation, and deterministic stability via lexicographical tie-breaking — is developed in full in the companion paper [[Vega26Hallelujah]](#references). We refer the reader to that work for the complete proofs; reproducing them here is outside the scope of the present compendium.

**Implication:** The proven properties established in [[Vega26Hallelujah]](#references) ensure that the optimal weighted vertex cover on the reduced graph $G'$ provides near-optimal guidance for the cover on $G$, forming the rigorous backbone of the $S_2$ component in Hvala and a principled foundation *beyond* pure heuristic intuition.

> **Working Assumption (Closing the $\sqrt{2}$ gap via the ensemble).** We adopt, as a stronger assumption, that the component-wise $\arg\min$ of the five Hvala heuristics inherits the lower-bound and stability properties of the degree-1 vertex cover core (established in [[Vega26Hallelujah]](#references)) on *every* connected component, and that the non-covered regime of one heuristic is always absorbed by at least one of the other four. Under this assumption, the Strong Hypothesis follows.

**Gap:** The properties established in [[Vega26Hallelujah]](#references), and their extension in this Working Assumption, do not yet constitute a complete worst-case proof of $\rho < \sqrt{2}$. Closing this gap is explicitly left to forthcoming work (Section 7).

---

## 7. Forthcoming Work

The companion paper [[Vega26Hallelujah,Vega25HallelujahPreprint]](#references) develops the *single* degree-1 weighted vertex cover heuristic in full theoretical detail. It serves as the rigorous backbone of the present paper's core component.

Building on that companion paper, we plan a sequel that:

1. Proves (or refutes) the Working Assumptions of Sections 4.2 and 6.5, turning the Strong Hypothesis into a theorem — or producing a constructive counterexample.
2. Extends the orthogonality analysis of Section 4 to a formal taxonomy of graph regimes, with quantitative worst-case bounds per heuristic.
3. Integrates a sixth and seventh heuristic targeting the residual regime (dense irregular graphs) identified by the Gemini–Vega stress tests.

> **Suggestion (Path to P = NP).** We suggest that the most promising route to converting the Strong Hypothesis into a theorem is (i) sharpening the structural taxonomy, (ii) re-using the lower-bound preservation already proven in [[Vega26Hallelujah]](#references), and (iii) invoking a component-wise $\arg\min$ bound. We offer this as a research direction, not a result.

---

## 8. Conclusion

This work demonstrates that the Hvala algorithm — a compendium of five mutually-compensating heuristics, built around the degree-1 weighted vertex cover reduction whose rigorous analysis appears in the companion paper [[Vega26Hallelujah,Vega25HallelujahPreprint]](#references) — achieves exceptional empirical performance on Minimum Vertex Cover, with no tested instance exceeding ratio 1.071 across 201+ diverse graphs.

We have been explicit about the epistemic status of every claim, using dedicated environments: results are tagged as Strong Hypothesis, Working Assumption, Empirical Confirmation, or Suggestion. In particular:

- **The Strong Hypothesis** states that Hvala achieves $\rho < \sqrt{2}$ on *every* graph, which would imply **P = NP** under SETH. This is clearly flagged as an unproven hypothesis.
- **The Working Assumptions** of Sections 4.2 and 6.5 are explicit stronger assumptions under which the Strong Hypothesis becomes a theorem.
- **The Empirical Confirmations** of Section 6 record what the experiments on 201+ instances actually show.
- **The Suggestions** of Sections 4.2 and 7 are invitations to further work.

Whether the Strong Hypothesis ultimately holds (solving P versus NP) or fails (revealing the limits of empirical validation and the importance of worst-case analysis), we believe the structured epistemic presentation — and the explicit focus on *how* the five heuristics cover each other's worst cases — advances our collective understanding of ensemble approximation. We invite vigorous scrutiny, attempted refutation, and independent validation from the theoretical computer science community.

**Algorithm Availability:** The Hvala algorithm is publicly available for independent verification:
- **PyPI:** https://pypi.org/project/hvala
- **Installation:** `pip install hvala`
- **Usage:** `from hvala.algorithm import find_vertex_cover`
- **Source Code:** Available for inspection and verification.

**Companion Paper:** The degree-1 weighted vertex cover reduction is fully analysed in [[Vega26Hallelujah]](#references); a preprint is available at [[Vega25HallelujahPreprint]](#references).

---

## Acknowledgment

The author is sincerely grateful to Iris, Marilin, Sonia, Yoselin, Arelis, Anissa, Liuva, Yudit, Gretel, Gema, and Blaquier, as well as Israel, Arderi, Juan Carlos, Yamil, Alejandro, Aroldo, Yary, Reinaldo, Alex, Emmanuel, and Michael for their constant support. Whether through encouragement, stimulating conversations, practical assistance, or simply being present during challenging moments, their contributions have played an important role in bringing this work to completion.

---

## References

**[karp2009reducibility]** Karp, Richard M. (2009). *Reducibility Among Combinatorial Problems.* In 50 Years of Integer Programming 1958–2008 (pp. 219–241). Springer, Berlin, Heidelberg. DOI: [10.1007/978-3-540-68279-0_8](https://doi.org/10.1007/978-3-540-68279-0_8).

**[papadimitriou1998combinatorial]** Papadimitriou, Christos H. and Steiglitz, Kenneth (1998). *Combinatorial Optimization: Algorithms and Complexity.* Courier Corporation, Mineola, New York.

**[karakostas2009better]** Karakostas, George (2009). *A Better Approximation Ratio for the Vertex Cover Problem.* ACM Transactions on Algorithms, 5(4), 1–8. DOI: [10.1145/1597036.1597045](https://doi.org/10.1145/1597036.1597045).

**[karpinski1996approximating]** Karpinski, Marek and Zelikovsky, Alexander (1996). *Approximating Dense Cases of Covering Problems.* DIMACS Series in Discrete Mathematics and Theoretical Computer Science, 26, 147–164. Providence, Rhode Island.

**[dinur2005hardness]** Dinur, Irit and Safra, Samuel (2005). *On the Hardness of Approximating Minimum Vertex Cover.* Annals of Mathematics, 162, 439–485. DOI: [10.4007/annals.2005.162.439](https://doi.org/10.4007/annals.2005.162.439).

**[khot2017independent]** Khot, Subhash and Minzer, Dor and Safra, Muli (2017). *On Independent Sets, 2-to-2 Games, and Grassmann Graphs.* Proceedings of the 49th Annual ACM SIGACT Symposium on Theory of Computing, 576–589. Montreal, Canada. DOI: [10.1145/3055399.3055432](https://doi.org/10.1145/3055399.3055432).

**[dinur2018towards]** Dinur, Irit and Khot, Subhash and Kindler, Guy and Minzer, Dor and Safra, Muli (2018). *Towards a proof of the 2-to-1 games conjecture?* Proceedings of the 50th Annual ACM SIGACT Symposium on Theory of Computing, 376–389. Los Angeles, California. DOI: [10.1145/3188745.3188804](https://doi.org/10.1145/3188745.3188804).

**[khot2018pseudorandom]** Khot, Subhash and Minzer, Dor and Safra, Muli (2018). *Pseudorandom Sets in Grassmann Graph Have Near-Perfect Expansion.* 2018 IEEE 59th Annual Symposium on Foundations of Computer Science, 592–601. Paris, France. DOI: [10.1109/FOCS.2018.00062](https://doi.org/10.1109/FOCS.2018.00062).

**[khot2008vertex]** Khot, Subhash and Regev, Oded (2008). *Vertex Cover Might Be Hard to Approximate to Within* $2-\epsilon$*.* Journal of Computer and System Sciences, 74(3), 335–349. DOI: [10.1016/j.jcss.2007.06.019](https://doi.org/10.1016/j.jcss.2007.06.019).

**[khot2002unique]** Khot, Subhash (2002). *On the Power of Unique 2-Prover 1-Round Games.* Proceedings of the 34th Annual ACM Symposium on Theory of Computing, 767–775. Montreal, Canada. DOI: [10.1145/509907.510017](https://doi.org/10.1145/509907.510017).

**[cai2017finding]** Cai, Shaowei and Lin, Jinkun and Luo, Chuan (2017). *Finding a Small Vertex Cover in Massive Sparse Graphs.* Journal of Artificial Intelligence Research, 59, 463–494. DOI: [10.1613/jair.5443](https://doi.org/10.1613/jair.5443).

**[zhang2023tivc]** Zhang, Yu and Wang, Shengzhi and Liu, Chanjuan and Zhu, Enqiang (2023). *TIVC: An Efficient Local Search Algorithm for Minimum Vertex Cover in Large Graphs.* Sensors, 23(18), 7831. DOI: [10.3390/s23187831](https://doi.org/10.3390/s23187831).

**[harris2024faster]** Harris, David G. and Narayanaswamy, N. S. (2024). *A Faster Algorithm for Vertex Cover Parameterized by Solution Size.* 41st International Symposium on Theoretical Aspects of Computer Science, 40:1–40:18. Clermont-Ferrand, France. DOI: [10.4230/LIPIcs.STACS.2024.40](https://doi.org/10.4230/LIPIcs.STACS.2024.40).

**[bar1985local]** Bar-Yehuda, R. and Even, S. (1985). *A Local-Ratio Theorem for Approximating the Weighted Vertex Cover Problem.* Annals of Discrete Mathematics, 25, 27–46.

**[luo2019local]** Luo, Chuan and Hoos, Holger H. and Cai, Shaowei and Lin, Qingwei and Zhang, Hongyu and Zhang, Dongmei (2019). *Local search with efficient automatic configuration for minimum vertex cover.* Proceedings of the 28th International Joint Conference on Artificial Intelligence, 1297–1304. Macao, China.

**[banharnsakun2023new]** Banharnsakun, Anan (2023). *A New Approach for Solving the Minimum Vertex Cover Problem Using Artificial Bee Colony Algorithm.* Decision Analytics Journal, 6, 100175. DOI: [10.1016/j.dajour.2023.100175](https://doi.org/10.1016/j.dajour.2023.100175).

**[RA15]** Rossi, Ryan and Ahmed, Nesreen (2015). *The Network Data Repository with Interactive Graph Analytics and Visualization.* Proceedings of the AAAI Conference on Artificial Intelligence, 29(1). DOI: [10.1609/aaai.v29i1.9277](https://doi.org/10.1609/aaai.v29i1.9277).

**[Vega26Hallelujah]** Vega, Frank (2026). *An Approximate Solution to the Minimum Vertex Cover Problem: The Hallelujah Algorithm.* International Journal of Parallel, Emergent and Distributed Systems. Taylor & Francis. DOI: [10.1080/17445760.2026.2660724](https://doi.org/10.1080/17445760.2026.2660724). Accepted for publication.

**[Vega25HallelujahPreprint]** Vega, Frank (2025). *An Approximate Solution to the Minimum Vertex Cover Problem: The Hallelujah Algorithm (Preprint, v2).* Available at: [https://www.preprints.org/manuscript/202510.2392/v2](https://www.preprints.org/manuscript/202510.2392/v2). Public preprint corresponding to the accepted IJPEDS article.

**[Vega25Hvala]** Vega, Frank (2025). *The Hvala Algorithm.* Available at: [https://dev.to/frank_vega_987689489099bf/the-hvala-algorithm-5395](https://dev.to/frank_vega_987689489099bf/the-hvala-algorithm-5395).

**[Vega25Resistire]** Vega, Frank (2025). *The Resistire Experiment.* Available at: [https://dev.to/frank_vega_987689489099bf/the-resistire-experiment-632](https://dev.to/frank_vega_987689489099bf/the-resistire-experiment-632).

**[Vega25Creo]** Vega, Frank (2025). *The Creo Experiment.* Available at: [https://dev.to/frank_vega_987689489099bf/the-creo-experiment-2i1b](https://dev.to/frank_vega_987689489099bf/the-creo-experiment-2i1b).

**[Vega25Gemini]** Vega, Frank (2025). *The Gemini-Vega Validation.* Available at: [https://dev.to/frank_vega_987689489099bf/the-gemini-vega-validation-27i2](https://dev.to/frank_vega_987689489099bf/the-gemini-vega-validation-27i2).

**[NPBench]** Roars. *NP-Complete Benchmark Instances.* Available at: [https://roars.dev/npbench/](https://roars.dev/npbench/). Vertex cover benchmark collection.

**[LargeGraphs]** Cai, Shaowei. *Large Graphs Collection for Vertex Cover Benchmarking.* Available at: [https://lcs.ios.ac.cn/~caisw/graphs.html](https://lcs.ios.ac.cn/~caisw/graphs.html). Network Data Repository collection.

---

**MSC (2020):** 05C69 (Covering and packing), 68Q25 (Analysis of algorithms and problem complexity), 90C27 (Combinatorial optimization), 68W25 (Approximation algorithms)

---

**Documentation**  
Available as PDF at *[An Approximate Solution to the Minimum Vertex Cover Problem: The Hvala Algorithm](https://www.preprints.org/manuscript/202506.0875/v12)*.

The Hvala algorithm is available as a Python package: [https://pypi.org/project/hvala](https://pypi.org/project/hvala)  
Source code and full experimental data are provided in the supplementary materials.
