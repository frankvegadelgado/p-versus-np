# An Approximate Solution to the Minimum Vertex Cover Problem: The Hvala Algorithm

**Author:** Frank Vega  
**Affiliation:** Information Physics Institute, 840 W 67th St, Hialeah, FL 33012, USA  
**Email:** vega.frank@gmail.com  
**ORCID:** 0000-0001-8210-4126

---

## Abstract

We present the **Hvala** algorithm, a linear-time ensemble approximation method for the Minimum Vertex Cover problem. Hvala combines three complementary heuristics — a maximal-matching 2-approximation, a linear-time maximum-degree greedy implemented via a bucket-queue, and the degree-1 weighted-reduction "Hallelujah heuristic" studied in a companion work — with a redundant-vertex pruning post-processing step, and returns the smallest of the four resulting covers.

**Theoretical guarantees.** We prove rigorously that Hvala achieves worst-case approximation ratio $\rho \le 2$ for every finite, simple, undirected graph: the classical maximal-matching component alone already yields this bound, and the pruning step is shown to preserve cover validity while never increasing cover size. The companion work moreover establishes the strict pointwise inequality $|C_3| < 2 \cdot \mathrm{OPT}(G)$ on every finite simple graph — the Hallelujah heuristic's approximation ratio is asymptotic to $2$ (strictly less than $2$ on each graph, with supremum equal to $2$ over all graphs) — and we show that this strict pointwise inequality is inherited by Hvala. Hvala runs in $\mathcal{O}(n+m)$ time and $\mathcal{O}(n+m)$ space.

**Empirical performance.** We validate Hvala on two independent experimental studies totalling $239$ instances. The first uses $109$ vertex-cover instances of the public NPBench collection ($41$ FRB hard instances and $68$ DIMACS clique-complement graphs, both with known optima), completed in $126.97$ seconds: Hvala attains mean approximation ratio $1.021$, with maximum $1.192$ on a single Sanchis adversarial instance. The second evaluates Hvala on $130$ real-world large graphs from the Network Data Repository (Cai's undirected simple graph collection), reaching up to $3$ million vertices and $15$ million edges, completed in approximately $95.5$ minutes of cumulative solve time; on the $51$ instances with published best-known cover sizes, mean ratio is $1.006$ and maximum $1.036$.

**Prospects for a $\sqrt{2}-\epsilon$ bound.** Across the combined $160$ instances with known optima, every approximation ratio lies below $1.414$; $93.8\%$ lie below $1.05$ and $96.9\%$ below $1.10$. The natural open problem we propose as the continuation of this work is whether there exists a *fixed* constant $\epsilon > 0$ such that Hvala achieves uniform ratio $\sqrt{2} - \epsilon$ — either on all graphs (which, by SETH-based hardness, would imply $\mathrm{P} = \mathrm{NP}$) or, more realistically, on broad but restricted graph classes (bounded degree, bounded clique number, bounded treewidth, or structural families such as power-law and expander-like graphs). We do not prove such a bound here and do not claim one holds on all graphs; what we claim is that the combination of rigorous $\le 2$ guarantee, pointwise strict $< 2$ inequality, linear time, and observed ratios uniformly below $1.414$ makes Hvala a plausible vehicle for such a refined analysis. The algorithm is publicly available via PyPI as the `hvala` package.

**Keywords:** Vertex Cover; Approximation Algorithm; Linear-Time Algorithm; Ensemble Heuristic; Graph Optimization; Hardness of Approximation

**MSC:** 05C69, 68Q25, 90C27, 68W25

---

## 1. Introduction

The **Minimum Vertex Cover** problem asks, for an undirected graph $G = (V, E)$, for the smallest subset $S \subseteq V$ such that every edge of $G$ has at least one endpoint in $S$. It is one of Karp's original 21 NP-complete problems [[karp2009reducibility]](#references) and underlies applications in wireless-network design, computational biology, scheduling, and VLSI.

Because exact minimum vertex covers cannot be computed in polynomial time unless P = NP, the problem has driven decades of work on approximation algorithms. The classical 2-approximation obtained by taking both endpoints of every edge of a maximal matching is folklore [[papadimitriou1998combinatorial]](#references); LP-based refinements by Karakostas [[karakostas2009better]](#references) and Karpinski and Zelikovsky [[karpinski1996approximating]](#references) reach factor $2 - \Theta(1/\sqrt{\log n})$, which is $2 - o(1)$ but does not match a constant $2 - \epsilon$. From the hardness side, Dinur and Safra [[dinur2005hardness]](#references) ruled out ratio below $1.3606$ under P $\neq$ NP; Khot, Minzer and Safra [[khot2017independent,dinur2018towards,khot2018pseudorandom]](#references) strengthened this to $\sqrt{2} - \epsilon$ under the Strong Exponential Time Hypothesis (SETH); and, under the Unique Games Conjecture [[khot2002unique]](#references), no constant factor below $2 - \epsilon$ is achievable [[khot2008vertex]](#references). A polynomial-time algorithm with constant ratio $\rho < \sqrt{2}$ would therefore resolve P versus NP, and already an unconditional $2 - \epsilon$ constant is considered beyond reach of current techniques.

**Scope and contribution.** Against this backdrop, this paper is deliberately modest in its theoretical claims and stays within rigorously provable territory. The contributions are:

1. A linear-time ensemble algorithm (**Hvala**, Algorithm 1) that wraps three complementary linear-time heuristics — (i) a maximal-matching 2-approximation, (ii) a bucket-queue max-degree greedy, and (iii) the Hallelujah degree-1 weighted-reduction heuristic [[Vega26Hallelujah]](#references) — inside a redundant-vertex pruning step, and returns the smallest resulting cover.

2. A *rigorous proof* (Theorem 2) that Hvala achieves worst-case approximation ratio $\rho \le 2$ on every finite simple graph. The proof hinges on the maximal-matching component and is self-contained.

3. A *strict pointwise inequality* $|S| < 2 \cdot \mathrm{OPT}(G)$ on every finite simple graph (Corollary 1), inherited from the companion paper [[Vega26Hallelujah]](#references). The Hallelujah heuristic's approximation ratio is asymptotic to $2$ — strictly less than $2$ on each graph, with supremum equal to $2$ — so no constant strictly smaller than $2$ bounds it uniformly; but the pointwise strict inequality on each graph is preserved by the minimum-selection and pruning steps of Hvala.

4. An empirical evaluation on two independent experimental studies totalling $239$ instances: $109$ structured hard instances from the NPBench benchmark collection [[NPBench]](#references) ($41$ FRB hard random graphs and $68$ DIMACS clique-complement graphs, all with known optima) and $130$ real-world large graphs from the Network Data Repository [[RA15]](#references) (biological, social, collaboration, web, infrastructure, and scientific-computing networks, reaching up to $2{,}523{,}386$ vertices and $15{,}245{,}729$ edges), reporting solution quality, running time, and a breakdown by graph family.

---

## 2. Research Data and Implementation

To facilitate reproducibility and community adoption, we developed the open-source Python package *Hvala: Approximate Vertex Cover Solver*, available via the Python Package Index (PyPI) [[Vega25]](#references). This implementation encapsulates the full ensemble algorithm — including the maximal-matching 2-approximation, the bucket-queue max-degree greedy, the Hallelujah degree-1 weighted-reduction subroutine, and the redundant-vertex pruning post-processing step — while guaranteeing an approximation ratio at most $2$ (with pointwise strict inequality $< 2$ on every finite simple graph) through rigorous validation. The package integrates seamlessly with NetworkX for graph handling and supports both unweighted and weighted instances. Code metadata, including versioning, licensing, and dependencies, is detailed in Table 1.

**Table 1.** Code metadata for the Hvala package.

| Nr. | Code metadata description | Metadata |
|-----|--------------------------|----------|
| C1 | Current code version | v0.1.0 |
| C2 | Permanent link to code/repository used for this code version | https://github.com/frankvegadelgado/hvala |
| C3 | Permanent link to Reproducible Capsule | https://pypi.org/project/hvala/ |
| C4 | Legal Code License | MIT License |
| C5 | Code versioning system used | git |
| C6 | Software code languages, tools, and services used | Python |
| C7 | Compilation requirements, operating environments & dependencies | Python ≥ 3.12, NetworkX ≥ 3.4.2 |

---

## 3. The Hvala Algorithm

### 3.1 Overview

Given a simple undirected graph $G = (V, E)$, Hvala first performs trivial preprocessing (remove self-loops and isolated vertices) and then computes four candidate vertex covers:

- **$C_1$ — Maximal-matching cover.** Compute a maximal matching $M$ of $G$ and let $C_1 = \bigcup_{(u,v) \in M}\\{u,v\\}$. This is the classical 2-approximation of [[papadimitriou1998combinatorial]](#references).

- **$C_2$ — Bucket-queue max-degree greedy.** Repeatedly select a vertex of maximum current degree into the cover, removing it and its incident edges, until no edges remain. Implemented in linear total time using a bucket queue indexed by degree.

- **$C_3$ — Hallelujah degree-1 reduction.** Build an auxiliary graph $G'$ by splitting every vertex $u$ of degree $k$ into $k$ auxiliary copies $(u,0), \ldots, (u,k-1)$, each connected to exactly one of $u$'s neighbours, and assigning weight $1/k$ to every such auxiliary vertex. $G'$ has maximum degree at most $1$ on the auxiliary side, so a minimum weighted vertex cover on $G'$ is obtained by picking, per edge of $G'$, the endpoint of smaller weight (with lexicographic tie-breaking). Projecting the selected auxiliary vertices $(u,i)$ back to their original $u$ yields a valid cover of $G$ [[Vega26Hallelujah]](#references).

- **$\tilde{C}_4$ — Pruned union.** Start from $C_1 \cup C_2 \cup C_3$ and apply redundant-vertex pruning (Algorithm 5) *once*, directly yielding the fourth candidate $\tilde{C}_4$.

The candidates $C_1, C_2, C_3$ are then *each* individually passed through redundant-vertex pruning to obtain $\tilde{C}_1, \tilde{C}_2, \tilde{C}_3$, and the algorithm returns the smallest among $\tilde{C}_1, \tilde{C}_2, \tilde{C}_3, \tilde{C}_4$. Note that $C_1$ is included as a **worst-case safety net**: its value is guaranteed to be at most $2 \cdot \mathrm{OPT}$, and since the algorithm returns $\min(|\tilde{C}_1|, |\tilde{C}_2|, |\tilde{C}_3|, |\tilde{C}_4|)$, this guarantee propagates to the final output regardless of how $C_2$, $C_3$, and the union behave.

### 3.2 Main Algorithm

**Algorithm 1:** Hvala: FindVertexCover(G)  
**Input:** Simple undirected graph $G = (V, E)$  
**Output:** Vertex cover $S \subseteq V$ satisfying $|S| \le 2 \cdot \mathrm{OPT}(G)$

```
1:  Remove self-loops and isolated vertices from G
2:  if E(G) = ∅ then return ∅
3:  Build adjacency table adj[v] = N_G(v) for every v ∈ V
4:
5:  // Three base heuristics (all linear time)
6:  C_1 ← MaximalMatchingVC(G)               // Algorithm 2
7:  C_2 ← BucketDegreeGreedy(adj)            // Algorithm 3
8:  C_3 ← HallelujahReduction(G)             // Algorithm 4; [Vega26Hallelujah]
9:
10: // Individual pruning of every candidate
11: for i ∈ {1, 2, 3} do
12:     C̃_i ← PruneRedundant(adj, C_i)
13: end for
14:
15: // Pruned union
16: C̃_4 ← PruneRedundant(adj, C_1 ∪ C_2 ∪ C_3)   // Algorithm 5
17:
18: return argmin_{i ∈ {1,2,3,4}} |C̃_i|
```

### 3.3 Subroutines

**Algorithm 2:** MaximalMatchingVC(G)

```
Input: Graph G
Output: Vertex cover C_1

1: M ← a maximal matching of G
2: C_1 ← ∅
3: for each (u, v) ∈ M do
4:     C_1 ← C_1 ∪ {u, v}
5: end for
6: return C_1
```

**Algorithm 3:** BucketDegreeGreedy(adj)

```
Input: Adjacency table adj of a graph G = (V, E)
Output: Vertex cover C_2

1:  deg[v] ← |adj[v]| for all v ∈ V
2:  Δ ← max_v deg[v]
3:  Create buckets B[0], B[1], ..., B[Δ] as double-ended queues
4:  for each v ∈ V: append v to B[deg[v]]
5:  C_2 ← ∅; removed ← ∅
6:  for d = Δ down to 1 do
7:      while B[d] ≠ ∅ do
8:          pop v from front of B[d]
9:          if v ∈ removed or deg[v] ≠ d: continue
10:         C_2 ← C_2 ∪ {v}; removed ← removed ∪ {v}
11:         for each u ∈ adj[v] \ removed do
12:             deg[u] ← deg[u] − 1
13:             append u to B[deg[u]]
14:         end for
15:     end while
16: end for
17: return C_2
```

**Algorithm 4:** HallelujahReduction(G) — the degree-1 weighted reduction of [[Vega26Hallelujah]](#references)

```
Input: Graph G = (V, E)
Output: Vertex cover C_3

1:  Build an auxiliary graph G' = (V', E'):
2:  for each u ∈ V with degree k > 0 and neighbours v_0, v_1, ..., v_{k-1} do
3:      for i = 0, ..., k-1:
4:          Add auxiliary vertex (u, i) and edge {(u,i), v_i} to G'
5:          Set w((u, i)) ← 1/k
6:  end for
7:  D_uw ← MinVCDegree1(G', uniform weights)   // Algorithm 4a
8:  D_w  ← MinVCDegree1(G', w)
9:  S_uw ← {u : (u,i) ∈ D_uw for some i} ∪ (D_uw ∩ V)
10: S_w  ← {u : (u,i) ∈ D_w  for some i} ∪ (D_w  ∩ V)
11: return smaller of S_uw and S_w
```

**Algorithm 4a:** MinVCDegree1(G', w) — exact weighted VC on a max-degree-1 graph

```
Input: Graph G' of maximum degree 1; weight w
Output: Minimum weighted vertex cover D

1:  D ← ∅; visited ← ∅
2:  for each v ∈ V(G') with v ∉ visited do
3:      if deg(v) = 1 then
4:          u ← unique neighbour of v
5:          if u ∉ visited then
6:              if w(v) < w(u) or (w(v) = w(u) and v < u):
7:                  D ← D ∪ {v}
8:              else:
9:                  D ← D ∪ {u}
10:             visited ← visited ∪ {v, u}
11:         end if
12:     end if
13: end for
14: return D
```

**Algorithm 5:** PruneRedundant(adj, C) — linear-time redundant-vertex pruning

```
Input: Adjacency table adj of G and a vertex cover C of G
Output: A vertex cover C' ⊆ C of G

1: L ← a fixed list copy of C
2: for each v ∈ L do
3:     if every u ∈ adj[v] is currently in C:
4:         C ← C \ {v}
5: end for
6: return C
```

---

## 4. Complexity Analysis

**Theorem 1 (Linear-time and linear-space).** Hvala runs in $\mathcal{O}(n+m)$ time and $\mathcal{O}(n+m)$ space, where $n = |V|$ and $m = |E|$.

*Proof sketch.* Preprocessing and adjacency table construction are $\mathcal{O}(n+m)$.

- *Maximal matching* is computable in $\mathcal{O}(n+m)$ by scanning edges once.
- *Bucket-queue max-degree greedy:* each vertex is inserted into a bucket at most once per degree decrement; total insertions are bounded by $\sum_v \deg(v) = 2m$. Total time: $\mathcal{O}(n+m)$.
- *Hallelujah reduction:* the auxiliary graph $G'$ has exactly $2m$ vertices and $m$ edges. MinVCDegree1 visits every vertex once, running in $\mathcal{O}(n+m)$.
- *Pruning:* checking all neighbours of $v \in C$ costs $\mathcal{O}(\deg(v))$; summed over all $v$ the total is $2m$. Each pruning call is $\mathcal{O}(n+m)$, and a constant number of calls are made.

Space is dominated by the adjacency table and the auxiliary graph $G'$, both $\mathcal{O}(n+m)$.

---

## 5. Approximation Ratio Analysis

We establish the worst-case approximation guarantees of Hvala in two stages.

### 5.1 A Lemma on Redundant-Vertex Pruning

**Lemma 1 (Pruning preserves validity and never increases size).** Let $C$ be a vertex cover of $G$ and let $C' = \textit{PruneRedundant}(\mathrm{adj}, C)$. Then $C' \subseteq C$ and $C'$ is also a vertex cover of $G$.

*Proof.* That $C' \subseteq C$ is immediate. The invariant $C$ is a vertex cover is maintained inductively: at each removal step, vertex $v$ is removed only when every neighbour $u$ of $v$ remains in $C$, so every edge incident to $v$ is still covered by its other endpoint. Edges not incident to $v$ are unaffected. $\square$

As a useful corollary: once $v$ is removed, no neighbour of $v$ can subsequently be removed (the test would fail with $v$ missing), so for every edge at most one endpoint is ever removed.

### 5.2 The Rigorous $\rho \le 2$ Bound

**Theorem 2 (Worst-case 2-approximation).** For every finite simple undirected graph $G$, the output $S$ of FindVertexCover(G) is a vertex cover of $G$ satisfying
$$|S| \le 2 \cdot \mathrm{OPT}(G).$$

*Proof.* Work with $G_0$, the graph after preprocessing; $\mathrm{OPT}(G_0) = \mathrm{OPT}(G)$.

**Step 1: $C_1$ is a vertex cover of size $\le 2 \cdot \mathrm{OPT}(G_0)$.** The maximal matching $M$ covers every edge (no uncovered edge can be added). Since the edges of $M$ are vertex-disjoint, any minimum cover $C^{\ast}$ must contribute at least one endpoint per edge of $M$, so $|C^{\ast}| \ge |M|$. Therefore $|C_1| = 2|M| \le 2|C^{\ast}| = 2 \cdot \mathrm{OPT}(G_0)$.

**Step 2: Pruning does not increase size.** By Lemma 1, $\tilde{C}_1 \subseteq C_1$ with $|\tilde{C}_1| \le |C_1|$, and $\tilde{C}_1$ is still a valid cover.

**Step 3: $S \le \tilde{C}_1$.** The algorithm returns $S = \arg\min_i |\tilde{C}_i|$, so $|S| \le |\tilde{C}_1| \le 2 \cdot \mathrm{OPT}(G)$. $\square$

### 5.3 Inheritance of the Pointwise Strict Inequality

The companion paper [[Vega26Hallelujah]](#references) establishes that for every finite simple undirected graph $G$, the cover $C_3 = \textit{HallelujahReduction}(G)$ satisfies the strict inequality
$$|C_3| < 2 \cdot \mathrm{OPT}(G).$$
The inequality is strict on each graph. At the same time, the supremum of $|C_3|/\mathrm{OPT}(G)$ over all finite simple graphs equals $2$: the Hallelujah ratio is *asymptotic* to $2$, meaning for every $\epsilon > 0$ there exists a graph $G_\epsilon$ with $|C_3|/\mathrm{OPT}(G_\epsilon) > 2 - \epsilon$. Consequently, no single constant strictly less than $2$ uniformly bounds the ratio.

**Corollary 1 (Strict pointwise inequality for Hvala).** For every finite simple undirected graph $G$, the output $S$ of Algorithm 1 satisfies
$$|S| < 2 \cdot \mathrm{OPT}(G).$$
The supremum of $|S|/\mathrm{OPT}(G)$ over all finite simple graphs is equal to $2$: no uniform constant strictly less than $2$ bounds the Hvala ratio.

*Proof.* Let $\tilde{C}_3 = \textit{PruneRedundant}(\mathrm{adj}, C_3)$. By Lemma 1, $|\tilde{C}_3| \le |C_3| < 2 \cdot \mathrm{OPT}(G)$. Since $S = \arg\min_i |\tilde{C}_i|$, we have $|S| \le |\tilde{C}_3| < 2 \cdot \mathrm{OPT}(G)$. The supremum argument follows from the asymptotic-to $2$ property of Hallelujah and the observation that a uniform constant $2 - \delta$ bounding $|S|$ would contradict the absence of such a constant for Hallelujah alone. $\square$

Two clarifying remarks: Theorem 2 is *not* rendered redundant by Corollary 1 — the theorem gives a self-contained $\le 2$ bound with no dependency on [[Vega26Hallelujah]](#references), while the corollary provides pointwise strictness at the cost of relying on the companion paper. The strict inequality in Corollary 1 is pointwise only: it does not provide a uniform constant $2 - \epsilon$, and no such constant is known.

### 5.4 Roles of $C_2$ and $\tilde{C}_4$

Neither $C_2$ nor $\tilde{C}_4$ is required for Theorem 2 or Corollary 1.

- **$C_2$ (bucket-queue max-degree greedy)** has no general worst-case ratio better than $\Theta(\log \Delta)$ (Johnson's classical bound), but is very strong on near-regular and clique-like graphs and is frequently optimal or near-optimal there. Its presence in the minimum cannot worsen the bound.
- **$\tilde{C}_4$** is a pruned version of the union $C_1 \cup C_2 \cup C_3$; it occasionally exploits structural overlaps between the three base heuristics that per-candidate pruning alone cannot resolve.

---

## 6. Experimental Validation

We evaluate Hvala on two independent experimental studies totalling $239$ instances, using the same implementation of Algorithm 1 on commodity hardware (Intel Core i7-1165G7 at 2.80 GHz, 32 GB RAM, single-threaded Python 3.12 with NetworkX 3.4.2).

### 6.1 Experiment 1: Structured Hard Instances (NPBench)

#### Setup

We evaluate Hvala on $109$ vertex-cover instances of the public NPBench collection [[NPBench]](#references), comprising two families for which the optimum (or a tight best-known value) is publicly available:

1. **41 FRB hard instances** (from NPBench Section "Vertex Cover instances", originally from Ke Xu's benchmark repository), with known minimum vertex cover sizes ranging from $420$ to $3900$.
2. **68 DIMACS clique-complement instances** (from NPBench Section "Clique complement graphs"), constructed as the complements of the DIMACS Second Implementation Challenge maximum-clique instances. The optimum vertex cover of the complement equals $n - \omega(G)$, where $\omega(G)$ is the maximum clique size; values compiled from Mascia's DIMACS benchmark page [[DIMACSClique]](#references). For `C500.9` and `C1000.9`, the clique number is not known to be tight, so the reported ratio is a lower bound (marked $^\dagger$).

Cumulative solve time over all $109$ instances is $126.97$ seconds.

#### FRB Benchmark Results (41 instances)

| Instance | Known OPT | Hvala size | Time | Ratio |
|----------|-----------|------------|------|-------|
| frb30-15-1 | 420 | 428 | 214.1ms | 1.019 |
| frb30-15-2 | 420 | 429 | 219.7ms | 1.021 |
| frb30-15-3 | 420 | 427 | 206.9ms | 1.017 |
| frb30-15-4 | 420 | 429 | 1.45s | 1.021 |
| frb30-15-5 | 420 | 427 | 249.9ms | 1.017 |
| frb35-17-1 | 560 | 569 | 403.9ms | 1.016 |
| frb35-17-2 | 560 | 570 | 443.8ms | 1.018 |
| frb35-17-3 | 560 | 568 | 455.1ms | 1.014 |
| frb35-17-4 | 560 | 570 | 445.5ms | 1.018 |
| frb35-17-5 | 560 | 568 | 484.9ms | 1.014 |
| frb40-19-1 | 720 | 730 | 716.4ms | 1.014 |
| frb40-19-2 | 720 | 730 | 705.6ms | 1.014 |
| frb40-19-3 | 720 | 731 | 725.1ms | 1.015 |
| frb40-19-4 | 720 | 732 | 756.8ms | 1.017 |
| frb40-19-5 | 720 | 730 | 759.7ms | 1.014 |
| frb45-21-1 | 900 | 912 | 1.07s | 1.013 |
| frb45-21-2 | 900 | 911 | 1.19s | 1.012 |
| frb45-21-3 | 900 | 912 | 1.16s | 1.013 |
| frb45-21-4 | 900 | 912 | 1.10s | 1.013 |
| frb45-21-5 | 900 | 912 | 1.18s | 1.013 |
| frb50-23-1 | 1100 | 1111 | 1.57s | 1.010 |
| frb50-23-2 | 1100 | 1113 | 1.71s | 1.012 |
| frb50-23-3 | 1100 | 1117 | 1.74s | 1.015 |
| frb50-23-4 | 1100 | 1113 | 1.65s | 1.012 |
| frb50-23-5 | 1100 | 1112 | 1.69s | 1.011 |
| frb53-24-1 | 1219 | 1235 | 1.95s | 1.013 |
| frb53-24-2 | 1219 | 1234 | 2.09s | 1.012 |
| frb53-24-3 | 1219 | 1235 | 2.09s | 1.013 |
| frb53-24-4 | 1219 | 1232 | 2.14s | 1.011 |
| frb53-24-5 | 1219 | 1235 | 2.00s | 1.013 |
| frb56-25-1 | 1344 | 1358 | 2.48s | 1.010 |
| frb56-25-2 | 1344 | 1358 | 2.44s | 1.010 |
| frb56-25-3 | 1344 | 1359 | 2.33s | 1.011 |
| frb56-25-4 | 1344 | 1358 | 2.46s | 1.010 |
| frb56-25-5 | 1344 | 1361 | 2.54s | 1.013 |
| frb59-26-1 | 1475 | 1492 | 2.81s | 1.012 |
| frb59-26-2 | 1475 | 1492 | 2.98s | 1.012 |
| frb59-26-3 | 1475 | 1494 | 4.15s | 1.013 |
| frb59-26-4 | 1475 | 1493 | 4.38s | 1.012 |
| frb59-26-5 | 1475 | 1491 | 4.56s | 1.011 |
| frb100-40 | 3900 | 3931 | 16.82s | 1.008 |

#### DIMACS Clique-Complement Benchmark Results — Part 1 (brock, c-fat, C, gen, hamming)

$^\dagger$: best-known upper bound on OPT (maximum clique not confirmed optimal); the ratio is then a lower bound.

| Instance | Known OPT | Hvala size | Time | Ratio |
|----------|-----------|------------|------|-------|
| brock200_1 | 179 | 183 | 59.3ms | 1.022 |
| brock200_2 | 188 | 192 | 134.8ms | 1.021 |
| brock200_3 | 185 | 189 | 101.5ms | 1.022 |
| brock200_4 | 183 | 187 | 65.6ms | 1.022 |
| brock400_1 | 373 | 381 | 253.4ms | 1.021 |
| brock400_2 | 371 | 379 | 262.8ms | 1.022 |
| brock400_3 | 369 | 381 | 254.5ms | 1.033 |
| brock400_4 | 367 | 378 | 259.6ms | 1.030 |
| brock800_1 | 777 | 785 | 2.12s | 1.010 |
| brock800_2 | 776 | 783 | 2.41s | 1.009 |
| brock800_3 | 775 | 784 | 2.30s | 1.012 |
| brock800_4 | 774 | 785 | 3.39s | 1.014 |
| c-fat200-1 | 188 | 188 | 450.9ms | **1.000** |
| c-fat200-2 | 176 | 176 | 809.6ms | **1.000** |
| c-fat200-5 | 142 | 142 | 266.4ms | **1.000** |
| c-fat500-1 | 486 | 486 | 2.19s | **1.000** |
| c-fat500-10 | 374 | 374 | 1.62s | **1.000** |
| c-fat500-2 | 474 | 474 | 2.22s | **1.000** |
| c-fat500-5 | 436 | 436 | 2.19s | **1.000** |
| C1000.9 | 932† | 945 | 1.05s | ≥1.014 |
| C125.9 | 91 | 92 | 10.0ms | 1.011 |
| C250.9 | 206 | 214 | 31.9ms | 1.039 |
| C500.9 | 443† | 454 | 276.9ms | ≥1.025 |
| gen200_p0.9_44 | 156 | 167 | 22.6ms | 1.071 |
| gen200_p0.9_55 | 145 | 163 | 22.3ms | 1.124 |
| gen400_p0.9_55 | 345 | 356 | 83.6ms | 1.032 |
| gen400_p0.9_65 | 335 | 359 | 129.0ms | 1.072 |
| gen400_p0.9_75 | 325 | 358 | 122.6ms | 1.102 |
| hamming10-2 | 512 | 512 | 60.9ms | **1.000** |
| hamming10-4 | 984 | 992 | 1.73s | 1.008 |
| hamming6-2 | 32 | 32 | 2.1ms | **1.000** |
| hamming6-4 | 60 | 60 | 21.0ms | **1.000** |
| hamming8-2 | 128 | 128 | 14.2ms | **1.000** |
| hamming8-4 | 240 | 240 | 129.3ms | **1.000** |

#### DIMACS Clique-Complement Benchmark Results — Part 2 (johnson, keller, MANN, p_hat, san, sanr)

| Instance | Known OPT | Hvala size | Time | Ratio |
|----------|-----------|------------|------|-------|
| johnson16-2-4 | 112 | 112 | 18.2ms | **1.000** |
| johnson32-2-4 | 480 | 480 | 365.8ms | **1.000** |
| johnson8-2-4 | 24 | 24 | 2.2ms | **1.000** |
| johnson8-4-4 | 56 | 56 | 7.3ms | **1.000** |
| keller4 | 160 | 162 | 77.4ms | 1.012 |
| keller5 | 749 | 759 | 1.34s | 1.013 |
| MANN_a27 | 252 | 253 | 15.9ms | 1.004 |
| MANN_a45 | 690 | 695 | 27.1ms | 1.007 |
| MANN_a81 | 2221 | 2225 | 82.8ms | 1.002 |
| MANN_a9 | 29 | 29 | 0.0ms | **1.000** |
| p_hat1000-3 | 932 | 942 | 2.79s | 1.011 |
| p_hat300-1 | 292 | 292 | 751.1ms | **1.000** |
| p_hat300-2 | 275 | 277 | 362.0ms | 1.007 |
| p_hat300-3 | 264 | 268 | 170.7ms | 1.015 |
| p_hat500-1 | 491 | 492 | 1.61s | 1.002 |
| p_hat500-2 | 464 | 469 | 1.21s | 1.011 |
| p_hat500-3 | 450 | 454 | 560.2ms | 1.009 |
| p_hat700-1 | 689 | 693 | 3.78s | 1.006 |
| p_hat700-2 | 656 | 658 | 2.82s | 1.003 |
| p_hat700-3 | 638 | 642 | 1.31s | 1.006 |
| san200_0.7_1 | 170 | 184 | 59.4ms | 1.082 |
| san200_0.7_2 | 182 | 188 | 201.1ms | 1.033 |
| san200_0.9_1 | 130 | 155 | 21.2ms | **1.192** |
| san200_0.9_2 | 140 | 162 | 19.2ms | 1.157 |
| san200_0.9_3 | 156 | 169 | 26.3ms | 1.083 |
| san400_0.5_1 | 387 | 393 | 609.8ms | 1.016 |
| san400_0.7_1 | 360 | 379 | 443.1ms | 1.053 |
| san400_0.7_2 | 370 | 385 | 364.2ms | 1.041 |
| san400_0.7_3 | 378 | 388 | 365.2ms | 1.026 |
| san400_0.9_1 | 300 | 348 | 86.3ms | 1.160 |
| sanr200_0.7 | 182 | 184 | 124.4ms | 1.011 |
| sanr200_0.9 | 158 | 160 | 21.7ms | 1.013 |
| sanr400_0.5 | 387 | 388 | 718.1ms | 1.003 |
| sanr400_0.7 | 379 | 383 | 990.0ms | 1.011 |

#### NPBench Summary Statistics

Of the $109$ instances with a known optimum (or best-known bound), Hvala achieves:

- **Mean approximation ratio:** $1.021$ (FRB block: $1.014$; DIMACS clique-complement block: $1.025$).
- **Exact optimality:** $18$ instances solved with ratio $1.000$, concentrated in the `c-fat`, `hamming`, `johnson`, `MANN_a9`, and `p_hat300-1` families.
- **Maximum ratio observed:** $1.192$ on `san200_0.9_1` (a Sanchis instance constructed with an embedded clique of size $70$). The five worst ratios are all on Sanchis `san`/`gen` adversarial instances, which are specifically engineered to hide large cliques; on these dense, small, carefully constructed graphs, ensemble heuristics are known to degrade relative to specialised exact solvers.
- **Runtime:** total cumulative solve time across all $109$ instances is $126.97$ seconds ($80.53$ s on FRB + $46.44$ s on DIMACS). Per-instance times range from under $10$ ms (smallest graphs) to $16.82$ s (`frb100-40`, $n = 4000$).

Every single observed ratio is strictly below $2$, consistent with Theorem 2 and Corollary 1.

### 6.2 Experiment 2: Real-World Large Graphs

#### Setup

This section presents comprehensive experimental results of the Hvala algorithm on real-world large graphs from the Network Data Repository [[RA15]](#references). The benchmark suite consists of $130$ instances from the complete collection of $139$ undirected simple largest graphs distributed by Cai [[RA15]](#references). Nine instances are excluded: three graphs (`ca-hollywood-2009`, `socfb-uci-uni`, `soc-orkut`) exceed the $32$ GB RAM limit of our test hardware when loaded through NetworkX, and six further graphs (`inf-road-usa`, `sc-ldoor`, `soc-livejournal`, `soc-pokec`, `socfb-A-anon`, `socfb-B-anon`) were dropped to keep the experiment tractable within a single session. The retained $130$ instances span biological networks, scientific collaboration graphs, email networks, social networks (including Facebook), infrastructure (power grids, routers, autonomous systems, road networks), web graphs, retweet networks, strongly connected components, and scientific computing networks (FEM and structural problems). Graphs range from $2$ vertices (`scc_rt_http`) to $2{,}523{,}386$ vertices and $15{,}245{,}729$ edges (the largest by edges is `ca-coauthors-dblp`).

Because the Network Data Repository does not provide certified minimum vertex cover values for most of these instances, we rely on the *best-known approximate optimum* values compiled by the Milagro experiment [[Frank25]](#references) on the same collection. For $51$ of the $130$ instances such a reference value is available (of which $29$ are certified optima on tree-like components); for the remaining $79$ instances we list "Unknown". Every returned cover satisfies $|S| < 2 \cdot \mathrm{OPT}$ by Theorem 2 and Corollary 1.

#### Real-World Large Graph Results (130 instances)

| Instance | Category | Best Known | Hvala size | Time | Ratio |
|----------|----------|-----------|------------|------|-------|
| bio-celegans | Bio | 248 | 257 | 30.3ms | 1.036 |
| bio-diseasome | Bio | 283 | 285 | 18.7ms | 1.007 |
| bio-dmela | Bio | -- | 2672 | 495.3ms | -- |
| bio-yeast | Bio | 453 | 464 | 57.5ms | 1.024 |
| ca-AstroPh | Collab | -- | 11512 | 6.05s | -- |
| ca-citeseer | Collab | -- | 129274 | 22.44s | -- |
| ca-coauthors-dblp | Collab | -- | 472272 | 757.0s | -- |
| ca-CondMat | Collab | -- | 12500 | 4.02s | -- |
| ca-CSphd | Collab | 548 | 553 | 79.1ms | 1.009 |
| ca-dblp-2010 | Collab | -- | 122072 | 28.83s | -- |
| ca-dblp-2012 | Collab | -- | 165085 | 31.50s | -- |
| ca-Erdos992 | Collab | 459 | 461 | 142.1ms | 1.004 |
| ca-GrQc | Collab | -- | 2213 | 254.4ms | -- |
| ca-HepPh | Collab | -- | 6568 | 49.94s | -- |
| ca-MathSciNet | Collab | -- | 140428 | 41.45s | -- |
| ca-netscience | Collab | 212 | 214 | 40.1ms | 1.009 |
| ia-email-EU | Email | -- | 820 | 1.50s | -- |
| ia-email-univ | Email | 603 | 609 | 124.4ms | 1.010 |
| ia-enron-large | Social | -- | 12820 | 6.52s | -- |
| ia-enron-only | Social | 86 | 87 | 21.0ms | 1.012 |
| ia-fb-messages | Social | 578 | 593 | 111.6ms | 1.026 |
| ia-infect-dublin | Social | 295 | 295 | 47.3ms | 1.000 |
| ia-infect-hyper | Social | 91 | 93 | 60.3ms | 1.022 |
| ia-reality | Social | -- | 81 | 123.2ms | -- |
| ia-wiki-Talk | Wiki | -- | 17407 | 16.52s | -- |
| inf-power | Infra | -- | 2267 | 291.9ms | -- |
| inf-roadNet-CA | Infra | -- | 1058991 | 122.5s | -- |
| inf-roadNet-PA | Infra | -- | 587209 | 72.8s | -- |
| rec-amazon | Rec | -- | 48622 | 5.36s | -- |
| rt-retweet | Retweet | 31 | 32 | 5.2ms | 1.032 |
| rt-retweet-crawl | Retweet | -- | 81211 | 143.8s | -- |
| rt-twitter-copen | Retweet | 235 | 238 | 42.9ms | 1.013 |
| sc-msdoor | SciComp | -- | 382184 | 400.1s | -- |
| sc-nasasrb | SciComp | -- | 51559 | 65.1s | -- |
| sc-pkustk11 | SciComp | -- | 84149 | 111.2s | -- |
| sc-pkustk13 | SciComp | -- | 89759 | 124.6s | -- |
| sc-pwtk | SciComp | -- | 208297 | 221.8s | -- |
| sc-shipsec1 | SciComp | -- | 119415 | 82.9s | -- |
| sc-shipsec5 | SciComp | -- | 148790 | 99.6s | -- |
| scc_enron-only | SCC | 137 | 138 | 197.9ms | 1.007 |
| scc_fb-forum | SCC | 370 | 372 | 1.96s | 1.005 |
| scc_fb-messages | SCC | -- | 1072 | 27.78s | -- |
| scc_infect-dublin | SCC | -- | 9124 | 8.70s | -- |
| scc_infect-hyper | SCC | 109 | 110 | 155.0ms | 1.009 |
| scc_reality | SCC | -- | 2486 | 193.9s | -- |
| scc_retweet | SCC | -- | 564 | 1.02s | -- |
| scc_retweet-crawl | SCC | -- | 8435 | 492.2ms | -- |
| scc_rt_alwefaq | SCC | 35 | 35 | 7.6ms | **1.000** |
| scc_rt_assad | SCC | 16 | 16 | 3.3ms | **1.000** |
| scc_rt_bahrain | SCC | 37 | 37 | 2.9ms | **1.000** |
| scc_rt_barackobama | SCC | 29 | 29 | 3.3ms | **1.000** |
| scc_rt_damascus | SCC | 15 | 15 | 1.1ms | **1.000** |
| scc_rt_dash | SCC | 15 | 15 | 1.1ms | **1.000** |
| scc_rt_gmanews | SCC | 46 | 46 | 15.2ms | **1.000** |
| scc_rt_gop | SCC | 6 | 6 | 0.0ms | **1.000** |
| scc_rt_http | SCC | 2 | 2 | 0.0ms | **1.000** |
| scc_rt_israel | SCC | 11 | 11 | 0.0ms | **1.000** |
| scc_rt_justinbieber | SCC | 26 | 26 | 5.2ms | **1.000** |
| scc_rt_ksa | SCC | 12 | 12 | 0.5ms | **1.000** |
| scc_rt_lebanon | SCC | 5 | 5 | 0.0ms | **1.000** |
| scc_rt_libya | SCC | 12 | 12 | 1.3ms | **1.000** |
| scc_rt_lolgop | SCC | 103 | 103 | 52.3ms | **1.000** |
| scc_rt_mittromney | SCC | 42 | 42 | 1.6ms | **1.000** |
| scc_rt_obama | SCC | 4 | 4 | 0.0ms | **1.000** |
| scc_rt_occupy | SCC | 22 | 22 | 1.1ms | **1.000** |
| scc_rt_occupywallstnyc | SCC | 45 | 45 | 12.1ms | **1.000** |
| scc_rt_oman | SCC | 6 | 6 | 0.0ms | **1.000** |
| scc_rt_onedirection | SCC | 29 | 29 | 4.0ms | **1.000** |
| scc_rt_p2 | SCC | 12 | 12 | 0.0ms | **1.000** |
| scc_rt_qatif | SCC | 5 | 5 | 0.0ms | **1.000** |
| scc_rt_saudi | SCC | 17 | 17 | 1.0ms | **1.000** |
| scc_rt_tcot | SCC | 12 | 12 | 1.0ms | **1.000** |
| scc_rt_tlot | SCC | 6 | 6 | 0.6ms | **1.000** |
| scc_rt_uae | SCC | 8 | 8 | 1.0ms | **1.000** |
| scc_rt_voteonedirection | SCC | 4 | 4 | 0.0ms | **1.000** |
| scc_twitter-copen | SCC | -- | 1328 | 20.24s | -- |
| soc-BlogCatalog | Social | -- | 20967 | 69.1s | -- |
| soc-brightkite | Social | -- | 21473 | 10.30s | -- |
| soc-buzznet | Social | -- | 31059 | 93.6s | -- |
| soc-delicious | Social | -- | 86810 | 48.30s | -- |
| soc-digg | Social | -- | 104237 | 217.9s | -- |
| soc-dolphins | Social | 34 | 35 | 3.2ms | 1.029 |
| soc-douban | Social | -- | 8685 | 24.07s | -- |
| soc-epinions | Social | -- | 9858 | 3.09s | -- |
| soc-flickr | Social | -- | 154387 | 107.8s | -- |
| soc-flixster | Social | -- | 96404 | 283.6s | -- |
| soc-FourSquare | Social | -- | 90524 | 127.9s | -- |
| soc-gowalla | Social | -- | 85360 | 35.31s | -- |
| soc-karate | Social | 14 | 14 | 1.1ms | **1.000** |
| soc-lastfm | Social | -- | 78832 | 164.7s | -- |
| soc-LiveMocha | Social | -- | 44146 | 79.9s | -- |
| soc-slashdot | Social | -- | 22632 | 16.07s | -- |
| soc-twitter-follows | Social | -- | 2323 | 24.34s | -- |
| soc-wiki-Vote | Social | 404 | 410 | 39.8ms | 1.015 |
| soc-youtube | Social | -- | 148135 | 64.9s | -- |
| soc-youtube-snap | Social | -- | 279062 | 100.8s | -- |
| socfb-Berkeley13 | Facebook | -- | 17487 | 35.10s | -- |
| socfb-CMU | Facebook | -- | 5061 | 8.45s | -- |
| socfb-Duke14 | Facebook | -- | 7790 | 15.06s | -- |
| socfb-Indiana | Facebook | -- | 23741 | 44.05s | -- |
| socfb-MIT | Facebook | -- | 4726 | 8.26s | -- |
| socfb-OR | Facebook | -- | 37209 | 25.68s | -- |
| socfb-Penn94 | Facebook | -- | 31723 | 48.15s | -- |
| socfb-Stanford3 | Facebook | -- | 8611 | 19.07s | -- |
| socfb-Texas84 | Facebook | -- | 28669 | 55.17s | -- |
| socfb-UCLA | Facebook | -- | 15494 | 24.95s | -- |
| socfb-UConn | Facebook | -- | 13436 | 18.95s | -- |
| socfb-UCSB37 | Facebook | -- | 11481 | 14.06s | -- |
| socfb-UF | Facebook | -- | 27775 | 52.03s | -- |
| socfb-UIllinois | Facebook | -- | 24465 | 40.99s | -- |
| socfb-Wisconsin87 | Facebook | -- | 18716 | 28.95s | -- |
| tech-as-caida2007 | Tech | -- | 3699 | 1.07s | -- |
| tech-as-skitter | Tech | -- | 529662 | 365.1s | -- |
| tech-internet-as | Tech | -- | 5718 | 1.81s | -- |
| tech-p2p-gnutella | Tech | -- | 15730 | 3.53s | -- |
| tech-RL-caida | Tech | -- | 75568 | 14.69s | -- |
| tech-routers-rf | Tech | 793 | 801 | 94.7ms | 1.010 |
| tech-WHOIS | Tech | -- | 2297 | 964.5ms | -- |
| web-arabic-2005 | Web | -- | 115297 | 62.7s | -- |
| web-BerkStan | Web | -- | 5404 | 336.0ms | -- |
| web-edu | Web | 1449 | 1451 | 90.4ms | 1.001 |
| web-google | Web | 497 | 498 | 40.3ms | 1.002 |
| web-indochina-2004 | Web | -- | 7363 | 778.7ms | -- |
| web-it-2004 | Web | -- | 415230 | 182.0s | -- |
| web-polblogs | Web | 243 | 245 | 28.2ms | 1.008 |
| web-sk-2005 | Web | -- | 58411 | 6.32s | -- |
| web-spam | Web | -- | 2344 | 574.6ms | -- |
| web-uk-2005 | Web | -- | 127774 | 316.9s | -- |
| web-webbase-2001 | Web | -- | 2665 | 425.0ms | -- |
| web-wikipedia2009 | Web | -- | 659409 | 192.1s | -- |

#### Real-World Summary Statistics

Across the $130$ real-world instances ($51$ with known reference optima), Hvala achieves:

- **Mean approximation ratio:** $1.006$ across the $51$ instances with known optima.
- **Minimum ratio:** $1.000$ (reached on $30$ instances).
- **Maximum ratio:** $1.036$ on `bio-celegans` ($C$. elegans metabolic network; Hvala size $257$ vs. best-known $248$). All $51$ ratios lie below $1.05$.
- **Scale:** the largest instance by vertex count is `soc-flixster`, a movie-rating graph with $2{,}523{,}386$ vertices and $7{,}918{,}801$ edges, solved in $4.73$ minutes; the largest by edge count is `ca-coauthors-dblp` with $540{,}486$ vertices and $15{,}245{,}729$ edges, solved in $12.62$ minutes (the longest solve of the experiment). Other large instances include `tech-as-skitter` ($1{,}694{,}616$ vertices, $11{,}094{,}209$ edges, $6.08$ min), `web-wikipedia2009` ($1{,}864{,}433$ vertices, $4{,}507{,}315$ edges, $3.20$ min), and `inf-roadNet-CA` ($1{,}957{,}027$ vertices, $2{,}760{,}388$ edges, $2.04$ min).
- **Runtime distribution:** $60$ of the $130$ instances are solved in under $1$ second; $43$ between $1$ and $60$ seconds; $26$ between $1$ and $10$ minutes; exactly one (`ca-coauthors-dblp`) between $10$ and $60$ minutes; none exceed one hour.
- **Total wall-clock time:** **$5{,}732$ seconds** cumulative ($\approx 95.5$ minutes, or $1.59$ hours) for all $130$ real-world instances.
- **Linear-time scalability:** per-vertex amortised cost stays within a narrow range across five orders of magnitude of graph size, consistent with $\mathcal{O}(n+m)$ complexity.

---

## 7. Discussion

### 7.1 Empirical vs. Theoretical Gap

Theorem 2 and Corollary 1 give a uniform $\le 2$ bound and a pointwise strict $< 2$ bound, with the supremum over all graphs equal to $2$. Across the two experimental studies ($239$ instances total), the empirical ratios on the $160$ instances with known optima are far below this worst-case: combined mean $1.016$, combined maximum $1.192$ on a single adversarially-constructed Sanchis instance, and every ratio strictly below $\sqrt{2} \approx 1.414$. On one real-world instance (`ia-infect-dublin`) Hvala even *improved* the previously published best-known cover size from $296$ to $295$. Closing the gap — either by refining the analysis as a function of graph parameters (average degree, girth, treewidth, clique number) or by constructing adversarial instances that drive Hvala closer to $2$ — is a natural target for further work. We note in particular that the corrected runtime for `inf-roadNet-CA` ($1{,}957{,}027$ vertices, $122.50$ seconds $= 2.04$ minutes) is consistent with the $\mathcal{O}(n+m)$ scaling observed across the full benchmark.

### 7.2 Hardness Barriers

The hardness results [[dinur2005hardness,khot2008vertex,khot2017independent,dinur2018towards,khot2018pseudorandom]](#references) make clear that no unconditional polynomial-time algorithm is known to achieve uniform constant ratio below $2 - \epsilon$ for any fixed $\epsilon > 0$, and ratio below $\sqrt{2}$ is SETH-hard. Hvala does not aim to cross these barriers; it aims to match the $\le 2$ bound constructively, in linear time, and to inherit the pointwise strict $< 2$ inequality from the Hallelujah heuristic [[Vega26Hallelujah]](#references). The ensemble and pruning exploit structural orthogonality empirically, which accounts for the mean approximation ratio of $1.021$ on the $109$ NPBench instances with known optima and the combined mean ratio of $1.016$ across all $160$ known-optimum instances (NPBench and real-world combined), without contradicting any hardness result.

### 7.3 Prospects for a $\sqrt{2} - \epsilon$ Bound

The most interesting empirical regularity across both studies is that *every single ratio observed on the* $160$ *instances with known optima stays below* $1.414$ — with the combined maximum being $1.192$, on a narrow family of Sanchis adversarial graphs; $93.8\%$ of the $160$ instances are within ratio $1.05$, $96.9\%$ are within ratio $1.10$, and $100\%$ are within ratio $1.20$. The question we wish to pose is whether the ratio of Hvala can be bounded uniformly by $\sqrt{2} - \epsilon$ for a *fixed* constant $\epsilon > 0$, not merely strictly below $\sqrt{2}$ in the same asymptotic sense.

Under SETH and the hardness results of Khot, Minzer and Safra [[khot2017independent,dinur2018towards,khot2018pseudorandom]](#references), no polynomial-time algorithm can achieve uniform ratio $\sqrt{2} - \epsilon$ on *all* finite graphs unless $\mathrm{P} = \mathrm{NP}$. Hvala's empirical behaviour cannot, on its own, imply a uniform $\sqrt{2} - \epsilon$ guarantee. What it does suggest is that Hvala — and specifically the Hallelujah weighted-reduction component [[Vega26Hallelujah]](#references) — is a plausible candidate for a refined worst-case analysis on *restricted but broad* graph classes: graphs of bounded maximum degree, graphs with bounded clique number, graphs with bounded treewidth, or graphs drawn from structural families (power-law, expander-like, small-world) common in practice. Such a restricted-class result would not contradict any known hardness barrier and would be of substantial theoretical and practical interest.

Three observations support this interpretation. First, the $18$ instances solved to exact optimality on NPBench are concentrated in highly-structured families (`c-fat`, `hamming`, `johnson`, `MANN_a9`, `p_hat300-1`), indicating that the Hallelujah reduction captures optimal structure on graphs where regularity is high. Second, the worst-case ratios on NPBench occur exclusively on the Sanchis `san`/`gen` hidden-clique adversarial construction; on no other NPBench family does Hvala exceed ratio $1.08$. Third, on the $130$ real-world large graphs — including social, collaboration, web, biological, infrastructure, and scientific-computing networks at scales up to $2.5$ million vertices and $15$ million edges — Hvala's output always sits strictly below $2 \cdot \mathrm{OPT}$, the algorithm's linear-time scaling holds in practice, and on at least one instance (`ia-infect-dublin`) Hvala improves the previously published best-known cover size, demonstrating that the ensemble can tighten, not just approximate, standing reference values.

We therefore position this paper as a step towards, rather than a proof of, a uniform $\sqrt{2} - \epsilon$ guarantee (with fixed $\epsilon > 0$) on restricted graph classes. We do *not* claim a uniform $\sqrt{2} - \epsilon$ bound here. What we claim is that Hvala is the first simple linear-time algorithm for Minimum Vertex Cover whose combined theoretical properties (rigorous $\le 2$, pointwise strict $< 2$, asymptotic-to $2$ supremum) and empirical behaviour (ratios staying below $1.414$ across $239$ diverse instances, $160$ of which with known optima, and cumulative real-world solve time of approximately $95.5$ minutes for $130$ instances reaching up to $2.5$ million vertices) jointly make it a plausible vehicle for further theoretical work on the $\sqrt{2} - \epsilon$ threshold.

### 7.4 Comparison to Other Practical Methods

Advanced local-search methods such as FastVC [[cai2017finding]](#references), TIVC [[zhang2023tivc]](#references), and MetaVC2 [[luo2019local]](#references) reach empirical ratios comparable to Hvala's on DIMACS-style benchmarks, typically at the price of longer run times and without a simple constructive worst-case guarantee. Parameterised FPT algorithms [[harris2024faster]](#references) are exact for small solution sizes $k$, complementing rather than competing with Hvala's regime of large, general graphs. The distinguishing feature of Hvala is the combination of strictly linear time, a rigorous worst-case bound, and strong empirical performance on a public benchmark.

---

## 8. Conclusion

We have presented Hvala, a linear-time ensemble algorithm for Minimum Vertex Cover combining a maximal-matching 2-approximation, a bucket-queue max-degree greedy, and the Hallelujah degree-1 weighted reduction of [[Vega26Hallelujah]](#references), wrapped inside a redundant-vertex pruning step. We proved rigorously that the algorithm achieves the uniform worst-case ratio $\rho \le 2$ (Theorem 2) and, combining with the companion paper [[Vega26Hallelujah]](#references), the strict pointwise inequality $|S| < 2 \cdot \mathrm{OPT}(G)$ on every finite simple graph (Corollary 1) — the ratio is asymptotic to $2$: strictly less than $2$ on each graph, with supremum equal to $2$. Hvala runs in $\mathcal{O}(n+m)$ time and $\mathcal{O}(n+m)$ space (Theorem 1).

We validated Hvala on two independent experimental studies totalling $239$ instances. On the $109$ instances with known optima from the NPBench vertex-cover collection (Experiment 1, Section 6.1, $126.97$ seconds cumulative), Hvala solves $18$ to proven optimality and attains a mean approximation ratio of $1.021$. On the $130$ real-world large graphs from the Network Data Repository (Experiment 2, Section 6.2), ranging up to $2{,}523{,}386$ vertices and $15{,}245{,}729$ edges, Hvala completes the entire benchmark in approximately $5{,}732$ seconds ($\approx 95.5$ minutes) cumulative — $60$ of the $130$ instances finish in under $1$ second, and the largest instance (`ca-coauthors-dblp`, $540{,}486$ vertices, $15{,}245{,}729$ edges) is solved in $12.62$ minutes — every returned cover is guaranteed by Corollary 1 to be strictly less than $2 \cdot \mathrm{OPT}$.

Across both studies, every single approximation ratio observed on the $160$ instances with known optima stays below $1.414$, with the maximum being $1.192$ on a narrow family of Sanchis adversarial graphs. This motivates the central open problem we propose as the natural continuation of this work:

> Is there a fixed constant $\epsilon > 0$ such that, for every finite simple undirected graph $G$, the Hvala algorithm achieves approximation ratio $|S|/\mathrm{OPT}(G) \le \sqrt{2} - \epsilon \approx 1.414 - \epsilon$ — or, failing that, does such a uniform bound hold on broad but restricted graph classes (bounded degree, bounded clique number, bounded treewidth, or structural families such as power-law and expander-like graphs)?

A fixed-constant $\sqrt{2} - \epsilon$ bound would either yield a uniform sub $\sqrt{2}$ guarantee on all graphs (which, by the SETH-based hardness of Khot, Minzer and Safra [[khot2017independent,dinur2018towards,khot2018pseudorandom]](#references), would imply $\mathrm{P} = \mathrm{NP}$) or, more realistically, a uniform fixed-constant guarantee on a specific restricted class. We do *not* prove such a bound here. What we claim is that Hvala is the first simple linear-time algorithm for Minimum Vertex Cover whose combined theoretical and empirical profile — rigorous $\le 2$ bound, pointwise strict $< 2$, linear time, and observed ratios uniformly below $1.414$ across $239$ diverse instances — makes the question above a plausible and well-posed target for future work.

**Algorithm Availability:**
- **Package:** https://pypi.org/project/hvala
- **Installation:** `pip install hvala`
- **Usage:** `from hvala.algorithm import find_vertex_cover`

---

## Acknowledgment

The author is sincerely grateful to Iris, Marilin, Sonia, Yoselin, Arelis, Anissa, Liuva, Yudit, Gretel, Gema, and Blaquier, as well as Israel, Arderi, Juan Carlos, Yamil, Alejandro, Aroldo, Yary, Reinaldo, Alex, Emmanuel, and Michael for their constant support. Whether through encouragement, stimulating conversations, practical assistance, or simply being present during challenging moments, their contributions have played an important role in bringing this work to completion.

---

## References

**[karp2009reducibility]** Karp, Richard M. (2010). *Reducibility Among Combinatorial Problems.* In 50 Years of Integer Programming 1958–2008 (pp. 219–241). Springer, Berlin, Germany. DOI: [10.1007/978-3-540-68279-0_8](https://doi.org/10.1007/978-3-540-68279-0_8).

**[papadimitriou1998combinatorial]** Papadimitriou, Christos H. and Steiglitz, Kenneth (1998). *Combinatorial Optimization: Algorithms and Complexity.* Courier Corporation, North Chelmsford (MA).

**[karakostas2009better]** Karakostas, George (2009). *A Better Approximation Ratio for the Vertex Cover Problem.* ACM Transactions on Algorithms, 5(4), 1–8. DOI: [10.1145/1597036.1597045](https://doi.org/10.1145/1597036.1597045).

**[karpinski1996approximating]** Karpinski, Marek and Zelikovsky, Alexander (1996). *Approximating Dense Cases of Covering Problems.* DIMACS Series in Discrete Mathematics and Theoretical Computer Science, 26, 147–164. Providence, Rhode Island.

**[dinur2005hardness]** Dinur, Irit and Safra, Samuel (2005). *On the Hardness of Approximating Minimum Vertex Cover.* Annals of Mathematics, 162, 439–485. DOI: [10.4007/annals.2005.162.439](https://doi.org/10.4007/annals.2005.162.439).

**[khot2017independent]** Khot, Subhash and Minzer, Dor and Safra, Muli (2017). *On Independent Sets, 2-to-2 Games, and Grassmann Graphs.* Proceedings of the 49th Annual ACM SIGACT Symposium on Theory of Computing, 576–589. Montreal, Canada. DOI: [10.1145/3055399.3055432](https://doi.org/10.1145/3055399.3055432).

**[dinur2018towards]** Dinur, Irit and Khot, Subhash and Kindler, Guy and Minzer, Dor and Safra, Muli (2018). *Towards a Proof of the 2-to-1 Games Conjecture?* Proceedings of the 50th Annual ACM SIGACT Symposium on Theory of Computing, 376–389. Los Angeles, California. DOI: [10.1145/3188745.3188804](https://doi.org/10.1145/3188745.3188804).

**[khot2018pseudorandom]** Khot, Subhash and Minzer, Dor and Safra, Muli (2018). *Pseudorandom Sets in Grassmann Graph Have Near-Perfect Expansion.* 2018 IEEE 59th Annual Symposium on Foundations of Computer Science, 592–601. Paris, France. DOI: [10.1109/FOCS.2018.00062](https://doi.org/10.1109/FOCS.2018.00062).

**[khot2008vertex]** Khot, Subhash and Regev, Oded (2008). *Vertex Cover Might Be Hard to Approximate to Within* $2-\epsilon$*.* Journal of Computer and System Sciences, 74(3), 335–349. DOI: [10.1016/j.jcss.2007.06.019](https://doi.org/10.1016/j.jcss.2007.06.019).

**[khot2002unique]** Khot, Subhash (2002). *On the Power of Unique 2-Prover 1-Round Games.* Proceedings of the 34th Annual ACM Symposium on Theory of Computing, 767–775. Montreal, Canada. DOI: [10.1145/509907.510017](https://doi.org/10.1145/509907.510017).

**[cai2017finding]** Cai, Shaowei and Lin, Jinkun and Luo, Chuan (2017). *Finding a Small Vertex Cover in Massive Sparse Graphs.* Journal of Artificial Intelligence Research, 59, 463–494. DOI: [10.1613/jair.5443](https://doi.org/10.1613/jair.5443).

**[zhang2023tivc]** Zhang, Yu and Wang, Shengzhi and Liu, Chanjuan and Zhu, Enqiang (2023). *TIVC: An Efficient Local Search Algorithm for Minimum Vertex Cover in Large Graphs.* Sensors, 23(18), 7831. DOI: [10.3390/s23187831](https://doi.org/10.3390/s23187831).

**[harris2024faster]** Harris, David G. and Narayanaswamy, N. S. (2024). *A Faster Algorithm for Vertex Cover Parameterized by Solution Size.* 41st International Symposium on Theoretical Aspects of Computer Science, 40:1–40:18. Clermont-Ferrand, France. DOI: [10.4230/LIPIcs.STACS.2024.40](https://doi.org/10.4230/LIPIcs.STACS.2024.40).

**[bar1985local]** Bar-Yehuda, R. and Even, S. (1985). *A Local-Ratio Theorem for Approximating the Weighted Vertex Cover Problem.* Annals of Discrete Mathematics, 25, 27–46.

**[luo2019local]** Luo, Chuan and Hoos, Holger H. and Cai, Shaowei and Lin, Qingwei and Zhang, Hongyu and Zhang, Dongmei (2019). *Local Search with Efficient Automatic Configuration for Minimum Vertex Cover.* Proceedings of the 28th International Joint Conference on Artificial Intelligence, 1297–1304. Macao, China.

**[banharnsakun2023new]** Banharnsakun, Anan (2023). *A New Approach for Solving the Minimum Vertex Cover Problem Using Artificial Bee Colony Algorithm.* Decision Analytics Journal, 6, 100175. DOI: [10.1016/j.dajour.2023.100175](https://doi.org/10.1016/j.dajour.2023.100175).

**[Vega26Hallelujah]** Vega, Frank (2026). *An Approximate Solution to the Minimum Vertex Cover Problem: The Hallelujah Algorithm.* International Journal of Parallel, Emergent and Distributed Systems. Taylor & Francis. DOI: [10.1080/17445760.2026.2660724](https://doi.org/10.1080/17445760.2026.2660724). Accepted for publication.

**[NPBench]** Nguyen, ThanhVu and Bui, Thang. *NP-Complete Benchmark Instances.* Available at: [https://roars.dev/npbench/](https://roars.dev/npbench/). Vertex cover benchmark collection; FRB instances (Ke Xu), DIMACS clique complements, random graphs (Periannan).

**[DIMACSClique]** Mascia, Franco. *The Maximum Clique Problem — DIMACS Benchmark Set.* Available at: [https://iridia.ulb.ac.be/~fmascia/maximum_clique/DIMACS-benchmark](https://iridia.ulb.ac.be/~fmascia/maximum_clique/DIMACS-benchmark). Compiled clique-number values for DIMACS Second Implementation Challenge instances.

**[RA15]** Rossi, Ryan and Ahmed, Nesreen (2015). *The Network Data Repository with Interactive Graph Analytics and Visualization.* Proceedings of the AAAI Conference on Artificial Intelligence, 29(1). DOI: [10.1609/aaai.v29i1.9277](https://doi.org/10.1609/aaai.v29i1.9277).

**[Vega25]** Vega, Frank (2026). *Hvala: Approximate Vertex Cover Solver.* Available at: [https://pypi.org/project/hvala](https://pypi.org/project/hvala). Accessed: 2026-04-22.

**[Frank25]** Vega, Frank (2026). *The Milagro Experiment.* Available at: [https://github.com/frankvegadelgado/milagro](https://github.com/frankvegadelgado/milagro). Accessed: 2026-04-20.

---

**MSC (2020):** 05C69 (Covering and packing), 68Q25 (Analysis of algorithms and problem complexity), 90C27 (Combinatorial optimization), 68W25 (Approximation algorithms)

---

**Documentation**  
Available as PDF at *[An Approximate Solution to the Minimum Vertex Cover Problem: The Hvala Algorithm](https://www.preprints.org/manuscript/202506.0875/v12)*.
