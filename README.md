# 📘 Final Year Project — Obtaining Coherent Configurations on Non-Distance Regular Graphs

This project explores the structure of **coherent configurations** derived from graph operations applied to known strongly regular graphs. It focuses on identifying algebraic closures through matrix operations and 2-dimensional Weisfeiler–Leman refinement, with applications in combinatorial design and spectral graph theory.

---

## 📌 Project Title
**Obtaining Coherent Configurations on Non-Distance Regular Graphs**

---

## 🎓 Academic Info
- **Institution**: Nanyang Technological University (NTU), Singapore
- **Degree**: BSc. Mathematical and Computer Sciences (Honours with Highest Distinction)
- **Grade**: A
- **Tools Used**: SageMath, Python, NumPy, C++ (for WL-algorithm), LaTeX

---

## 🧠 Problem Statement

Given a strongly regular graph (SRG), can we obtain a minimal coherent configuration by applying:
- Vertex deletion
- Seidel switching
- Edge switching (especially on block graphs from orthogonal arrays)?

The goal was to observe **how far such graphs are from being association schemes**, and whether their **coherent closure** stabilizes into identifiable algebraic patterns.

---

## 🔍 Key Concepts

- **Coherent Configuration**: A partition of the edge set (including self-loops) satisfying closure under composition and adjacency.
- **Weisfeiler-Leman Algorithm (2-WL)**: A refinement method for graph colorings that helps detect non-isomorphism and builds matrix closures.
- **Wielandt’s Principle**: A tool to confirm minimality of the generated algebra from adjacency matrices.

---

## 🔬 Methods

1. Constructed initial graphs from known SRGs:
   - Rook graphs $R(n,k)$
   - Triangular graphs $T(n)$

2. Applied operations:
   - Vertex deletion
   - Seidel switching
   - Block edge-switching using OA(2, n)

3. Computed adjacency matrices and applied closure under multiplication using SageMath.

4. Verified coherent rank using:
   - Manual combinatorial tracking
   - Implementation of WL-refinement algorithm (k-dim) in C++

5. Validated closure with explicit applications of **Wielandt’s Principle**.

---

## 🧩 Key Observations

- For block graphs from **orthogonal arrays OA(2, n)**:
  - Edge switching consistently yielded a **coherent rank of 15**
  - Type matrices had a **stable block diagonal structure**:
    - Diagonal 1: \( n-1 \) blocks of type \( I \) and \( J - I \)
    - Diagonal 2: fibres with \( J - I_{n-1} \) and 0

- Switching on rook graphs disrupted symmetry, but algebraic closure retained regularity under WL refinement.

---

## 📊 Deliverables

- Full SageMath scripts for closure computation
- Graph visualizations and adjacency matrix logs
- PDF write-up (LaTeX) including matrix patterns and theory
- Experimental results on type matrix decomposition

---

## 📘 Learnings

- Deepened understanding of combinatorial algebra and its connection to matrix representation theory
- Gained practical experience using SageMath for algebraic graph theory
- Bridged graph operations with coherent configuration theory

---

## 📁 Repository Status

Currently unpublished. For access or collaboration, please contact:

📧 yujiho54@gmail.com  
🔗 [LinkedIn](https://linkedin.com/in/jingrui19)

---
