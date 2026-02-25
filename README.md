# 🧬 PhD Thesis Repository

## **Transcriptional dynamics of the flower–fruit transition in *Vanilla planifolia* Andrews (Orchidaceae):

A molecular phenology and co-expression network approach**

---

## 📌 Overview

This repository contains the analytical workflow developed for the doctoral thesis focused on understanding the **molecular regulation of the flower-to-fruit transition** in *Vanilla planifolia*.

Each folder corresponds to a **step of the research pipeline**.

---

## 🎯 Objective

To describe the transcriptional dynamics associated with the molecular phenology of the flower-to-fruit transition in *V. planifolia* by integrating:

* Differential gene expression
* Co-expression network analysis (WGCNA)
* Comparative orthology

to identify candidate regulatory genes.

---

## 🧭 Workflow

### **1. Molecular phenology (RNA-seq & DEG)**

* Ovary transcriptomes from four stages: Pre-Pol, Pol, Post-Pol, Fer
* De novo assembly, annotation, DEG analysis
* Functional enrichment and stage-specific signatures
* Labeling of TFs and epigenetic regulators

---

### **2. Comparative orthology**

* Integration of transcriptomes/genomes from:
  *Arabidopsis*, *Solanum*, *Vitis*, and *Vanilla*
* Identification of conserved orthogroups
* Creation of the **C-Orth** label

---

### **3. Co-expression network (WGCNA)**

* Network construction from normalized counts
* Module detection and functional enrichment
* Identification of hub genes and intramodular cores

---

### **4. Integration & prioritization**

* Combination of structural (DEG, Hub) and functional labels
* Intersection analysis (UpSet plots)
* Prioritization of candidate regulatory genes

---

## 🗂 Repository Structure

Each directory represents one step:

```
01_Molecular_Phenology/
02_Comparative_Orthology/
03_WGCNA_Network/
04_Integration_Prioritization/
```

---

## 👩‍🔬 Author
Olga Andrea Hernández Miranda
PhD candidate in Biologicas Sciences in Experimental Biology

---
<img width="1498" height="524" alt="image" src="https://github.com/user-attachments/assets/be3a1f62-8439-40fe-81ad-e71b245ed461" />
