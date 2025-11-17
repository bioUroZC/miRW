# **miRW: A Multi-Omics Random Walk Framework for Independent Construction of Sample-Specific Protein–Protein Interaction Networks in Cancer**

This repository provides a complete implementation of the *miRW* framework, designed to compute sample-specific protein-protein interaction networks using Random Walk with Restart (RWR).
The workflow supports two modes:

* **S1:** Seed setting S₁ (prior knowledge + expression): proteins whose expression levels exceeded the 90th percentile in the sample or belonged to a predefined set of 1,972 proteins associated with 11 cancer-related cellular states (invasion, metastasis, proliferation, angiogenesis, apoptosis, cell cycle, differentiation, DNA damage, DNA repair, hypoxia, and inflammation)
* 
* **S2:** Seed setting S₂ (expression only): proteins whose expression exceeded the 90th percentile in the sample.

Both produce importance measures for each sample and network link.

## 🔧 **1. Installation**

Clone the repository:

```bash
git clone https://github.com/bioUroZC/miRW.git
cd miRW
```

Install required Python packages:

```bash
pip install -r requirements.txt
```

Required packages include:

```
numpy
pandas
networkx
tqdm
scipy
```

---

## 📁 **2. Input Data Files**

The pipeline requires three input files:

### **2.1 CellMarkers.csv**

A list of cell-type–specific seed nodes.

| Column   | Description           |
| -------- | --------------------- |
| Marker   | Gene/miRNA identifier |
| CellType | Cell type category    |

---

### **2.2 Links.csv**

Network edges.

| Column | Description             |
| ------ | ----------------------- |
| Node1  | Endpoint of interaction |
| Node2  | Endpoint of interaction |

---

### **2.3 exprSet.csv**

Gene or miRNA expression matrix.

* **Rows:** genes/miRNAs
* **Columns:** samples
* **Values:** expression levels

Example:

```
          Sample1   Sample2   Sample3 ...
GeneA       3.2       2.1       1.5
GeneB       4.8       7.2       6.1
...
```

---

## ▶️ **3. Running the Pipeline**

The complete pipeline can be launched using:

```bash
python run_miRW.py \
  --seed data_example/CellMarkers.csv \
  --links data_example/Links.csv \
  --expr data_example/exprSet.csv \
  --outdir results/
```

### Arguments

| Argument   | Description                              |
| ---------- | ---------------------------------------- |
| `--seed`   | Path to CellMarkers.csv                  |
| `--links`  | Path to Links.csv                        |
| `--expr`   | Path to exprSet.csv                      |
| `--outdir` | Output directory (created automatically) |

---

## 📊 **4. Output Files**

The program generates four output files:

### **4.1 miRWImpP.csv**

Importance scores (*Imp*) from **S1 (seed-based RWR)**.

### **4.2 miRWFlowP.csv**

Flow scores from **S1**.

---

### **4.3 miRWImpE.csv**

Importance scores from **S2 (seed-free RWR)**.

### **4.4 miRWFlowE.csv**

Flow scores from **S2**.

Each output contains:

| Column    | Description                        |
| --------- | ---------------------------------- |
| Sample    | Sample identifier                  |
| link      | Combined interaction (Node1_Node2) |
| miRW-Imp  | Importance score                   |
| miRW-Flow | Flow score                         |

Example:

```
Sample,link,miRW-Imp
TCGA_01A,ARF5_RAB11FIP3,2.18618
TCGA_01A,RAB11A_RAB11FIP3,2.51365
...
```

---

## 📐 **5. Code Structure**

```
miRW-RWR/
│
├── miRWfunc1.py      # Functions for S1 (seed-based RWR)
├── miRWfunc2.py      # Functions for S2 (seed-free RWR)
├── run_miRW.py       # Main pipeline script
│
├── data_example/     # Example input data
│
├── results/          # Outputs (created after running)
├── requirements.txt
└── README.md
```

---

## 🧪 **6. Example Workflow**

To reproduce an example run:

```bash
python run_miRW.py \
  --seed data_example/CellMarkers.csv \
  --links data_example/Links.csv \
  --expr data_example/exprSet.csv \
  --outdir demo_output/
```

---

## 📚 **7. Citation**

If you use this tool in your research, please cite:

```

```

---

## 📄 **8. License**

This project is released under the **MIT License**, allowing free academic and commercial use.

---

## ❓ **9. Questions / Issues**

If you encounter problems or have questions, feel free to open an issue on GitHub.

