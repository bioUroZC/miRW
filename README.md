# **miRW: A Multi-Omics Random Walk Framework for Independent Construction of Sample-Specific Interaction Networks**

This repository provides a complete implementation of the **miRW** framework, designed to compute sample-specific edge importance scores from gene or miRNA expression data using **Random Walk with Restart (RWR)**.

The workflow supports two computation modes:

* **S1 — Seed-based RWR (optional)**
  Uses predefined biological marker genes (from `CellMarkers.csv`) as seed nodes.

* **S2 — Seed-free RWR (required)**
  Computes edge importance using expression data only.

If `CellMarkers.csv` is not provided, the pipeline automatically runs **S2 only**.


# 🔧 **1. Installation**

Clone the repository:

```bash
git clone https://github.com/bioUroZC/miRW.git
cd miRW/miRWsingleLayer
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

# 📁 **2. Input Data Files**

Example input files are stored in:

```
miRWsingleLayer/data_example/
```

The pipeline uses up to **three** CSV files.

---

## **2.1 CellMarkers.csv (optional)**

Cell-type–specific seed nodes for mode S1.

| Column   | Description           |
| -------- | --------------------- |
| Marker   | Gene/miRNA identifier |
| CellType | Cell-type category    |

If the file is missing, **S1 is skipped automatically**.

---

## **2.2 Links.csv (required)**

Interaction network edges.

| Column | Description             |
| ------ | ----------------------- |
| Node1  | Endpoint of interaction |
| Node2  | Endpoint of interaction |

---

## **2.3 exprSet.csv (required)**

Expression matrix.

* **Rows:** genes or miRNAs
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

# ▶️ **3. Running the Pipeline**

The main workflow is implemented in:

```
miRWsingleLayer/miRW.py
```

### **Run S1 + S2 (when CellMarkers.csv exists):**

```bash
python miRW.py \
  --seed data_example/CellMarkers.csv \
  --links data_example/Links.csv \
  --expr data_example/exprSet.csv \
  --outdir results/
```

### **Run S2 only (no seed file):**

```bash
python miRW.py \
  --links data_example/Links.csv \
  --expr data_example/exprSet.csv \
  --outdir results/
```

### Arguments

| Argument   | Description                              |
| ---------- | ---------------------------------------- |
| `--seed`   | Path to CellMarkers.csv (optional)       |
| `--links`  | Path to Links.csv (required)             |
| `--expr`   | Path to exprSet.csv (required)           |
| `--outdir` | Output directory (created automatically) |

---

# 📊 **4. Output Files**

Results are saved in the folder specified by `--outdir`.

---

## **4.1 S1 Outputs (only if seed file exists)**

| File            | Description                           |
| --------------- | ------------------------------------- |
| `miRWImpP.csv`  | Importance scores from seed-based RWR |
| `miRWFlowP.csv` | Flow scores from seed-based RWR       |

---

## **4.2 S2 Outputs (always generated)**

| File            | Description                          |
| --------------- | ------------------------------------ |
| `miRWImpE.csv`  | Importance scores from seed-free RWR |
| `miRWFlowE.csv` | Flow scores from seed-free RWR       |

---

### Output Format

| Column    | Description                    |
| --------- | ------------------------------ |
| Sample    | Sample identifier              |
| link      | Sorted node pair (NodeA_NodeB) |
| miRW-Imp  | Importance score               |
| miRW-Flow | Flow score                     |

Example:

```
Sample,link,miRW-Imp
TCGA_01A,ARF5_RAB11FIP3,2.18618
TCGA_01A,RAB11A_RAB11FIP3,2.51365
```

---

# 📐 **5. Code Structure**

```
miRW/
│
├── miRWsingleLayer/
│   ├── data_example/
│   │   ├── CellMarkers.csv      # Optional
│   │   ├── Links.csv            # Required
│   │   └── exprSet.csv          # Required
│   │
│   ├── miRW.py                  # Main pipeline script
│   ├── miRWfunc1.py             # Seed-based RWR functions (S1)
│   ├── miRWfunc2.py             # Seed-free RWR functions (S2)
│   ├── requirements.txt
│   └── ReadMe.txt
│
└── README.md
```

---

# 🧪 **6. Example Workflow**

Example run that performs both S1 and S2:

```bash
python miRW.py \
  --seed data_example/CellMarkers.csv \
  --links data_example/Links.csv \
  --expr data_example/exprSet.csv \
  --outdir demo_output/
```

---

# 📚 **7. Citation**

If you use **miRW** in your research, please cite:

*(Add your publication or manuscript reference here)*

---

# 📄 **8. License**

This project is released under the **MIT License**, allowing both academic and commercial use.

---

# ❓ **9. Questions / Issues**

If you encounter any issues or have feature requests, please open an Issue on GitHub:

👉 [https://github.com/bioUroZC/miRW/issues](https://github.com/bioUroZC/miRW/issues)

---

