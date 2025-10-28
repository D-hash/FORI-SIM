# FORI-SIM: An ILP-Based Algorithm for Graph Similarity Search

This repository contains the code for **FORI-SIM**, the algorithm presented in  
**"Accelerating Graph Similarity Search through Integer Linear Programming"**, submitted to the IEEE ICDE 2026 conference.

---

## 🛠️ Installation Requirements

### Required Software

- **GUROBI**
- **GEDLIB**  
  https://github.com/dbblumenthal/gedlib
- **LIBLSAP** *(header-only library)*  
  https://forge.greyc.fr/projects/liblsap/  
  Add the following to your CMake options:  
LIBLSAP_ROOT=~/liblsap

 ⚠️ **Note:** You must modify the following line in `liblsap/cpp/include/lsap.hh`:  
Change  
> ```cpp
> static const DataType _zero = 0;
> ```  
to  
> ```cpp
> static constexpr DataType _zero = 0;
> ```

---

### 🔧 Configuration Variables

#### CMake Options

-DGUROBI_DIR=/path/to/gurobi\<version>\/\<OSdependent\>/

#### Environment Variables
LIBLSAP_ROOT=/path/to/liblsap GEDLIB_ROOT=/path/to/gedlib GUROBI_HOME="/path/to/gurobi\<version\>/\<OSdependent\>/"


---

## 📊 Replicating Our Results

Our results can be replicated by running the scripts *verification_uniform.sh* and *verification_non_uniform.sh* for uniform and non-uniform edit cost cases, respectively.

---

## 🔍 RelateExternald Implementations

### FORI-LP and FORI

To compute lower bounds or optimal GED solutions, we use the implementation from:  
https://github.com/meffertj/FORI-GED

### BM Lower Bound

We use the BM lower bound implementation from GEDLIB:  
https://github.com/dbblumenthal/gedlib

### LS and A*-BMAO

We adapted the LS implementation and used the original A*-BMAO code from:  
https://github.com/LijunChang/Graph_Edit_Distance
