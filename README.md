# enzyASM
---------------------

> [!WARNING]
What started as enzyASM has now evolved into a more functional Python package called [QMzyme](https://qmzyme.readthedocs.io/). enzyASM is not currently supported and interested users are encouraged to use QMzyme.

---

A toolkit for enzyme active site modeling (ASM).
See requirements.txt for dependencies and suggested packages

### Current Functionality
- Creates a reduced active site model from a protein PDB code or file, suitable for truncated model generate, such as the quantum chemical cluster approach.
- Couples with AQME (https://github.com/jvalegre/aqme) to automate QM input file generation, QM output parsing and error checking.
- Scripy generation to visualize model generation in PyMOL.

### Ongoing Improvements
This code is in an early stage, and there are many improvements I hope to make. I am also open to suggestions via Pull Requests if there is anything you think would be beneficial to include.
- Semi empirical QM conformational sampling via CREST-xTB (using AQME to automate that process)
- Adding additional truncation scheme methods
- Automating QM residue selection schemes: code will generate input files necessary to perform selection calculations (including but not limited to Fukui Shift Analysis and Charge Shift Analysis) and automate result analysis to generate a model with a validated QM region.
- Broadening applications to QM/MM modeling
- Automating ligand docking and shape complementarity evaluations
- Automating reaction coordinate scanning 
- Automating residue mutatgenesis screening
- Automating ligand candidate screening
- Integration with molecular mechanics software
