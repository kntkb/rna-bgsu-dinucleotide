# Create dinucleotide RNA dataset for QCArchive submission (RNA-dinucleotide-single-points)


Download json/csv files from [RNA BGSU](https://www.bgsu.edu/research/rna.html)
------
[Internal loop motifs](http://rna.bgsu.edu/rna3dhub/motifs/release/il/3.61)
- InternalLoopMotifAtlasRelease3.61.json

[Hairpin loop motifs](http://rna.bgsu.edu/rna3dhub/motifs/release/hl/3.61)
- HairpinLoopMotifAtlasRelease3.61.json

[Representative RNA structures](http://rna.bgsu.edu/rna3dhub/nrlist/release/3.245)
- nrlist_3.245_4.0A.csv



Extract structure data using [RNA BGSU APIs](https://www.bgsu.edu/research/rna/APIs.html)
------

- `download_motifs.ipynb`  
    - Download internal/hairpin loop motifs
    - Note that these are not properly formatted in CIF format where html/css codes are also included in the file.
    
- `download_junctions_from_nrlist.ipynb`  
    - Download junction loops 
    - Since junction loops are not provided as part of the [RNA 3D Motif Atlas](http://rna.bgsu.edu/rna3dhub/motifs) at the time of analysis, junction loops were extracted from [Representative RNA structures (version 3.245, resolution<4.0)](http://rna.bgsu.edu/rna3dhub/nrlist/release/3.233).



Check downloaded cif/pdb files
------

- `check_motif.ipynb`  
    - Check for missing atoms, inconsistent residue numbers, residue names. Following files failed to pass the test. Filename added with "warning" or "missing_atoms" and ignored for further process. 
        - `HL_44730.2.cif, HL_51090.1.cif, HL_01181.4.cif, HL_70505.1.cif, HL_35188.1.cif, HL_49873.1.cif, HL_19239.1.cif, HL_48810.1.cif, IL_46464.1.cif, IL_51971.1.cif, IL_77341.1.cif, IL_39900.1.cif, IL_16160.1.cif, IL_50112.1.cif, IL_57188.3.cif, IL_89836.1.cif, IL_22427.1.cif, IL_80617.1.cif, IL_64414.1.cif, J3_7RQB_030.cif, J3_7RQ8_031.cif, J3_7RQ8_003.cif, J3_7RQB_003.cif`



Split loops into three consecutive bases and perform clustering for each base set
------

- `split_motif.ipynb`  
    - Hairpins, internal, and junction loops were split into dinucleotides (two connected bases).
    - Atoms `P`, `OP1`, and `OP2` were deleted from 5' base.

- `cluster_motifs.ipynb`
    - Each triple base set was clustered with bottom-up (agglomerative) hierarchical clustering with average linking.
    - Euclidean distances between set of pre-defined atom pairs were used as input features. Also tried torsion angles as input features for clustering but internal distances were better.
    - Distance threhold of "1" was used for clustering. This was defined so that the number of cluster centroids from each triple base set were ~4000 structures. The maximum RMSD measured from the centroid structure for each cluster was also monitored to ensure structures are well-clustered. Maximum RMSD from all clusters are below ~2.4 Angstroms.



Edit 5' base residues and add TER from base pair catalogs and base triple database structures
------

- `edit_base.ipynb`  
    - Atoms `P`, `OP1`, and `OP2` were deleted from 5' base.
    - TER was added to insure bases were disconnected



Add hydrogen atoms and minimize structure
------

- `script/check_files.py`
    - All files were checked if they could be loaded with PDBFixer
    - Warning found for `BP_tWS_AG.pdb` and `BP_cWS_AC.pdb`. TER was inserted in the wrong position. Manually fixed.

- `script/openmm_implicit_minimizer.py`
    - Added hydrogen prior to minimization
    - Implicit solvent applied (GBN2)
    - Heavy atom restraint applied (30 kcal/mol/A^2)
    - 10 step NVT to check nothing fuzzy with the input structure (***removed from latest version***)
        - Errors raised for some Triple base structures downloaded from [RNA Base Triple Database](http://rna.bgsu.edu/triples/triples.php)
        - It turns out that some triple base structures (e.g. Triple_tWH_cSS_CAG.pdb/Triple_tWH_cSS_CAG.cif) have overlapping atoms... These errors were raised for modeled structures. No problems detected for experimental structures (i.e. "exemplar").   



Test creating dataset prior to QCArchive submission
------

- `qca-dataset-submission_TEST/generate-dataset.ipynb`
    - [Isomorphic AssertionError](https://github.com/choderalab/rna_bgsu/issues/1) issue
        - motivated to check all strtuctures and remove invalid structures (`check_structure.py`)
- `qca-dataset-submission_TEST/check_structures.py`
    - check pdb to rdkit conversion, stererocenters, stereobonds (input pdb: structures by exported by `script/openmm_implicit_minimizer.py`)
    - add suffix to pdb if warnings and/or problems were found and move problematic strucutres to a different directory
    >cd minimized
    >mkdir error
    >mv *.pdb.\* error
    - **pdb structures that passed the check and conveted to rdkit mol object are saved as pickle**
    - `rdmols_basepairCATALOG.pkl`, `rdmols_loopMOTIFS.pkl`, `rdmols_triplebase_exemplar.pkl` are used for QCArchive submission
    - `rdmols_triplebase.pkl` was excluded from the dataset because these are modeled structures and invalid structures were found for several cases.
- `qca-dataset-submission_TEST/convert_pdb2sdf.sh`
    - convert pdb structures that passed `qca-dataset-submission_TEST/check_structures.py` to sdf using schrodinger software (pdbconvert: pdb -> mae / canvasConvert: mae -> sdf)
    - compress pdb/sdf for basepairCATALOG, loopMOTIFS, triplebaseDB, triplebaseDB_exemplar, respectively
        - files are stored in `minimized` dircetory
