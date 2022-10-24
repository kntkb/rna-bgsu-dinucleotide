# Create dinucleotide RNA dataset for QCArchive submission

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
        HL: `HL_44730.2.cif, HL_51090.1.cif, HL_01181.4.cif, HL_70505.1.cif, HL_35188.1.cif, HL_49873.1.cif, HL_19239.1.cif, HL_48810.1.cif`  
        IL: `IL_46464.1.cif, IL_51971.1.cif, IL_16160.1.cif, IL_50053.1.cif, IL_50112.1.cif, IL_57188.3.cif, IL_89836.1.cif, IL_22427.1.cif, IL_80617.1.cif, IL_77034.1.cif`  
        J3: `J3_6ZVK_022.cif, J3_7RQB_030.cif, J3_5VPP_042.cif, J3_2UXC_006.cif, J3_1N33_006.cif, J3_5E81_021.cif, J3_5IBB_015.cif, J3_6CZR_010.cif, J3_1N32_006.cif, J3_6ORD_043.cif, J3_4GCW_002.cif, J3_7OIF_022.cif, J3_5VPO_043.cif, J3_4GCW_001.cif, J3_5OT7_010.cif, J3_6GSK_040.cif, J3_6NTA_015.cif, J3_6GSL_014.cif, J3_4TUE_016.cif, J3_5IB8_015.cif, J3_7U2I_003.cif, J3_4V4J_003.cif, J3_4V5K_012.cif, J3_7ACJ_021.cif, J3_2UXD_006.cif, J3_4V5G_039.cif, J3_5VPP_015.cif, J3_7U2I_030.cif, J3_4V5K_034.cif, J3_7OT5_023.cif, J3_6Q9A_022.cif, J3_1HNW_006.cif, J3_6NTA_008.cif, J3_6GSL_034.cif, J3_5IBB_041.cif, J3_4DR7_006.cif, J3_4V5G_015.cif, J3_6ORD_015.cif, J3_6CFJ_031.cif, J3_7AC7_022.cif, J3_4B3R_006.cif, J3_7NWG_002.cif, J3_5VPO_016.cif, J3_7RQB_003.cif, J3_6O8W_009.cif, J3_5T7V_001.cif, J3_6NTA_042.cif, J3_6GSK_015.cif, J3_6GSL_041.cif, J3_4TUE_043.cif, J3_5IB8_041.cif, J3_6Q97_017.cif`



Split loops into dinucleotides (two consectutive nucleotides) and perform clustering
------

- `split_motif.ipynb`  
    - Hairpins, internal, and junction loops were split into dinucleotides (two connected bases).
    - Atoms `P`, `OP1`, and `OP2` were deleted from 5' base.

- `cluster_motifs.ipynb`
    - Each dinucleotide was clustered with bottom-up (agglomerative) hierarchical clustering with average linking.
    - Euclidean distances between set of pre-defined atom pairs were used as input features.




Add hydrogen atoms and minimize structure
------

- `script/openmm_implicit_minimizer.py`
    - Added hydrogen prior to minimization
    - Implicit solvent applied (GBN2)
    - Heavy atom restraint applied (30 kcal/mol/A^2)
    - all minimized files are stored in `minimized/pdb`

