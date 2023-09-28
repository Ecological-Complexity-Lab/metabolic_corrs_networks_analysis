# metabolic_corrs_networks_analysis
sharon's project- creating networks of correlations between metabolites in different bacteria strains (before/after evolution) and analyzing them to catch the effect of evolution 

Workflow:
in the file "LCMS_output_analysis"
1. LCMS output comes in the shape of chromatograms, in raw files (located in a hard-drive (number 2) in Shimon's lab). Used compound discoverer to convert the data into a csv file called "pos_compound_results_no_filter.xlsx". 5652 metabolic features were found. The workfloe of CD includes filtration, by amounts in the blank samples (filtered if the ratio os above above 2.5) and intensity (only above 1e6). after filtration 2801 metabolic features remained in the file "Compounds_after_filtrations_pos.xlsx" in the "compounds" sheet.  
2. consalidation of annotated metabolites- manually grouped rows of metabolites that were annotated by CD to un-annotated metabolites with similar m/z, RT, formula. if the metabolite fits- checked the peak with Xcalibur to make sure it's the same metabolite. consalidated 481 into 180. 
3. consalidation of un-annotated metabolites- used martin's code- by RT, changed it so the grouping would happen only if the m/z diiference fits to a known adduct. 2565 metabolic featueres (both annotated and un-annotated) into 2118. 
4. finished with a file called "unique_19_each_585_new.xlsx"

5. 
