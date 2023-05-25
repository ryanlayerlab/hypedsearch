Instructions:

    1. Run comet using files in the "input_spectra" folder and the "UniProt_mouse.fasta" file found in "hypedsearch/data/database" - done :)
    2. Name the output files to (filename).txt - done :)
    3. Put output files in the raw_comet_results folder - done :)
    4. Check the settings in BMEM_config.yaml - done :)
    5. Run "building_hypedsearch_database.py" - done :)
    6. Run "BMEM_run.py"
    7. Run "building_hybrid_database.py"
    8. The hybrid database will be found in "hypedsearch_hybrid_database". Take this and run Comet using files in the "input_spectra" folder and the new hybrid database
    9. Name the output_files to hybrid_(filename).txt
    10. Put these output files in hybrid_comet_results
    11. Run "comparing_hypedseach_and_comet.py"
    12. Results will be in the "hypedsearch_comet_comparison" folder