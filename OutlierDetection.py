import HomologyAnalysis;

positives = ['Scutellaria-baicalensis.fasta']
negatives = [
    'Oryza-sativa_rna.fna',
    'Triticum-urartu.fsa_nt',
    'Arabidopsis-thaliana_rna.fna',
    'Citrus-clementina.fsa_nt',
    'Phoenix-dactylifera.fsa_nt',
    'Vitis-vinifera.fsa_nt'
    ]

speciesDictionary = {}
for sequence_file in positives + negatives:
    hmm_location = "Methyltransf_2.hmm"
    hammer_location = "hmmsearch"
    hammer_result_file = "search_results" + sequence_file + ".txt"
    analysis = HomologyAnalysis.HomologyAnalysis()
    analysis.load_sequence_into_memory(sequence_file)
    analysis.convert_to_protein()
    analysis.hammer_analysis(hmm_location, hammer_location, hammer_result_file)
    resultDictionary = analysis.produce_final_fasta(hammer_result_file)
    speciesDictionary[sequence_file] = resultDictionary


print(speciesDictionary)
