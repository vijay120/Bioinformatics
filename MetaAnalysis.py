from tqdm import tqdm
from HomologyAnalysis import HomologyAnalysis
import os
from Bio.Emboss.Applications import NeedleCommandline

"""
Run my initial analysis
"""
def analyze(sequence_file, hmm_location, hammer_location, hammer_result_file, output_fasta):
	analysis = HomologyAnalysis()
	analysis.load_sequence_into_memory(sequence_file)
	analysis.convert_to_protein()
	analysis.hammer_analysis(hmm_location, hammer_location, hammer_result_file)
	analysis.produce_final_fasta(hammer_result_file, output_fasta)


MAKE_FASTAS = True
ALIGN_CONTROLS = True
PAIRWISE = True
GRAPH = True

if __name__ == "__main__":
	hmm_location = "Methyltransf_2.hmm"
	hmmr_search = "hmmsearch"
	hmmr_align = "hmmalign"
	hmmr_build = "hmmbuild"
	hammer_result_file = "search_results.txt"
	output_fasta = "output.fasta"

	"""
	Create fasta file for each containign the proteins close to our HMM file.
	"""
	if MAKE_FASTAS:
		fna_files = [
			"Scutellaria-baicalensis.fasta",
			"Arabidopsis-thaliana_rna.fna",
			"Citrus-clementina.fsa_nt",
			"Oryza-sativa_rna.fna",
			"Phoenix-dactylifera.fsa_nt",
			"Triticum-urartu.fsa_nt",
			"Vitis-vinifera.fsa_nt"
		]
		for fna in tqdm(fna_files):
			analyze(fna,
			        hmm_location,
			        hmmr_search,
			        fna.split(".")[0] + hammer_result_file,
			        fna.split(".")[0] + output_fasta)

	"""
	Put negative control files into one pool, align with each other and generate a new HMM based on that alignment.
	"""
	if ALIGN_CONTROLS:
		combine_negative_controls = ["Arabidopsis-thaliana_rnaoutput.fasta",
		                             "Citrus-clementinaoutput.fasta",
		                             "Oryza-sativa_rnaoutput.fasta",
		                             "Phoenix-dactyliferaoutput.fasta",
		                             "Triticum-urartuoutput.fasta",
		                             "Vitis-viniferaoutput.fasta"]


		"""
		Create one negative controls file
		"""
		with open("Negative_Controls.fasta", "w") as nc:
			for f in combine_negative_controls:
				with open(f, "r") as cf:
					for line in cf:
						if ">" in line:
							nc.write(line.replace("|", "_"))
						else:
							nc.write(line)
					nc.write("\n")

		"""
		Use hmmer to align the negative control file to the HMM of the methyltransferase
		"""
		command = hmmr_align \
		          + " -o negative_controls.align " \
		          + hmm_location \
		          + " Negative_Controls.fasta"
		os.system(command)
	#
	# 	"""
	# 	Build a new HMM based on the aligned sequences
	# 	"""
	# 	command = hmmr_build \
	# 	          + " out.hmm " \
	# 	          + "negative_controls.align"
	# 	os.system(command)
	#
	# 	"""
	# 	Use that new HMM to search the outlier species
	# 	"""
	# 	command = hmmr_search + " " + "out.hmm" + " " + "Scutellaria-baicalensisoutput.fasta" + " > " + "Part2_Result.txt"
	# 	os.system(command)
	#
	# """
	# Analysis of pairwise scoring between negative control file and each protein in the outlier species.
	# """
	# if PAIRWISE:
	# 	with open("Scutellaria-baicalensisoutput.fasta", "r") as outlier:
	# 		sequence = True
	# 		sequence_name = None
	# 		current_file = None
	# 		for line in outlier:
	# 			if ">" in line:
	# 				if current_file is not None:
	# 					current_file.close()
	# 					needle_cline = NeedleCommandline(asequence="pairwise_results/" + sequence_name + ".fasta",
	# 			                                 bsequence="Negative_Controls.fasta",
	# 			                                 gapopen=10,
	# 			                                 gapextend=0.5,
	# 			                                 outfile= "pairwise_results/" + sequence_name + ".txt")
	# 					needle_cline()
	#
	# 				sequence_name = line.replace("|", "_")[1:len(line)].strip()
	# 				current_file = open("pairwise_results/" + sequence_name + ".fasta", "w")
	# 				current_file.write(line.replace("|", "_"))
	#
	# 			elif sequence_name is not None:
	# 				current_file.write(line)
	#
	# """
	# Initially was going to graph if needed, but the sorting of the averages was telling enough.
	# """
	# if GRAPH:
	# 	pairwise_results = [f for f in os.listdir("pairwise_results") if f.endswith(".txt")]
	# 	averages = []
	# 	for result in pairwise_results:
	# 		scores = []
	# 		with open("pairwise_results/" + result, 'r') as r:
	# 			for line in r:
	# 				if "Score" in line:
	# 					scores.append(float(line.split("# Score: ")[1]))
	# 		scores.sort()
	# 		averages.append((reduce(lambda x, y: x + y, scores) / len(scores), result))
	# 	averages.sort()
	# 	print averages
	#
	#
	#
