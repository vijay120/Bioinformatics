import os


class HomologyAnalysis:
	CODON_MAP = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
	             "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
	             "TAT": "Y", "TAC": "Y", "TAA": "STOP", "TAG": "STOP",
	             "TGT": "C", "TGC": "C", "TGA": "STOP", "TGG": "W",
	             "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
	             "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
	             "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
	             "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
	             "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
	             "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
	             "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
	             "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
	             "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
	             "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
	             "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
	             "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
	             }

	def __init__(self):
		self.mrna_transcripts = {}
		self.protein_transcripts = {}
		self.temp_fasta_location = "protein_sequences.fasta"

	def load_sequence_into_memory(self, seq_file):
		"""
		Reads in a fna file and parses it for mRNA sequences, key is the header file.
		:param sequence_file: FNA file that is loded into memory as a bunch of mRNA sequences.
		:return:
		"""
		with open(seq_file, 'r') as mrna_file:
			current_key = None
			for line in mrna_file:
				if line[0] == ">":
					current_key = line
					self.mrna_transcripts[current_key] = ""
				else:
					self.mrna_transcripts[current_key] += line.strip()

	def convert_to_protein(self):
		"""
		Converts the mRNA sequences to the expected protein frames.
		:return:
		"""
		with open("protein_sequences.fasta", 'w') as protein_file:
			for key in self.mrna_transcripts.keys():
				protein = self._convert_mrna_to_protein(self.mrna_transcripts[key])
				if protein is not None:
					self.protein_transcripts[key.split(" ")[0][1:len(key)]] = protein
					protein_file.write(key.split(" ")[0] + "\n")
					protein_file.write(protein + "\n")
					protein_file.write("\n")
				else:
					pass

	def hammer_analysis(self, hmm_file, hammer_search_binary, result_file):
		"""
		Uses the hammer package to create a file.

		:param hmm_file: Where the hmm file that has been made or downloded is found.
		:param hammer_search_binary: Where the hammer binary is found hmmsearch
		:param result_file: Output of hammer
		:return:
		"""
		command = hammer_search_binary + " " + hmm_file + " " + self.temp_fasta_location + " > " + result_file
		os.system(command)

	def produce_final_fasta(self, hammer_result_file, final_fasta_file):
		"""
		Parses hammer output and puts the output into a FASTA file.

		:param hammer_result_file: The result file of hammer
		:param final_fasta_file: The output fasta file
		:return:
		"""
		matching_sequences = []
		with open(hammer_result_file, "r") as results:
			collect = False
			skip_one = False
			for line in results:

				""" End of interesting stuff """
				if "inclusion threshold" in line:
					continue

				if collect and skip_one:
					write_line = line.strip()
					if len(write_line) > 0:
						if "NM_" in line:
							sequence_name = "NM_" + str(line.split("NM_")[1].split(" ")[0])
							matching_sequences.append(sequence_name)
						elif "gi|" in line:
							sequence_name = "gi|" + str(line.split("gi|")[1].split(" ")[0])
							matching_sequences.append(sequence_name)
						elif "TR" in line:
							sequence_name = "TR" + str(line.split("TR")[1].split(" ")[0])
							matching_sequences.append(sequence_name)
					else:
						break
				else:
					skip_one = True

				""" Beginning of matched sequences """
				if "E-value" in line:
					collect = True

		with open(final_fasta_file, "w") as result:
			for m in matching_sequences:
				result.write(">" + m + "\n")

				# Add stop codon.
				try:
					self.protein_transcripts[m] += "*"
				except KeyError:
					print "Error"
					print m
					print self.protein_transcripts
					exit(1)

				# Each line has 80 characters
				lines = [self.protein_transcripts[m][i:i + 80] for i in range(0, len(self.protein_transcripts[m]), 80)]

				for line in lines:
					result.write(line + "\n")
				result.write("\n")

	def _get_sequence_frame(self, sequence, frame_index=0):
		"""
		Gets the starting ATG file starting at a position in the frame.

		:param sequence: mRNA sequence
		:param frame_index: Where to start looking
		:return: The index of the A in the start codon
		"""
		for index in range(frame_index, len(sequence) - 2):
			"""
			Looks for a start codon.  Short circuits so if an A isn't found it just keeps going.
			"""
			if sequence[index] == "A" and sequence[index + 1] == "T" and sequence[index + 2] == "G":
				return index
		return None

	def _convert_mrna_to_protein(self, sequence):
		"""
		Looks for protein sequences larger than 100 that are framed.

		:param sequence: mRNA sequence
		:return: Protein sequence
		"""
		frame_index = self._get_sequence_frame(sequence)

		while frame_index is not None:
			current_sequence = self._convert_frame_to_protein(frame_index, sequence)
			"""
			Filter out smaller subsequences here
			"""
			if current_sequence is None:
				return None

			if len(current_sequence) > 100:
				return current_sequence
			else:
				frame_index = self._get_sequence_frame(sequence, frame_index + 1)

		return None

	def _convert_frame_to_protein(self, frame_index, sequence):
		"""
		Takes a starting of a frame and makes a protein until it hits the STOP codon.
		:param frame_index: Where to start in the sequence
		:param sequence: mRNA sequence
		:return: Created protein sequence
		"""
		protein_sequence = ""
		for frame_pos in range(frame_index, len(sequence) - 2, 3):
			codon = sequence[frame_pos] + sequence[frame_pos + 1] + sequence[frame_pos + 2]
			if HomologyAnalysis.CODON_MAP.get(codon) is None:
				""" FNA supports variable sequences, so we will too """
				protein_sequence += "-"

			elif HomologyAnalysis.CODON_MAP[codon] == "STOP":
				return protein_sequence
			else:
				protein_sequence += HomologyAnalysis.CODON_MAP[codon]


if __name__ == "__main__":
	sequence_file = "Arabidopsis-thaliana_rna.fna"
	hmm_location = "Methyltransf_2.hmm"
	hammer_location = "hmmsearch"
	hammer_result_file = "search_results.txt"
	output_fasta = "output.fasta"

	analysis = HomologyAnalysis()
	"""
	Load mRNA sequence into file
	"""
	analysis.load_sequence_into_memory(sequence_file)

	"""
	Convert to protein
	"""
	analysis.convert_to_protein()

	"""
	User hammer to get hmm matches
	"""
	analysis.hammer_analysis(hmm_location, hammer_location, hammer_result_file)

	"""
	Parse hammer output and output resulting file.
	"""
	analysis.produce_final_fasta(hammer_result_file, output_fasta)
