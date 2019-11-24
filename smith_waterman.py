#!/usr/local/bin/python3
import copy
import os.path
import sys

"""
This program uses the  Smith-Waterman algorithm to compute the optimal local alignment(s) of two sequences.
To execute, simply run the script from the command line and follow the intructions.

"""

class Sequence:
	"""
	A class that allows input sequences to be stored as objects with name, 
	status and sequence strings as parameters
	"""
	def __init__(self, name, status, sequence_str):
		self.name = name
		self.status = status
		self.sequence_str = sequence_str
		

class Alignment:
	"""
	A class that stores temprorary information about the alignment
	while it is being computed and finalizes and stores
	the final alignment string after it is done.
	"""
	def __init__(self, q_str, s_str, g_str, end, score, q_name, s_name, i = None, j = None):
		"""
		Initializes alignment.
		Parameters:
		q_str (str): string for the first sequence
		s_str (str): string for the second sequence
		end (list): list that contain the endpoints of the alignment, 
			i.e. the indices at which the alignment end in the two strings
		score (float): the raw alignment score
		i (int): current index of sequence 1
		j (int): current index of sequence 2
		q_name (str): name of sequence 1
		s_name (str): name of sequence 2 
		"""
		self.q_str = q_str
		self.s_str = s_str
		self.g_str = g_str
		self.q_name = q_name
		self.s_name = s_name
		self.end = end
		self.score = score
		self.done = False
		if i == None and j == None:
			self.i = end[0]
			self.j = end[1]
		else:
			self.i = i
			self.j = j
	
	def set_starting_points(self, i, j):
		"""
		Sets the starting points, i.e. the indices at which the alignment starts. To be used when the alignment is done.
		Parameters:
		i (int): index of sequence 1
		j (int): index of sequence 2
		"""
		self.q_start = i + 1
		self.s_start = j + 1
		
	def finalize(self):
		"""
		Finalizes and formats the string representing the final alignment 
		and a header providing information about the alignment.
		No parameters.
		"""
		final = ''
		i = self.q_start
		j = self.s_start
		k = 0
		
		while k < len(self.q_str):
			final += 'Query:' + self.findspace(i) + str(i) + ' ' + self.q_str[k:k+60] + ' ' + str(min([self.end[0], i + 59-self.q_str[k:k+60].count('-')])) + '\n           ' + self.g_str[k:k+60] + '\n' + 'Sbjct:' + self.findspace(j) + str(j) + ' ' + self.s_str[k:k+60] + ' ' + str(min([self.end[1], j + 59 -self.s_str[k:k+60].count('-')])) + '\n\n'
			i = i + 60-self.q_str[k:k+60].count('-')
			j = j + 60-self.s_str[k:k+60].count('-')
			k = k + 60
			
		self.final = final
		self.header = '******* Query: ' + self.q_name + ' *******\n******* Subject: ' + self.s_name + ' *******\n******* ALIGNMENT NUMBER: ' + str(self.nr) + ', SCORE: ' + str(self.score) + ' *******\n'
	
	def findspace(self, n):
		"""
		Returns the correct number of spaces to make the formatted alignment appear correctly.
		"""
		if n < 10:
			return '   '
		elif n < 100:
			return '  '
		elif n < 1000:
			return ' '
		else:
			return ''
				 
	def set_number(self, n):
		self.nr = n


def read_scoring_matrix(filename):
	"""
	Reads the scoring matrix from file.
	Parameters:
	filename (str): the files name/path

	Returns:
	output_dict (dict): a nested dictionary formatted so that output_dict['A']['B']
	gives the score for a (mis)match between amino acids A and B
	"""
	infile = open(filename, 'r')
	output_dict = {}
	for line in infile:
		if line[0] != '#': # ignore comments at the beginning of the file
			if line[0] == ' ':
				sequence = line.split() # store the sequence in which the amino acids appear in the table
			else:
				array = line.split()
				current = array[0] # the current amino acid to be added as key to the dictionary
				for i in range(len(array)):
					if i>0:
						if current not in output_dict:
							output_dict[current] = {sequence[i-1]: float(array[i])}
						else:
							output_dict[current][sequence[i-1]] = float(array[i])
	infile.close()
	return output_dict

def read_sequence(filename, status, input_name):
	"""
	Reads a Fasta-formatted file and returns the sequence as a Sequence instance.
	
	Parameters:
	filename (str): the files name/path
	status (str): Query of Sbjct
	
	Returns:
	sequence (Sequence instance)
	"""
	infile = open(filename, 'r')
	sequence_str = ''
	for line in infile:
		if line[0] == '>':
			if input_name == '': # uses the header in the fasta file as name if the user has not provided one
				name = line[1:].rstrip() # Removes the > sign and stores the name/information about the sequence
			else:
				name = input_name
		else:
			sequence_str = sequence_str + line.rstrip()
	sequence = Sequence(name, status, sequence_str)
	infile.close()
	return sequence
	


def g(gap_opening, gap_extension, gap_length):
	"""
	Calculates and return the gap penalty, which is given by gap_opening + gap_extension*gap_length
	
	Parameters:
	gap_opening (float or int): Gap opening penalty
	gap_extension (float or int): Gap extension penalty
	gap_length (float or int): Gap length
	
	Returns:
	gap_opening + gap_extension * gap_length
	"""
	return gap_opening + gap_extension*gap_length


def compute_alignments(query_filename, subject_filename, scoring_matrix_filename, gap_opening, gap_extension, find_all, query_name, subject_name):
	query = read_sequence(query_filename, 'Query', query_name)
	subject = read_sequence(subject_filename, 'Sbjct', subject_name)
	scoring_dict = read_scoring_matrix(scoring_matrix_filename)
	numbers_matrix = [[0 for j in range(len(subject.sequence_str)+1)] for i in range(len(query.sequence_str)+1)] # initialize number matrix with only zeros
	direction_matrix = [[ ['s'] for j in range(len(subject.sequence_str)+1)] for i in range(len(query.sequence_str)+1)] # initialize direction matrix with only s's
	max_value = 0
	max_positions = [[0,0]]
	number_count = 1
	
	# Compute numbers_matrix and direction_matrix:
	for i in range(1, len(query.sequence_str) +1):
		for j in range(1, len(subject.sequence_str)+1):
			pair = numbers_matrix[i-1][j-1] + scoring_dict[query.sequence_str[i-1]][subject.sequence_str[j-1]] # note that index i-1 and j-1 in query.sequence_str and subject.sequence_str, respectively, corresponds to index [i][j] in the matrix. This is beacuse the matrix has an extra row and an extra collumn corresponding to initial inserts/deletions
			insert = 0
			insert_length = 0
			for k in range(j): # Calculate the total score for insertions of all lengths ending at current position, and store the score and length for the gap with the highest score
				new_ivalue = numbers_matrix[i][k] - g(gap_opening, gap_extension, j-k)
				if new_ivalue > insert:
					insert = new_ivalue
					insert_length = j-k
			delete = 0
			delete_length = 0
			for k in range(i): # Calculate the total score for deletions of all lengths ending at current position, and store the score and length for the gap with the highest score
				new_dvalue = numbers_matrix[k][j] - g(gap_opening, gap_extension, i-k)
				if new_dvalue > delete:
					delete = new_dvalue
					delete_length = i-k
			current_value = max([0, pair, delete, insert]) # Find the maximum of the match/mismatch, deletion, insertion scores and 0.
			if current_value > max_value:
				max_value = current_value
				max_positions = [[i,j]]
			elif current_value == max_value and find_all: # if the program is set to find all optimalt local alignments, the program stores all positions (if more than one) having max value. if not, the program stores only one position
				max_positions.append([i,j])
			cell_value = max([0, pair, delete, insert])
			direction = []
			if cell_value == 0:
				direction.append('s')
			else:
				if cell_value == pair:
					direction.append('p')
				if cell_value == delete:
					direction.append('d' + str(delete_length)) # Deletion "arrows" are stored in the format 'di', where i is an integer.
				if cell_value == insert:
					direction.append('i' + str(insert_length))
			if not find_all: # if the program is not set to find all alignments, it only store one arrow per cell.
				direction = direction[:1]
			numbers_matrix[i][j] = cell_value
			direction_matrix[i][j] = direction
	
	# Find alignments
	alignments = []
	for max_position in max_positions: # Initialize alignments from the starting position(s) (only one maximum position was stored earlier in the program if find_all is False)
		alignments.append(Alignment('', '', '', max_position, max_value, query.name, subject.name))
	found_alignments = []
	now_computing = len(alignments) # This variable keeps track of how many alignments are currently being computed (if the user chooses to compute all alignments)
	while now_computing > 0:
		new_alignments = []
		for a in alignments:
			if not a.done:
				alignmentstart_found = False
				i = a.i
				j = a.j
				current_directions = direction_matrix[i][j]
				for current_direction in current_directions:
					if current_direction == 'p': # if match/mismatch
						q_add = query.sequence_str[i-1] # note that index i-1 and j-1 in query.sequence_str and subject.sequence_str, respectively, corresponds to index [i][j] in the matrix 
						s_add = subject.sequence_str[j-1]
						if q_add == s_add:
							g_add = q_add
						else:
							if scoring_dict[q_add][s_add]>0:
								g_add = '+'
							else:
								g_add = ' '
						i = i-1
						j = j-1
					elif current_direction[0] == 'i': # if insertion
						q_add = ''
						g_add = '' 
						s_add = ''
						length = int(current_direction[1:]) # unpack insertion length
						for l in range(length):
							q_add = '-' + q_add
							g_add = ' ' + g_add
							s_add = subject.sequence_str[j-1-l] + s_add
						j = j-length
					elif current_direction[0] == 'd': # if deletion
						q_add = ''
						g_add = '' 
						s_add = ''
						length = int(current_direction[1:]) # unpack deletion length
						for l in range(length):
							s_add = '-' + s_add
							g_add = ' ' + g_add
							q_add = query.sequence_str[i-1-l] + q_add
						i = i-length
					else: # 's' found, indicating startpoint, i.e. that the alignment computation is done for this alignment
						alignmentstart_found = True
						a.done = True
						a.set_starting_points(i,j)
						a.set_number(number_count)
						a.finalize()
						found_alignments.append(a) # done alignments are stored in the list found_alingments, which will be returned later.
						number_count += 1
						now_computing -= 1
					if not alignmentstart_found: # initialize new alignment that "inherit" the information hitherto computed
						new_alignments.append(Alignment(q_add + a.q_str, s_add + a.s_str, g_add + a.g_str, a.end, a.score, query.name, subject.name, i, j))
						now_computing += 1
				
				alignments = copy.deepcopy(new_alignments) # copy the new list of alignments to the 'alignment' pointer
				now_computing = len(alignments)
	return found_alignments


if __name__ == '__main__': # The script from here to the end is just interaction with the user for input/output when the script is run from terminal
	exit_loop = False
	print('\nWelcome to the sequence alignment program!\nYou will be asked to provide file names/paths for two sequence files in Fasta format, and the file name/path for one scoring matrix in the following format:\n\n# initial lines of comments \n# preceded by "#"\n# if needed\n   A  B  C\nA  1 -1 -2\nB -1  1  0\nC -2  0  1\n\nYou can at any time type exit to exit the program.\n')
	while not exit_loop:
		exit == False
		query_filename = ''
		while not os.path.isfile(query_filename):
			query_filename = input('Input file name/path for sequence 1: ')
			if query_filename == 'exit':
				exit_loop = True
				break
			if not os.path.isfile(query_filename):
				print('Error! File ' +	query_filename + ' not found!')
		if exit_loop:
			break
		
		query_name = input('\nInput short name for sequence 1 (or press enter without writing anything if you want to use the header from the Fasta file) : ')
		if query_name == 'exit':
			break

		subject_filename = ''
		while not os.path.isfile(subject_filename):
			subject_filename = input('\nInput file name/path for sequence 2: ')
			if subject_filename == 'exit':
				exit_loop = True
				break
			if not os.path.isfile(subject_filename):
				print('Error! File ' +	subject_filename + ' not found!')
		if exit_loop:
			break
		
		subject_name = input('\nInput short name for sequence 2 (or press enter without writing anything if you want to use the header from the Fasta file) : ')
		if subject_name == 'exit':
			break
		
		
		matrix_filename = ''
		while not os.path.isfile(matrix_filename):
			matrix_filename = input('\nInput file name/path for scoring matrix: ')
			if matrix_filename == 'exit':
				exit_loop = True
				break
			if not os.path.isfile(matrix_filename):
				print('Error! File ' +	matrix_filename + ' not found!')
		if exit_loop:
			break
	
		usedefault = ''
		while not (usedefault == 'y' or usedefault == 'n'):
			usedefault = input('\n\nDefault settings are: \nGap opening penalty: 11\nGap extension penalty: 1\nFind all optimal alignments (as opposed to only one).\n\nDo you want to use these settings?\nType y if you want to use default setings or type n if you want to change these: ')
			if usedefault == 'exit':
				exit_loop = True
				break
			if not (usedefault == 'y' or usedefault == 'n'):
				print('Please type y or n')
	
		if exit_loop:
			break
	
		if usedefault == 'n':	 
			input_is_number = False
			while not input_is_number:
				gap_opening = input('Input gap opening penalty: ')
				if gap_opening == 'exit':
					exit_loop = True
					break
				try:
					gap_opening = float(gap_opening)
				except ValueError:
					print('Error! Gap opening penalty must be a number!')
				else:
					input_is_number = True
			if exit_loop:
				break
		
			input_is_number = False
			while not input_is_number:
				gap_extension = input('Input gap extension penalty: ')
				if gap_extension == 'exit':
					exit_loop = True
					break
				try:
					gap_extension = float(gap_extension)
				except ValueError:
					print('Error! Gap opening penalty must be a number!')
				else:
					input_is_number = True
			if exit_loop:
				break
		
			find_all_in = ''
			while not (find_all_in == 'y' or find_all_in == 'n'):
				find_all_in = input('\n\nDo you want to find all optimal alignments or only one?\nType y if you want to find all, or n if you want to find only one: ')
				if find_all_in == 'exit':
					exit_loop = True
					break
				if not (find_all_in == 'y' or find_all_in == 'n'):
					print('Please type y or n\n')
			if find_all_in == 'n':
				find_all = False
			else:
				find_all = True
				
		if usedefault == 'y':
			gap_opening = 11.0
			gap_extension = 1.0
			find_all = True
	
		if exit_loop:
			break
		print('\nComputing alignment(s)... Please wait...\n\n')
		
		try:
			alignments = compute_alignments(query_filename, subject_filename, matrix_filename, gap_opening, gap_extension, find_all, query_name, subject_name)
		except:
			print('\n*****ERROR*****\nAn error occured. Please make sure that all your files are in a valid format and that all letters in your sequences are represented in the scoring matrix and try again.\n')
			sys.exit(1)
		output_string = ''
		for a in alignments:
			output_string += a.header
			output_string += a.final
			output_string += '\n'
		print(output_string)
		save_to_file = ''
		while not (save_to_file == 'y' or save_to_file == 'n'):
			save_to_file = input('Do you want to save the alignment(s) to file? y or n: ')
			if save_to_file == 'exit':
				exit_loop == True
				break
			if not (save_to_file == 'y' or save_to_file == 'n'):
				print('Please type y or n\n')
		if save_to_file == 'y':
			filename = ''
			while len(filename)<1:
				filename = input('Please type file name: ')
				if filename == 'exit':
					exit_loop = True
					break
				if os.path.isfile(filename):
					cont = ''
					while not (cont == 'y' or cont == 'n'):
						cont = input('File ' + filename + ' already exsts. Are you sure you want to overwrite it? y or n: ')
						if cont == 'exit':
							exit_loop = True
							break
						if not (cont == 'y' or cont == 'n'):
							print('Please type y or n')
					if cont == 'y':
						outfile = open(filename, 'w')
						outfile.write(output_string)
						outfile.close()
					elif cont == 'n':
						filename == ''
				else:
					outfile = open(filename, 'w')
					outfile.write(output_string)
					outfile.close()
		go_again = ''
		while not (go_again == 'y' or go_again == 'n'):
			go_again = input('\nDo you want to do another sequence alignment?\nType y to do another, or n to exit: ')
			if go_again == 'exit':
				exit_loop = True
				break
			if not (go_again == 'y' or go_again == 'n'):
				print('Please type y or n')
		if go_again == 'n':
			exit_loop = True
			break
		
		
		
		
		
		
		
		
				
				
	
