#! /usr/bin/python3.9
import brain
import brain_util as bu
import numpy as np
#import pptree
import json
import copy

from collections import namedtuple
from collections import defaultdict
from enum import Enum








# BrainAreas
LEX = "LEX"
DET = "DET"
#SUBJ = "SUBJ"
OBJ = "OBJ"
#VERB = "VERB"
PREP = "PREP"
PREP_P = "PREP_P"
ADJ = "ADJ"
ADVERB = "ADVERB"
VERB_PLUR = "VERB_PLUR"
SUBJ_PLUR = "SUBJ_PLUR"
VERB_SING = "VERB_SING"
SUBJ_SING = "SUBJ_SING"



# Unique to Russian
NOM = "NOM"
ACC = "ACC"
DAT = "DAT"

# Fixed area stats for explicit areas
LEX_SIZE = 21

# Actions
DISINHIBIT = "DISINHIBIT"
INHIBIT = "INHIBIT"
# Skip firing in this round, just activate the word in LEX/DET/other word areas.
# All other rules for these lexical items should be in PRE_RULES.
ACTIVATE_ONLY = "ACTIVATE_ONLY"
CLEAR_DET = "CLEAR_DET"

#AREAS = [LEX, LEX_SING, DET, SUBJ, OBJ, VERB, ADJ, ADVERB, PREP, PREP_P]
AREAS = [LEX, DET, SUBJ_PLUR, OBJ, VERB_PLUR, ADJ, ADVERB, PREP, PREP_P, VERB_SING, SUBJ_SING]

EXPLICIT_AREAS = [LEX]
RECURRENT_AREAS = [OBJ, ADJ, ADVERB, PREP, PREP_P, VERB_PLUR, SUBJ_PLUR, VERB_SING, SUBJ_SING]



AreaRule = namedtuple('AreaRule', ['action', 'area', 'index'])
FiberRule = namedtuple('FiberRule', ['action', 'area1', 'area2', 'index'])
FiringRule = namedtuple('FiringRule', ['action'])
OtherRule = namedtuple('OtherRule', ['action'])

def generic_noun(index, form):

	if form == 'plural':
		VERB =  VERB_PLUR
		SUBJ = SUBJ_PLUR

	if form == 'singular':
		VERB = VERB_SING 
		SUBJ = SUBJ_SING 
	

	return {
		"index": index,
		"PRE_RULES": [
		FiberRule(DISINHIBIT, LEX, SUBJ, 0), 
		FiberRule(DISINHIBIT, LEX, OBJ, 0),
		FiberRule(DISINHIBIT, LEX, PREP_P, 0),
		FiberRule(DISINHIBIT, DET, SUBJ, 0),
		FiberRule(DISINHIBIT, DET, OBJ, 0),
		FiberRule(DISINHIBIT, DET, PREP_P, 0),
		FiberRule(DISINHIBIT, ADJ, SUBJ, 0),
		FiberRule(DISINHIBIT, ADJ, OBJ, 0),
		FiberRule(DISINHIBIT, ADJ, PREP_P, 0),
		FiberRule(DISINHIBIT, VERB, OBJ, 0),
		FiberRule(DISINHIBIT, PREP_P, PREP, 0),
		FiberRule(DISINHIBIT, PREP_P, SUBJ, 0),
		FiberRule(DISINHIBIT, PREP_P, OBJ, 0),
		],
		"POST_RULES": [
		AreaRule(INHIBIT, DET, 0),
		AreaRule(INHIBIT, ADJ, 0),
		AreaRule(INHIBIT, PREP_P, 0),
		AreaRule(INHIBIT, PREP, 0),
		FiberRule(INHIBIT, LEX, SUBJ, 0),
		FiberRule(INHIBIT, LEX, OBJ, 0),
		FiberRule(INHIBIT, LEX, PREP_P, 0),
		FiberRule(INHIBIT, ADJ, SUBJ, 0),
		FiberRule(INHIBIT, ADJ, OBJ, 0),
		FiberRule(INHIBIT, ADJ, PREP_P, 0),
		FiberRule(INHIBIT, DET, SUBJ, 0),
		FiberRule(INHIBIT, DET, OBJ, 0),
		FiberRule(INHIBIT, DET, PREP_P, 0),
		FiberRule(INHIBIT, VERB, OBJ, 0),
		FiberRule(INHIBIT, PREP_P, PREP, 0),
		FiberRule(INHIBIT, PREP_P, VERB, 0),
		FiberRule(DISINHIBIT, LEX, SUBJ, 1),
		FiberRule(DISINHIBIT, LEX, OBJ, 1),
		FiberRule(DISINHIBIT, DET, SUBJ, 1),
		FiberRule(DISINHIBIT, DET, OBJ, 1),
		FiberRule(DISINHIBIT, ADJ, SUBJ, 1),
		FiberRule(DISINHIBIT, ADJ, OBJ, 1),
		FiberRule(INHIBIT, PREP_P, SUBJ, 0),
		FiberRule(INHIBIT, PREP_P, OBJ, 0),
		FiberRule(INHIBIT, VERB, ADJ, 0),
		]
	}

def generic_trans_verb(index, form):

	if form == 'singular':
		VERB = VERB_SING 
		SUBJ = SUBJ_SING 

	if form == 'plural':
		VERB = VERB_PLUR 
		SUBJ = SUBJ_PLUR 

	return {
		"index": index,
		"PRE_RULES": [
		FiberRule(DISINHIBIT, LEX, VERB, 0),
		FiberRule(DISINHIBIT, VERB, SUBJ, 0),
		FiberRule(DISINHIBIT, VERB, ADVERB, 0),
		AreaRule(DISINHIBIT, ADVERB, 1),
		],
		"POST_RULES": [
		FiberRule(INHIBIT, LEX, VERB, 0),
		AreaRule(DISINHIBIT, OBJ, 0),
		AreaRule(INHIBIT, SUBJ, 0),
		AreaRule(INHIBIT, ADVERB, 0),
		FiberRule(DISINHIBIT, PREP_P, VERB, 0),
		]
	}	

	'''''''''''''''''


	if form == 'plural':
		
		return {
		"index": index,
		"PRE_RULES": [
		FiberRule(DISINHIBIT, LEX, VERB_SING, 0),
		FiberRule(DISINHIBIT, LEX, VERB_PLUR, 0),
		FiberRule(DISINHIBIT, VERB_SING, SUBJ_SING, 0),
		FiberRule(DISINHIBIT, VERB_PLUR, SUBJ_PLUR, 0),
		FiberRule(DISINHIBIT, VERB_SING, ADVERB, 0),
		FiberRule(DISINHIBIT, VERB_PLUR, ADVERB, 0),
		AreaRule(DISINHIBIT, ADVERB, 1),
		],
		"POST_RULES": [
		FiberRule(DISINHIBIT, LEX, VERB_SING, 0),
		FiberRule(DISINHIBIT, LEX, VERB_PLUR, 0),
		AreaRule(DISINHIBIT, OBJ, 0),
		AreaRule(INHIBIT, SUBJ_SING, 0),
		AreaRule(INHIBIT, SUBJ_PLUR, 0),
		AreaRule(INHIBIT, ADVERB, 0),
		FiberRule(DISINHIBIT, PREP_P, VERB_SING, 0),
		FiberRule(DISINHIBIT, PREP_P, VERB_PLUR, 0)
		]
		}

	'''''''''''
	

def generic_intrans_verb(index, form):

	if form == 'plural':
		VERB =  VERB_PLUR
		SUBJ = SUBJ_PLUR

	if form == 'singular':
		VERB = VERB_SING 
		SUBJ = SUBJ_SING 
	

	return {
		"index": index,
		"PRE_RULES": [
		FiberRule(DISINHIBIT, LEX, VERB, 0),
		FiberRule(DISINHIBIT, VERB, SUBJ, 0),
		FiberRule(DISINHIBIT, VERB, ADVERB, 0),
		AreaRule(DISINHIBIT, ADVERB, 1),
		],
		"POST_RULES": [
		FiberRule(INHIBIT, LEX, VERB, 0),
		AreaRule(INHIBIT, SUBJ, 0),
		AreaRule(INHIBIT, ADVERB, 0),
		FiberRule(DISINHIBIT, PREP_P, VERB, 0),
		]
	}

def generic_copula(index, form):

	if form == 'plural':
		VERB =  VERB_PLUR
		SUBJ = SUBJ_PLUR

	if form == 'singular':
		VERB = VERB_SING 
		SUBJ = SUBJ_SING 
	

	return {
		"index": index,
		"PRE_RULES": [
		FiberRule(DISINHIBIT, LEX, VERB, 0),
		FiberRule(DISINHIBIT, VERB, SUBJ, 0),
		],
		"POST_RULES": [
		FiberRule(INHIBIT, LEX, VERB, 0),
		AreaRule(DISINHIBIT, OBJ, 0),
		AreaRule(INHIBIT, SUBJ, 0),
		FiberRule(DISINHIBIT, ADJ, VERB, 0)
		]
	}

def generic_adverb(index, form):


	return {
		"index": index,
		"PRE_RULES": [
		AreaRule(DISINHIBIT, ADVERB, 0),
		FiberRule(DISINHIBIT, LEX, ADVERB, 0)
		],
		"POST_RULES": [
		FiberRule(INHIBIT, LEX, ADVERB, 0),
		AreaRule(INHIBIT, ADVERB, 1),
		]

	}

def generic_determinant(index, form):


	return {
		"index": index,
		"PRE_RULES": [
		AreaRule(DISINHIBIT, DET, 0),
		FiberRule(DISINHIBIT, LEX, DET, 0)
		],
		"POST_RULES": [
		FiberRule(INHIBIT, LEX, DET, 0),
		FiberRule(INHIBIT, VERB_PLUR, ADJ, 0),
		FiberRule(INHIBIT, VERB_SING, ADJ, 0)
		]
	}

def generic_adjective(index, form):
	if form == 'plural':
		VERB =  VERB_PLUR


	if form == 'singular':
		VERB = VERB_SING 

	

	return {
		"index": index,
		"PRE_RULES": [
		AreaRule(DISINHIBIT, ADJ, 0),
		FiberRule(DISINHIBIT, LEX, ADJ, 0)
		],
		"POST_RULES": [
		FiberRule(INHIBIT, LEX, ADJ, 0),
		FiberRule(INHIBIT, VERB, ADJ, 0),
		]

	}

def generic_preposition(index, form):

	return {
		"index": index,
		"PRE_RULES": [
			AreaRule(DISINHIBIT, PREP, 0),
			FiberRule(DISINHIBIT, LEX, PREP, 0),
		],
		"POST_RULES": [
			FiberRule(INHIBIT, LEX, PREP, 0),
			AreaRule(DISINHIBIT, PREP_P, 0),
			FiberRule(INHIBIT, LEX, SUBJ_SING, 1),
			FiberRule(INHIBIT, LEX, SUBJ_PLUR, 1),
			FiberRule(INHIBIT, LEX, OBJ, 1),
			FiberRule(INHIBIT, DET, SUBJ_SING, 1),
			FiberRule(INHIBIT, DET, SUBJ_PLUR, 1),
			FiberRule(INHIBIT, DET, OBJ, 1),
			FiberRule(INHIBIT, ADJ, SUBJ_SING, 1),
			FiberRule(INHIBIT, ADJ, SUBJ_PLUR, 1),
			FiberRule(INHIBIT, ADJ, OBJ, 1),
		]
	}

LEXEME_DICT = {
	"the" : generic_determinant(0, 'singular'),
	"a": generic_determinant(1, 'singular'),
	"dogs" : generic_noun(2, 'plural'),
	"cat" : generic_noun(3, 'singular'),
	"mice" : generic_noun(4, 'singular'),
	"people" : generic_noun(5, 'plural'),
	"chase" : generic_trans_verb(6, 'plural'),
	"loves" : generic_trans_verb(7, 'singular'),
	"bites" : generic_trans_verb(8, 'singular'),
	"of" : generic_preposition(9, 'singular'),
	"big": generic_adjective(10, 'singular'),
	"bad": generic_adjective(11, 'singular'),
	"run": generic_intrans_verb(12, 'singular'),
	"fly": generic_intrans_verb(13, 'singular'),
	"quickly": generic_adverb(14, 'singular'),
	"in": generic_preposition(15, 'singular'),
	"are": generic_copula(16, 'singular'),
	"man": generic_noun(17, 'singular'),
	"woman": generic_noun(18, 'singular'),
	"saw": generic_trans_verb(19, 'plural'),
	"sees": generic_trans_verb(20, 'singular')
}



ENGLISH_READOUT_RULES = {
	VERB_SING: [LEX,  SUBJ_SING, OBJ, PREP_P, ADVERB, ADJ],
	SUBJ_SING: [LEX, DET, ADJ, PREP_P],
	VERB_PLUR: [LEX,  SUBJ_PLUR, OBJ, PREP_P, ADVERB, ADJ],
	SUBJ_PLUR: [LEX, DET, ADJ, PREP_P],
	OBJ: [LEX, DET, ADJ, PREP_P],
	PREP_P: [LEX,  PREP, ADJ, DET],
	PREP: [LEX],
	ADJ: [LEX],
	DET: [LEX],
	ADVERB: [LEX],
	LEX: [],
}


class ParserBrain(brain.Brain):
	def __init__(self, p, lexeme_dict={}, all_areas=[], recurrent_areas=[], initial_areas=[], readout_rules={}):
		brain.Brain.__init__(self, p)
		self.lexeme_dict = lexeme_dict
		self.all_areas = all_areas
		self.recurrent_areas = recurrent_areas
		self.initial_areas = initial_areas

		self.fiber_states = defaultdict()
		self.area_states = defaultdict(set)
		self.activated_fibers = defaultdict(set)
		self.readout_rules = readout_rules
		self.initialize_states()

	def initialize_states(self):
		for from_area in self.all_areas:
			self.fiber_states[from_area] = defaultdict(set)
			for to_area in self.all_areas:
				self.fiber_states[from_area][to_area].add(0)

		for area in self.all_areas:
			self.area_states[area].add(0)

		for area in self.initial_areas:
			self.area_states[area].discard(0)

	def applyFiberRule(self, rule):
		if rule.action == INHIBIT:
			self.fiber_states[rule.area1][rule.area2].add(rule.index)
			self.fiber_states[rule.area2][rule.area1].add(rule.index)
		elif rule.action == DISINHIBIT:
			self.fiber_states[rule.area1][rule.area2].discard(rule.index)
			self.fiber_states[rule.area2][rule.area1].discard(rule.index)

	def applyAreaRule(self, rule):
		if rule.action == INHIBIT:
			self.area_states[rule.area].add(rule.index)
		elif rule.action == DISINHIBIT:
			self.area_states[rule.area].discard(rule.index)

	def applyRule(self, rule):
		if isinstance(rule, FiberRule):
			self.applyFiberRule(rule)
			return True
		if isinstance(rule, AreaRule):
			self.applyAreaRule(rule)
			return True
		return False

	def parse_project(self):
		project_map = self.getProjectMap()
		self.remember_fibers(project_map)
		self.project({}, project_map)

	# For fiber-activation readout, remember all fibers that were ever fired.
	def remember_fibers(self, project_map):
		for from_area, to_areas in project_map.items():
			self.activated_fibers[from_area].update(to_areas)

	def recurrent(self, area):
		return (area in self.recurrent_areas)

	# TODO: Remove brain from ProjectMap somehow
	# perhaps replace Parser state with ParserBrain:Brain, better design
	def getProjectMap(self):
		proj_map = defaultdict(set)
		for area1 in self.all_areas:
			if len(self.area_states[area1]) == 0:
				for area2 in self.all_areas:
					if area1 == LEX and area2 == LEX:
						continue
					if len(self.area_states[area2]) == 0:
						if len(self.fiber_states[area1][area2]) == 0:
							if self.areas[area1].winners:
								proj_map[area1].add(area2)
							if self.areas[area2].winners:
								proj_map[area2].add(area2)
		return proj_map

	def activateWord(self, area_name, word):
		area = self.areas[area_name]
		k = area.k

		assembly_start = self.lexeme_dict[word]["index"]*k
		area.winners = list(range(assembly_start, assembly_start+k))
		area.fix_assembly()

	def activateIndex(self, area_name, index):
		area = self.areas[area_name]
		k = area.k
		assembly_start = index*k
		area.winners = list(range(assembly_start, assembly_start+k))
		area.fix_assembly()

	def interpretAssemblyAsString(self, area_name):
		return self.getWord(area_name, 0.7)

	def getWord(self, area_name, min_overlap=0.7):
		if not self.areas[area_name].winners:
			raise Exception("Cannot get word because no assembly in " + area_name)
		winners = set(self.areas[area_name].winners)
		area_k = self.areas[area_name].k
		threshold = min_overlap * area_k
		for word, lexeme in self.lexeme_dict.items():
			word_index = lexeme["index"]
			word_assembly_start = word_index * area_k
			word_assembly = set(range(word_assembly_start, word_assembly_start + area_k))
			if len((winners & word_assembly)) >= threshold:
				return word
		return None

	def getActivatedFibers(self):
		# Prune activated_fibers pased on the readout_rules
		pruned_activated_fibers = defaultdict(set)
		for from_area, to_areas in self.activated_fibers.items():
			for to_area in to_areas:
				if to_area in self.readout_rules[from_area]:
					pruned_activated_fibers[from_area].add(to_area)

		return pruned_activated_fibers


class EnglishParserBrain(ParserBrain):
	def __init__(self, p, non_LEX_n=10000, non_LEX_k=100, LEX_k=21, 
		default_beta=0.2, LEX_beta=1.0, recurrent_beta=0.05, interarea_beta=0.5, verbose=False):
		ParserBrain.__init__(self, p, 
			lexeme_dict=LEXEME_DICT,
			all_areas=AREAS, 
			recurrent_areas=RECURRENT_AREAS, 
			initial_areas=[LEX, SUBJ_SING, VERB_SING, SUBJ_PLUR, VERB_PLUR],
			readout_rules=ENGLISH_READOUT_RULES)
		self.verbose = verbose

		## singular
		LEX_n = LEX_SIZE * LEX_k
		self.add_explicit_area(LEX, LEX_n, LEX_k, default_beta)


		DET_k = LEX_k
		self.add_area(SUBJ_SING, non_LEX_n, non_LEX_k, default_beta)
		self.add_area(SUBJ_PLUR, non_LEX_n, non_LEX_k, default_beta)
		#self.add_area(SUBJ, non_LEX_n, non_LEX_k, default_beta)
		self.add_area(OBJ, non_LEX_n, non_LEX_k, default_beta)
		#self.add_area(VERB, non_LEX_n, non_LEX_k, default_beta)
		self.add_area(VERB_SING, non_LEX_n, non_LEX_k, default_beta)
		self.add_area(VERB_PLUR, non_LEX_n, non_LEX_k, default_beta)
		self.add_area(ADJ, non_LEX_n, non_LEX_k, default_beta)
		self.add_area(PREP, non_LEX_n, non_LEX_k, default_beta)
		self.add_area(PREP_P, non_LEX_n, non_LEX_k, default_beta)
		self.add_area(DET, non_LEX_n, DET_k, default_beta)
		self.add_area(ADVERB, non_LEX_n, non_LEX_k, default_beta)





		# LEX: all areas -> * strong, * -> * can be strong
		# non LEX: other areas -> * (?), LEX -> * strong, * -> * weak
		# DET? Should it be different?
		custom_plasticities = defaultdict(list)
		for area in RECURRENT_AREAS:
			custom_plasticities[LEX].append((area, LEX_beta))
			custom_plasticities[area].append((LEX, LEX_beta))
			custom_plasticities[area].append((area, recurrent_beta))
			for other_area in RECURRENT_AREAS:
				if other_area == area:
					continue
				custom_plasticities[area].append((other_area, interarea_beta))

		self.update_plasticities(area_update_map=custom_plasticities)

	def getProjectMap(self):
		proj_map = ParserBrain.getProjectMap(self)
		# "War of fibers"
		if LEX in proj_map and len(proj_map[LEX]) > 4:  # because LEX->LEX
			raise Exception("Got that LEX projecting into many areas: " + str(proj_map[LEX]))
		return proj_map


	def getWord(self, area_name, min_overlap=0.7):
		word = ParserBrain.getWord(self, area_name, min_overlap)
		if word:
			return word
		if not word and area_name == DET:
			winners = set(self.areas[area_name].winners)
			area_k = self.areas[area_name].k
			threshold = min_overlap * area_k
			nodet_index = DET_SIZE - 1
			nodet_assembly_start = nodet_index * area_k
			nodet_assembly = set(range(nodet_assembly_start, nodet_assembly_start + area_k))
			if len((winners & nodet_assembly)) > threshold:
				return "<null-det>"
		# If nothing matched, at least we can see that in the parse output.
		return "<NON-WORD>"



class ParserDebugger():
	def __init__(self, brain, all_areas, explicit_areas):
		self.b = brain
		self.all_areas = all_areas
		self.explicit_areas = explicit_areas

	def run(self):
		command = input("DEBUGGER: ENTER to continue, 'P' for PEAK \n")
		while command:
			if command == "P":
				self.peak()
				return
			elif command:
				print("DEBUGGER: Command not recognized...")
				command = input("DEBUGGER: ENTER to continue, 'P' for PEAK \n")
			else:
				return

	def peak(self):
		remove_map = defaultdict(int)
		# Temporarily set beta to 0
		self.b.no_plasticity = True
		self.b.save_winners = True

		for area in self.all_areas:
			self.b.areas[area].unfix_assembly()
		while True:
			test_proj_map_string = input("DEBUGGER: enter projection map, eg. {\"VERB\": [\"LEX\"]}, or ENTER to quit\n")
			if not test_proj_map_string:
				break
			test_proj_map = json.loads(test_proj_map_string)
			# Important: save winners to later "remove" this test project round 
			to_area_set = set()
			for _, to_area_list in test_proj_map.items():
				for to_area in to_area_list:
					to_area_set.add(to_area)
					if not self.b.areas[to_area].saved_winners:
						self.b.areas[to_area].saved_winners.append(self.b.areas[to_area].winners)

			for to_area in to_area_set:
				remove_map[to_area] += 1

			self.b.project({}, test_proj_map)
			for area in self.explicit_areas:
				if area in to_area_set:
					area_word = self.b.interpretAssemblyAsString(area)
					print("DEBUGGER: in explicit area " + area + ", got: " + area_word)

			print_assemblies = input("DEBUGGER: print assemblies in areas? Eg. 'LEX,VERB' or ENTER to cont\n")
			if not print_assemblies:
				continue
			for print_area in print_assemblies.split(","):
				print("DEBUGGER: Printing assembly in area " + print_area)
				print(str(self.b.areas[print_area].winners))
				if print_area in self.explicit_areas:
					word = self.b.interpretAssemblyAsString(print_area)
					print("DEBUGGER: in explicit area got assembly = " + word)

		# Restore assemblies (winners) and w values to before test projections
		for area, num_test_projects in remove_map.items():
			self.b.areas[area].winners = self.b.areas[area].saved_winners[0]
			self.b.areas[area].w = self.b.areas[area].saved_w[-num_test_projects - 1]
			self.b.areas[area].saved_w = self.b.areas[area].saved_w[:(-num_test_projects)]
		self.b.no_plasticity = False
		self.b.save_winners = False
		for area in self.all_areas:
			self.b.areas[area].saved_winners = []

	

# strengthen the assembly representing this word in LEX
# possibly useful way to simulate long-term potentiated word assemblies 
# so that they are easily completed.
def potentiate_word_in_LEX(b, word, rounds=21):
	b.activateWord(LEX, word)
	for _ in range(21):
		b.project({}, {LEX: [LEX]})

# "dogs chase cats" experiment, what should happen?
# simplifying assumption 1: after every project round, freeze assemblies
# exp version 1: area not fired into until LEX fires into it 
# exp version 2: project between all disinhibited fibers/areas, forming some "ghosts"

# "dogs": open fibers LEX<->SUBJ and LEX<->OBJ but only SUBJ disinhibited
# results in "dogs" assembly in LEX<->SUBJ (reciprocal until stable, LEX frozen)
# in version 2 would also have SUBJ<->VERB, so LEX<->SUBJ<->VERB overall

# "chase": opens fibers LEX<->VERB and VERB<->OBJ, inhibit SUBJ, disi
# results in "chase" assembly in LEX<->VERB
# in version 2 would also havee VERB<->OBJ

# "cats": 


# Readout types
class ReadoutMethod(Enum):
	FIXED_MAP_READOUT = 1
	FIBER_READOUT = 2
	NATURAL_READOUT = 3





def parse(sentence="mice loves cat", language="English", p=0.1, LEX_k=21, LEX_k_p = 5,
	project_rounds=21, verbose=False, debug=False, readout_method=ReadoutMethod.FIBER_READOUT):

	if language == "English":
		b = EnglishParserBrain(p, LEX_k=LEX_k, verbose=verbose)
		lexeme_dict = LEXEME_DICT
		all_areas = AREAS
		explicit_areas = EXPLICIT_AREAS
		readout_rules = ENGLISH_READOUT_RULES


	parseHelper(b, sentence, p, LEX_k, project_rounds, verbose, debug, 
		lexeme_dict,all_areas, explicit_areas, readout_method, readout_rules)


def parseHelper(b, sentence, p, LEX_k, project_rounds, verbose, debug, 
	lexeme_dict, all_areas, explicit_areas, readout_method, readout_rules):
	debugger = ParserDebugger(b, all_areas, explicit_areas)

	sentence = sentence.split(" ")

	extreme_debug = False

	for word in sentence:
		if word not in lexeme_dict:
			print(word, 'is not in dictionnary, moving to the next word')
			continue


		lexeme = lexeme_dict[word]
		b.activateWord(LEX, word)
		if verbose:
			print("Activated word: " + word)
			print(b.areas[LEX].winners)

		for rule in lexeme["PRE_RULES"]:
			b.applyRule(rule)

		proj_map = b.getProjectMap()
		for area in proj_map:
			if area not in proj_map[LEX]:
				b.areas[area].fix_assembly()
				if verbose:
					print("FIXED assembly bc not LEX->this area in: " + area)
			elif area != LEX:
				b.areas[area].unfix_assembly()
				b.areas[area].winners = []
				if verbose:
					print("ERASED assembly because LEX->this area in " + area)

		proj_map = b.getProjectMap()
		if verbose:
			print("Got proj_map = ")
			print(proj_map)

		for i in range(project_rounds):
			b.parse_project()
			if verbose:
				proj_map = b.getProjectMap()
				print("Got proj_map = ")
				print(proj_map)
			if extreme_debug and word == "a":
				print("Starting debugger after round " + str(i) + "for word" + word)
				debugger.run()

		#if verbose:
		#	print("Done projecting for this round")
		#	for area_name in all_areas:
		#		print("Post proj stats for " + area_name)
		#		print("w=" + str(b.areas[area_name].w))
		#		print("num_first_winners=" + str(b.areas[area_name].num_first_winners))

		for rule in lexeme["POST_RULES"]:
			b.applyRule(rule)

		if debug:
			print("Starting debugger after the word " + word)
			debugger.run()
			

	# Readout
	# For all readout methods, unfix assemblies and remove plasticity.
	b.no_plasticity = True
	for area in all_areas:
		b.areas[area].unfix_assembly()

	dependencies = []
	def read_out(area, mapping):
		to_areas = mapping[area]
		b.project({}, {area: to_areas})
		this_word = b.getWord(LEX)

		for to_area in to_areas:
			if to_area == LEX:
				continue
			b.project({}, {to_area: [LEX]})
			other_word = b.getWord(LEX)

			if other_word ==  '<NON-WORD>':
				print('syntactic error')


			dependencies.append([this_word, other_word, to_area])

		for to_area in to_areas:
			if to_area != LEX:
				read_out(to_area, mapping)


	def treeify(parsed_dict, parent):
		for key, values in parsed_dict.items():
			key_node = pptree.Node(key, parent)
			if isinstance(values, str):
				_ = pptree.Node(values, key_node)
			else:
				treeify(values, key_node)

	

	if readout_method == ReadoutMethod.FIXED_MAP_READOUT:
		# Try "reading out" the parse.
		# To do so, start with final assembly in VERB
		# project VERB->SUBJ,OBJ,LEX

		parsed = {VERB_SING: read_out(VERB_PLUR, readout_rules)}

		print("Final parse dict: ")
		print(parsed)

		root = pptree.Node(VERB_SING)
		treeify(parsed[VERB_SING], root)



	if readout_method == ReadoutMethod.FIXED_MAP_READOUT:
		# Try "reading out" the parse.
		# To do so, start with final assembly in VERB
		# project VERB->SUBJ,OBJ,LEX

		parsed = {VERB_PLUR: read_out(VERB_PLUR, readout_rules)}

		print("Final parse dict: ")
		print(parsed)

		root = pptree.Node(VERB_PLUR)
		treeify(parsed[VERB_PLUR], root)




		

	if readout_method == ReadoutMethod.FIBER_READOUT:
		activated_fibers = b.getActivatedFibers()
		if verbose:
			print("Got activated fibers for readout:")
			print(activated_fibers)

		read_out(VERB_SING, activated_fibers)
		print("Got dependencies: ")
		print(dependencies)

	if readout_method == ReadoutMethod.FIBER_READOUT:
		activated_fibers = b.getActivatedFibers()
		if verbose:
			print("Got activated fibers for readout:")
			print(activated_fibers)

		read_out(VERB_PLUR, activated_fibers)
		print("Got dependencies: ")
		print(dependencies)

		# root = pptree.Node(VERB)
		#treeify(parsed[VERB], root)

	# pptree.print_tree(root)


def main():
    parse()

if __name__ == "__main__":
    main()


# TODOs
# BRAIN
# fix brain.py to work when no-assembly areas are projected in 

# PARSER IMPLEMENTATION
# Factor out debugger of parse
# Factor out read-out, possibly other aspects of parse
# consider areas where only A->B needed not A<->B, easy to fix
# for example, SUBJ/OBJ->DET, etc?

# PARSER CONCEPTUAL
# 1) NATURAL READ OUT: 
	# "Fiber-activation read out": Remember fibers that were activated
	# "Lexical-item read out": Get word from V, see rules (not sufficient but recovers basic structure)

# 2) PREP area: of, others
# "brand of toys", to merge brand<->of<->toys, look for activated noun areas
# for example if OBJ is the only one, we're done
# if multiple, recency? (first instance of lookahead/memory!)

# 3) Intransitive verbs (in particular wrt read out)

# RESEARCH IDEAS
# 1) Russian experiment (free word order)
# 2) Grammaticality, detect some sort of error for non-grammatical


