#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
# import numpy as np
# import copy
# import random
# import platform
# import subprocess
# import sys
# import os
# import pickle
# from multiprocessing import Pool

import modules_1
################################################################################
class SelectSet:
	def __init__(self, nw_cond, target_cond, rnd_cond):
		self.nw_cond = nw_cond

		self.nw_model = nw_cond[0]
		self.strand = nw_cond[1]
		self.n_strand = nw_cond[2]
		self.n_segments = nw_cond[3]
		self.n_cell = nw_cond[4]
		self.n_sc = nw_cond[5]
		self.l_bond = nw_cond[6]
		self.c_n = nw_cond[7]

		self.target_cond = target_cond

		self.restart = rnd_cond[0] 
		self.cond_top = rnd_cond[1]
		self.n_hist = rnd_cond[2]

	############################
	def select_set(self):
		# ネットワークを設定
		if self.nw_model == "Regular_NW":
			calcd_data_dic = self.regnw_setup()
		else:
			calcd_data_dic = self.rndnw_setup()

		return calcd_data_dic
	
	#######################################
	def regnw_setup(self):
		nwsetup = modules_1.RegNWSetup.NWSetup(self.nw_cond, self.target_cond)
		calcd_data_dic = nwsetup.calc_all()
		return calcd_data_dic

	def rndnw_setup(self):
		make8 = modules_1.RndNWSetup.Make8(self.nw_cond)
		base_top_list = make8.make_8chain_dic()

		if self.restart == '':
			# トポロジーの異なるネットワークを探索して、任意の多重度のネットワークポリマーの代数的連結性の分布関数を策定
			mod = modules_1.RndNWSetup.ModifyTop(base_top_list, self.nw_cond, self.cond_top, self.target_cond, self.n_hist)
			candidate_list, target_dir = mod.top_search()
		else:
			sel = modules_1.RndNWSetup.Select(self.restart, self.n_hist, self.target_cond)
			candidate_list, target_dir = sel.top_select()

		sel = modules_1.RndNWSetup.Select(self.restart, self.n_hist, self.target_cond)
		top_dic_list = sel.nw_search(candidate_list, target_dir)

		###########################################
		# ターゲットとなるネットワーク全体の辞書を設定。
		setup = modules_1.RndNWSetup.SetUp(top_dic_list, base_top_list, self.n_segments, self.n_sc)
		calcd_data_dic = setup.make_data_dic()

		return calcd_data_dic
