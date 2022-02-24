#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################
import os
import sys
import codecs
from UDFManager import UDFManager
#######################################################
#
def setupcondition():
	# check 'calc_condition.udf' and make it.
	findudf()
	# Read udf and setup initial conditions.
	basic_cond, nw_cond, sim_cond, rnd_cond, target_cond = read_and_setcondition()

	return basic_cond, nw_cond, sim_cond, rnd_cond, target_cond
	

###########################################
# check 'calc_condition.udf' and make it.
def findudf():
	if not os.path.isfile('./calc_condition.udf'):
		print()
		print('In this directory, no "calc_condition.udf" is found !')
		print('New one will be generated.')
		print('Please, modify and save it !\n')
		makenewudf()
		input('Press ENTER to continue...')
	return
###########################################
# make new udf when not found.
def makenewudf():
	contents = '''
	\\begin{def}
		CalcCond:{
			Cognac_ver:select{"cognac112"} "使用する Cognac のバージョン",
			Cores: int "計算に使用するコア数を指定"
			} "計算の条件を設定"
		TargetCond:{
			Model:{TargetModel:select{"Regular_NW", "Random_NW"} "ネットワークのモデルを選択",
				Regular_NW:{chains:select{"3_Chain_S", "3_Chain_D", "4_Chain", "6_Chain", "8_Chain"} "分岐の数と種類を選択"
					} "規則構造での条件を入力",
				Random_NW:{chains:select{"3_Chain", "4_Chain", "5_Chain", "6_Chain", "7_Chain"} "分岐の数と種類を選択",
					Calc_Topolpgy:select{"Calc", "Read"} "ランダムネットワークの「計算を行うか、読み込むか」を選択",
						Calc:{pre_sampling:int "プレサンプリング数", pre_try:int "プレサンプリング時の再トライ数", sampling:int "サンプリング数", try:int "サンプリング時の再トライ数", n_parallel:int "並行計算のCPU数"} "ランダムサーチ計算する場合の条件を設定",
						Read:{dir_name:string} "過去の計算結果のディレクトリを記入",
					N_histgram:int "ヒストグラムの分割数"
					} "ランダム構造での条件を入力"
				} "シミュレーションの条件を設定"
			NetWork:{N_Segments: int "ストランド中のセグメント数", 
					N_Subchain: int "各セグメントの側鎖の数", 
					N_UnitCells: int "一辺あたりのユニットセルの数"
				} "ネットワークの条件を設定"
			Multiplisity:{Set_or_Calc:select{"Set", "Calc"} "多重度を設定するかどうかのフラッグ",
					Set:{Multiplicity: int} "多重度を設定",
					Calc:{TargetDensity:float} "多重度を自動設定した場合の密度を設定 \\n設定した密度になるように多重度を設定"
				} "多重度設定に関する設定"
			Shrinkage:{Shrink:select{"Yes", "No"} "ストランドを自然長から圧縮するかどうかのフラッグ \\n非圧縮時には、多重度に応じて密度が変化",
				Yes:{Control:select{"Density", "Shrink"} "圧縮する場合に、密度コントロールにするか、圧縮率を決めるかを設定", 
				Density:{target_density: float} "目標とする密度を設定", 
				Shrinkage:{value: float} "ストランドの圧縮比率を設定"
				}
				} "ストランドを自然長から圧縮するかどうかを設定"
			TopologyType:{
				Type:select{"Entangled", "NO_Entangled"} "ネットワーク・トポロジーを選択",
					Entangled:{Step_rfc[]: float "Slow Push Off での rfc 条件"} "密度、末端間距離を設定値に合わせるように多重度を自動設定。\\n絡み合いが入るように初期化",
					NO_Entangled:{ExpansionRatio: float "NPT 計算での初期膨張率", StepPress[]: float "NPT 計算での圧力変化"} "密度、末端間距離を設定値に合わせるように多重度を自動設定。\\n絡み合いが入らないようにNPTで縮める。"
				} "ネットワーク・トポロジーを選択",
			} "計算ターゲットの条件を設定"
		SimulationCond:{
			Equilib_Condition:{
					repeat: int "平衡化計算の繰り返し数",
					Time:{delta_T: double, Total_Steps: int, Output_Interval_Steps: int} "平衡化計算の時間条件を入力"
				} "平衡化計算の時間条件を入力",
			l_bond: float "シミュレーションでのボンドの自然長"
			} "シミュレーションの条件を設定"
	\end{def}

	\\begin{data}
		CalcCond:{"cognac112",1}
TargetCond:{
	{"Random_NW", {"4_Chain"}{"4_Chain","Read",{100,100,100,100,1}{"4_chains_3_cells_100_trials_100_sampling"}50}}
	{20, 0, 3}
	{"Set", {1}{0.85}}
	{"Yes", {"Density", {0.85}{0.5}}}
	{"Entangled",
		{[1.0730000,1.0000000,0.9000000,0.8000000]},
		{2.0000000, [0.2000000,0.5000000,1.0000000,2.0000000,3.0000000,4.5000000]}
		}
	}
SimulationCond:{
	{4,{1.00000000000000e-02,200000,1000}}
	0.97
	}

\end{data}
	'''
	###
	with codecs.open('./calc_condition.udf', 'w', 'utf_8') as f:
		f.write(contents)
	return

###########################################
# Read udf and setup initial conditions
def read_and_setcondition():
	dic={'y':True,'yes':True,'q':False,'quit':False}
	while True:
		# read udf
		basic_cond, nw_cond, sim_cond, rnd_cond = readconditionudf()
		# select
		condsetup = CondSetup(nw_cond, sim_cond)
		target_cond = condsetup.calc_conditions()
		print('Change UDF: type [r]eload')
		print('Quit input process: type [q]uit')
		inp = input('Condition is OK ==> [y]es >> ').lower()
		if inp in dic:
			inp = dic[inp]
			break
		print('##### \nRead Condition UDF again \n#####\n\n')
	if inp:
		print("\n\nSetting UP progress !!")
		return basic_cond, nw_cond, sim_cond, rnd_cond, target_cond
	else:
		sys.exit("##### \nQuit !!")

####################################
# Read condition udf
def readconditionudf():
	u = UDFManager('calc_condition.udf')
	u.jump(-1)
	##################
	# 使用するCognacのバージョン
	ver_cognac = u.get('CalcCond.Cognac_ver')
	# 計算に使用するコア数
	core = u.get('CalcCond.Cores')
	# ベースとするUDFの名前
	base_udf = "base_uin.udf"
	blank_udf = ver_cognac + '.udf'
	###########
	basic_cond = [ver_cognac, blank_udf, base_udf, core]
	#######################################################
	## 計算ターゲット

	###################
	## Networkモデルの設定
	nw_model = u.get('TargetCond.Model.TargetModel')
	###################
	## Networkモデルの設定
	restart = ''
	cond_top = []
	n_hist = 0
	#
	if nw_model == "Regular_NW":
		strand_type = u.get('TargetCond.Model.Regular_NW.chains')
	elif nw_model == "Random_NW":
		strand_type = u.get('TargetCond.Model.Random_NW.chains')
		calc = u.get('TargetCond.Model.Random_NW.Calc_Topolpgy')
		n_hist = u.get('TargetCond.Model.Random_NW.N_histgram')
		if calc == 'Read':
			restart = u.get('TargetCond.Model.Random_NW.Read.dir_name')
			if not os.path.exists(os.path.join(restart, 'init.pickle')):
				exit("##########\ntarget directory does not exists.")
			elif n_strand != int(restart.split('_')[0]):
				sys.exit("##########\nnumber of strands: selected n_strand is different from original Calculation.")
			elif n_cell != int(restart.split('_')[2]):
				sys.exit("##########\nnumber of cells: selected n_cell is different from original Calculation.")
		elif calc == 'Calc':
			cond_top = u.get('TargetCond.Model.Random_NW.Calc')
	################
	if strand_type == "3_Chain" or strand_type == "3_Chain_S" or strand_type == "3_Chain_D":
		n_strand = 3
	elif strand_type == "4_Chain":
		n_strand = 4
	elif strand_type == "5_Chain":
		n_strand = 5
	elif strand_type == "6_Chain":
		n_strand = 6
	elif strand_type == "7_Chain":
		n_strand = 7
	elif strand_type == "8_Chain":
		n_strand = 8
	###################
	## ポリマー鎖の設定
	n_segments = u.get('TargetCond.NetWork.N_Segments')
	n_sc = u.get('TargetCond.NetWork.N_Subchain')
	n_cell = u.get('TargetCond.NetWork.N_UnitCells')
	#####
	if u.get('TargetCond.Multiplisity.Set_or_Calc') == 'Set':
		multi = u.get('TargetCond.Multiplisity.Set.Multiplicity')
		density = 0
	elif u.get('TargetCond.Multiplisity.Set_or_Calc') == 'Calc':
		multi = 0
		density = u.get('TargetCond.Multiplisity.Calc.Density')
	#####
	if u.get('TargetCond.Shrinkage.Shrink') == 'Yes':
		if u.get('TargetCond.Shrinkage.Yes.Control') == 'Density':
			density = u.get('TargetCond.Shrinkage.Yes.Density.target_density')
			shrinkage = 0.
		elif u.get('TargetCond.Shrinkage.Yes.Control') == 'Shrink':
			shrinkage = u.get('TargetCond.Shrinkage.Yes.Shrinkage.value')
			density = 0
	elif u.get('TargetCond.Shrinkage.Shrink') == 'No':
		shrinkage = 0
		# density = -1
	#####
	topology = u.get('TargetCond.Type.Topology')
	if topology == 'Entangled':
		step_rfc = u.get('TargetCond.Type.Entangled.Step_rfc[]')
		expand = 1.0
		step_press = []
	elif topology == 'NO_Entangled':
		step_rfc = []
		expand = u.get('TargetCond.Type.NO_Entangled.ExpansionRatio')
		step_press = u.get('TargetCond.Type.NO_Entangled.StepPress[]')
	##########
	## シミュレーションの条件
	equilib_repeat = u.get('SimulationCond.Equilib_Condition.repeat')
	equilib_time = u.get('SimulationCond.Equilib_Condition.repeat')
	#####
	l_bond = u.get('SimulationCond.l_bond')
	#####
	if n_segments <= 10:
		c_n = 1.5
	elif n_segments <= 20:
		c_n = 1.65
	elif n_segments <= 40:
		c_n = 1.7
	else:
		c_n = 1.75
	#########################################################################################
	nw_cond = [nw_model, strand_type, n_strand, n_segments, n_cell, n_sc, l_bond, c_n]
	sim_cond = [topology, multi, density, shrinkage, expand, step_press, step_rfc, equilib_repeat, equilib_time]
	rnd_cond = [restart, cond_top, n_hist]

	return basic_cond, nw_cond, sim_cond, rnd_cond

######################################
##### Setup Calculation COnditions
######################################
class CondSetup:
	def __init__(self, nw_cond, sim_cond):
		self.nw_model = nw_cond[0]
		self.strand = nw_cond[1]
		self.n_strand = nw_cond[2]
		self.n_segments = nw_cond[3]
		self.n_cell = nw_cond[4]
		self.n_sc = nw_cond[5]
		self.l_bond = nw_cond[6]
		self.c_n = nw_cond[7]
		#
		self.topology = sim_cond[0]
		self.multi = sim_cond[1]
		self.density = sim_cond[2]
		self.shrinkage = sim_cond[3]
		self.expand = sim_cond[4]
		self.step_press = sim_cond[5]
		self.rfc = sim_cond[6]
		self.equilib_repeat = sim_cond[7]
		self.equilib_time = sim_cond[8]
		
	############################################
	##### ネットワークポリマーの諸量を計算 ######
	############################################
	#-----ネットワークポリマーの諸量を計算
	def calc_conditions(self):
		## 計算システムの諸量を計算して、出力
		e2e, n_chains, n_beads_unit, org_unitcell = self.set_length()
		target_cond = self.init_calc(e2e, n_chains, n_beads_unit, org_unitcell)
		return target_cond

	#####################
	#
	def set_length(self):
		e2e = self.l_bond*((self.n_segments + 1)*self.c_n)**0.5					# 理想鎖状態での末端間距離

		if self.nw_model == "Regular_NW":
			if self.strand == "3_Chain_S":
				n_chains = 12						        					# サブチェインの本数
				n_beads_unit = 8 + self.n_segments*(1 + self.n_sc)*n_chains		# ユニットセル当たりの粒子数
				org_unitcell = (2*2**0.5)*e2e				        			# 理想鎖状態でのユニットセル長
			elif self.strand == "3_Chain_D":
				n_chains = 24						       
				n_beads_unit = 16 + self.n_segments*(1 + self.n_sc)*n_chains	
				org_unitcell = (2*2**0.5)*e2e		
			elif self.strand == "4_Chain":
				n_chains = 16						      
				n_beads_unit = 8 + self.n_segments*(1 + self.n_sc)*n_chains	
				org_unitcell = (4*3**0.5)*e2e/3			 
			elif self.strand == "6_Chain":
				n_chains = 3						  
				n_beads_unit = 1 + self.n_segments*(1 + self.n_sc)*n_chains		
				org_unitcell = e2e						      
			elif self.strand == "8_Chain":
				n_chains = 8						   
				n_beads_unit = 2 + self.n_segments*(1 + self.n_sc)*n_chains   
				org_unitcell = (2*3**0.5)*e2e/3	

		elif self.nw_model == "Random_NW":
			n_chains = self.n_strand
			n_beads_unit = 2 + self.n_segments*(1 + self.n_sc)*n_chains
			org_unitcell = (2*3**0.5)*e2e/3	

		return e2e, n_chains, n_beads_unit, org_unitcell

	###############################################################
	def init_calc(self, e2e, n_chains, n_beads_unit, org_unitcell):
		calcd_density = 0
		if self.multi == 0:
			if self.shrinkage == 0.:
				# org_system = org_unitcell*self.n_cell							# e2e から決めたシステムサイズ
				calcd_multi = round(self.density*org_unitcell**3/n_beads_unit)	# 密度を設定値とした場合に必要な多重度
				calcd_density = n_beads_unit*calcd_multi/org_unitcell**3		# 上記の多重度での密度
				err_dens = round((calcd_density/self.density - 1)*100, 2) 		# 設定密度との誤差(%)
				if abs(err_dens) > 1:
					print(u"\n##### \n圧縮後の密度が、設定したい密度と 1% 以上違います。\nそれでも計算しますか？\n#####\n")
				single_net_atom = int(n_beads_unit*self.n_cell**3.)	    		# 一つのネットワーク中の粒子数
				total_net_atom = int(calcd_multi*single_net_atom)    			# 全システム中のネットワーク粒子数
				mod_unitcell = (calcd_multi*n_beads_unit/self.density)**(1/3)					# その際のシステムサイズ
				self.shrinkage = mod_unitcell/org_unitcell								# 収縮比
				system = mod_unitcell*self.n_cell
				vol = system**3									    			# システム体積
				print(self.n_cell)
				self.multi = calcd_multi
				# self.density = calcd_density
				
				mod_e2e = self.shrinkage*e2e											# 収縮後の末端間距離
				unit_cell = mod_unitcell
			elif self.shrinkage > 0.:
				sys.exit(u"\n##### \n多重度を自動計算にした場合、収縮条件は選択できません\n#####\n")
		elif self.multi != 0:
			if self.shrinkage == 0.:
				if self.density == 0.:
					err_dens = 0.
					system = org_unitcell*self.n_cell						# e2e から決めたシステムサイズ
					self.density = n_beads_unit*self.multi/org_unitcell**3	# 多重度での密度
					single_net_atom = int(n_beads_unit*self.n_cell**3.)	    # 一つのネットワーク中の粒子数
					total_net_atom = int(self.multi*single_net_atom)    	# 全システム中のネットワーク粒子数
					vol = system**3									    	# システム体積
					mod_e2e = e2e											# 収縮後の末端間距離
					unit_cell = org_unitcell
				elif self.density > 0.:
					sys.exit(u"\n##### \nSomething wrong!!\n#####\n")
			elif self.shrinkage != 0.:
				if self.density == 0.:
					err_dens = 0.
					mod_unitcell = org_unitcell*self.shrinkage
					system = mod_unitcell*self.n_cell						# e2e から決めたシステムサイズ
					self.density = n_beads_unit*self.multi/mod_unitcell**3	# 多重度での密度
					single_net_atom = int(n_beads_unit*self.n_cell**3.)	    # 一つのネットワーク中の粒子数
					total_net_atom = int(self.multi*single_net_atom)    	# 全システム中のネットワーク粒子数
					vol = system**3									    	# システム体積
					mod_e2e = self.shrinkage*e2e		
					unit_cell = mod_unitcell									# 収縮後の末端間距離
				elif self.density > 0.:
					err_dens = 0.
					mod_unitcell = (n_beads_unit*self.multi/self.density)**(1/3)
					self.shrinkage = mod_unitcell/org_unitcell
					single_net_atom = int(n_beads_unit*self.n_cell**3.)	    # 一つのネットワーク中の粒子数
					total_net_atom = int(self.multi*single_net_atom)    	# 全システム中のネットワーク粒子数
					system = mod_unitcell*self.n_cell						# e2e から決めたシステムサイズ
					vol = system**3									    	# システム体積
					mod_e2e = self.shrinkage*e2e											# 収縮後の末端間距離
					unit_cell = mod_unitcell
		else:
			sys.exit("Something Wrong!!")
		nu = n_chains/vol
		#
		text = "#########################################" + "\n"
		text += "ネットワークトポロジー\t\t" + str(self.nw_model) + "\n"
		text += "ネットワークモデル\t\t" + str(self.strand) + "\n"
		text += "#########################################" + "\n"
		text += "ストランド中のセグメント数\t" + str(self.n_segments) + "\n"
		text += "特性比:\t\t\t\t" + str(round(self.c_n, 2)) + "\n"
		text += "初期の末端間距離:\t\t" + str(round(e2e, 4)) + "\n"
		text += "#########################################" + "\n"
		text += "当初の単位ユニット:\t\t" + str(round(org_unitcell, 4)) + "\n"
		text += "一辺当たりの単位ユニット数\t" + str(self.n_cell) + "\n"
		# text += "当初のシステムサイズ:\t\t" + str(round(org_system, 4)) + "\n"
		text += "#########################################" + "\n"
		text += "ネットワークタイプ\t\t" + str(self.topology) + "\n"
		text += "#########################################" + "\n"
		text += "設定密度:\t\t\t" + str(self.density) + "\n"
		text += "多重度:\t\t\t\t" + str(self.multi) + "\n"
		text += "NW の全セグメント数:\t\t" + str(total_net_atom) + "\n"
		text += "収縮比:\t\t\t\t" + str(round(self.shrinkage, 4)) + "\n"
		text += "任意の多重度での密度:\t\t" + str(round(calcd_density, 4)) + "\n"
		text += "密度の誤差:\t\t\t" + str(round(err_dens, 3)) + "%" + "\n"
		# text += "収縮後の末端間距離:\t\t" + str(round(mod_e2e, 4)) + "\n"
		text += "システムサイズ:\t\t\t" + str(round(system, 4)) + "\n"
		text += "#########################################" + "\n"
		if self.topology == 'Entangled':
			text += "Slow Push Off での rfc 条件:\t" + ', '.join(map(str, self.rfc)) + "\n"
		else:
			text += "NPT 計算時の初期膨張率:\t\t" + str(self.expand) + "\n"
			text += "ステップ圧力:\t" + ', '.join(map(str, self.step_press)) + "\n"
		text += "#########################################" + "\n"
		text += "ストランドの数密度:\t\t" + str(round(nu, 5)) + "\n"
		text += "#########################################" + "\n"
		print(text)

		# if (topology == "Entangled" or topology == "NPT") and flag == 1:
		# 	print(u"\n##### \n圧縮後の密度が、設定したい密度と 1% 以上違います。\nそれでも計算しますか？\n#####\n")
		#
		with open("calc_conditions.txt", 'w') as f:
			f.write(text)
		#
		target_cond = [self.multi, system, unit_cell, total_net_atom, nu, self.nw_model, e2e, mod_e2e, self.shrinkage, err_dens]

		return target_cond
	
