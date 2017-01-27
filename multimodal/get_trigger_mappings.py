if __name__ == '__main__':
 	
	cue_modality_list = ['visual','tactile']
	cue_location_list = ['left','right']
	target_location_list = ['left','right']
	target_modality_list = ['visual','tactile']
	target_type_list = ['catch','target']

	########
	# Import libraries
	########
	import u3 #for labjack
	import time

	########
	# Initialize the labjack
	########
	labjack = u3.U3()
	labjack.configU3()
	labjack.getFeedback(u3.LED(State = False))

	def get_trials():
		trials = []
		for cue_modality in cue_modality_list:
			for cue_location in cue_location_list:
				for target_location in target_location_list:
					for target_modality in target_modality_list:
						for target_type in target_type_list:
							trials.append([cue_modality,cue_location,target_location,target_modality,target_type])
		return trials

	trial_list = get_trials()
	
	while len(trial_list)>0:
		cue_modality , cue_location , target_location, target_modality , target_type = trial_list.pop()

		if cue_location == 'left':
			if cue_modality == 'visual':
				labjack_to_eeg_cue_int = 10
				if target_location=='left':
					if target_modality=='visual':
						labjack_to_eeg_target_int = 20
					else:
						labjack_to_eeg_target_int = 21
				else:
					if target_modality=='visual':
						labjack_to_eeg_target_int = 22
					else:
						labjack_to_eeg_target_int = 23
			else:
				labjack_to_eeg_cue_int = 11
				if target_location=='left':
					if target_modality=='visual':
						labjack_to_eeg_target_int = 24
					else:
						labjack_to_eeg_target_int = 25
				else:
					if target_modality=='visual':
						labjack_to_eeg_target_int = 26
					else:
						labjack_to_eeg_target_int = 27
		else:
			if cue_modality == 'visual':
				labjack_to_eeg_cue_int = 12
				if target_location=='left':
					if target_modality=='visual':
						labjack_to_eeg_target_int = 28
					else:
						labjack_to_eeg_target_int = 29
				else:
					if target_modality=='visual':
						labjack_to_eeg_target_int = 30
					else:
						labjack_to_eeg_target_int = 31
			else:
				labjack_to_eeg_cue_int = 13
				if target_location=='left':
					if target_modality=='visual':
						labjack_to_eeg_target_int = 32
					else:
						labjack_to_eeg_target_int = 33
				else:
					if target_modality=='visual':
						labjack_to_eeg_target_int = 34
					else:
						labjack_to_eeg_target_int = 35
		if target_type=='catch':
			labjack_to_eeg_target_int = 36

		print [labjack_to_eeg_cue_int,labjack_to_eeg_target_int]
		#make sure all the labjack outputs are off
		labjack.getFeedback(u3.PortStateWrite(State = [0,0,0]))
		time.sleep(.2)
		#send trial info to labjack
		labjack.getFeedback(u3.PortStateWrite(State = [1,0,0]))
		time.sleep(.2)
		labjack.getFeedback(u3.PortStateWrite(State = [0,0,0]))
		time.sleep(.2)
		labjack.getFeedback(u3.PortStateWrite(State = [labjack_to_eeg_cue_int,0,0]))
		time.sleep(.2)
		labjack.getFeedback(u3.PortStateWrite(State = [0,0,0]))
		time.sleep(.2)
		labjack.getFeedback(u3.PortStateWrite(State = [labjack_to_eeg_target_int,0,0]))
		time.sleep(.2)

	labjack.getFeedback(u3.PortStateWrite(State = [0,0,0]))
	time.sleep(.2)
	labjack.getFeedback(u3.PortStateWrite(State = [2,0,0]))
	time.sleep(.2)
	labjack.getFeedback(u3.PortStateWrite(State = [0,0,0]))
	time.sleep(.2)
	labjack.getFeedback(u3.PortStateWrite(State = [99,0,0]))
	time.sleep(.2)
	labjack.getFeedback(u3.PortStateWrite(State = [0,0,0]))
	time.sleep(.2)
	labjack.getFeedback(u3.PortStateWrite(State = [40,0,0]))
	time.sleep(.2)
	labjack.getFeedback(u3.PortStateWrite(State = [0,0,0]))
	time.sleep(.2)
	labjack.getFeedback(u3.PortStateWrite(State = [42,0,0]))
	time.sleep(.2)
	labjack.getFeedback(u3.PortStateWrite(State = [0,0,0]))
	time.sleep(.2)
	labjack.getFeedback(u3.PortStateWrite(State = [43,0,0]))
	time.sleep(.2)
	labjack.getFeedback(u3.PortStateWrite(State = [0,0,0]))
	time.sleep(.2)
	labjack.getFeedback(u3.PortStateWrite(State = [2,0,0]))
	time.sleep(.2)
	labjack.getFeedback(u3.PortStateWrite(State = [0,0,0]))
	time.sleep(.2)
	labjack.getFeedback(u3.PortStateWrite(State = [99,0,0]))
