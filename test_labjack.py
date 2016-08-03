# test labjack 
import u3
import time

labjack = u3.U3()
labjack.configU3()
labjack.getFeedback(u3.LED(State = False))

print("""
CURRENT STIMULUS MARKERS:
	trigger_led_num = 8
	left_led_num = 11
	right_led_num = 9
	left_tact_num = 15
	right_tact_num = 13\n""")

stim_num_str = raw_input("What is the number of the tactamp stimulus you want to turn on along with the trigger led?\n>>> ")
stim_num = int(stim_num_str)

trigger_led_num = 8
labjack_to_tactamp_target_on_bits = [0,0,0,0,0,0,0,0]
labjack_to_tactamp_target_on_bits[trigger_led_num-8] = 1
labjack_to_tactamp_target_on_bits[stim_num-8] = 1

labjack_to_tactamp_target_on_bits_int = int(''.join(map(str,labjack_to_tactamp_target_on_bits[::-1])),2)
print("\nInteger Bit Value: "+str(labjack_to_tactamp_target_on_bits_int))

print("\nTesting Stimulus "+stim_num_str+"...")
for i in range(0,10):
	print("Cycle "+str(i))
	labjack.getFeedback(u3.PortStateWrite(State = [0,0,0]))
	time.sleep(5)
	labjack.getFeedback(u3.PortStateWrite(State = [0,labjack_to_tactamp_target_on_bits_int,0]))
	time.sleep(5)

labjack.getFeedback(u3.PortStateWrite(State = [0,0,0]))
labjack.close()
