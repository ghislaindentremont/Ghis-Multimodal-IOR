
# small script to prepend necessary information at beginning of vmrk file 

import os

for file in os.listdir():
    if file.endswith(".vmrk"):
        with open(file, "r+") as f:
		     old = f.read() # read everything in the file
		     f.seek(0) # rewind
		     f.write("Brain Vision Data Exchange Marker File, Version 1.0\n\n" + old) # write the new line before


