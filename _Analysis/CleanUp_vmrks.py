import os

# filedir = input("What is the vmrk directory?\n>>> ")
# filedir = "/Volumes/Seagate Backup Plus Drive/Experiments/m'ultimodal_ior/_Data/forMNE/new_data"
filedir = "/Users/ghislaindentremont/Documents/Multimodal_IOR/Ghis_EEG"

# filename = input("What is the vmrk file name?\n>>> ")
participant = "e45"
filename = "multimodal_ior_%s.vmrk" % participant

os.chdir(filedir)  # change to file's folder

# rename original it hasn't already been 'cleaned up'
already_done = os.path.isfile("original_vmrks/original_" + filename)

if not already_done:
    os.rename(filename, "original_vmrks/original_" + filename)

    newFile = open(filename, "w")  # open file in which error-free vmrk info will be put in
    oldFile = open("original_vmrks/original_" + filename, 'r')
    done = False
    dataStarted = False

    def get_number_string(lyne):
        if line == '':
            number_str = "1"  # to get out of bottom while loop when file complete
        else:
            number_str = lyne.replace(' ', '').split(',')[1].strip().strip('S').strip('R')  # this takes the number
        return(number_str)

    while not done:
        if not dataStarted:
            line = oldFile.readline()
            if line[0:3] == "Mk2":
                dataStarted = True
                number_string = get_number_string(line)
            else:
                newFile.write(line)
        if dataStarted:
            if line == '':
                done = True
            else:
                if number_string == "1":
                    lines = line
                    line = oldFile.readline()
                    number_string = get_number_string(line)
                    while number_string != "1":
                        if number_string in ["99","98","97","43","42"]:
                            break
                        else:
                            lines += line
                        line = oldFile.readline()
                        number_string = get_number_string(line)
                    if number_string == "1":
                        newFile.write(lines)
                    else:
                        while number_string != "1":
                            line = oldFile.readline()
                            number_string = get_number_string(line)
                else:
                    print("ERROR")

    newFile.close()
    oldFile.close()
else:
    print("I HAVE ALREADY 'CLEANED UP' THIS FILE")