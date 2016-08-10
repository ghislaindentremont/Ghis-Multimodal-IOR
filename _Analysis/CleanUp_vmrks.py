import fnmatch
import os

filedir = input("What is the vmrk directory?\n>>> ")
filename = input("What is the vmrk file name?\n>>> ")

os.chdir(filedir)  # change to file's folder
newFile = open("clean_" + filename, "w")  # open file in which error-free vmrk info will be put in
oldFile = open(filename, 'r')
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
                    if number_string in ["99","98"]:
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
