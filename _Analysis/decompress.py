# taken from: http://pythonhug.blogspot.ca/2012/09/decompress-file-using-gzip.html

# used to decompress eeg files in my case

import gzip
import os

for file in os.listdir():
	if file.endswith(".eeg.gz"):
		
		sourceFile = file
		
		if (os.path.exists(sourceFile)):
			pass
		else:
			print ("Sorry, cant find the file! Please check if the filename and path are correct!\n")
			exit(0)
		destFile = file[:-3]

		zipFile = gzip.open(sourceFile,"rb")

		unCompressedFile = open(destFile,"wb")

		decoded = zipFile.read()

		unCompressedFile.write(decoded)

		zipFile.close()

		unCompressedFile.close()