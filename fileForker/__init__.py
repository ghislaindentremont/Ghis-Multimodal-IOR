import billiard
billiard.forking_enable(0)

########
# Define a class that spawns a new process
########
class childClass:
	def __init__(self,childFile):
		self.childFile = childFile
		self.initDict = {}
		self.qTo = billiard.Queue()
		self.qFrom = billiard.Queue()
		self.started = False
	def f(self,childFile,qTo,qFrom,initDict):
		execfile(childFile)
		import sys
		sys.exit()
	def start(self):
		if self.started:
			print 'Oops! Already started this child.'
		else:
			self.process = billiard.Process( target=self.f , args=(self.childFile,self.qTo,self.qFrom,self.initDict,) )
			self.process.start()
			self.started = True
	def isAlive(self):
		return self.process.is_alive()
	def stop(self,killAfter=None):
		if not self.started:
			print 'Oops! Not started yet!'
		else:
			self.qTo.put('quit')
			self.process.join(timeout=killAfter)
			if self.process.is_alive():
				self.process.terminate()
		return None
