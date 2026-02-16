import os, sys

class ProgressCounter(object):
	def __init__(self, **kwargs):
		super(ProgressCounter, self).__init__()
		self.iter = 0
	def __call__(self, file):
		if str(file).find('area') != -1:
			self.iter = self.iter + 1
			print "================= {} ======================".format(file)
			print "Iteration number: {}".format(self.iter)
		if str(file).find('log') != -1:
			areafile = "%s.%s.area" % (str(file).split('.')[0], str(file).split('.')[1])
			if os.path.isfile(areafile):
				total=self.file_len(areafile)-1
				done=self.search_pattern(areafile, '-1')
				progress = 100*done/total
				print "Mesh quality progress: %d%%" % (progress)
				if progress == 100:
					print "Exiting..."
					print "Reached 100% of progress!"
					os._exit(0)

	# get total lines in file
	def file_len(self, fname):
		with open(fname) as f:
			data = f.readlines()
		return len(data)
	# count pattern in file
	def search_pattern(self, filename, pattern):
		i=0
		with open(filename, 'r') as f:
			content = f.read().splitlines()
		for line in content:
			if line.find('-1') != -1:
				i=i+1
		return i