from multiprocessing import Process, Queue

class SFFProcess():
	def __init__(self, target, args, processes=8):
		self.queue = Queue()
		self.p_list = [Process(target=target, args=args+(d, self.queue, processes)) for d in range(processes)]

	def start(self):
		for p in self.p_list:
			p.start()

	def terminate(self):
		for p in self.p_list:
			if p.is_alive():
				p.terminate()

	def join(self):
		for p in self.p_list:
			p.join()

	def close(self):
		self.join()
		for p in self.p_list:
			p.close()

def sffprocess(target, args, processes=8):
	return SFFProcess(target, args, processes)
