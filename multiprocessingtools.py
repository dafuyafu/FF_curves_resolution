'''
	https://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic
	recently accessed: 20200923
'''

import multiprocessing
import multiprocessing.pool

class NoDaemonProcess(multiprocessing.Process):
	def _get_daemon(self):
		return False
	def _set_daemon(self, value):
		pass
	daemon = property(_get_daemon, _set_daemon)

class SFFPool(multiprocessing.pool.Pool):
	Process = NoDaemonProcess

# def test():
#     print("Creating 5 (non-daemon) workers and jobs in main process.")
#     pool = MyPool(5)

#     result = pool.map(work, [randint(1, 5) for x in range(5)])

#     pool.close()
#     pool.join()
#     print(result)

# if __name__ == '__main__':
#     test()

class StopEval(Exception):
	pass