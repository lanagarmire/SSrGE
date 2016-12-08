""" generic class to perform multi-processing """

from multiprocessing import Process
from multiprocessing import Queue

from time import sleep


class MPI():
    """generic multiprocessing class"""
    def __init__(self, input_list,
                 ProcessClass,
                 nb_processes=1,
                 verbose=True):
        """
        id_list: list    should be a list of input
        ProcessClass: class with process method and id attribute
        """
        self.input_queue = Queue()
        self.processes = []
        self.verbose = verbose
        self.nb_processes = nb_processes
        self.ProcessClass = ProcessClass

        for inpt in input_list:
            self.input_queue.put(inpt)

        for i in range(nb_processes):
            self.processes.append(
                MultiprocessingInstance(
                    input_queue=self.input_queue,
                    ProcessClass=ProcessClass,
                    id=i)
            )

    def _run(self):
        for p in self.processes:
            p.start()

        while self.input_queue.qsize():
            for p in self.processes:
                if p.exitcode:
                    raise KeyboardInterrupt
            sleep(1)

    def run(self):
        if self.verbose:
            rep = raw_input(
                'launching {0} processes with class {1} continue? (Y/n)'\
                .format(self.nb_processes, self.ProcessClass))
            if rep != 'Y':
                return

        try:
            self._run()

        except KeyboardInterrupt:
            for p in self.processes:
                p.terminate()

class MultiprocessingInstance(Process):
    """
    generic multiprocessing class
    """
    def __init__(self, input_queue, ProcessClass, id):
        """
        input_queue: Multiprocessing.Queue
        ProcessClass: class with process method and id attribute
        """
        self.input_queue = input_queue
        self.id = id
        self.process_instance = ProcessClass(id=id)
        Process.__init__(self)

    def run(self):
        while self.input_queue.qsize():
            try:
                sample = self.input_queue.get(True, 0.2)
            except Exception as e:
                print "exception:{0}".format(e)
                continue
            else:
                print "sample for sample {0} with id {1}"\
                    .format(sample, self.id)
                self.process_instance.process(sample)
