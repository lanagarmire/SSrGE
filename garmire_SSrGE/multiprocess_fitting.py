from multiprocessing import Queue
from multiprocessing import Process

from contextlib import contextmanager
import signal

import numpy as np

from time import sleep
from sys import stdout

from garmire_SSrGE.config import MIN_OBS_FOR_REGRESS
from garmire_SSrGE.config import TIME_LIMIT

from collections import Counter

from numpy import hstack
from scipy.sparse import hstack as shstack
from scipy.sparse import issparse

import warnings

from time import time


class TimeoutException(Exception): pass


@contextmanager
def time_limit(seconds):
    warnings.catch_warnings()
    warnings.simplefilter("ignore")

    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")

    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)

    try:
        yield
    finally:
        signal.alarm(0)


def debug():
    """
    #### DEBUG ####
    **** Test function ****
    """
    from garmire_SSrGE.examples import create_example_matrix_v1
    from sklearn.linear_model import Lasso


    X, Y, W = create_example_matrix_v1()

    multi_test = BatchFitting(I_mat=X,
                              O_mat=Y,
                              model=Lasso,
                              model_params={'alpha': 0.01},
                              nb_processes=1,
                              only_nonzero=False,
                              min_obs_for_regress=0,
                              cis_model=None)
    g_index, coefs, intercepts = multi_test.run()

    return g_index, coefs, intercepts

class MultiProcessFitting(Process):
    def __init__(self,
                 input_queue,
                 output_queue,
                 model,
                 model_params,
                 matrix,
                 process_id,
                 time_limit=TIME_LIMIT,
                 min_obs_for_regress=MIN_OBS_FOR_REGRESS,
                 only_nonzero=False,
                 cis_model=None):
        """ """
        Process.__init__(self)
        self.input_queue = input_queue
        self.output_queue = output_queue
        self.model = model
        self.model_params = model_params
        self.matrix = matrix
        self.process_id = process_id
        self.only_nonzero = only_nonzero
        self.cis_model = cis_model
        self.time_limit = time_limit
        self.min_obs_for_regress = min_obs_for_regress

    def run(self):
        """ """
        model = self.model(**self.model_params)

        while not self.input_queue.empty():
            try:
                gene_i, y, data = self.input_queue.get(True, 0.1)
            except Exception as e:
                continue

            if self.only_nonzero:
                matrix, y = self._clean_matrix(y)
            else:
                matrix = self.matrix

            if self.cis_model:
                matrix = self._matrix_to_cis_model(matrix, gene_i)

            if y.shape[0] > self.min_obs_for_regress and \
               not isinstance(matrix, type(None)):
                try:
                    with time_limit(self.time_limit):
                        model.fit(X=matrix, y=y, **data)
                except Exception as e:
                    intercept = np.nan
                    coefs = np.empty(self.matrix.shape[1])
                    coefs[:] = np.nan
                    print('\n exception found for linear model:{0}\n skipping'\
                        .format(e))
                else:
                    if self.cis_model:
                        coefs = np.zeros(self.matrix.shape[1])
                        coefs[self.cis_model[gene_i]] = model.coef_
                    else:
                        coefs = model.coef_

                    intercept = model.intercept_

            else:
                intercept = np.nan
                coefs = np.empty(self.matrix.shape[1])
                coefs[:] = np.nan

            coefs = Counter({i:np.abs(coefs[i])
                             for i in np.nonzero(np.nan_to_num(coefs))[0]})

            while True:
                try:
                    self.output_queue.put((gene_i, coefs, intercept), timeout=0.1)
                except Exception as e:
                    continue
                else:
                    break

    def _matrix_to_cis_model(self, matrix, gene_i):
        """ """
        if not self.cis_model[gene_i]:
            return None
        return matrix.T[self.cis_model[gene_i]].T

    def _clean_matrix(self, y):
        """ """
        index = np.nonzero(y)[0]
        return self.matrix[index], y[index]


class BatchFitting():
    """ """
    def __init__(
            self,
            I_mat,
            O_mat,
            model,
            model_params,
            nb_processes=1,
            time_limit=TIME_LIMIT,
            min_obs_for_regress=MIN_OBS_FOR_REGRESS,
            add_y_index=False,
            only_nonzero=False,
            cis_model=None,
            ):
        self.I_mat = I_mat
        self.O_mat = O_mat
        self.model = model
        self.model_params = model_params
        self.nb_processes = nb_processes
        self.add_y_index = add_y_index
        self.only_nonzero = only_nonzero
        self.cis_model = cis_model
        self.time_limit = time_limit
        self.min_obs_for_regress = min_obs_for_regress

    def run(self):
        """ run batch fitting """

        res = self._run()

        if isinstance(res, Exception):
            raise res

        return res

    def _kill_processes(self):
        """ """
        for process in self.processes_list:
            process.terminate()

    def _get_qsize(self, output_queue):
        """ """
        while True:
            try:
                with time_limit(self.time_limit):
                    out_qsize = output_queue.qsize()

                break
            except Exception as e:
                print('exception was found for qsize:', e)
                continue

        return out_qsize

    def _run(self):
        """custom unordered multiprocessing"""
        input_queue = Queue()
        output_queue = Queue()
        self.processes_list = []
        res_list = []
        qsize = self.O_mat.shape[0]
        i = 0

        for y in self.O_mat:
            data = {}
            if self.add_y_index:
                data['y_index'] = i
            input_queue.put((i, y, data))
            i += 1

        for i in range(self.nb_processes):
            self.processes_list.append(
                MultiProcessFitting(
                    input_queue=input_queue,
                    output_queue=output_queue,
                    model=self.model,
                    model_params=self.model_params,
                    matrix=self.I_mat,
                    process_id=i,
                    time_limit=self.time_limit,
                    only_nonzero=self.only_nonzero,
                    min_obs_for_regress=self.min_obs_for_regress,
                    cis_model=self.cis_model)
            )

        for process in self.processes_list:
            process.start()

        terminate = False

        j = 0
        prog = ['/', '-', '\\', '|']

        start_time = time()

        while True:
            for process in self.processes_list:
                if process.exitcode:
                    print('error with process with id: {0} terminating'\
                        .format(process.process_id))
                    terminate = True
                    break

            if terminate:
                break

            out_qsize = self._get_qsize(output_queue)

            stdout.write('\r{0} / {1} models done {2} (time: {3})'\
                         .format(out_qsize, qsize, prog[j], time() - time() - start_time))
            stdout.flush()

            j += 1

            if j == 4:
                j = 0

            if out_qsize >= qsize:
                break

            sleep(0.5)

        if terminate:
            print 'one of the process raised an exception'\
                '\n killing process...'
            self._kill_processes()

            return Exception('process not finished correctly!')

        print '\n'

        for i in range(qsize):
            res_list.append(output_queue.get())
            stdout.write('\r{0} / {1} results loaded'\
                         .format(i + 1, qsize))
            stdout.flush()

        del output_queue
        del input_queue

        self._kill_processes()

        return zip(*res_list)


if __name__ == "__main__":
    debug()
