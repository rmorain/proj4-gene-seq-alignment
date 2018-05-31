#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))




import time



class GeneSequencing:
    def __init__( self ):
        pass

    def align_all( self, sequences, banded, align_length ):
        #print(banded)
        #print(align_length)
        results = []
        for i in range(len(sequences)):
            jresults = []
            for j in range(len(sequences)):
                s = {'align_cost':i+j,
                     'seqi_first100':'abc-easy  DEBUG:(seq{}, {} chars,align_len={}{})'.format(i+1,
                         len(sequences[i]), align_length, ',BANDED' if banded else ''),
                     'seqj_first100':'as-123--  DEBUG:(seq{}, {} chars,align_len={}{})'.format(j+1,
                         len(sequences[j]), align_length, ',BANDED' if banded else '')}
                jresults.append(s)
            results.append(jresults)
        return results


