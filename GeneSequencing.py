#!/usr/bin/python3

from which_pyqt import PYQT_VER
from enum import Enum
import math

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))




import time



class GeneSequencing:


    def __init__( self ):
        # Cost of each alignment action
        self.SUB_COST = 1
        self.INDEL_COST = 5
        self.MATCH_COST = -3
        self.BANDED_D = 3
        pass

    # def align_all( self, sequences, banded, align_length ):
    #     #print(banded)
    #     #print(align_length)
    #     results = []
    #     for i in range(len(sequences)):
    #         jresults = []
    #         for j in range(len(sequences)):
    #             s = {'align_cost':i+j,
    #                  'seqi_first100':'abc-easy  DEBUG:(seq{}, {} chars,align_len={}{})'.format(i+1,
    #                      len(sequences[i]), align_length, ',BANDED' if banded else ''),
    #                  'seqj_first100':'as-123--  DEBUG:(seq{}, {} chars,align_len={}{})'.format(j+1,
    #                      len(sequences[j]), align_length, ',BANDED' if banded else '')}
    #             jresults.append(s)
    #         results.append(jresults)
    #     return results


    def align_all(self, sequences, banded, align_length):
        """
        Align all the sequences
        :param sequences: List of strings
        :param banded: Boolean - indicates if we should use the banded input
        :param align_length: The k value; Integer. How many characters the algorithm should align
        :return: A list of jresults
                Each jresult is a dictionary that has an 'align_cost': int, 'seqi_first100': 'string'
                and 'seqj_first100': 'string'
        """

        # Create results list
        results = []
        count = 1
        # Align each sequence with every other sequence
        for i in range(len(sequences)):
            jresults = []
            for j in range(len(sequences)):
                # Check if we are going to index two different sequences
                if i <= j:
                    # align_seq() returns a prepared dictionary of align cost and character alignment strings
                    # print("Count: {}" .format(count))
                    count += 1
                    jresults.append(self.align_sequence(sequences[i][:align_length], sequences[j][:align_length], banded))
                else:
                    jresults.append({'align_cost': -1, 'seqi_first100': 'nothing', 'seqj_first100': 'nothing'})
            results.append(jresults)
        # Return the results
        # print("done!")
        return results


    def align_sequence(self, s1, s2, banded):
        """
        Do an unrestricted alignment of the sequences
        :param s1: A string. Represents one sequence
        :param s2: A string. Represents a different sequence
        :param align_length: INT How many characters we align
        :return: A dictionary with the parameters:
                    'align_cost': int, 'seqi_first100': 'string', 'seqj_first100': 'string'
        """
        # Set up the grid
        grid = self.get_grid(s1, s2)
        # Fill in the grid
        self.fill_grid(s1, s2, grid, banded)
        # Minimum cost is in the bottom right corner
        align_cost = grid[len(s1)][len(s2)]['align_cost']
        buffers = self.extract_first100(s1, s2, grid)
        # print('Sequence 1: {}' .format(buffers['ibuffer']))
        # print('Sequence 2: {}' .format(buffers['jbuffer']))
        return {'align_cost': align_cost, 'seqi_first100': buffers['ibuffer'], 'seqj_first100': buffers['jbuffer']}


    def get_grid(self, s1, s2):
        """
        Initialize a grid to use for the sequence alignment
        :param s1: A string. Represents one sequence
        :param s2: A string. Represents a different sequence
        :param align_length: INT How many characters we align
        :return: 2D list
        """
        # Create a list
        grid = []
        for i in range(len(s1)+1):
            grid.append([])     # Add a row to the grid
            for j in range(len(s2)+1):
                grid[i].append({'align_cost': None, 'arrow': None})    # Add each cell
        return grid

    def fill_grid(self, s1, s2, grid, banded):

        # Start with a zero in the first cell
        grid[0][0]['align_cost'] = 0
        if not banded:
            # Fill in the first row
            for j in range(1, len(s2)+1):
                grid[0][j]['align_cost'] = self.INDEL_COST*j
            # Fill in the first col
            for i in range(1, len(s1)+1):
                grid[i][0]['align_cost'] = self.INDEL_COST*i
        else:
            # Fill in the first row
            for j in range(1, self.BANDED_D + 1):
                grid[0][j]['align_cost'] = self.INDEL_COST * j
            # Fill in the first col
            for i in range(1, self.BANDED_D + 1):
                grid[i][0]['align_cost'] = self.INDEL_COST * i


        # Fill in the rest of the grid
        # For each row
        for i in range(1, len(s1)+1):
            # For each cell in the row
            for j in range(1, len(s2)+1):
                if banded and abs(i - j) > self.BANDED_D:
                    pass
                else:
                    # Get the candidate cost of each neighboring cell
                    if grid[i-1][j]['align_cost'] is not None:
                        top_cost = grid[i-1][j]['align_cost'] + self.INDEL_COST
                    else:
                        top_cost = math.inf
                    if grid[i][j-1]['align_cost'] is not None:
                        left_cost = grid[i][j-1]['align_cost'] + self.INDEL_COST
                    else:
                        left_cost = math.inf
                    top_left_cost = grid[i-1][j-1]['align_cost']

                    # Check if we have a match or a sub
                    if s1[i-1] == s2[j-1]:
                        top_left_cost += self.MATCH_COST
                    else:
                        top_left_cost += self.SUB_COST

                    # Assign the minimum to the grid
                    costs = [top_cost, left_cost, top_left_cost]
                    min_cost = math.inf
                    min_index = -1
                    for x in range(3):
                        if costs[x] < min_cost:
                            min_cost = costs[x]
                            min_index = x

                    grid[i][j] = {'align_cost': min_cost, 'arrow': min_index}

    def extract_first100(self, s1, s2, grid):
        # Start at the bottom corner
        i = len(grid) - 1
        j = len(grid[0]) - 1
        ibuffer = []
        jbuffer = []
        while True:
            if i == 0 or j == 0:
                # If the last characters are aligned just add them to the buffer
                while i != j:
                    if i > j:
                        ibuffer.insert(0, s1[i-1])
                        jbuffer.insert(0, '-')
                        i -= 1
                    else:
                        ibuffer.insert(0, '-')
                        jbuffer.insert(0, s2[i-1])
                        j -= 1
                break

            # Check if we should add a character
            if grid[i][j]['arrow'] == 2:
                ibuffer.insert(0, s1[i-1])
                jbuffer.insert(0, s2[j-1])
                i -= 1
                j -= 1
            # Otherwise add an indel
            else:
                ibuffer.insert(0, '-')
                jbuffer.insert(0, '-')
                if grid[i][j]['arrow'] == 0:
                    i -= 1
                else:
                    j -= 1

        return {'ibuffer': ''.join(ibuffer[:100]), 'jbuffer': ''.join(jbuffer[:100])}

