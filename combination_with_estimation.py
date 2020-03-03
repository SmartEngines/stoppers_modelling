from combination import Cell, cell_dist, merge_cells, levmetric_ocr
from treap import Treap

class TreapBasedSequenceStructure:
    '''
    Data structure for sequence of y_1jk, y_2jk, ..., y_njk based on balanced binary search trees
    '''
    def __init__(self):
        self.treap = Treap()

    def get_sum(self):
        '''
        Returns a sum of inserted elements (subtree_sum of the root node of a tree)
        '''
        return 0.0 if self.treap.root == -1 else self.treap.nodes[self.treap.root].subtree_sum

    def insert(self, element):
        '''
        Inserts a new element (insertion order is not important)
        '''
        self.treap.add_element(element)

    def get_modelling_sum(Y, n):
        '''
        Returns approximation of the sum of distances required for next combined result modelling
        '''
        ret = 0.0
        for cell in Y:
            for c in cell.keys():
                s = cell[c].get_sum()
                lc, ls = cell[c].treap.get_lower(cell[c].treap.root, s / n)
                ret += s * lc - n * ls
        # quickly computed approximation of GLD
        ret /= (n * (n + 1))
        # post-summation normalization: a rough approximation of normalized GLD
        return 2 * ret / (ret + 2 * len(Y)) 

class ListBasedSequenceStructure:
    def __init__(self):
        self.elements = []  # a list of inserted elements
        self.sum = 0.0      # accumulator for sum of elements

    def get_sum(self):
        '''
        Returns a sum of inserted elements
        '''
        return self.sum

    def insert(self, element):
        '''
        Inserts a new element (order IS IMPORTANT due to the implementation of get_modelling_sum)
        '''
        self.elements.append(element)
        self.sum += element

    def get_modelling_sum(Y, n):
        '''
        Returns approximation of the sum of distances required for next combined result modelling
        '''
        d = [0.0 for _ in range(n)] # separate GLD for each combined sample
        for cell in Y:
            for c in cell.keys():
                s = cell[c].get_sum()
                for i, e in enumerate(cell[c].elements):
                    d[i] += abs(s - n * e)
        for i in range(n):
            d[i] /= (2 * n * (n + 1))
            # normalizing each GLD separately
            d[i] = 2.0 * d[i] / (d[i] + 2 * len(Y))
        return sum(d)


class AlignmentWithEstimation:
    '''
    Represents dynamic alignment of a sequence of input text string recognition results
    with estimation of expected distance to the next combined result.

    Implementation assumes that all input strings will have character classification results
    with the same alphabet. 

    Implementation also assums that all input recognition results have the same weight
    '''
    def __init__(self, empty_weight, SequenceStructure):
        self.base = None                    # aligned and combined result
        self.base_samples = 0               # number of processed samples
        self.es   = '@'                     # empty character label
        self.ec   = Cell({self.es: 1.0})    # empty cell template
        self.ew   = empty_weight            # weight of empty character (for string result generation)

        self.Y = None                       # structure for estimation computation
        self.alphabet = None                # alphabet of character classification results

        self.SequenceStructure = SequenceStructure # class for storing sequences y_1jk, y_2jk, ..., y_njk
        
    def add_string(self, arg_s):
        '''
        Adds a new text string recognition result
        '''
        # input can be either a simple string or a sequence of Cells
        s = [Cell({c: 1.0}) for c in arg_s] if isinstance(arg_s, str) else [c.clone() for c in arg_s]
        
        empty_input_string = False # whether an empty string was passed as an input
        if len(s) == 0:
            # empty sequence is normalized as a sequence of a single empty Cell
            s = [self.ec.clone()]
            # populating alphabet if it is initialized
            if self.alphabet is not None:
                for c in self.alphabet:
                    if c != self.es:
                        s[0].vars[c] = 0.0
            empty_input_string = True

        if self.base is None:
            # the first input becomes the new base
            self.base = [c.normalized() for c in s]
            self.base_samples = 1
            
            self.Y = []
            for cell in self.base:
                y = {}
                for c, v in cell.vars.items():
                    y[c] = self.SequenceStructure()
                    y[c].insert(v)
                self.Y.append(y)
            if not empty_input_string:
                # initializing alphabet
                self.alphabet = ''.join(list(s[0].vars.keys()))
            return

        # if this is a first non-empty input string the stored matrix Y has to be updated
        if not empty_input_string and self.alphabet is None:
            # initializing alphabet
            self.alphabet = ''.join(list(s[0].vars.keys()))
            for i in range(len(self.Y)):
                for c in self.alphabet:
                    if c != self.es and c not in self.Y[i].keys():
                        self.Y[i][c] = self.SequenceStructure()
                        for _ in range(self.base_samples):
                            self.Y[i][c].insert(0.0)
           
        # Alignment DP matrix direction labels
        PATH_UNDEFINED = -1
        PATH_MATCHED = 0
        PATH_UNMATCHED_BASE = 1 # 'deletion'
        PATH_UNMATCHED_S = 2 # 'insertion'
                    
        # levenshtein DP matrix: <len base> x <len s>
        dp = [[] for i in range(len(self.base) + 1)]
        path = [[] for i in range(len(self.base) + 1)]
        for i in range(len(self.base) + 1):
            dp[i] = [0.0 for j in range(len(s) + 1)]
            path[i] = [PATH_UNDEFINED for j in range(len(s) + 1)]
            
        # dp[i][j] = levenshtein(base[preflen i], s[preflen j])
        dp[0][0] = 0.0
        path[0][0] = PATH_UNDEFINED
        for s_preflen in range(1, len(s) + 1):
            dp[0][s_preflen] = cell_dist(s[s_preflen - 1], self.ec) + dp[0][s_preflen - 1]
            path[0][s_preflen] = PATH_UNMATCHED_S
        for b_preflen in range(1, len(self.base) + 1):
            dp[b_preflen][0] = cell_dist(self.base[b_preflen - 1], self.ec) + dp[b_preflen - 1][0]
            path[b_preflen][0] = PATH_UNMATCHED_BASE
        
        for b_preflen in range(1, len(self.base) + 1):
            for s_preflen in range(1, len(s) + 1):
                # cell from base is aligned with empty cell
                pen_unmatched_base = cell_dist(self.base[b_preflen - 1], self.ec) + dp[b_preflen - 1][s_preflen]
                # cell from input is aligned with empty cell
                pen_unmatched_s    = cell_dist(s[s_preflen - 1], self.ec) + dp[b_preflen][s_preflen - 1]
                # cells from base and input are aligned together
                pen_matched        = cell_dist(self.base[b_preflen - 1], s[s_preflen - 1]) + dp[b_preflen - 1][s_preflen - 1]
                
                dp[b_preflen][s_preflen] = min(pen_unmatched_base, pen_unmatched_s, pen_matched)
                # determining the path label to trace the path after DP
                if pen_matched < min(pen_unmatched_base, pen_unmatched_s):
                    path[b_preflen][s_preflen] = PATH_MATCHED
                elif pen_unmatched_base < pen_unmatched_s:
                    path[b_preflen][s_preflen] = PATH_UNMATCHED_BASE
                else:
                    path[b_preflen][s_preflen] = PATH_UNMATCHED_S
        
        # tracing path: building new base
        new_base = []
        new_Y = []
        
        cur_b_preflen = len(self.base)
        cur_s_preflen = len(s)
        while cur_b_preflen > 0 or cur_s_preflen > 0:
            if path[cur_b_preflen][cur_s_preflen] == PATH_MATCHED:
                new_base.append(merge_cells(self.base[cur_b_preflen - 1], s[cur_s_preflen - 1], self.base_samples, 1))
                # updating Y
                new_Y.append(self.Y[cur_b_preflen - 1])
                for c, v in s[cur_s_preflen - 1].vars.items():
                    new_Y[-1][c].insert(v)
                cur_b_preflen -= 1
                cur_s_preflen -= 1                
            elif path[cur_b_preflen][cur_s_preflen] == PATH_UNMATCHED_BASE:
                new_base.append(merge_cells(self.base[cur_b_preflen - 1], self.ec, self.base_samples, 1))
                # updating Y
                new_Y.append(self.Y[cur_b_preflen - 1])
                if self.alphabet is None:
                    new_Y[-1][self.es].insert(1.0)
                else:
                    for c in self.alphabet:
                        if c != self.es:
                            new_Y[-1][c].insert(0.0)
                    new_Y[-1][self.es].insert(1.0)
                cur_b_preflen -= 1
            elif path[cur_b_preflen][cur_s_preflen] == PATH_UNMATCHED_S:
                new_base.append(merge_cells(self.ec, s[cur_s_preflen - 1], self.base_samples, 1))
                # updatind Y
                y = {}
                for c, v in s[cur_s_preflen - 1].vars.items():
                    y[c] = self.SequenceStructure()
                if self.alphabet is None:
                    y[self.es].insert(1.0)
                else:
                    for c in self.alphabet:
                        if c != self.es:
                            for i in range(self.base_samples):
                                y[c].insert(0.0)
                    y[self.es].insert(1.0)
                for c, v in s[cur_s_preflen - 1].vars.items():
                    y[c].insert(v)
                new_Y.append(y)
                cur_s_preflen -= 1
        new_base.reverse()
        new_Y.reverse()

        # setting a new base
        self.base = new_base
        self.base_samples += 1
        self.Y = new_Y

    def get_string_result(self):
        '''
        Returns the string representation of a base, to compare with ground truth
        '''        
        ret = ''
        for cell in self.base:
            best_nonempty = cell.best_key_not_from(self.es)
            best_nonempty_score = 0.0 if best_nonempty == '' else cell.vars[best_nonempty]
            empty_score         = 0.0 if self.es not in cell.vars.keys() else cell.vars[self.es]
            if best_nonempty_score > empty_score * self.ew:
                ret += best_nonempty
        return ret      

    def get_modelling_sum(self):
        '''
         Returns approximation of the sum of distances required for next combined result modelling
        '''
        return self.SequenceStructure.get_modelling_sum(self.Y, self.base_samples)