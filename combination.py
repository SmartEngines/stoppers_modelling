# This module contains implementation of the text string recognition results
# combination algorithm which accounts for per-character alternatives

class Cell:
    '''
    Object of Cell class represents a single character classification result
    '''
    def __init__(self, varmap):
        '''
        self.vars: {<class label>: <membership estimation>}
        '''
        self.vars = varmap
    def __str__(self):
        return str(self.vars)
    def __repr__(self):
        return str(self.vars)
    def clone(self):
        return Cell({k: v for k, v in self.vars.items()})    
    
    def normalize(self):
        '''
        Normalizes the current Cell (sum of membership estimations brought to 1)
        '''
        s = 0.0
        for k, v in self.vars.items():
            s += v
        for k in self.vars.keys():
            self.vars[k] /= s
    
    def normalized(self):
        '''
        Returns a normalized clone of the current Cell
        '''
        ret = self.clone()
        ret.normalize()
        return ret

    def best_key(self):
        '''
        Returns class label with highest membership estimation
        '''
        ret = ''
        best_value = -1.0
        for k, v in sorted(self.vars.items()):
            if v > best_value:
                best_value = v
                ret = k
        return ret

    def best_key_not_from(self, exclusion_mask):
        '''
        Returns class label with highest membership estimation, excluding
        character labels from 'exclusion_mask'
        '''
        ret = ''
        best_value = -1.0
        for k, v in sorted(self.vars.items()):
            if k in exclusion_mask:
                continue
            if v > best_value:
                best_value = v
                ret = k
        return ret
    

def cell_dist(a, b):
    '''
    Implements a scaled taxicab metric for Cell objects
    '''
    na = a.normalized()
    nb = b.normalized()

    ret = 0.0
    for a_key in na.vars.keys():
        if a_key in nb.vars.keys():
            ret += abs(na.vars[a_key] - nb.vars[a_key])
        else:
            ret += na.vars[a_key]
    for b_key in nb.vars.keys():
        if b_key not in na.vars.keys():
            ret += nb.vars[b_key]

    return ret / 2.0    


def merge_cells(a, b, wa, wb):
    '''
    Combines Cells a and b with respective weights wa and wb
    Output membership estimations for each class are weighted average of corresponding
    ones from a and b
    '''
    na = a.normalized()
    nb = b.normalized()
    
    ret = na.clone()
    for b_key in nb.vars.keys():
        if b_key in ret.vars.keys(): # keys in both cells
            ret.vars[b_key] = (wa * ret.vars[b_key] + wb * nb.vars[b_key]) / (wa + wb)
        else: # only b-keys
            ret.vars[b_key] = (wb * nb.vars[b_key]) / (wa + wb)
    for a_key in na.vars.keys():
        if a_key not in nb.vars.keys(): # only a-keys
            ret.vars[a_key] = (wa * ret.vars[a_key]) / (wa + wb)
    
    return ret
    
def levmetric_ocr(ocr_string1, ocr_string2):
    '''
    Computes a normalized Generalized Levenshtein Distance for two text string
    recognition results with per-character alternatives
    '''
    if (len(ocr_string1) == 0) and (len(ocr_string2) == 0):
        return 0.0;

    # levenshtein DP matrix: <len base> x <len s>
    dp = [[] for i in range(len(ocr_string1) + 1)]
    for i in range(len(ocr_string1) + 1):
        dp[i] = [0.0 for j in range(len(ocr_string2) + 1)]

    # dp[i][j] is a prefix-wise distance levenshtein(a[preflen i], b[preflen j])
    dp[0][0] = 0.0
    for s_preflen in range(1, len(ocr_string2) + 1):
        dp[0][s_preflen] = cell_dist(ocr_string2[s_preflen - 1], Cell({'@': 1.0})) + dp[0][s_preflen - 1]
    for b_preflen in range(1, len(ocr_string1) + 1):
        dp[b_preflen][0] = cell_dist(ocr_string1[b_preflen - 1], Cell({'@': 1.0})) + dp[b_preflen - 1][0]
    
    for b_preflen in range(1, len(ocr_string1) + 1):
        for s_preflen in range(1, len(ocr_string2) + 1):
            # cell from ocr_string1 is aligned with empty cell
            pen_unmatched_base = cell_dist(ocr_string1[b_preflen - 1], Cell({'@': 1.0})) + dp[b_preflen - 1][s_preflen]
            # cell from ocr_string2 is aligned with empty cell
            pen_unmatched_s    = cell_dist(ocr_string2[s_preflen - 1], Cell({'@': 1.0})) + dp[b_preflen][s_preflen - 1]
            # cells of ocr_string1 and ocr_string2 are aligned together
            pen_matched        = cell_dist(ocr_string1[b_preflen - 1], ocr_string2[s_preflen - 1]) + dp[b_preflen - 1][s_preflen - 1]
            
            dp[b_preflen][s_preflen] = min(pen_unmatched_base, pen_unmatched_s, pen_matched)

    # Generalized Levenshtein Distance value
    levenshtein = dp[len(ocr_string1)][len(ocr_string2)]
    # returning normalized value
    return 2.0 * levenshtein / (len(ocr_string1) + len(ocr_string2) + levenshtein)


class Alignment:
    '''
    Represents dynamic alignment of a sequence of input text string recognition results
    '''
    def __init__(self, empty_weight):
        '''
        Ctor accepts a single parameter (empty_weight) which only affects the generation
        of the final string representation of the combined result (and does not affect the
        combination algorithm itself)
        '''
        self.base = None                    # aligned and combined result
        self.base_weight = 0.0              # sum of weights of combined inputs
        self.es   = '@'                     # empty character label
        self.ec   = Cell({self.es: 1.0})    # empty cell template
        self.ew   = empty_weight            # weight of empty character (for string result generation)
        
    def clone(self):
        '''
        Clones the full structure
        '''
        ret = Alignment(self.ew)
        ret.base = [elem.clone() for elem in self.base]
        ret.base_weight = self.base_weight
        return ret
            
    def add_string(self, arg_s, weight):
        '''
        Adds a new text string recognition result with weight 'weight'
        '''
        # input can be either a simple string or a sequence of Cells
        s = [Cell({c: 1.0}) for c in arg_s] if isinstance(arg_s, str) else [c.clone() for c in arg_s]
        if len(s) == 0:
            # empty sequence is normalized as a sequence of a single empty Cell
            s = [self.ec.clone()]
        if self.base is None:
            # the first input becomes the new base
            self.base = [c.normalized() for c in s]
            self.base_weight = weight
            return
           
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
        
        cur_b_preflen = len(self.base)
        cur_s_preflen = len(s)
        while cur_b_preflen > 0 or cur_s_preflen > 0:
            if path[cur_b_preflen][cur_s_preflen] == PATH_MATCHED:
                new_base.append(merge_cells(self.base[cur_b_preflen - 1], s[cur_s_preflen - 1], self.base_weight, weight))
                cur_b_preflen -= 1
                cur_s_preflen -= 1                
            elif path[cur_b_preflen][cur_s_preflen] == PATH_UNMATCHED_BASE:
                new_base.append(merge_cells(self.base[cur_b_preflen - 1], self.ec, self.base_weight, weight))
                cur_b_preflen -= 1
            elif path[cur_b_preflen][cur_s_preflen] == PATH_UNMATCHED_S:
                new_base.append(merge_cells(self.ec, s[cur_s_preflen - 1], self.base_weight, weight))
                cur_s_preflen -= 1
        new_base.reverse()
        
        # setting a new base
        self.base = new_base
        self.base_weight += weight

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