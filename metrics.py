# This python module contains implementations for basic functions 
# required for the experiments

def normalize_string(string):
    '''
    Normalizes a recognition result string by converting it to upper case and 
    replacing the letter O to the digit 0 for ease of recognition results comparison
    '''
    return string.upper().replace('O', '0')

def levenshtein(a, b):
    '''
    Computes a common Levenshtein distance between two strings X and Y
    Complexity O(|X|*|Y|)
    '''
    if min(len(a), len(b)) == 0:
        return max(len(a), len(b))
    
    dp = [[] for i in range(len(a) + 1)]
    for i in range(len(a) + 1):
        dp[i] = [0 for j in range(len(b) + 1)]
        
    # dp[i][j] is a prefix-wise distance levenshtein(a[preflen i], b[preflen j])
    for b_preflen in range(len(b) + 1):
        dp[0][b_preflen] = b_preflen
    for a_preflen in range(1, len(a) + 1):
        dp[a_preflen][0] = a_preflen
        for b_preflen in range(1, len(b) + 1):
            # insertion or deletion (penalty 1)
            dp[a_preflen][b_preflen] = 1 + min(dp[a_preflen - 1][b_preflen], dp[a_preflen][b_preflen - 1])
            # substitution penalty
            subs_penalty = 0 if (a[a_preflen - 1] == b[b_preflen - 1]) else 1
            # substitution
            dp[a_preflen][b_preflen] = min(dp[a_preflen][b_preflen], subs_penalty + dp[a_preflen - 1][b_preflen - 1])
    
    return dp[len(a)][len(b)]

def end_to_end(result, ideal):
    '''
    Computes an 'end-to-end' metric between strings (0 for equal strings only)
    '''
    return 0 if normalize_string(result) == normalize_string(ideal) else 1

def d_levenshtein(a, b):
    '''
    Computes a Levenshtein distance for strings (with prior strings normalization)
    '''
    return 1.0 * levenshtein(normalize_string(a), normalize_string(b))

def levmetric(a, b):
    '''
    Computes a Normalized Levenshtein distance for strings (with prior normalization)
    '''
    if len(a) == 0 and len(b) == 0:
        return 0.0
    l = d_levenshtein(a, b)
    return (2.0 * l) / (len(a) + len(b) + l)