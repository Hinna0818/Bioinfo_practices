import numpy as np

def cos_similarity(vec1, vec2):
    '''
    compute cosine-similarity between two input vectors
    '''
    if len(vec1) != len(vec2):
        return None
    
    if not isinstance(vec1, np.ndarray):
        vec1 = np.array(vec1)
    
    if not isinstance(vec2, np.ndarray):
        vec2 = np.array(vec2)
    
    dot_product = np.dot(vec1, vec2)
    len1 = np.sqrt(np.sum(vec1 ** 2))
    len2 = np.sqrt(np.sum(vec2 ** 2))
    res = dot_product / (len1 * len2)

    return np.round(res, 4)

if __name__ == "__main__":
    vec1 = [5,1,3,4,2]
    vec2 = [4,2,3,5,1]
    output = cos_similarity(vec1, vec2)
    print(output)