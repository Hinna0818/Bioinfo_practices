from two_sequence_align import matrix_global, traceback_global

class twoseq_GlobalAligner():
    def __init__(self, seq1, seq2, rules = [4, 3, 4]):
        self.seq1 = seq1
        self.seq2 = seq2
        self.rules = rules
        self.score_matrix = None
        self.aligned_seq1 = None
        self.aligned_seq2 = None
    
    def align(self):
        """
        运行全局比对，生成得分矩阵并进行回溯
        """
        self.score_matrix = matrix_global(self.seq1, self.seq2, self.rules)
        self.aligned_seq1, self.aligned_seq2 = traceback_global(self.seq1, self.seq2, self.score_matrix, self.rules)

    def print_matrix(self):
        """
        打印得分矩阵
        """
        if self.score_matrix is None:
            print("Please run .align() first.")
        else:
            print("Score Matrix:")
            print(self.score_matrix)

    def print_alignment(self):
        """
        打印最终比对结果
        """
        if self.aligned_seq1 is None or self.aligned_seq2 is None:
            print("Please run .align() first.")
        else:
            print("\nAlignment Result:")
            print("Seq1:", ''.join(self.aligned_seq1))
            print("      ", ''.join(self._match_line()))
            print("Seq2:", ''.join(self.aligned_seq2))

    def _match_line(self):
        """
        生成中间的匹配符号行
        """
        match_line = []
        for a, b in zip(self.aligned_seq1, self.aligned_seq2):
            if a == b:
                match_line.append("|")
            elif a == '-' or b == '-':
                match_line.append(" ")
            else:
                match_line.append(".")
        return match_line

## test
if __name__ == "__main__":
    a1 = ['a', 't', 't', 'c', 'c', 'a', 'a', 'g']
    a2 = ['t', 't', 'c', 'g', 'a', 'g', 't']

    aligner = twoseq_GlobalAligner(a1, a2, rules=[4, 3, 4])
    aligner.align()
    aligner.print_matrix()
    aligner.print_alignment()
