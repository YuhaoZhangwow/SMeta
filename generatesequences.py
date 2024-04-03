from random import choice, randint, random

def generate_sequence(length, base_sequence=None, similarity=0.8):
    """生成一个长度为length的随机核酸序列，可以选择一个基序列和相似度"""
    nucleotides = ['A', 'T', 'C', 'G']
    if base_sequence:
        # 生成相似序列
        sequence = ''.join([b if random() < similarity else choice(nucleotides) for b in base_sequence])
        # 如果生成的序列比所需的短，用随机碱基填充到指定长度
        sequence += ''.join(choice(nucleotides) for _ in range(length - len(sequence)))
    else:
        # 生成完全随机序列
        sequence = ''.join(choice(nucleotides) for _ in range(length))
    return sequence

# 生成10条相似的序列，每条长度约1000
similar_sequences = [generate_sequence(1000) for _ in range(10)]

# 生成另一个fasta文件的序列，其中5条与相似序列相似，5条完全不同，长度约200
other_sequences = [generate_sequence(200, base_sequence=choice(similar_sequences), similarity=0.9) for _ in range(5)] + \
                  [generate_sequence(200) for _ in range(5)]

# 准备生成fasta格式的文件内容
def to_fasta(sequences, similar=True):
    """将序列列表转换为fasta格式的字符串"""
    fasta_content = ""
    for i, seq in enumerate(sequences):
        header = f">seq{'_similar' if similar else ''}_{i+1}\n"
        fasta_content += header + '\n'.join([seq[j:j+80] for j in range(0, len(seq), 80)]) + "\n"
    return fasta_content

# 生成两个fasta文件的内容
similar_fasta_content = to_fasta(similar_sequences)
other_fasta_content = to_fasta(other_sequences, similar=False)

# 保存到文件
with open('/mnt/data/similar_sequences.fasta', 'w') as f:
    f.write(similar_fasta_content)

with open('/mnt/data/other_sequences.fasta', 'w') as f:
    f.write(other_fasta_content)

# 返回文件路径
similar_fasta_file = '/mnt/data/similar_sequences.fasta'
other_fasta_file = '/mnt/data/other_sequences.fasta'
similar_fasta_file, other_fasta_file
