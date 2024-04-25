from Sequence import Sequence
from Dotplot import Dotplot
from SequenceAlignment import SequenceAlignment

#seq1 = Sequence("MVHLTDAEKSAVSCLWAKVNPDEVGGEALGRLLVVYPWTQRYFDSFGDLSSASAIMGNPKVKAHGKKVITAFNEGLKNLDNLKGTFASLSELHCDKLHVDPENFRLLGNAIVIVLGHHLGKDFTPAAQAAFQKVVAGVATALAHKYH", "human", "P02089")
#seq2 = Sequence("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH", "mouse", "P68871")
#seq1 = Sequence("MVHHH")
#seq2 = Sequence("MVHLT")
#dp = Dotplot.from_sequences(seq1, seq2)

#print(dp)
#output_seq1, output_seq2 = dp.needleman_wunsch(10,-20,-20)
seq1 ="AATCG"
seq2 ="AACG"
dna = SequenceAlignment(seq1, seq2)
output_seq1, output_seq2 = dna.needleman_wunsch(10,-20,-20)
print(output_seq1)
print("------------")
print(output_seq2)