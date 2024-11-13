from math import log, sqrt, exp, isnan, factorial
import random

genetic_code = 1 # default code type
translate_table = [
 "FFLLSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "1-Standard Code",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS!!VVVVAAAADDEEGGGG", "2-Vertebrate Mitochondrial Code",
 "FFLLSSSSYY!!CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "3-Yeast Mitochondrial Code",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "4-Mold Mitochondrial Code",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG", "5-Invertebrate Mitochondrial Code",
 "FFLLSSSSYYQQCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "6-Ciliate, Dasycladacean and Hexamita Code",
 "", "7-",
 "", "8-",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "9-Echinoderm and Flatworm Mitochondrial Code",
 "FFLLSSSSYY!!CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "10-Euplotid Nuclear Code",
 "FFLLSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "11-Bacterial and Plant Plastid Code",
 "FFLLSSSSYY!!CC!WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "12-Alternative Yeast Nuclear Code",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG", "13-Ascidian Mitochondrial Code",
 "FFLLSSSSYYY!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "14-Alternative Flatworm Mitochondrial Code",
 "FFLLSSSSYY!QCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "15-Blepharisma Nuclear Code",
 "FFLLSSSSYY!LCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "16-Chlorophycean Mitochondrial Code",
 "", "17-",
 "", "18-",
 "", "19-",
 "", "20-",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "21-Trematode Mitochondrial Code",
 "FFLLSS!SYY!LCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "22-Scenedesmus obliquus mitochondrial Code",
 "FF!LSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "23-Thraustochytrium Mitochondrial Code",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "24-Rhabdopleuridae Mitochondrial Code",
 "FFLLSSSSYY!!CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "25-Candidate Division SR1 and Gracilibacteria Code",
 "FFLLSSSSYY!!CC!WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "26-Pachysolen tannophilus Nuclear Code",
 "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "27-Karyorelict Nuclear Code",
 "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "28-Condylostoma Nuclear Code",
 "FFLLSSSSYYYYCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "29-Mesodinium Nuclear Code",
 "FFLLSSSSYYEECC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "30-Peritrich Nuclear Code",
 "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "31-Blastocrithidia Nuclear Code",
 "", "32-",
 "FFLLSSSSYYY!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "33-Cephalodiscidae Mitochondrial UAA-Tyr Code"
]


seq_name = "seq"       # 序列名称
length = 0          # 序列长度
GC = [0.0] * 4      # GC 内容
NA = -1
NUMBER_OF_RATES = 6
SMALLVALUE = 1e-6  # Small value to avoid floating-point precision issues

class Base:
    def __init__(self):
        # 内部变量
        self.Si = [0.0] * 5
        self.Vi = [0.0] * 5
        self.L = [0.0] * 5
        self.KAPPA = [1] * NUMBER_OF_RATES
        
        self.SEKa = self.SEKs = self.AICc = self.lnL = self.AkaikeWeight = NA
        self.Ka = self.Ks = self.Sd = self.Nd = self.S = self.N = self.snp = self.t = self.kappa = 0.0
        self.model = ""
        self.name = ""
        self.S = 0
        self.N = 0
        self.Sd = 0
        self.Nd = 0
        self.snp = 0
        self.Ka = 0
        self.Ks = 0
        self.SEKa = 0
        self.SEKs = 0
        self.Ps = 0
        self.Pn = 0
        self.t = 0
        self.lnL = 0
        self.AICc = 0
        self.AkaikeWeight = 0
        self.model = ""
        # NUMBER_OF_RATES = 6
        self.DNASIZE = 4
        self.XSIZE = self.DNASIZE * self.DNASIZE
        self.KAPPA = [0] * NUMBER_OF_RATES

    def writeFile(self, result, filename):
        try:
            with open(filename, 'a') as out:
                out.write(result)
        except IOError as e:
            print(f"Cannot write to file {filename}: {e}")
    
    def parse_output(self):
        result = ""

        # Sequence name
        result = self.add_string(result, seq_name)
        # Method name
        result = self.add_string(result, self.name)

        # Ka
        if self.Ka < SMALLVALUE:
            tmp = "NA"
        else:
            tmp = str(self.Ka)
        result = self.add_string(result, tmp)
        
        # Ks
        if self.Ks < SMALLVALUE:
            tmp = "NA"
        else:
            tmp = str(self.Ks)
        result = self.add_string(result, tmp)
        
        # Ka/Ks
        if self.Ks < SMALLVALUE or isnan(self.Ks) or isnan(self.Ka):
            tmp = "NA"
        else:
            tmp = str(self.Ka / self.Ks)
        result = self.add_string(result, tmp)

        # Fisher's test: p_value
        if self.Sd < SMALLVALUE or self.Nd < SMALLVALUE or self.S < SMALLVALUE or self.N < SMALLVALUE:
            tmp = "NA"
        else:
            tmp = str(self.fisher(self.Sd, self.Nd, self.S - self.Sd, self.N - self.Nd))
        result = self.add_string(result, tmp)

        # Length of compared pairwise sequences
        result = self.add_string(result, str(length))
        
        # Synonymous(S) sites
        if self.S < SMALLVALUE:
            tmp = "NA"
        else:
            tmp = str(self.S)
        result = self.add_string(result, tmp)

        # Nonsynonymous(N) sites
        if self.N < SMALLVALUE:
            tmp = "NA"
        else:
            tmp = str(self.N)
        result = self.add_string(result, tmp)

        # L[0], L[2], L[4] only for Prof.Li's series(LWL85, LPB93...)
        if self.L[0] < SMALLVALUE and self.L[2] < SMALLVALUE and self.L[4] < SMALLVALUE:
            tmp = "NA"        
        else:        
            tmp = f"{str(self.L[0])}:{str(self.L[2])}:{str(self.L[4])}"
        result = self.add_string(result, tmp)

        # Substitutions
        result = self.add_string(result, str(self.snp))

        # Synonymous(Sd) Substitutions(Nd)
        if self.Sd > SMALLVALUE:
            tmp = str(self.Sd)        
        else:
            tmp = "NA"
        result = self.add_string(result, tmp)
        
        # Nonsynonymous Substitutions(Nd)
        if self.Nd > SMALLVALUE:
            tmp = str(self.Nd)        
        else:
            tmp = "NA"
        result = self.add_string(result, tmp)

        # Si for Li's series' methods(LWL85, LPB93...)
        if self.Si[0] != 0.0 or self.Si[2] != 0.0 or self.Si[4] != 0.0: # Si[0], Si[2], Si[4]
            tmp  = f"{str(self.Si[0])}:{str(self.Si[2])}:{str(self.Si[4])}"
        else:
            tmp = "NA"
        result = self.add_string(result, tmp)

        # Vi for Li's series' methods(LWL85, LPB93...)
        if self.Vi[0] != 0.0 or self.Vi[2] != 0.0 or self.Vi[4] != 0.0: # Vi[0], Vi[2], Vi[4]
            tmp  = f"{str(self.Vi[0])}:{str(self.Vi[2])}:{str(self.Vi[4])}"
        else:
            tmp = "NA"
        result = self.add_string(result, tmp)

        # Divergence time or distance t = (S*Ks+N*Ka)/(S+N)
        if self.t < SMALLVALUE:
            tmp = "NA"
        else:
            tmp = str(self.t)
        result = self.add_string(result, tmp)

        # Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)
        tmp = ":".join(str(k) for k in self.KAPPA[:-1])
        tmp += f":{str(self.KAPPA[-1])}"
        result = self.add_string(result, tmp)

        # GC Content
        tmp = f"{str(GC[0])}({str(GC[1])}:{str(GC[2])}:{str(GC[3])})"
        result = self.add_string(result, tmp)
        
        # Maximum Likelihood Value
        if isnan(self.lnL):
            tmp = "NA"
        else:
            tmp = str(self.lnL)
        result = self.add_string(result, tmp)
        
        # AICc
        if isnan(self.AICc):
            tmp = "NA"
        else:
            tmp = str(self.AICc)
        result = self.add_string(result, tmp)

        # Akaike weight in model selection
        if isnan(self.AkaikeWeight):
            tmp = "NA"
        else:
            tmp = str(self.AkaikeWeight)
        result = self.add_string(result, tmp)
        
        # Selected Model according to AICc
        if self.model == "" or len(self.model) == 0:
            tmp = "NA"
        else:
            tmp = self.model
        result = self.add_string(result, tmp, "\n")

        # Uncomment and modify if Standard Errors are needed
        # Standard Errors
        # if math.isnan(SEKa):
        #     tmp = "NA"
        # else:
        #     tmp = str(SEKa)
        # result = self.add_string(result, tmp, "\t")
        # 
        # if math.isnan(SEKs):
        #     tmp = "NA"
        # else:
        #     tmp = str(SEKs)
        # result = self.add_string(result, tmp, "\n")
        
        return result
    
    
    def parse_input(self, input_file):
        '''
        convert msa output 2 seq dict
        '''
        
        
    
    def add_string(self, result, string, flag="\t"):
        # Format string for outputing into file
        result += string
        result += flag
        
        return result
    
    def get_random(self):
        # Generate a random integer
        random.randint(0, 2**31 - 1)
    
    def convert_char(self, ch):
        # Convert a char T,C,A,G into a digit 0,1,2,3, respectively
        ret = -1
        if ch == 'T' or ch == 'U':
            ret = 0
        elif ch == 'C':
            ret = 1
        elif ch == 'A':
            ret = 2
        elif ch == 'G':
            ret = 3
        return ret
    
    def convert_int(self, i):
        # Convert a digit 0,1,2,3 into a char T,C,A,G, respectively
        ch = '-'
        if i == 0:
            ch = 'T'
        elif i == 1:
            ch = 'C'
        elif i == 2:
            ch = 'A'
        elif i == 3:
            ch = 'G'
        return ch
    
    def string_to_upper(self, string:str):
        # Convert a string to uppercase
        return string.upper()
    
    def getID(self, codon):
        # Return the codon's id from codon table
        return (self.convert_char(codon[0]) * self.XSIZE + self.convert_char(codon[1]) * self.DNASIZE + self.convert_char(codon[2]))
    
    def get_amino_acid(self, codon):
        if isinstance(codon, int):
            return translate_table[2*(genetic_code-1)][codon]
        if isinstance(codon, list):print(codon)
        return translate_table[2*(genetic_code-1)][self.getID(codon)]
    
    def get_num_nonsense(self, genetic_code):
        # Get the number of stop codon in a given genetic code table
        CODON = 64
        num = 0
        for i in range(CODON):
            if self.get_amino_acid(i) == '!':
                num += 1
        return num
    
    
    def get_codon(self, IDcodon):
        # Return a codon according to the id
        codon_li = ['T'] * 3
        if 0 <= IDcodon < 64:
            codon_li[0] = self.convert_int(IDcodon // 16)
            codon_li[1] = self.convert_int((IDcodon % 16) // 4)
            codon_li[2] = self.convert_int(IDcodon % 4)
            codon = "".join(codon_li)
        return codon
    
    def sum_array(self, arr, size):
        # Sum array's elements
        return sum(arr[:size])
    
    def init_array(self, x, n, value=0):
        # Init value to array
        x[:n] = [value] * n
    
    def scale_array(self, scale, x, n):
        # Elements in array are mutipled by scale
        for i in range(n):
            scale[i] *= x
    
    def norm(self, x, n):
        # Sqrt of the sum of the elements' square
        return sqrt(sum(x[i] * x[i] for i in range(n)))
    
    def copy_array(self, from_list, to_list, n):
        # Copy array's values one by one: to[] = from[]
        to_list[:n] = from_list[:n]
    
    def innerp(self, x, y, n):
        # Sum of 'n' products multiplied by two elements x[], y[]
        return sum(x[i] * y[i] for i in range(n))
    
    def init_identity_matrix(self, a, n):
        # Set x[i,j]=0 when x!=j and x[i,j]=1 when x=j
        for i in range(n):
            for j in range(n):
                a[i][j] = 1.0 if i == j else 0.0
    
    def fisher(self, sd, s, nd, n):
        # Compute p-value by Fisher exact test to justify the validity of ka/ks

        matrix = [sd, s, nd, n]
        R = [matrix[0] + matrix[2], matrix[1] + matrix[3]]
        C = [matrix[0] + matrix[1], matrix[2] + matrix[3]]
        sum_val = R[0] + R[1]

        numerator = self.factorial(R[0])
        numerator += self.factorial(R[1])
        numerator += self.factorial(C[0])
        numerator += self.factorial(C[1])

        fac_sum = self.factorial(sum_val)
        denominator = fac_sum
        for i in range(4):
            denominator += self.factorial(matrix[i])

        prob_current = exp(numerator - denominator)

        prob_total = 0.0
        for i in range(int(R[0] + SMALLVALUE)):
            matrix[0] = i
            matrix[1] = C[0] - i
            matrix[2] = R[0] - i
            matrix[3] = R[1] - C[0] + i
            if all(x > SMALLVALUE for x in matrix):
                denominator = fac_sum
                for j in range(4):
                    denominator += self.factorial(matrix[j])
                temp = exp(numerator - denominator)
                if temp <= prob_current:
                    prob_total += temp
        # print(f"prob_total is {prob_total}")
        return prob_total
    
    def factorial(self, n):
        # factorial
        temp = 1.0
        if n > 0:
            n += 1
            x = 0
            x += 0.1659470187408462e-06 / (n + 7)
            x += 0.9934937113930748e-05 / (n + 6)
            x -= 0.1385710331296526 / (n + 5)
            x += 12.50734324009056 / (n + 4)
            x -= 176.6150291498386 / (n + 3)
            x += 771.3234287757674 / (n + 2)
            x -= 1259.139216722289 / (n + 1)
            x += 676.5203681218835 / n
            x += 0.9999999999995183
            temp = log(x) - 5.58106146679532777 - n + (n - 0.5) * log(n + 6.5)
        return temp


class NG86(Base):
    def __init__(self):
        super().__init__()
        self.GAMMA = 0
        self.Ps = 0.0
        self.Pn = 0.0
        self.name = "NG"

    def get_codon_site(self, codon)->None:
        
        temp = ""
        syn = 0.0
        stop = 0
        if (self.get_amino_acid(codon) == '!'):return
        for i in range(0, 3, 2):
            for j in range(4):
                temp = list(codon)
                if (j != self.convert_char(temp[i])):
                    temp[i] = self.convert_int(j)
                    if (self.get_amino_acid(''.join(temp)) == '!'):
                        stop += 1
                    else:
                        if (self.get_amino_acid(''.join(temp)) == self.get_amino_acid(codon)):
                            syn += 1
        # print(f"syn is {syn}")
        # print(f"stop is {stop}")
        # print(stop)
        self.S += (syn/3.0)
        
        self.N += (3-stop/3.0-syn/3.0)

    def get_codon_difference(self, codon1, codon2)->None:
        # Count codon's differences
        diff = [-1] * 3
        num = 0
        stop = 0
        path = 1
        sd_temp = 0.0
        nd_temp = 0.0
        temp1 = ""
        temp2 = ""
        if self.get_amino_acid(codon1) == '!' or self.get_amino_acid(codon2) == '!':
            return

        for i in range(len(codon1)):
            if codon1[i] != codon2[i]:
                diff[num] = i
                num += 1

        if num == 0:
            return

        self.snp += num

        path = factorial(num)
        if num == 1:
            if self.get_amino_acid(codon1) == self.get_amino_acid(codon2):
                sd_temp += 1
            else:
                nd_temp += 1

        if num == 2:
            for i in range(num):
                for j in range(num):
                    if i != j:
                        temp1 = list(codon1)
                        temp1[diff[i]] = codon2[diff[i]]
                        if self.get_amino_acid(''.join(temp1)) != '!':
                            if self.get_amino_acid(''.join(temp1)) == self.get_amino_acid(codon1):
                                sd_temp += 1
                            else:
                                nd_temp += 1
                            temp2 = list(temp1)
                            temp2[diff[j]] = codon2[diff[j]]
                            if self.get_amino_acid(''.join(temp2)) == self.get_amino_acid(''.join(temp1)):
                                sd_temp += 1
                            else:
                                nd_temp += 1
                        else:
                            stop += 1

        if num == 3:
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        if i != j and i != k and j != k:
                            temp1 = list(codon1)
                            temp1[diff[i]] = codon2[diff[i]]
                            temp2 = list(temp1)
                            temp2[diff[j]] = codon2[diff[j]]
                            if self.get_amino_acid(''.join(temp1)) != '!' and self.get_amino_acid(''.join(temp2)) != '!':
                                if self.get_amino_acid(''.join(temp1)) == self.get_amino_acid(codon1):
                                    sd_temp += 1
                                else:
                                    nd_temp += 1
                                if self.get_amino_acid(''.join(temp2)) == self.get_amino_acid(''.join(temp1)):
                                    sd_temp += 1
                                else:
                                    nd_temp += 1
                                if self.get_amino_acid("".join(temp2)) == self.get_amino_acid(codon2):
                                    sd_temp += 1
                                else:
                                    nd_temp += 1
                            else:
                                stop += 1
        

        if path == stop:
            if num == 2:
                self.Sd += 0.5
                self.Nd += 1.5
            else:
                self.Sd += 1.0
                self.Nd += 2.0
        else:
            self.Sd += sd_temp / (path - stop)
            self.Nd += nd_temp / (path - stop)

    def preprocess(self, seq1, seq2):
        # Preprocess
        for i in range(0, len(seq1), 3):
            self.get_codon_site(seq1[i:i+3])
            self.get_codon_site(seq2[i:i+3])
            self.get_codon_difference(seq1[i:i+3], seq2[i:i+3])
        self.S /= 2.0
        self.N /= 2.0
        y = len(seq1) / (self.S + self.N)
        self.S *= y
        self.N *= y

    def kaks_formula(self, p):
        NA = -1
        # Jukes and Cantor's one-parameter formula
        d = 1 - (4 * p) / 3
        if d < 0.0:
            d = NA
        else:
            if self.GAMMA == 6 or self.GAMMA == -1:
                self.name = "GNG"
            if self.GAMMA == 6:
                d = pow(d, -1.0 / 0.6) - 1
                if d < 0.0:
                    d = NA
                else:
                    d = (3 * d * 0.6) / 4.0
            else:
                d = log(d)
                if d > 0.0:
                    d = NA
                else:
                    d = (-3.0) * d / 4.0
        return d
    
    def run(self, seq1, seq2)->None:
        # Main function of calculating kaks
        self.preprocess(seq1, seq2)
        # self.Sd = 86.333333
        # self.Nd = 99.666667
        # self.S = 192.310304
        # self.N = 677.689696
        self.Ks = self.kaks_formula(self.Sd / self.S)
        self.Ka = self.kaks_formula(self.Nd / self.N)
        self.t = (self.S * self.Ks + self.N * self.Ka) / (self.S + self.N)
        print(self.parse_output())

class NONE(NG86):

    def __init__(self):
        super().__init__()
        self.name = "NONE"

    def run(self, seq1, seq2):
        self.preprocess(seq1, seq2)
        self.Ks = self.Sd / self.S
        self.Ka = self.Nd / self.N
        self.t = (self.S * self.Ks + self.N * self.Ka) / (self.S + self.N)
        print(self.Ka)
        print(self.Ks)
        return self.parse_output()


if __name__ == "__main__":
    ins = NG86()
    ins.run("ATGGACATTGAAGCATATTTTGAAAGAATTGGCTATAAGAACTCTAGGAACAAATTGGACTTGGAAACATTAACTGACATTCTTGAGCACCAGATCCGGGCTGTTCCCTTTGAGAACCTTAACATGCATTGTGGGCAAGCCATGGAGTTGGGCTTAGAGGCTATTTTTGATCACATTGTAAGAAGAAACCGGGGTGGGTGGTGTCTCCAGGTCAATCAACTTCTGTACTGGGCTCTGACCACAATCGGTTTTCAGACCACAATGTTAGGAGGGTATTTTTACATCCCTCCAGTTAACAAATACAGCACTGGCATGGTTCACCTTCTCCTGCAGGTGACCATTGACGGCAGGAATTACATTGTCGATGCTGGGTCTGGAAGCTCCTCCCAGATGTGGCAGCCTCTAGAATTAATTTCTGGGAAGGATCAGCCTCAGGTGCCTTGCATTTTCTGCTTGACAGAAGAGAGAGGAATCTGGTACCTGGACCAAATCAGGAGAGAGCAGTATATTACAAACAAAGAATTTCTTAATTCTCATCTCCTGCCAAAGAAGAAACACCAAAAAATATACTTATTTACGCTTGAACCTCGAACAATTGAAGATTTTGAGTCTATGAATACATACCTGCAGACGTCTCCAACATCTTCATTTATAACCACATCATTTTGTTCCTTGCAGACCCCAGAAGGGGTTTACTGTTTGGTGGGCTTCATCCTCACCTATAGAAAATTCAATTATAAAGACAATACAGATCTGGTCGAGTTTAAAACTCTCACTGAGGAAGAGGTTGAAGAAGTGCTGAAAAATATATTTAAGATTTCCTTGGGGAGAAATCTCGTGCCCAAACCTGGTGATGGATCCCTTACTATT","ATGGACATCGAAGCATACTTTGAAAGGATTGGTTACAAGAACTCAGTGAATAAATTGGACTTAGCCACATTAACTGAAGTTCTTCAGCACCAGATGCGAGCAGTTCCTTTTGAGAATCTTAACATGCATTGTGGAGAAGCCATGCATCTGGATTTACAGGACATTTTTGACCACATAGTAAGGAAGAAGAGAGGTGGATGGTGTCTCCAGGTTAATCATCTGCTGTACTGGGCTCTGACCAAAATGGGCTTTGAAACCACAATGTTGGGAGGATATGTTTACATAACTCCAGTCAGCAAATATAGCAGTGAAATGGTCCACCTTCTAGTACAGGTGACCATCAGTGACAGGAAGTACATTGTGGATTCCGCCTATGGAGGCTCCTACCAGATGTGGGAGCCTCTGGAATTAACATCTGGGAAGGATCAGCCTCAGGTGCCTGCCATCTTCCTTTTGACAGAGGAGAATGGAACCTGGTACTTGGACCAAATCAGAAGAGAGCAGTATGTTCCAAATGAAGAATTTGTTAACTCAGACCTCCTTGAAAAGAACAAATATCGAAAAATCTACTCCTTTACTCTTGAGCCCCGAGTTATCGAGGATTTTGAATATGTGAATAGCTATCTTCAGACATCGCCAGCATCTGTGTTTGTAAGCACATCGTTCTGTTCCTTGCAGACCTCGGAAGGGGTTCACTGTTTAGTGGGCTCCACCTTTACAAGTAGGAGATTCAGCTATAAGGACGATGTAGATCTGGTTGAGTTTAAATATGTGAATGAGGAAGAAATAGAAGATGTACTGAAAACCGCATTTGGCATTTCTTTGGAGAGAAAGTTTGTGCCCAAACATGGTGAACTAGTTTTTACTATT")
    # ins = NONE()
    # ins.run("ATGGACATTGAAGCATATTTTGAAAGAATTGGCTATAAGAACTCTAGGAACAAATTGGACTTGGAAACATTAACTGACATTCTTGAGCACCAGATCCGGGCTGTTCCCTTTGAGAACCTTAACATGCATTGTGGGCAAGCCATGGAGTTGGGCTTAGAGGCTATTTTTGATCACATTGTAAGAAGAAACCGGGGTGGGTGGTGTCTCCAGGTCAATCAACTTCTGTACTGGGCTCTGACCACAATCGGTTTTCAGACCACAATGTTAGGAGGGTATTTTTACATCCCTCCAGTTAACAAATACAGCACTGGCATGGTTCACCTTCTCCTGCAGGTGACCATTGACGGCAGGAATTACATTGTCGATGCTGGGTCTGGAAGCTCCTCCCAGATGTGGCAGCCTCTAGAATTAATTTCTGGGAAGGATCAGCCTCAGGTGCCTTGCATTTTCTGCTTGACAGAAGAGAGAGGAATCTGGTACCTGGACCAAATCAGGAGAGAGCAGTATATTACAAACAAAGAATTTCTTAATTCTCATCTCCTGCCAAAGAAGAAACACCAAAAAATATACTTATTTACGCTTGAACCTCGAACAATTGAAGATTTTGAGTCTATGAATACATACCTGCAGACGTCTCCAACATCTTCATTTATAACCACATCATTTTGTTCCTTGCAGACCCCAGAAGGGGTTTACTGTTTGGTGGGCTTCATCCTCACCTATAGAAAATTCAATTATAAAGACAATACAGATCTGGTCGAGTTTAAAACTCTCACTGAGGAAGAGGTTGAAGAAGTGCTGAAAAATATATTTAAGATTTCCTTGGGGAGAAATCTCGTGCCCAAACCTGGTGATGGATCCCTTACTATT","ATGGACATCGAAGCATACTTTGAAAGGATTGGTTACAAGAACTCAGTGAATAAATTGGACTTAGCCACATTAACTGAAGTTCTTCAGCACCAGATGCGAGCAGTTCCTTTTGAGAATCTTAACATGCATTGTGGAGAAGCCATGCATCTGGATTTACAGGACATTTTTGACCACATAGTAAGGAAGAAGAGAGGTGGATGGTGTCTCCAGGTTAATCATCTGCTGTACTGGGCTCTGACCAAAATGGGCTTTGAAACCACAATGTTGGGAGGATATGTTTACATAACTCCAGTCAGCAAATATAGCAGTGAAATGGTCCACCTTCTAGTACAGGTGACCATCAGTGACAGGAAGTACATTGTGGATTCCGCCTATGGAGGCTCCTACCAGATGTGGGAGCCTCTGGAATTAACATCTGGGAAGGATCAGCCTCAGGTGCCTGCCATCTTCCTTTTGACAGAGGAGAATGGAACCTGGTACTTGGACCAAATCAGAAGAGAGCAGTATGTTCCAAATGAAGAATTTGTTAACTCAGACCTCCTTGAAAAGAACAAATATCGAAAAATCTACTCCTTTACTCTTGAGCCCCGAGTTATCGAGGATTTTGAATATGTGAATAGCTATCTTCAGACATCGCCAGCATCTGTGTTTGTAAGCACATCGTTCTGTTCCTTGCAGACCTCGGAAGGGGTTCACTGTTTAGTGGGCTCCACCTTTACAAGTAGGAGATTCAGCTATAAGGACGATGTAGATCTGGTTGAGTTTAAATATGTGAATGAGGAAGAAATAGAAGATGTACTGAAAACCGCATTTGGCATTTCTTTGGAGAGAAAGTTTGTGCCCAAACATGGTGAACTAGTTTTTACTATT")
