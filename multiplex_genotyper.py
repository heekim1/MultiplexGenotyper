
import sys
import unittest

class Sample(object):
    def __init__(self,id,genotype):
        self.id = id
        self.genotype = genotype
    def __eq__(self, other):
        return self.id == other.id
    def __hash__(self):
        return hash(self.id)


class MultiplexGenotyper():
    '''Data is list of "genotype,sample_id_1,sample_id_2,...,sample_id_n".
    Samples starting with "NORM" shall be mapped to the list of "mapped_norm", while smaples
    starting with "MUT" shall be initiated to the list of "unmapped_mut". The samples in "unmapped_mut"
    shall be mapped to the list of "mapped_mut". If the list of "unmapped_mut" is [], all of samples
    are mapped to either "NORM" or "MUT". If any sample maps to both "NORM" and "MUT", it retunrs
    "INCONSISTENT". If the list of "unmapped_mut" is not empty, it returns "NONUNIQUE".
    '''
    def __init__(self,data):
        self.inconsistent = False
        self.nonunique = False
        self.mapped_norm = set()
        self.mapped_mut = set()
        self.unmapped_mut = []
        self.data = data

    def _set_norm_and_unmapped_mut(self):
        for line in self.data:
            cols = line.split(",")
            if cols[0] == 'NORM':
                for sample_id in cols[1:]:
                    self.mapped_norm.add(Sample(int(sample_id),'NORM'))
            else:
                self.unmapped_mut.append(sorted(cols[1:],key=lambda x: int(x)))

    def _get_mapped_mut(self):
        '''smaples in unmmaped_mut are evaluated in loopinh until any of samples has inconsistancy or 
           any of nonunique samples is detected.
        '''
        while not self.inconsistent:
            new_unmapped_mut = []
            for sample_list in sorted(self.unmapped_mut,key=len):
                if len(sample_list) == 1:
                    sample = Sample(int(sample_list[0]),'MUT')
                    self.mapped_mut.add(sample)
                    # check INCONSISTENT
                    if sample in self.mapped_norm:
                        self.inconsistent = True
                        break
                else:
                    # check INCONSISTENT. if sample exists in both normal and mutant, 
                    # sample is inconsistent.
                    inconsistent_sample_list = filter(lambda x: Sample(int(x),'NORM') in self.mapped_norm and \
                                                    Sample(int(x),'MUT') in self.mapped_mut, sample_list)
                    if len(inconsistent_sample_list):
                        self.inconsistent = True
                        break

                    # remove sample if sample exists in mapped_norm or mapped_mut
                    reduce_sample_list = filter(lambda x: Sample(int(x),'NORM') not in self.mapped_norm and \
                                                Sample(int(x),'MUT') not in self.mapped_mut, sample_list)

                    # check if sample(s) in unmapped_mut has been resolved 
                    if len(reduce_sample_list) > 0:
                        new_unmapped_mut.append(reduce_sample_list)
            
            # exit loop if unmapped_mut == [] or no change from the previous unmapped_mut.
            if self.unmapped_mut == [] or self.unmapped_mut == new_unmapped_mut:
                if len(self.unmapped_mut)  > 0:
                    self.nonunique = True
                break
            self.unmapped_mut = new_unmapped_mut

    def _union_sort_print(self):
        for obj in sorted(self.mapped_norm.union(self.mapped_mut),key=lambda x: x.id):
            print "%d,%s" % (obj.id,obj.genotype)


    def _get_report(self):
        if self.inconsistent:
            print "INCONSISTENT"
        elif self.nonunique:
            print "NONUNIQUE"
        else:
            print "MUT COUNT: ",len(self.mapped_mut)
            print "NORM COUNT: ",len(self.mapped_norm)
            self._union_sort_print()
            
    def run(self):
        self._set_norm_and_unmapped_mut()
        self._get_mapped_mut()
        self._get_report()


def get_cases():
    with open(sys.argv[1],'r') if len(sys.argv) > 1 else sys.stdin as fp:
        test_case = []
        for line in fp:
            line = line.strip()
            if line == '':
                yield test_case
                test_case = []
            else:
                cols = line.split(",")
                test_case.append(line)

def main():
    for test_case in get_cases():
        mp = MultiplexGenotyper(test_case)
        mp.run()
        print ""

class TestMultiplexGenotyper(unittest.TestCase):
    def setUp(self):
        self.test_case = ['NORM,11,44,118,11', 'MUT,49,83', 'NORM,122,43,98', 
                          'MUT,83', 'NORM,146,122,121', 'MUT,49', 'NORM,11,100,50']
        self.mg = MultiplexGenotyper(self.test_case)
    def test_set_norm_and_unmapped_mut(self):
        self.mg._set_norm_and_unmapped_mut()
        assert(len(self.mg.mapped_norm) == 10)
    def test_get_mapped_mut(self):
        self.mg._set_norm_and_unmapped_mut()
        self.mg._get_mapped_mut()
        print "len(self.mg.mapped_mut) : ", len(self.mg.mapped_mut)
        assert(len(self.mg.mapped_mut) == 2)
    def test_run(self):
        self.mg._set_norm_and_unmapped_mut()
        self.mg._get_mapped_mut()
        assert(self.mg.inconsistent == False and self.mg.nonunique == False)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMultiplexGenotyper)
    #unittest.TextTestRunner(verbosity=2).run(suite)
    main()