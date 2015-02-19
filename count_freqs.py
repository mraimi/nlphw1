#! /usr/bin/python

__author__="Daniel Bauer <bauer@cs.columbia.edu>"
__date__ ="$Sep 12, 2011"
__rare__ ="_RARE_"
__STOP__ = ("STOP")
__STAR__ = ("*") 
__tag_set__ = ("I-LOC", "B-LOC", "I-ORG", "B-ORG", "I-PER", "O", "I-MISC", "B-MISC");

import sys
from collections import defaultdict
import math

"""
Count n-gram frequencies in a CoNLL NER data file and write counts to
stdout. 
"""

def simple_conll_corpus_iterator(corpus_file):
    """
    Get an iterator object over the corpus file. The elements of the
    iterator contain (word, ne_tag) tuples. Blank lines, indicating
    sentence boundaries return (None, None).
    """
    l = corpus_file.readline()
    while l:
        line = l.strip()
        if line: # Nonempty line
            # Extract information from line.
            # Each line has the format
            # word pos_tag phrase_tag ne_tag
            fields = line.split(" ")
            ne_tag = fields[-1]
            #phrase_tag = fields[-2] #Unused
            #pos_tag = fields[-3] #Unused
            word = " ".join(fields[:-1])
            yield word, ne_tag
        else: # Empty line
            yield (None, None)                        
        l = corpus_file.readline()

def sentence_iterator(corpus_iterator):
    """
    Return an iterator object that yields one sentence at a time.
    Sentences are represented as lists of (word, ne_tag) tuples.
    """
    current_sentence = [] #Buffer for the current sentence
    for l in corpus_iterator:        
            if l==(None, None):
                if current_sentence:  #Reached the end of a sentence
                    yield current_sentence
                    current_sentence = [] #Reset buffer
                else: # Got empty input stream
                    sys.stderr.write("WARNING: Got empty input file/stream.\n")
                    raise StopIteration
            else:
                current_sentence.append(l) #Add token to the buffer

    if current_sentence: # If the last line was blank, we're done
        yield current_sentence  #Otherwise when there is no more token
                                # in the stream return the last sentence.

def get_ngrams(sent_iterator, n):
    """
    Get a generator that returns n-grams over the entire corpus,
    respecting sentence boundaries and inserting boundary tokens.
    Sent_iterator is a generator object whose elements are lists
    of tokens.
    """
    for sent in sent_iterator:
         #Add boundary symbols to the sentence
         w_boundary = (n-1) * [(None, "*")]
         w_boundary.extend(sent)
         w_boundary.append((None, "STOP"))
         #Then extract n-grams
         ngrams = (tuple(w_boundary[i:i+n]) for i in xrange(len(w_boundary)-n+1))
         for n_gram in ngrams: #Return one n-gram at a time
            yield n_gram        


class Hmm(object):
    """
    Stores counts for n-grams and emissions. 
    """

    def __init__(self, n=3):
        assert n>=2, "Expecting n>=2."
        self.n = n
        self.emission_counts = defaultdict(int)
        self.ngram_counts = [defaultdict(int) for i in xrange(self.n)]
        self.all_states = set()
        self.emiss_prob_log = defaultdict(float)
        self.emiss_prob = dict()
        self.emiss_high_prob = dict()
        self.simple_counts = defaultdict(int)
        self.tri_ests = dict()
        

    def train(self, corpus_file):
        """
        Count n-gram frequencies and emission probabilities from a corpus file.
        """
        ngram_iterator = \
            get_ngrams(sentence_iterator(simple_conll_corpus_iterator(corpus_file)), self.n)

        for ngram in ngram_iterator:
            #Sanity check: n-gram we get from the corpus stream needs to have the right length
            assert len(ngram) == self.n, "ngram in stream is %i, expected %i" % (len(ngram, self.n))

            tagsonly = tuple([ne_tag for word, ne_tag in ngram]) #retrieve only the tags            
            for i in xrange(2, self.n+1): #Count NE-tag 2-grams..n-grams
                self.ngram_counts[i-1][tagsonly[-i:]] += 1
            
            if ngram[-1][0] is not None: # If this is not the last word in a sentence
                self.ngram_counts[0][tagsonly[-1:]] += 1 # count 1-gram
                self.emission_counts[ngram[-1]] += 1 # and emission frequencies

            # Need to count a single n-1-gram of sentence start symbols per sentence
            if ngram[-2][0] is None: # this is the first n-gram in a sentence
                self.ngram_counts[self.n - 2][tuple((self.n - 1) * ["*"])] += 1

    def write_counts(self, output, printngrams=[1,2,3]):
        """
        Writes counts to the output file object.
        Format:

        """
        # First write counts for emissions
        for word, ne_tag in self.emission_counts:            
            output.write("%i WORDTAG %s %s\n" % (self.emission_counts[(word, ne_tag)], ne_tag, word))
            #print str(self.emission_counts[(word, ne_tag)])+ " " + ne_tag + " " + word
            print self.ngram_counts
        # Then write counts for all ngrams
        for n in printngrams:            
            for ngram in self.ngram_counts[n-1]:
                ngramstr = " ".join(ngram)
                output.write("%i %i-GRAM %s\n" %(self.ngram_counts[n-1][ngram], n, ngramstr))

    def read_counts(self, corpusfile):

        self.n = 3
        self.emission_counts = defaultdict(int)
        self.ngram_counts = [defaultdict(int) for i in xrange(self.n)]
        self.all_states = set()

        for line in corpusfile:
            parts = line.strip().split(" ")
            count = float(parts[0])
            if parts[1] == "WORDTAG":
                ne_tag = parts[2]
                word = parts[3]
                self.emission_counts[(word, ne_tag)] = count
                self.all_states.add(ne_tag)
            elif parts[1].endswith("GRAM"):
                n = int(parts[1].replace("-GRAM",""))
                ngram = tuple(parts[2:])
                self.ngram_counts[n-1][ngram] = count
    
    def emission_gen(self):
        for emiss_word, emiss_tag in self.emission_counts:
            l = math.log(self.emission_counts[(emiss_word, emiss_tag)]/ float(self.ngram_counts[0][emiss_tag,]))
            self.emiss_prob_log[(emiss_word, emiss_tag)] = l
            p = self.emission_counts[(emiss_word, emiss_tag)]/ float(self.ngram_counts[0][emiss_tag,])
            self.emiss_prob[(emiss_word, emiss_tag)] = p
        #print self.emiss_prob

    def rare_replace(self):

        for emiss_word, emiss_tag in self.emission_counts:
            self.simple_counts[emiss_word] += self.emission_counts[(emiss_word, emiss_tag)]
        
        keys = self.emission_counts.keys()
        for emiss_word, emiss_tag in keys:
            if self.simple_counts[emiss_word] < 5:
                self.emission_counts[(__rare__, emiss_tag)] += self.emission_counts[(emiss_word, emiss_tag)]
                del self.emission_counts[(emiss_word, emiss_tag)]

    def rare_entity_tag(self):

        # First construct a dictionary, emiss_high_prob, that only contains the 
        # highest probability for any given key and its corresponding tag

        test_file = open('ner_dev.dat', 'r')
        for emiss_word, emiss_tag in self.emiss_prob_log:
            if emiss_word not in self.emiss_high_prob or self.emiss_prob_log[(emiss_word, emiss_tag)]>(self.emiss_high_prob[emiss_word])[1]:
                self.emiss_high_prob[emiss_word] = (emiss_tag, self.emiss_prob_log[(emiss_word, emiss_tag)])
        
        # Read through the test file, look up the 
        l = test_file.readline()
        while l:
            if l == "\n":
                #print ""
                l = test_file.readline()
                continue
            l = l.strip()
            '''
            if l not in self.emiss_high_prob:
                print l + " " + (self.emiss_high_prob[__rare__])[0]+ " " + str((self.emiss_high_prob[__rare__])[1])
            else:
                print l + " " + (self.emiss_high_prob[l.strip()])[0] + " " + str((self.emiss_high_prob[l.strip()])[1])
            '''
            l = test_file.readline()

        test_file.close()

    def trigram_estimate(self, y1, y2, y3):
        numerator = float((self.ngram_counts[2])[(y1,y2,y3)])
        #print self.ngram_counts[2]
        denominator = float((self.ngram_counts[1])[(y1,y2)])
        if denominator == 0:
            print "zero denominator error. exiting.\n"
            sys.exit()
        return float(numerator/denominator)

    def trigram_file_est(self, fname):
        test_file = open(fname, 'r')
        l = test_file.readline()
        while l:
            #print l
            l = l.strip()
            words = l.split()
            if words[1] != "3-GRAM":
                l = test_file.readline()
                continue
            # This is assuming the format as seen in ner_counts.dat e.g.
            # 244 3-GRAM I-PER I-PER I-PER
            self.tri_ests[(words[2], words[3], words[4])] \
            = self.trigram_estimate(words[2], words[3], words[4])
            #print words[2] + " " + words[3] + " " + words[4] + " " \
            #+ str(math.log(self.tri_ests[(words[2], words[3], words[4])]))
            l = test_file.readline()
        test_file.close()
    
    def init_ksets(self, n, k_sets):
        for x in xrange(0, n+1):
            if x == n:
                k_sets[(x)] = [__STOP__]
            else:
                k_sets[(x)] = __tag_set__
            #print k_sets
        #print self.k_sets

    def check_emiss(self, word, tag):
        if (word,tag) not in self.emiss_prob:
            if (__rare__,tag) in self.emiss_prob:
                return self.emiss_prob[(__rare__,tag)]
            else:
                return 0;
        else:
            return self.emiss_prob[(word, tag)]

    def check_trigram(self, tag1, tag2, tag3):
        if (tag1,tag2,tag3) not in self.tri_ests:
            return 0
        else:
            return self.tri_ests[(tag1, tag2, tag3)]


    def viterbi(self):
        test_file = open('ner_dev.dat', 'r')
        word = test_file.readline()
        s = ""
        while word:
            pi = dict()
            bp = dict()
            ts = dict()
            k_sets = dict()
            
            if word != "\n":
                s += " " + word.strip()
                word = test_file.readline()
                continue
            word = test_file.readline()
            words = s.strip().split()
            s = ""
            length = len(words)
            #print "length: " + str(length)
            self.init_ksets(length, k_sets)
            #print k_sets
            pi[(0,"*","*")] = 1;
            if length > 0:
                for v in k_sets[(0)]:
                    # print "pi of 0 is " + str(pi[(0,"*","*")])
                    # print v
                    # print "trigram est " + str(self.check_trigram("*","*",v))
                    # print "word: " + words[0]

                    # if (words[0],v) not in self.emiss_prob:
                    #     emiss = self.emiss_prob[(__rare__,v)]
                    # else:
                    #     emiss = self.emiss_prob[(words[0],v)]
                    # print "emission prob " + str(self.check_emiss(words[0],v))
                    pi[(1,"*",v)] = pi[(0,"*","*")]*self.check_trigram("*","*", v)*self.check_emiss(words[0],v)
            if length > 1:
                for u in k_sets[(0)]:
                        for v in k_sets[(1)]:
                            pi[(2,u,v)] = pi[(1,"*",u)]*self.check_trigram("*",u,v)*self.check_emiss(words[1],v)
                            bp[(2,u,v)] = "*"
            if length > 2:    
                for k in xrange(3,length):
                    for u in k_sets[(k-1)]:
                        for v in k_sets[(k)]:
                            max_w = None
                            max_p = 0.0
                            for w in k_sets[(k-2)]:
                                if k == 2:
                                    w = "*"
                                # if (k-1, w, u) not in pi:
                                #     print str((k-1, w, u)) + " Not in PI!!!!"
                                #     sys.exit()
                                # print "prob is " + str(float(pi[(k-1,w,u)]*self.check_trigram(w,u,v)*self.check_emiss(words[k],v)))
                                # print "k is " + str(k)
                                # print "word is " + words[k]
                                # print "w is " + w
                                # print "pi of k -1 " + str(pi[(k-1,w,u)])
                                # print "trigram est " + str(self.check_trigram(w,u,v))
                                # print "emission " + str(self.check_emiss(words[k],v))
                                if pi[(k-1,w,u)]*self.check_trigram(w,u,v)*self.check_emiss(words[k],v) > max_p:
                                    # print "FOUND MAX FOR " + str((w, u, v))
                                    max_w = w
                                    max_p = pi[(k-1,w,u)]*self.check_trigram(w,u,v)*self.check_emiss(words[k],v)
                            pi[(k,u,v)] = max_p
                            #print " k,u,v " + str(k) + " " + u + " " + v
                            bp[(k,u,v)] = max_w
                # for tag in pi:
                #     if tag[0] == 2:
                #         print tag
                #         print pi[tag]

              #  print "size of pi is " + str(len(pi))
                max_u = None
                max_v = None
                max_p = 0.0    
                for u in k_sets[(length-2)]:
                    for v in k_sets[(length-1)]:
                        #print "value to compare: " + str(pi[(length),u,v]*self.tri_ests[(__STOP__,u,v)])
                        # print "pi_n u v " + str(pi[(length-1),u,v])
                        # print "trigram stop,u,v " + str(self.check_trigram(u,v,__STOP__))

                        if pi[(length-1),u,v]*self.check_trigram(u,v,__STOP__) > max_p:
                       #     print "FOUND A MAX"
                            max_u = u
                            max_v = v
                            #max_p = pi[(length-1),u,v]*self.check_trigram(__STOP__,u,v)
              #  print "max_u " + max_u
              #  print "max_v " + max_v
                ts[(length-2)] = max_u
                ts[(length-1)] = max_v

                #print ts[(length-1)]
                #sys.exit()
                #print bp
                for k in range(length-3,-1,-1):
                    #print str(k)
                    ts[(k)] = bp[(k+2), ts[(k+1)], ts[(k+2)]]


               
                for i in xrange (0,length):
                    if i == 0:
                        print words[i] + " " + str(ts[(i)]) + " " + str(pi[(i,__STAR__,__STAR__)])
                    elif i == 1:
                        print words[i] + " " + str(ts[(i)]) + " " + str(pi[(i,__STAR__,ts[(i)])])
                    else:
                        print words[i] + " " + str(ts[(i)]) + " " + str((pi[(i,ts[(i-1)],ts[(i)])]))

def usage():
    print """
    python count_freqs.py [input_file] > [output_file]
        Read in a named entity tagged training input file and produce counts.
    """

if __name__ == "__main__":

    if len(sys.argv)!=2: # Expect exactly one argument: the training data file
        usage()
        sys.exit(2)

    try:
        input = file(sys.argv[1],"r")
    except IOError:
        sys.stderr.write("ERROR: Cannot read inputfile %s.\n" % arg)
        sys.exit(1)
    
    # Initialize a trigram counter
    counter = Hmm(3)
    
    # Collect counts
    counter.train(input)
    

    #question 4
    # Replace rare words
    counter.rare_replace()
    
    # Collect emission probabilities
    counter.emission_gen()
    
    # Use rare probabilities to tag entities
    counter.rare_entity_tag()
    

    # question 5a
    # trigram estimation
    # counter.trigram_estimate('I-ORG', 'I-ORG', 'I-PER')
    
    # trigram estimation from file
    counter.trigram_file_est('ner_counts.dat')

    #question 5b
    counter.viterbi()

    # Write the counts
    # counter.write_counts(sys.stdout)
